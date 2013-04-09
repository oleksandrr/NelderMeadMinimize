(* ::Package:: *)

BeginPackage["NelderMeadMinimize`"];

ClearAll[NelderMeadMinimize];
SetAttributes[NelderMeadMinimize, HoldAll];
Options[NelderMeadMinimize] = {
    "InitialPoints" -> Automatic, "RandomSeed" -> Sequence[0, Method -> "MKL"],
    "ConvergenceTolerance" -> $MachineEpsilon, "HistoryLength" -> Automatic, MaxIterations -> Infinity,
    "ReflectRatio" -> 1, "ExpandRatio" -> 2, "ContractRatio" -> 1/2, "ShrinkRatio" -> 1/2,
    "Diagnostics" -> False, CompilationTarget :> $CompilationTarget
   } // Sort;

ClearAll[NelderMeadMinimize`Dump`CompiledNelderMead];
Options[NelderMeadMinimize`Dump`CompiledNelderMead] = {
    "HistoryLength" -> Automatic,
    "ReflectRatio" -> 1, "ExpandRatio" -> 2, "ContractRatio" -> 1/2, "ShrinkRatio" -> 1/2,
    "ReturnValues" -> "AugmentedOptimizedParameters",
    "Diagnostics" -> False, CompilationTarget :> $CompilationTarget
   } // Sort;

(* This value will be set on return if diagnostics are enabled *)
ClearAll[NelderMeadMinimize`Dump`OperationCounts];

(* Separable and nonseparable hyperellipsoid functions *)
ClearAll[NelderMeadMinimize`Dump`Hyperellipsoid];
ClearAll[NelderMeadMinimize`Dump`RotatedHyperellipsoid];

(* Rosenbrock's function *)
ClearAll[NelderMeadMinimize`Dump`Rosenbrock];

Begin["`Private`"];

(* Caller for the compiled Nelder-Mead code *)
NelderMeadMinimize[
   minimand_, vars : {varseq__Symbol},
   opts : OptionsPattern[NelderMeadMinimize]
  ] :=
  BlockRandom@Block[{varseq, dimension = Length[vars]},
    SeedRandom@OptionValue["RandomSeed"];
    With[{
      compiledMinimizer = NelderMeadMinimize`Dump`CompiledNelderMead[
        Switch[Head[minimand],
         Function | CompiledFunction, minimand,
         _, Function[vars, minimand]
        ], vars,
        Sort@FilterRules[
          {opts}~Join~FilterRules[Options[NelderMeadMinimize], Except[{opts}]],
          Options[NelderMeadMinimize`Dump`CompiledNelderMead]
         ], "ReturnValues" -> "AugmentedOptimizedParameters"
       ],
      initialPoints = Cases[
        OptionValue["InitialPoints"] // N,
        List@Repeated[_Real, {dimension}],
        {0, 1}, dimension + 1
       ]
     },
     With[{
       result = compiledMinimizer[
         initialPoints~Join~RandomReal[{0, 1}, {dimension + 1 - Length[initialPoints], dimension}],
         OptionValue["ConvergenceTolerance"], If[# === Infinity, -1, #] & @ OptionValue["MaxIterations"]
        ]
      },
      {First[result], Thread[vars -> Rest[result]]}
     ]
    ]
   ];

(* Deal with univariate case *)
NelderMeadMinimize[
   minimand_, var_Symbol,
   opts : OptionsPattern[NelderMeadMinimize]
  ] :=
  NelderMeadMinimize[minimand, {var}, opts];

(* Produces inlinable code for use inside Compile (where Apply is not supported directly) *)
ClearAll[apply];
SetAttributes[apply, HoldRest];
apply[f : _Function | _CompiledFunction, vars : {___Symbol}] :=
  Function[args, Block[vars, vars = args; f @@ vars]];

(* Version-specific compilable Ordering *)
ClearAll[ordering];
ordering = If[$VersionNumber >= 8, Ordering,
   (* The Mathematica 6 and 7 runtimes don't have support for Ordering *)
   Function[{lst}, IntegerPart@Last@Transpose@Sort@Transpose[{lst, Range@Length[lst]}]]
  ];

(* Inlinable utility function *)
ClearAll[meanAbsoluteDifferences];
meanAbsoluteDifferences = If[$VersionNumber >= 8,
   Function[{lst}, Mean@Abs@Differences[lst]],
   Function[{lst}, Plus @@ Abs[Most[lst] - Rest[lst]]/(Length[lst] - 1)]
  ];

(* Produces compiled code for the Nelder-Mead algorithm with the objective function inlined *)
NelderMeadMinimize`Dump`CompiledNelderMead[
   objectiveFunction : _Function | _CompiledFunction, vars : {__Symbol},
   opts : OptionsPattern[NelderMeadMinimize`Dump`CompiledNelderMead]
  ] :=
  NelderMeadMinimize`Dump`CompiledNelderMead[
    Verbatim[objectiveFunction], vars,
    opts
   ] =
   With[{
     (* Inlined option values *)
     historyLength = If[# === Automatic, 10 Length[vars], #] & @ OptionValue["HistoryLength"],
     reflectRatio = OptionValue["ReflectRatio"], expandRatio = OptionValue["ExpandRatio"],
     contractRatio = OptionValue["ContractRatio"], shrinkRatio = OptionValue["ShrinkRatio"],
     (* Other inlined values *)
     origin = ConstantArray[0., Length[vars]],
     infinity = $MaxMachineNumber,
     epsilon = $MachineEpsilon,
     (* Inlined functions *)
     f = apply[objectiveFunction, vars],
     ordering = ordering,
     meanAbsoluteDifferences = meanAbsoluteDifferences,
     (* Return values *)
     return = Switch[OptionValue["ReturnValues"],
       "OptimizedParameters", Function[Null, simplex[[best]]],
       "AugmentedOptimizedParameters", Function[Null, simplex[[best]]~Prepend~values[[best]]],
       "Simplex", Function[Null, simplex],
       "AugmentedSimplex", Function[Null, Transpose[Transpose[simplex]~Prepend~values]],
       _, Function[Null]
      ],
     (* Code transformations to be applied before compilation *)
     transformations = If[TrueQ@OptionValue["Diagnostics"],
       HoldPattern[h_[pre___, diagnostic[diags___], post___]] :> h[pre, diags, post],
       HoldPattern[h_[pre___, __diagnostic, post___]] :> h[pre, post]
      ],
     (* Options passed to Compile *)
     compileOpts = Sequence @@ If[$VersionNumber >= 8, {
         (* Mathematica 8 and above offer improved behaviour using these options *)
         CompilationOptions -> {"ExpressionOptimization" -> True, "InlineCompiledFunctions" -> True},
         RuntimeOptions -> {"CompareWithTolerance" -> False, "EvaluateSymbolically" -> False},
         CompilationTarget -> OptionValue[CompilationTarget]
        }, {
         (* Undocumented optimization level option for Mathematica 6 and 7 *)
         "CompileOptimizations" -> All
        }
       ]
    },
    Compile @@ (
      Hold[{{initialSimplex, _Real, 2}, {tolerance, _Real, 0}, {maxIterations, _Integer, 0}},
        Block[{
          (* Housekeeping *)
          history = Table[infinity, {historyLength}], iteration = maxIterations,
          (* Basic quantities *)
          simplex = initialSimplex, values = f /@ initialSimplex, order,
          (* Calculated points and function values *)
          centroid = origin,
          reflectedPoint = origin, reflectedValue = infinity,
          expandedPoint = origin, expandedValue = infinity,
          contractedPoint = origin, contractedValue = infinity,
          (* More readable indices into the simplex array *)
          best = 1, worst = -1, rest = Range[2, Length[initialSimplex]],
          diagnostic[
           (* Operation counts *)
           evaluations = Length[initialSimplex],
           reflections = 0, expansions = 0, contractions = 0, shrinkages = 0
          ]
         },
         While[
          (* Order simplex points by function value *)
          order = ordering[values];
          values = values[[order]]; simplex = simplex[[order]];
          (* Decrement and test iterator *)
          (iteration--) != 0,
          (* Check for convergence *)
          history[[1]] = values[[best]]; history = RotateLeft[history];
          If[meanAbsoluteDifferences[history] <= tolerance + epsilon meanAbsoluteDifferences[history],
           Break[]
          ];
          (* Find centroid of best (N - 1) points *)
          centroid = Mean@Most[simplex];
          (* Reflect *)
          reflectedPoint = centroid + reflectRatio (centroid - simplex[[worst]]);
          reflectedValue = f[reflectedPoint]; diagnostic[++evaluations];
          If[values[[best]] <= reflectedValue < values[[-2]],
           values[[worst]] = reflectedValue; simplex[[worst]] = reflectedPoint;
           diagnostic[++reflections]; Continue[]
          ];
          (* Expand *)
          If[reflectedValue < values[[best]],
           expandedPoint = centroid + expandRatio (reflectedPoint - centroid);
           expandedValue = f[expandedPoint]; diagnostic[++evaluations];
           If[expandedValue < reflectedValue,
            values[[worst]] = expandedValue; simplex[[worst]] = expandedPoint;
            diagnostic[++expansions]; Continue[],
            values[[worst]] = reflectedValue; simplex[[worst]] = reflectedPoint;
            diagnostic[++reflections]; Continue[]
           ]
          ];
          (* Contract *)
          If[reflectedValue < values[[worst]],
           (* Outside contraction *)
           contractedPoint = centroid + contractRatio (reflectedPoint - centroid);
           contractedValue = f[contractedPoint]; diagnostic[++evaluations];
           If[contractedValue <= reflectedValue,
            values[[worst]] = contractedValue; simplex[[worst]] = contractedPoint;
            diagnostic[++contractions]; Continue[]
           ],
           (* Inside contraction *)
           contractedPoint = centroid - contractRatio (centroid - simplex[[worst]]);
           contractedValue = f[contractedPoint]; diagnostic[++evaluations];
           If[contractedValue < values[[worst]],
            values[[worst]] = contractedValue; simplex[[worst]] = contractedPoint;
            diagnostic[++contractions]; Continue[]
           ]
          ];
          (* Shrink *)
          simplex[[rest]] = simplex[[best]] + shrinkRatio (simplex[[rest]] - simplex[[best]]);
          values[[rest]] = f /@ simplex[[rest]]; diagnostic[evaluations += Length[rest]];
          diagnostic[++shrinkages];
         ];
         diagnostic[
          NelderMeadMinimize`Dump`OperationCounts = {
            "Function evaluations" -> evaluations,
            "Reflections" -> reflections, "Expansions" -> expansions,
            "Contractions" -> contractions, "Shrinkages" -> shrinkages
           }
         ];
         return[]
        ], compileOpts
       ] //. transformations
     )
   ];

(* Unimodal convex separable function; minimum at {1, 1, ...} *)
NelderMeadMinimize`Dump`Hyperellipsoid =
  With[{x = {##}}, Plus @@ (Range@Length[x] (x - 1)^2)] &;

(* Unimodal convex nonseparable function; minimum at {1, 1, ...} *)
NelderMeadMinimize`Dump`RotatedHyperellipsoid =
  With[{x = {##}}, Plus @@ (FoldList[Plus, 0, x - 1]^2)] &;

(* Non-convex function; global minimum at {1, 1, ...} with secondary local minimum near {-1, 1, 1, ...} *)
NelderMeadMinimize`Dump`Rosenbrock =
  With[{x = Most[{##}], y = Rest[{##}]}, Plus @@ (100 (y - x^2)^2 + (1 - x)^2)] &;

End[];

EndPackage[];
