(* ::Package:: *)

BeginPackage["NelderMeadMinimize`"];

ClearAll[NelderMeadMinimize];
Attributes[NelderMeadMinimize] = {HoldAll};
Options[NelderMeadMinimize] = {
   "InitialPoints" -> Automatic, "RandomSeed" -> Sequence[0, Method -> "MKL"],
   "ConvergenceTolerance" -> $MachineEpsilon, MaxIterations -> Infinity,
   "HistoryLength" -> Automatic,
   "ReflectRatio" -> 1, "ExpandRatio" -> 2, "ContractRatio" -> 1/2, "ShrinkRatio" -> 1/2,
   CompilationTarget :> $CompilationTarget
  };

ClearAll[NelderMeadMinimize`Dump`CompiledNelderMead];
Options[NelderMeadMinimize`Dump`CompiledNelderMead] = {
   "HistoryLength" -> Automatic,
   "ReflectRatio" -> 1, "ExpandRatio" -> 2, "ContractRatio" -> 1/2, "ShrinkRatio" -> 1/2,
   CompilationTarget :> $CompilationTarget, "ReturnValues" -> "AugmentedOptimizedParameters"
  };

ClearAll[NelderMeadMinimize`Dump`Rosenbrock];

ClearAll[NelderMeadMinimize`Dump`Hyperellipsoid];

Begin["`Private`"];

(* Caller for the compiled Nelder-Mead code *)
NelderMeadMinimize[
   minimand_, vars : {__Symbol},
   opts : OptionsPattern[NelderMeadMinimize]
  ] :=
  BlockRandom@Block[vars,
    SeedRandom@OptionValue["RandomSeed"];
    Module[{
      dimensions = Length[vars], initialPoints,
      objectiveFunction = Switch[Head[minimand],
        Function | CompiledFunction, minimand,
        _, Function[vars, minimand]
       ],
      compiledMinimizer,
      maxIterations = If[# === Infinity, -1, #] & @ OptionValue["MaxIterations"]
     },
     compiledMinimizer = NelderMeadMinimize`Dump`CompiledNelderMead[
       objectiveFunction, vars,
       "HistoryLength" -> OptionValue["HistoryLength"],
       "ReflectRatio" -> OptionValue["ReflectRatio"], "ExpandRatio" -> OptionValue["ExpandRatio"],
       "ContractRatio" -> OptionValue["ContractRatio"], "ShrinkRatio" -> OptionValue["ShrinkRatio"],
       CompilationTarget -> OptionValue[CompilationTarget], "ReturnValues" -> "AugmentedOptimizedParameters"
      ];
     initialPoints = Cases[OptionValue["InitialPoints"] // N, {Repeated[_Real, {dimensions}]}];
     Block[{results},
      results = compiledMinimizer[
       initialPoints ~Join~ RandomReal[{0, 1}, {dimensions - Length[initialPoints] + 1, dimensions}],
       OptionValue["ConvergenceTolerance"], maxIterations
      ];
      {First[results], Thread[vars -> Rest[results]]}
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
apply[f : _Function | _CompiledFunction, vars : {__Symbol}] :=
  Function[arglist, Block[vars, vars = arglist; f @@ vars]];

(* Version-specific compilable Ordering *)
ClearAll[ordering];
ordering = If[$VersionNumber >= 8, Ordering,
   (* The Mathematica 6 and 7 runtimes don't have support for Ordering *)
   Function[{lst}, IntegerPart@Last@Transpose@Sort@Transpose[{lst, Range@Length[lst]}]]
  ];

(* Inlinable utility function *)
ClearAll[cumulativeAbsoluteDifferences];
cumulativeAbsoluteDifferences = If[$VersionNumber >= 8,
   Function[{lst}, Total@Abs@Differences[lst]/Length[lst]],
   Function[{lst}, Total@Abs[Most[lst] - Rest[lst]]/Length[lst]]
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
     totalAbsDiffs = cumulativeAbsoluteDifferences,
     (* Return values *)
     return = Switch[OptionValue["ReturnValues"],
       "OptimizedParameters", Function[Null, simplex[[best]]],
       "AugmentedOptimizedParameters", Function[Null, simplex[[best]] ~Prepend~ vals[[best]]],
       "Simplex", Function[Null, simplex],
       "AugmentedSimplex", Function[Null, Transpose[Transpose[simplex] ~Prepend~ vals]],
       All, Function[Null,
        results = {vals, simplex};(* This requires a call out of the VM *)
        {evaluations, reflections, expansions, contractions, shrinkages}
       ],
       _, Function[Null]
      ],
     (* Options to be passed to Compile *)
     compileopts = Sequence @@ If[$VersionNumber >= 8, {
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
    Compile[{{pts, _Real, 2}, {tol, _Real, 0}, {maxit, _Integer, 0}},
     Block[{
       (* Housekeeping *)
       history = Table[infinity, {historyLength}], iteration = maxit,
       (* Basic quantities *)
       simplex = pts, vals = f /@ pts, order,
       (* Calculated points and function values *)
       centroid = origin,
       reflectedPoint = origin, reflectedValue = infinity,
       expandedPoint = origin, expandedValue = infinity,
       contractedPoint = origin, contractedValue = infinity,
       (* More readable indices into the simplex array *)
       best = 1, worst = -1, rest = Range[2, Length[pts]],
       (* Operation counts (for debugging purposes) *)
       evaluations = Length[pts],
       reflections = 0, expansions = 0, contractions = 0, shrinkages = 0
      },
      While[
       (* Order simplex points by function value *)
       order = ordering[vals];
       vals = vals[[order]]; simplex = simplex[[order]];
       (* Decrement and test iterator *)
       (iteration--) != 0,
       (* Check for convergence *)
       history[[1]] = vals[[best]]; history = RotateLeft[history];
       If[totalAbsDiffs[history] <= tol + epsilon totalAbsDiffs[history],
        Break[]
       ];
       (* Find centroid of best (N - 1) points *)
       centroid = Mean@Most[simplex];
       (* Reflect *)
       reflectedPoint = centroid + reflectRatio (centroid - simplex[[worst]]);
       reflectedValue = f[reflectedPoint]; ++evaluations;
       If[vals[[best]] <= reflectedValue < vals[[-2]],
        vals[[worst]] = reflectedValue; simplex[[worst]] = reflectedPoint;
        ++reflections; Continue[]
       ];
       (* Expand *)
       If[reflectedValue < vals[[best]],
        expandedPoint = centroid + expandRatio (reflectedPoint - centroid);
        expandedValue = f[expandedPoint]; ++evaluations;
        If[expandedValue < reflectedValue,
         vals[[worst]] = expandedValue; simplex[[worst]] = expandedPoint;
         ++expansions; Continue[],
         vals[[worst]] = reflectedValue; simplex[[worst]] = reflectedPoint;
         ++reflections; Continue[]
        ];
       ];
       (* Contract *)
       If[reflectedValue < vals[[worst]],
        (* Outside contraction *)
        contractedPoint = centroid + contractRatio (reflectedPoint - centroid);
        contractedValue = f[contractedPoint]; ++evaluations;
        If[contractedValue <= reflectedValue,
         vals[[worst]] = contractedValue; simplex[[worst]] = contractedPoint;
         ++contractions; Continue[]
        ];,
        (* Inside contraction *)
        contractedPoint = centroid - contractRatio (centroid - simplex[[worst]]);
        contractedValue = f[contractedPoint]; ++evaluations;
        If[contractedValue < vals[[worst]],
         vals[[worst]] = contractedValue; simplex[[worst]] = contractedPoint;
         ++contractions; Continue[]
        ];
       ];
       (* Shrink *)
       simplex[[rest]] = simplex[[best]] + shrinkRatio (simplex[[rest]] - simplex[[best]]);
       vals[[rest]] = f /@ simplex[[rest]]; evaluations += Length[rest] - 1;
       ++shrinkages;
      ];
      return[]
     ], compileopts
    ]
   ];

(* Unimodal convex function; minimum at {1, 1, ...} *)
NelderMeadMinimize`Dump`Hyperellipsoid = 
  Block[{x = {##}}, Plus @@ (Range@Length[x] (x - 1)^2)] &;

(* Non-convex function; global minimum at {1, 1, ...} with secondary local minimum near {-1, 1, 1, ...} *)
NelderMeadMinimize`Dump`Rosenbrock = 
  Block[{x = Most[{##}], y = Rest[{##}]}, Plus @@ (100 (y - x^2)^2 + (1 - x)^2)] &;

End[];

EndPackage[];
