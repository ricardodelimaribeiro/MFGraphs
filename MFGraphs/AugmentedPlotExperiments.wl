(* ::Package:: *)

(* AugmentedPlotExperiments.wl
   Standalone comparison of three approaches to coloring the augmented
   infrastructure graph by u values, plus a 3D surface plot inspired by
   Cacace et al. (2017), Figure 3.

   Run cell-by-cell in the Mathematica front end, or evaluate the whole
   file with:
       Get["/path/to/MFGraphs/MFGraphs/AugmentedPlotExperiments.wl"]
*)

(* ------------------------------------------------------------------ *)
(* 1. LOAD PACKAGE                                                     *)
(* ------------------------------------------------------------------ *)

PrependTo[$Path, ParentDirectory @ DirectoryName @ $InputFileName];
Needs["MFGraphs`"];

(* ------------------------------------------------------------------ *)
(* 2. BUILD AND SOLVE A SCENARIO                                       *)
(*    Simple 3-vertex chain: 1 -> 2 -> 3                              *)
(*    entry at 1, exit at 3.  Swap in any other scenario below.       *)
(* ------------------------------------------------------------------ *)

scenario3v = makeScenario[<|
    "Model" -> <|
        "Vertices"  -> {1, 2, 3},
        "Adjacency" -> {{0,1,0},{0,0,1},{0,0,0}},
        "Entries"   -> {{1, 10}},
        "Exits"     -> {{3, 0}},
        "Switching" -> {}
    |>
|>];

sol3v = solveScenario[scenario3v];
sys3v = makeSystem[scenario3v, makeSymbolicUnknowns[scenario3v]];

(* Alias — swap these to experiment with other scenarios *)
s   = scenario3v;
sys = sys3v;
sol = sol3v;

(* ------------------------------------------------------------------ *)
(* 3. SHARED INFRASTRUCTURE                                            *)
(* ------------------------------------------------------------------ *)

augmented     = augmentAuxiliaryGraph[sys];
vertices      = augmented["Vertices"];
flowEdges     = augmented["FlowEdges"];
transEdges    = augmented["TransitionEdges"];
allAugEdges   = Join[flowEdges, transEdges];
edgeVariables = augmented["EdgeVariables"];
edgeKinds     = augmented["EdgeKinds"];
rules         = extractRules[sol];
auxEntryV     = systemData[sys, "AuxEntryVertices"];
auxExitV      = systemData[sys, "AuxExitVertices"];
exitCosts     = systemData[sys, "ExitCosts"];

(* Resolve u[a,b] with fallbacks for boundary vertices *)
getU[p_]          := (u @@ p) /. rules;
getEffectiveU[p_] := With[{uv = getU[p]},
    If[NumericQ[uv], uv,
        With[{cost = Lookup[exitCosts, p[[1]], Lookup[exitCosts, p[[2]], Missing[]]]},
            If[!MissingQ[cost], cost,
                If[MemberQ[auxEntryV, p[[1]]] || MemberQ[auxEntryV, p[[2]]],
                    With[{adj = getU[Reverse[p]]}, If[NumericQ[adj], adj, Missing[]]],
                    Missing[]]]]]];

uValues  = AssociationMap[getEffectiveU, vertices];
numericU = Select[uValues, NumericQ];
uMin     = Min[Values[numericU]];
uMax     = Max[Values[numericU]];

(* Blue (low) -> Red (high), matching the paper's colormap *)
uColorFn[t_] := Blend[{RGBColor[0.0,0.2,0.8], Cyan, Yellow, RGBColor[0.8,0.0,0.0]}, t];
uColorOf[v_] := If[NumericQ[uValues[v]],
    uColorFn[If[uMin == uMax, 0.5, Rescale[uValues[v], {uMin, uMax}]]],
    GrayLevel[0.75]];

(* Effective flow for thickness/opacity, reusing mfgAugmentedPlot logic *)
effectiveFlows = Association @ Map[
    With[{v1 = If[rules =!= {}, edgeVariables[#] /. rules /. _j -> 0, edgeVariables[#]]},
        # -> v1] &,
    allAugEdges];
numericJ = Cases[Values[effectiveFlows], _?NumericQ, Infinity];
maxJ     = Max[Append[numericJ, 1]];
flowThickness[e_] := AbsoluteThickness[
    If[NumericQ[effectiveFlows[e]], Rescale[effectiveFlows[e], {0, maxJ}, {1.5, 6}], 2.0]];
flowOpacity[e_] := Opacity[
    If[NumericQ[effectiveFlows[e]] && effectiveFlows[e] == 0, 0, 0.9]];

(* ------------------------------------------------------------------ *)
(* 4. APPROACH A — EdgeShapeFunction (stays inside Graph[])           *)
(*    Each directed edge is split into nSeg short Line segments,      *)
(*    each colored by the linearly interpolated u value.              *)
(*    Pros: keeps Graph interactivity.                                 *)
(*    Cons: Graph[] resamples pts[], so exact segment positions may   *)
(*          differ slightly from straight-line interpolation.         *)
(* ------------------------------------------------------------------ *)

nSeg = 30;  (* number of gradient segments per edge *)

plotA = Graph[vertices, allAugEdges,
    EdgeShapeFunction -> Function[{pts, e},
        With[{
            uS = uValues[e[[1]]],
            uT = uValues[e[[2]]],
            p1 = N[pts[[1]]],
            p2 = N[pts[[-1]]]
        },
            Join[
                (* Gradient segments *)
                Table[
                    With[{t0 = k/nSeg, t1 = (k+1)/nSeg},
                        {If[NumericQ[uS] && NumericQ[uT],
                            uColorFn[Rescale[(1-t0)*uS + t0*uT, {uMin, uMax}]],
                            GrayLevel[0.6]],
                         flowThickness[e], flowOpacity[e],
                         Line[{(1-t0)*p1 + t0*p2, (1-t1)*p1 + t1*p2}]}],
                    {k, 0, nSeg-1}],
                (* Arrowhead at 70% along the edge — inherits flow opacity *)
                With[{t = 0.7},
                    {{GrayLevel[0.3], flowOpacity[e], Arrowheads[{{0.025, 1}}],
                      Arrow[{(1-t)*p1 + t*p2, (1-(t+0.01))*p1 + (t+0.01)*p2}]}}]
            ]
        ]
    ],
    VertexStyle  -> (# -> uColorOf[#] & /@ vertices),
    VertexSize   -> 0.5,
    VertexLabels -> Map[# -> Placed[
        Style[StringReplace[ToString[#[[2]]], {"auxEntry" -> "in", "auxExit" -> "out"}],
              8, Black], Center] &, vertices],
    PlotLabel    -> Style["A: EdgeShapeFunction gradient", Bold, 12],
    ImageSize    -> 500
];

(* ------------------------------------------------------------------ *)
(* 5. APPROACH B — Explicit Graphics[]                                *)
(*    Extract coordinates from GraphEmbedding, then build Graphics    *)
(*    primitives directly.  Full control; no Graph interactivity.     *)
(*    Arrow direction drawn separately at low opacity.                *)
(* ------------------------------------------------------------------ *)

(* Derive a layout from a plain Graph render *)
layoutG = Graph[vertices, allAugEdges];
coords  = AssociationThread[vertices, GraphEmbedding[layoutG]];

gradientEdgePrims[e_, nSegB_:30] :=
    With[{
        v1 = e[[1]], v2 = e[[2]],
        p1 = N[coords[e[[1]]]], p2 = N[coords[e[[2]]]]
    },
        With[{uS = uValues[v1], uT = uValues[v2]},
            Table[
                With[{t0 = k/nSegB, t1 = (k+1)/nSegB},
                    {If[NumericQ[uS] && NumericQ[uT],
                        uColorFn[Rescale[(1-t0)*uS + t0*uT, {uMin, uMax}]],
                        GrayLevel[0.6]],
                     flowThickness[e], flowOpacity[e],
                     Line[{(1-t0)*p1 + t0*p2, (1-t1)*p1 + t1*p2}]}],
                {k, 0, nSegB-1}]]];

arrowPrims = Map[
    With[{p1 = N[coords[#[[1]]]], p2 = N[coords[#[[2]]]]},
        {GrayLevel[0.25], Opacity[0.5], Arrowheads[{{0.02, 0.65}}],
         Arrow[{p1, p2}]}] &,
    allAugEdges];

vertexPrimsB = Map[
    With[{p = coords[#]},
        {uColorOf[#], EdgeForm[{Thin, GrayLevel[0.3]}], Disk[p, 0.06]}] &,
    vertices];

labelPrimsB = Map[
    With[{p = coords[#]},
        Text[Style[StringReplace[ToString[#[[2]]], {"auxEntry" -> "in", "auxExit" -> "out"}],
                   7, Black, Bold], p]] &,
    vertices];

plotB = Graphics[
    {Flatten[Map[gradientEdgePrims, allAugEdges], 2],
     arrowPrims,
     vertexPrimsB,
     labelPrimsB},
    Frame     -> True,
    FrameTicks -> None,
    PlotLabel -> Style["B: Graphics[] gradient", Bold, 12],
    ImageSize -> 500,
    Background -> White
];

(* ------------------------------------------------------------------ *)
(* 6. APPROACH C — Graphics3D: network in x,y, u values as z          *)
(*    Inspired by Cacace et al. (2017), Figure 3 bottom panels.       *)
(*    - Network skeleton drawn in black at z = 0.                     *)
(*    - u value curves drawn in red at z = linear interpolation.      *)
(*    - Vertex spheres colored by u.                                  *)
(* ------------------------------------------------------------------ *)

nSeg3D = 60;

skeleton3D = Map[
    With[{p1 = N[coords[#[[1]]]], p2 = N[coords[#[[2]]]]},
        {Black, AbsoluteThickness[2],
         Line[{{p1[[1]], p1[[2]], 0}, {p2[[1]], p2[[2]], 0}}]}] &,
    allAugEdges];

uCurves3D = Map[
    With[{
        v1 = #[[1]], v2 = #[[2]],
        p1 = N[coords[#[[1]]]], p2 = N[coords[#[[2]]]]
    },
        With[{uS = uValues[v1], uT = uValues[v2]},
            If[NumericQ[uS] && NumericQ[uT],
                {Red, AbsoluteThickness[2.5],
                 Line[Table[
                     With[{t = k/nSeg3D},
                         {(1-t)*p1[[1]] + t*p2[[1]],
                          (1-t)*p1[[2]] + t*p2[[2]],
                          (1-t)*uS       + t*uT}],
                     {k, 0, nSeg3D}]]},
                {}]]] &,
    allAugEdges];

(* Vertical drop lines from u-curve to skeleton *)
dropLines3D = Map[
    With[{v = #, p = N[coords[#]]},
        If[NumericQ[uValues[v]],
            {GrayLevel[0.7], Dashed, AbsoluteThickness[0.5],
             Line[{{p[[1]], p[[2]], 0}, {p[[1]], p[[2]], uValues[v]}}]},
            {}]] &,
    vertices];

vertexSpheres3D = Map[
    With[{p = N[coords[#]]},
        {uColorOf[#], Sphere[{p[[1]], p[[2]], 0}, 0.04]}] &,
    vertices];

uSpheres3D = Map[
    With[{p = N[coords[#]]},
        If[NumericQ[uValues[#]],
            {uColorOf[#], Sphere[{p[[1]], p[[2]], uValues[#]}, 0.04]},
            {}]] &,
    vertices];

plotC = Graphics3D[
    {skeleton3D, uCurves3D, dropLines3D, vertexSpheres3D, uSpheres3D},
    Axes      -> True,
    AxesLabel -> {"x", "y", "u"},
    Boxed     -> False,
    PlotLabel -> Style["C: 3D surface (network \[Times] u)", Bold, 12],
    ImageSize -> 500,
    ViewPoint -> {2, -3, 1.5}
];

(* ------------------------------------------------------------------ *)
(* 7. COLORBAR (shared scale)                                         *)
(* ------------------------------------------------------------------ *)

colorBar = DensityPlot[y, {x, 0, 0.15}, {y, uMin, uMax},
    ColorFunction -> (uColorFn[Rescale[#, {uMin, uMax}]] &),
    ColorFunctionScaling -> False,
    FrameTicks -> {None, Automatic, None, None},
    PlotLabel  -> "u scale",
    ImageSize  -> {60, 300}];

(* ------------------------------------------------------------------ *)
(* 8. DISPLAY                                                          *)
(* ------------------------------------------------------------------ *)

Grid[{
    {plotA, plotB},
    {plotC, colorBar}
}, Spacings -> {2, 1}]
