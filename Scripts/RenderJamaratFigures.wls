#!/usr/bin/env wolframscript
(* Render the three Jamaratv9 richNetworkPlot figures from MFGraphs/Jamarat.wl
   as presentation-ready PDFs. Two versions per plot: with PlotLabel and without. *)

Needs["MFGraphs`"];

jamaratScenario = getExampleScenario[
    "Jamaratv9",
    {{1, 40}, {2, 45}},
    {{7, 3}, {8, 4}, {9, 2}}
];
jamaratSystem = makeSystem[jamaratScenario];

Print["Solving Jamaratv9 ..."];
{solveTime, jamaratSolRaw} = AbsoluteTiming[
    TimeConstrained[solveScenario[jamaratScenario], 1800, $TimedOut]
];
Print["Solve time: ", solveTime, "s"];

If[jamaratSolRaw === $TimedOut,
    Print["Solve timed out — aborting."];
    Exit[1]
];

jamaratSol = If[Head[jamaratSolRaw] =!= Association,
    <|"Rules" -> jamaratSolRaw, "Residual" -> True|>,
    jamaratSolRaw
];

outDir = FileNameJoin[{Directory[], "Results", "Jamaratv9_figures"}];
If[!DirectoryQ[outDir], CreateDirectory[outDir]];

(* Post-process the rendered Graph to bump font sizes and disc size for
   presentation legibility. The package hardcodes Style[..., 9] / Style[..., 8]
   for vertex/edge labels and VertexSize -> 0.58 inside the Graph. *)
beautify[Legended[g_Graph, lg_]] := Legended[beautify[g], lg];
beautify[g_Graph] := Module[{updated},
    updated = g /. {
        Style[expr_, 9, rest___]      :> Style[expr, 14, FontFamily -> "Helvetica", rest],
        Style[expr_, 8, rest___]      :> Style[expr, 13, FontFamily -> "Helvetica", rest],
        Style[expr_, 14, Bold]        :> Style[expr, 22, Bold, FontFamily -> "Helvetica"]
    };
    Graph[updated,
        VertexSize -> 0.85,
        VertexLabelStyle -> Directive[FontSize -> 14, FontFamily -> "Helvetica"]
    ]
];

exportPdf[stem_String, plot_] := Module[{path, decorated},
    path = FileNameJoin[{outDir, stem <> ".pdf"}];
    decorated = beautify[plot];
    Export[path, decorated, "PDF", ImageResolution -> 300];
    Print["Wrote ", path]
];

(* ---- Figure 1: augmented infrastructure (no solution) ---- *)
fig1Labeled = richNetworkPlot[
    jamaratScenario, jamaratSystem, <||>,
    PlotLabel -> "Jamarat augmented infrastructure",
    ImageSize -> 1400,
    ShowBoundaryValues -> False
];
fig1Plain = richNetworkPlot[
    jamaratScenario, jamaratSystem, <||>,
    PlotLabel -> None,
    ImageSize -> 1400,
    ShowBoundaryValues -> False
];
exportPdf["fig1_infrastructure_labeled", fig1Labeled];
exportPdf["fig1_infrastructure", fig1Plain];

(* ---- Figure 2: equilibrium ---- *)
fig2Labeled = richNetworkPlot[
    jamaratScenario, jamaratSystem, jamaratSol,
    PlotLabel -> "Jamarat: equilibrium",
    ImageSize -> 1400
];
fig2Plain = richNetworkPlot[
    jamaratScenario, jamaratSystem, jamaratSol,
    PlotLabel -> None,
    ImageSize -> 1400
];
exportPdf["fig2_equilibrium_labeled", fig2Labeled];
exportPdf["fig2_equilibrium", fig2Plain];

(* ---- Figure 3: equilibrium with value-color function ---- *)
fig3Labeled = richNetworkPlot[
    jamaratScenario, jamaratSystem, jamaratSol,
    PlotLabel -> "Jamarat equilibrium: value colors",
    ImageSize -> 1400,
    UseColorFunction -> True
];
fig3Plain = richNetworkPlot[
    jamaratScenario, jamaratSystem, jamaratSol,
    PlotLabel -> None,
    ImageSize -> 1400,
    UseColorFunction -> True
];
exportPdf["fig3_value_colors_labeled", fig3Labeled];
exportPdf["fig3_value_colors", fig3Plain];

Print["Done. ", Length[FileNames["*.pdf", outDir]], " PDFs in ", outDir];
