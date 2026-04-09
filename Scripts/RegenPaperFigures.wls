#!/usr/bin/env wolframscript
(* RegenPaperFigures.wls -- Regenerate the network diagrams for the paper
   "An Exact Solver for Stationary Mean-Field Games on Networks" by
   evaluating the Mathematica Graph expressions directly and exporting to PDF.

   The paper no longer uses the MLaTeX IMG-block pipeline, so this script is
   the authoritative source for the network figures in the paper.

   Outputs (to the paper's figures/ directory):
     example1.pdf   -- road merge   (original)
     example2.pdf   -- road fork    (original)
     example22.pdf  -- road fork    (augmented with auxiliary exit edges)
     Triangle.pdf   -- triangular network (with aux entry/exit edges)
     Jama.pdf       -- Jamarat bridge network

   Usage: wolframscript -file Scripts/RegenPaperFigures.wls
*)

figDir = "/Users/ribeirrd/Dropbox/Apps/Overleaf/Stationary solutions of a congested MFG on Networks/figures";
If[!DirectoryQ[figDir],
  Print["Figures directory not found: ", figDir];
  Exit[1]
];

(* --- example1: road merge, original graph --- *)
ex1 = Graph[{Labeled[1 \[UndirectedEdge] 2, Subscript["e", "12"]],
  Labeled[3 \[UndirectedEdge] 2, Subscript["e", "32"]],
  Labeled[2 \[UndirectedEdge] 4, Subscript["e", "24"]]},
  VertexLabels -> {1 -> Subscript["v", 1], 2 -> Subscript["v", 2],
    3 -> Subscript["v", 3], 4 -> Subscript["v", 4]},
  GraphLayout -> "RadialDrawing", VertexLabelStyle -> 14, EdgeLabelStyle -> 14];
Export[FileNameJoin[{figDir, "example1.pdf"}], ex1];
Print["Wrote example1.pdf"];

(* --- example2: road fork, original graph --- *)
ex2 = Graph[{Labeled[1 \[UndirectedEdge] 2, Subscript["e", "12"]],
  Labeled[2 \[UndirectedEdge] 3, Subscript["e", "23"]],
  Labeled[2 \[UndirectedEdge] 4, Subscript["e", "24"]]},
  VertexLabels -> {1 -> Subscript["v", 1], 2 -> Subscript["v", 2],
    3 -> Subscript["v", 3], 4 -> Subscript["v", 4]},
  GraphLayout -> "RadialDrawing", VertexLabelStyle -> 14, EdgeLabelStyle -> 14];
Export[FileNameJoin[{figDir, "example2.pdf"}], ex2];
Print["Wrote example2.pdf"];

(* --- example22: road fork, augmented graph --- *)
ex22 = Graph[{Labeled[1 \[UndirectedEdge] 2, Subscript["e", "12"]],
  Labeled[2 \[UndirectedEdge] 3, Subscript["e", "23"]],
  Labeled[2 \[UndirectedEdge] 4, Subscript["e", "24"]],
  Labeled[3 \[UndirectedEdge] 6, Subscript["e", Row[{3, OverTilde[3]}]]],
  Labeled[4 \[UndirectedEdge] 7, Subscript["e", Row[{4, OverTilde[4]}]]]},
  VertexLabels -> {1 -> Subscript["v", 1], 2 -> Subscript["v", 2],
    3 -> Subscript["v", 3], 4 -> Subscript["v", 4],
    6 -> Subscript["v", OverTilde[3]], 7 -> Subscript["v", OverTilde[4]]},
  EdgeStyle -> {3 \[UndirectedEdge] 6 -> Dashed, 4 \[UndirectedEdge] 7 -> Dashed},
  GraphLayout -> "RadialDrawing", VertexLabelStyle -> 14, EdgeLabelStyle -> 14];
Export[FileNameJoin[{figDir, "example22.pdf"}], ex22];
Print["Wrote example22.pdf"];

(* --- triangle: 3-vertex triangle with aux entry/exit --- *)
tri = Graph[{Labeled[1 \[UndirectedEdge] 2, Subscript["e", "12"]],
  Labeled[2 \[UndirectedEdge] 3, Subscript["e", "23"]],
  Labeled[1 \[UndirectedEdge] 3, Subscript["e", "13"]],
  Labeled[4 \[UndirectedEdge] 1, Subscript["e", Row[{OverTilde[1], 1}]]],
  Labeled[3 \[UndirectedEdge] 5, Subscript["e", Row[{3, OverTilde[3]}]]]},
  VertexLabels -> {1 -> Subscript["v", 1], 2 -> Subscript["v", 2],
    3 -> Subscript["v", 3], 4 -> Subscript["v", OverTilde[1]],
    5 -> Subscript["v", OverTilde[3]]},
  EdgeStyle -> {4 \[UndirectedEdge] 1 -> Dashed, 3 \[UndirectedEdge] 5 -> Dashed},
  GraphLayout -> "SpringEmbedding", VertexLabelStyle -> 14, EdgeLabelStyle -> 14];
Export[FileNameJoin[{figDir, "Triangle.pdf"}], tri];
Print["Wrote Triangle.pdf"];

(* --- Jama: Jamarat bridge --- *)
jama = Graph[{Labeled[1 \[UndirectedEdge] 2, Subscript["e", "12"]],
  Labeled[2 \[UndirectedEdge] 4, Subscript["e", "24"]],
  Labeled[1 \[UndirectedEdge] 3, Subscript["e", "13"]],
  Labeled[3 \[UndirectedEdge] 4, Subscript["e", "34"]],
  Labeled[3 \[UndirectedEdge] 5, Subscript["e", "35"]],
  Labeled[4 \[UndirectedEdge] 6, Subscript["e", "46"]],
  Labeled[5 \[UndirectedEdge] 6, Subscript["e", "56"]],
  Labeled[5 \[UndirectedEdge] 7, Subscript["e", "57"]],
  Labeled[6 \[UndirectedEdge] 7, Subscript["e", "67"]],
  Labeled[8 \[UndirectedEdge] 1, Subscript["e", Row[{OverTilde[1], 1}]]],
  Labeled[9 \[UndirectedEdge] 2, Subscript["e", Row[{OverTilde[2], 2}]]],
  Labeled[5 \[UndirectedEdge] 10, Subscript["e", Row[{5, OverTilde[5]}]]],
  Labeled[6 \[UndirectedEdge] 11, Subscript["e", Row[{6, OverTilde[6]}]]],
  Labeled[7 \[UndirectedEdge] 12, Subscript["e", Row[{7, OverTilde[7]}]]]},
  VertexLabels -> {1 -> Subscript["v", 1], 2 -> Subscript["v", 2],
    3 -> Subscript["v", 3], 4 -> Subscript["v", 4], 5 -> Subscript["v", 5],
    6 -> Subscript["v", 6], 7 -> Subscript["v", 7],
    8 -> Subscript["v", OverTilde[1]], 9 -> Subscript["v", OverTilde[2]],
    10 -> Subscript["v", OverTilde[5]], 11 -> Subscript["v", OverTilde[6]],
    12 -> Subscript["v", OverTilde[7]]},
  EdgeStyle -> {8 \[UndirectedEdge] 1 -> Dashed, 9 \[UndirectedEdge] 2 -> Dashed,
    5 \[UndirectedEdge] 10 -> Dashed, 6 \[UndirectedEdge] 11 -> Dashed,
    7 \[UndirectedEdge] 12 -> Dashed},
  VertexLabelStyle -> 14, EdgeLabelStyle -> 14];
Export[FileNameJoin[{figDir, "Jama.pdf"}], jama];
Print["Wrote Jama.pdf"];

Print["Done."];
