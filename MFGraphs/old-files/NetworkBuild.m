(* Wolfram Language package *)
DataToGraph::usage = 
"DataToGraph[<|\"Vertices List\" -> {1, 
  2}, \"Adjacency Matrix\" -> {{0, 1}, {0, 
   0}}, \"Entrance Vertices and Currents\" -> {{1, 
   I1}}, \"Exit Vertices and Terminal Costs\" -> {{2, 
   U1}}, \"Switching Costs\" -> {}|>] builds the graph with corresponding data.";

Begin["`Private`"]
DataToGraph[Data_?AssociationQ] :=
   
    Module[{BG, EntranceVertices, InwardVertices, InEdges, ExitVertices, OutwardVertices, OutEdges, AuxiliaryGraph, FG},
    BG = AdjacencyGraph[Data["Vertices List"], Data["Adjacency Matrix"], VertexLabels -> "Name"];
    EntranceVertices = First /@ Data["Entrance Vertices and Currents"];
    InwardVertices = AssociationThread[EntranceVertices, Unique["en"] & /@ EntranceVertices];
    (*Defines auxiliary vertices for the entrance vertices*)
    InEdges = MapThread[DirectedEdge, {InwardVertices /@ EntranceVertices, EntranceVertices}];
    ExitVertices = First /@ Data["Exit Vertices and Terminal Costs"];
    OutwardVertices = AssociationThread[ExitVertices, Unique["ex"] & /@ ExitVertices];
    OutEdges = MapThread[DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
    AuxiliaryGraph = Graph[Join[InEdges, OutEdges], VertexLabels -> "Name"];
    FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
    Association[{
    	"BG" -> BG, 
    	"EntranceVertices" -> EntranceVertices, 
    	"InwardVertices" -> InwardVertices,
    	"InEdges" -> InEdges, 
    	"ExitVertices" -> ExitVertices, 
    	"OutwardVertices" -> OutwardVertices, 
    	"OutEdges" -> OutEdges, 
    	"AuxiliaryGraph" -> AuxiliaryGraph, 
    	"FG" -> FG
    	 }]
    ]
End[]