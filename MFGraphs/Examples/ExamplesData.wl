(* ::Package:: *)

(* Wolfram Language package *)
(* Built-in test cases for MFGraphs *)

GetExampleData::usage =
"GetExampleData[n] returns an Association with keys \"Vertices List\", \"Adjacency Matrix\",
\"Entrance Vertices and Flows\", \"Exit Vertices and Terminal Costs\", and \"Switching Costs\"
from the test case test[n].
Supports test cases with 5 fields (basic), 6 fields (+cost function 'a'), or 7 fields (+alpha).";

DataG::usage = "DataG is a backward-compatibility alias for GetExampleData.";

GetExampleData::badfields =
"Unexpected number of fields (`1`) for test case `2`.";

(* Declare symbolic parameters in the MFGraphs` context so that user-level
   substitution rules like {I1 -> 100, U1 -> 0} match the symbols returned
   by GetExampleData. Without these declarations the symbols end up in
   MFGraphs`Private` and substitution silently fails. *)
I1::usage = "Symbolic parameter: input flow at entrance 1.";
I2::usage = "Symbolic parameter: input flow at entrance 2.";
I3::usage = "Symbolic parameter: input flow at entrance 3.";
U1::usage = "Symbolic parameter: terminal cost at exit 1.";
U2::usage = "Symbolic parameter: terminal cost at exit 2.";
U3::usage = "Symbolic parameter: terminal cost at exit 3.";
S1::usage  = "Symbolic parameter: switching cost 1.";
S2::usage  = "Symbolic parameter: switching cost 2.";
S3::usage  = "Symbolic parameter: switching cost 3.";
S4::usage  = "Symbolic parameter: switching cost 4.";
S5::usage  = "Symbolic parameter: switching cost 5.";
S6::usage  = "Symbolic parameter: switching cost 6.";
S7::usage  = "Symbolic parameter: switching cost 7.";
S8::usage  = "Symbolic parameter: switching cost 8.";
S9::usage  = "Symbolic parameter: switching cost 9.";
S10::usage = "Symbolic parameter: switching cost 10.";
S11::usage = "Symbolic parameter: switching cost 11.";
S12::usage = "Symbolic parameter: switching cost 12.";
S13::usage = "Symbolic parameter: switching cost 13.";
S14::usage = "Symbolic parameter: switching cost 14.";
S15::usage = "Symbolic parameter: switching cost 15.";
S16::usage = "Symbolic parameter: switching cost 16.";

Begin["`Private`"];

		test = Association[
    (* One vertex *)
    1 -> {
        (*VL=*){1},
        (*AM=*){{0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{1, U1}},
        (*SwitchingCostsData=*){}},

    (* 1 edge *)
    2 -> {
        (*VL=*){1, 2},
        (*AM=*){{0, 1}, {0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{2, U1}},
        (*SwitchingCostsData=*){}},

    (* 2 edges *)
    3 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){}},

    (* 3 edges *)
   4 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){}},

    (* 4 edges *)
    5 -> {
        (*VL=*){1, 2, 3, 4, 5},
        (*AM=*){{0, 1, 0, 0, 0}, {0, 0, 1, 0, 0},
        	    {0, 0, 0, 1, 0}, {0, 0, 0, 0 ,1}, {0, 0, 0, 0 ,0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{5, U1}},
        (*SwitchingCostsData=*){}},

   (* 9 edges *)
    6 -> {
        (*VL=*){1,2,3,4,5,6,7,8,9,10},
        (*AM=*){
        {0,1,0,0,0,0,0,0,0,0},
        {0,0,1,0,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0,0},
        {0,0,0,0,1,0,0,0,0,0},
        {0,0,0,0,0,1,0,0,0,0},
        {0,0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{10, U1}},
        (*SwitchingCostsData=*){}},

    (* Y 1-in 2-out 4 vertices no switching *)
    7 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}, {4, U2}},
        (*SwitchingCostsData=*){}},

    (* Y 1-in 2-out 4 vertices with switching *)
    8 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}, {4, U2}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {1, 2, 4, S2}, {3, 2, 1, S3}, {3, 2, 4, S4}, {4, 2, 1, S5}, {4, 2, 3, S6}}
        },

    (* Y 2-in 1-out 4 vertices no switching *)
    9 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}, {3, I2}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){}},

    (* Y 2-in 1-out 4 vertices with switching *)
    10 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}, {3, I2}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){{1, 2, 4, S1}, {1, 2, 3, S2}, {3, 2, 1, S3}, {3, 2, 4, S4}, {4, 2, 1, S5}, {4, 2, 3, S6}}
        },

    (* Attraction Problem *)
    11 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){
            {1, 2, 4, S1},
            {1, 2, 3, S2},
            {3, 2, 1, S3},
            {3, 2, 4, S4},
            {4, 2, 1, S5},
            {4, 2, 3, S6},
            {1, 3, 4, S7},
            {4, 3, 1, S8},
            {1, 3, 2, S9},
            {2, 3, 1, S10},
            {3, 4, 2, S11},
            {2, 4, 3, S12},
            {2, 3, 4, S13},
            {4, 3, 2, S14},
            {3, 1, 2, S15},
            {2, 1, 3, S16}}
        },

   (* Attraction without switching *)
    12 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){}},

    (* Attraction without edge in the middle *)
    13 -> {
    	(*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){{1,2,4,S1},{1,3,4,S2},{4,2,1,S3},{4,3,1,S4}}}/. {I1 -> 1, U1 -> 0, S1 -> 0, S2 -> 0, S3 -> 0, S4 -> 0},

   (* Triangle *)
    14 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {3, 2, 1, S2}}
        },

   (* Triangle with 2 entrances and 3 exits *)
    104 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
        (*DataIn=*){{1, I1},{2,I2}},
        (*FinalCosts=*){{1,0},{2,0},{3, 0}},
        (*SwitchingCostsData=*){}
        },

	(* Triangle with two exits *)
    "triangle with two exits" -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{2, U1},{3, U2}},
        (*SwitchingCostsData=*){}
        },

    (* Chain with two exits: 1 -> 2 -> 3, entry at 1, exits at 2 and 3 *)
    "chain with two exits" -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{2, U1}, {3, U2}},
        (*SwitchingCostsData=*){}
        },

    (* Numeric alias for chain with two exits *)
    105 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{2, U1}, {3, U2}},
        (*SwitchingCostsData=*){}
        },

    (* Y 2-in 1-out 3 vertices with switching *)
    15 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 0, 0}, {1, 0, 0}, {1, 0, 0}},
        (*DataIn=*){{2, I1}, {3, I2}},
        (*FinalCosts=*){{1, U1}},
        (*SwitchingCostsData=*){{2, 1, 3, S1}, {3, 1, 2,S2}}
        },

    (* Error in the switching costs *)
    16 -> {
        (*VL*){1, 2, 3},
        (*AM*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn*){{1, I1}},
        (*FinalCosts*){{3, U1}},
        (*SwitchingCostsData*){{1, 3, 2, S1}}
        },

    (* 2 edges entrance in the middle *)
    17 -> {
        (*VL=*) {1, 2, 3},
        (*AM*) {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn*) {{2, I1}},
        (*FinalCost*) {{1, U1}, {3,U2}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {3, 2, 1, S2}}
        },

     (* 2 edges *)
    18 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){ {1, 2, 3, S1}, {3, 2, 1, S2}}
        },

    (* Y 1-in 2-out 4 vertices with switching, 2 entrances *)
    19 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1},{4, I2}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {1, 2, 4, S2}, {3, 2, 1, S3}, {3, 2, 4, S4}, {4, 2, 1, S5}, {4, 2, 3, S6}}
        },

    20 -> {
        (*VL=*){1,2,3,4,5,6,7,8,9},
        (*AM=*){
        {0,0,1,0,0,0,0,0,0},
        {0,0,1,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0},
        {0,0,0,0,1,0,0,0,0},
        {0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,0,1,1,1},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{7, U1},{8, U2},{9, U3}},
        (*SwitchingCostsData=*){}},

    21 -> {
        (*VL=*){1,2,3,4,5,6,7,8,9,10,11,12},
        (*AM=*){
        {0,0,1,0,0,0,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0,0,0,0},
        {0,0,0,1,1,0,0,0,0,0,0,0},
        {0,0,0,0,0,1,0,0,0,0,0,0},
        {0,0,0,0,0,1,1,0,0,0,0,0},
        {0,0,0,0,0,0,0,1,0,0,0,0},
        {0,0,0,0,0,0,0,1,1,1,0,0},
        {0,0,0,0,0,0,0,0,1,0,1,0},
        {0,0,0,0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{10, U1},{11, U2},{12, U3}},
        (*SwitchingCostsData=*){}},

    22 -> {
        (*VL=*){1,2,3,4,5,6,7},
        (*AM=*){
        {0,1,1,0,0,0,0},
        {0,0,0,1,0,0,0},
        {0,0,0,1,1,0,0},
        {0,0,0,0,0,1,0},
        {0,0,0,0,0,1,1},
        {0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{5, U1},{6, U2},{7, U3}},
        (*SwitchingCostsData=*){}},

    "Jamaratv9"->
       {
        (*VL=*){1,2,3,4,5,6,7,8,9},
        (*AM=*){
        {0,1,1,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0},
        {0,0,0,1,1,0,0,0,0},
        {0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,1,1,0,0},
        {0,0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0}
        },
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{7, U1},{8, U2},{9, U3}},
        (*SwitchingCostsData=*){}},

    23 -> {
        (*VL=*){1,2,3,4,5,6},
        (*AM=*){
        {0,1,1,0,0,0},
        {0,0,1,0,0,0},
        {0,0,0,1,1,0},
        {0,0,0,0,1,1},
        {0,0,0,0,0,1},
        {0,0,0,0,0,0}},
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{4, U1},{5, U2},{6, U3}},
        (*SwitchingCostsData=*){}},

    27 -> {
        (*VL=*){1, 2},
        (*AM=*){{0, 1}, {1, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{2, U1}},
        (*SwitchingCostsData=*){}},

    "Braess split" -> {
        (*VL=*){1,2,3,4,5,6,7,8},
        (*AM=*){
        {0,1,1,0,0,0,0,0},
        {0,0,0,1,0,0,0,0},
        {0,0,0,0,1,0,0,0},
        {0,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->2},
        (*FinalCosts=*){{8, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,2,4,S1},{5,7,8,S2}}/.{S1->10,S2->10}},

    "Braess congest" -> {
        (*VL=*){1,2,3,4,5,6,7},
        (*AM=*){
        {0,1,1,0,0,0,0},
        {0,0,0,1,0,0,0},
        {0,0,0,1,0,0,0},
        {0,0,0,0,1,1,0},
        {0,0,0,0,0,0,1},
        {0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->2},
        (*FinalCosts=*){{7, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,2,4,S1},{4,6,7,S2}}/.{S1->10,S2->10}},

    "New Braess" -> {
        (*VL=*){1,2,3,4,5,6},
        (*AM=*){
        {0,1,0,1,0,0},
        {0,0,1,0,0,0},
        {0,0,0,0,0,1},
        {0,0,0,0,1,0},
        {0,0,0,0,0,1},
        {0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->4000},
        (*FinalCosts=*){{6, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,2,3,S1}, {3,2,1,S1},{4,5,6,S2}, {6,5,4,S2}}/.{S1->45,S2->45},
        (*a=*)Function[{j,edge},
        		Which[
        			edge === DirectedEdge[3,6] || edge === DirectedEdge[1,4],
        				j/100,
        			True, 0
        		]
        ]
        },

    "Big Braess split" -> {
        (*VL=*){1,2,3,4,5,6,7,8,9,10},
        (*AM=*){
        {0,1,1,0,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0,0},
        {0,0,0,0,0,1,0,0,0,0},
        {0,0,0,0,1,0,0,0,0,0},
        {0,0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->2},
        (*FinalCosts=*){{10, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,3,6,S1},{5,7,10,S2}}/.{S1->10,S2->10}},

    "Big Braess congest" -> {
        (*VL=*){1,2,3,4,5,6,7,8,9},
        (*AM=*){
        {0,1,1,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0},
        {0,0,0,0,1,0,0,0,0},
        {0,0,0,0,1,0,0,0,0},
        {0,0,0,0,0,1,1,0,0},
        {0,0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->2},
        (*FinalCosts=*){{9, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,3,5,S1},{5,7,9,S2}}/.{S1->10,S2->10}},

    "HRF Scenario 1" -> {
        (*VL=*){1,2,3,4,5,6,7,8,9,10},
        (*AM=*){
        {0,1,0,0,0,0,0,0,0,0},
        {0,0,1,1,0,0,0,0,0,0},
        {0,0,0,1,1,0,1,0,0,0},
        {0,0,0,0,1,1,1,0,0,0},
        {0,0,0,0,0,1,1,0,0,0},
        {0,0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,0,0,1,0,1},
        {0,0,0,0,0,0,0,0,0,0},
        {0,0,1,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}, {9, I2}},
        (*FinalCosts=*){{8, U1}, {10, U2}},
        (*SwitchingCostsData=*){}
        },

    "Paper example" -> {
        (*VL=*){1,2,3,4},
        (*AM=*){
        {0,1,0,0},
        {0,0,1,0},
        {0,0,0,1},
        {0,0,0,0}},
        (*DataIn=*){{1, I1},{3, I2}},
        (*FinalCosts=*){{2, U1},{4,U2}}/.{U1->0,U2->0},
        (*SwitchingCostsData=*){{1,2,3,S1},{3,2,1,S2},{2,3,4,S3},{4,3,2,S4}}},

    "Grid1020" -> {
    	(*VL=*)VertexList[GridGraph[{10,20}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{10,20}, DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{200, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        },

    "Grid0303" -> {
    	(*VL=*)VertexList[GridGraph[{3,3}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{3,3},DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{3*3, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        },

    "Grid0404" -> {
        (*VL=*)VertexList[GridGraph[{4,4}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{4,4},DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{4*4, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        },

    "Grid0505" -> {
        (*VL=*)VertexList[GridGraph[{5,5}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{5,5},DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{5*5, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        },

    "Grid0707" -> {
        (*VL=*)VertexList[GridGraph[{7,7}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{7,7},DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{7*7, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        },

    "Grid0710" -> {
        (*VL=*)VertexList[GridGraph[{7,10}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{7,10},DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{7*10, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        },

    "Grid1010" -> {
        (*VL=*)VertexList[GridGraph[{10,10}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{10,10},DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{10*10, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        }
];


(* --- GetExampleData: construct an Association from a test case --- *)

$ExampleDataFields5 = {
    "Vertices List",
    "Adjacency Matrix",
    "Entrance Vertices and Flows",
    "Exit Vertices and Terminal Costs",
    "Switching Costs"};

$ExampleDataFields6 = Join[$ExampleDataFields5, {"a"}];

$ExampleDataFields7 = Join[$ExampleDataFields6, {"alpha"}];

GetExampleData[n_] :=
    With[{data = test[n]},
        Switch[Length[data],
            5, AssociationThread[$ExampleDataFields5, data],
            6, AssociationThread[$ExampleDataFields6, data],
            7, AssociationThread[$ExampleDataFields7, data],
            _, (Message[GetExampleData::badfields, Length[data], n]; $Failed)
        ]
    ]

(* --- Backward compatibility alias --- *)
DataG = GetExampleData;

End[];
