(* Wolfram Language Test file *)

Test[
	MFGEquations = DataToEquations[DataG[26]];
    MFGEquations["criticalreduced"][[2]]//KeySort//Values
	,
	{}
    ,
    TestID -> "primeiro de uma longa lista"
]
