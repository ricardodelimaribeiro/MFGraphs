SetDirectory[ParentDirectory[ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]]];
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
example = GetExampleData[1];
Print["MFGraphs loaded: ", AssociationQ[example]];
Print["Example 1 keys count: ", Length[Keys[example]]];
Print["Example 1 has Vertices: ", KeyExistsQ[example, "Vertices"]];
