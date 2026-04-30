SetDirectory[ParentDirectory[ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]]];
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];

key = "Grid0404";
Print["Loading key: ", key];
data = GetExampleData[key];
Print["Loaded association: ", AssociationQ[data]];
Print["Keys: ", Keys[data]];

result = Check[
  DataToEquations[data],
  $Failed
];

Print["Result head: ", Head[result]];
If[result === $Failed,
  Print["DataToEquations returned $Failed"];
,
  Print["Association result: ", AssociationQ[result]];
  If[AssociationQ[result],
    Print["Result keys: ", Keys[result]]
  ];
];
