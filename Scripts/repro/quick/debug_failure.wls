Get["/Users/ribeirrd/Documents/GitHub/MFGraphs/MFGraphs/MFGraphs.wl"];

f = getExampleScenario["Inconsistent Y shortcut"];
sc = f[{{1, 100}}, {{3, 0}, {4, 0}}, Automatic];

If[FailureQ[sc],
  Print["Failure Message: ", sc[[2, "Message"]]];
  Print["Full Failure: ", InputForm[sc]];
];
