(* Exact, deliberately bounded paper-foundation collection. No cached solutions.
   Entries are prescribed rates; exits are terminal costs. All physical edges
   use the package's unit critical traversal law, including literature analogs. *)
$CuratedExactScenarios = <|
    "competing-exits" -> <|"Key" -> "chain with two exits",
        "Entries" -> {{1, 100}}, "Exits" -> {{2, 0}, {3, 10}},
        "Purpose" -> "Exit competition and an unused expensive exit; value indeterminacy."|>,
    "diamond" -> <|"Key" -> 12,
        "Entries" -> {{1, 80}}, "Exits" -> {{4, 0}},
        "Purpose" -> "Parallel routes, splitting and merging."|>,
    "triangle" -> <|"Key" -> "triangle with two exits",
        "Entries" -> {{1, 80}}, "Exits" -> {{2, 0}, {3, 0}},
        "Purpose" -> "A cycle with competing equal-cost exits and a dormant cross edge."|>,
    "braess-split" -> <|"Key" -> "Braess split",
        "Entries" -> {{1, 100}}, "Exits" -> {{8, 0}},
        "Purpose" -> "Parallel congested routes with canonical switching tolls; comparison topology."|>,
    "braess-congest" -> <|"Key" -> "Braess congest",
        "Entries" -> {{1, 100}}, "Exits" -> {{7, 0}},
        "Purpose" -> "Merge-split bottleneck and nonunique transition routing; no paradox claim without a controlled comparison."|>,
    "camilli-simple" -> <|"Key" -> "Camilli 2015 simple",
        "Entries" -> {{3, 50}, {4, 50}}, "Exits" -> {{1, 0}},
        "Purpose" -> "Verified Figure 1 topology, two sources and a cycle; stationary unit-edge analog only."|>,
    "grid3" -> <|"Key" -> "Grid0303",
        "Entries" -> {{1, 100}}, "Exits" -> {{9, 0}},
        "Purpose" -> "Structured cyclic network with transition-flow families."|>,
    "grid4" -> <|"Key" -> "Grid0404",
        "Entries" -> {{1, 100}}, "Exits" -> {{16, 0}},
        "Purpose" -> "Moderately sized structured showcase and branch-enumeration stress test."|>,
    "jamarat" -> <|"Key" -> "Jamaratv9",
        "Entries" -> {{1, 100}, {2, 100}}, "Exits" -> {{7, 0}, {8, 0}, {9, 0}},
        "Purpose" -> "Two entrances and three exits with many admissible transition routings."|>,
    "case23" -> <|"Key" -> 23,
        "Entries" -> {{1, 50}, {2, 50}}, "Exits" -> {{5, 0}, {6, 0}},
        "Purpose" -> "Hard merge-split case previously rescued by lexicographic branch-state enumeration."|>
|>;

curatedExactScenario[name_String] := With[{spec = $CuratedExactScenarios[name]},
    getExampleScenario[spec["Key"], spec["Entries"], spec["Exits"]]
];
