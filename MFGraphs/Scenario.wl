(*
   Scenario: typed scenario kernel for MFGraphs.

   Provides a typed wrapper around network-plus-parameters data so that
   scenarios can be constructed, validated, completed, and stored uniformly.

   Lifecycle:
     makeScenario[rawAssoc]  →  validates + completes + wraps in scenario[...] head
     validateScenario[s]     →  checks required keys; returns s or Failure[...]
     completeScenario[s]     →  fills Identity (hash), Benchmark defaults; returns s
     scenarioQ[x]            →  True iff x is a well-formed typed scenario

   Canonical top-level blocks inside scenario[<|...|>]:
     "Identity"    — name, version, contentHash
     "Model"       — raw network topology accepted by DataToEquations
     "Data"        — parameter substitution rules (e.g. {I1 -> 100, U1 -> 0})
     "Validation"  — consistency check results (populated after solving)
     "Benchmark"   — tier, timeout
     "Visualization" — optional plotting hints
     "Inheritance" — optional lineage/parent reference
*)

(* --- Public API declarations --- *)

scenario::usage =
"scenario[assoc] is the typed head for a MFGraphs scenario object. Use makeScenario \
to construct one; use ScenarioData to access keys.";

scenarioQ::usage =
"scenarioQ[x] returns True if x is a typed scenario[assoc_Association] object, \
False otherwise.";

makeScenario::usage =
"makeScenario[assoc] constructs a typed scenario from a raw association. The input \
must contain a \"Model\" key whose value is a network topology association accepted by \
DataToEquations (required keys: \"Vertices List\", \"Adjacency Matrix\", \
\"Entrance Vertices and Flows\", \"Exit Vertices and Terminal Costs\", \
\"Switching Costs\"). Optional keys: \"Data\" (parameter substitution rules), \
\"Identity\" (name, version), \"Benchmark\" (tier, timeout), \"Visualization\", \
\"Inheritance\". Returns a scenario[...] object on success or Failure[...] on error.";

validateScenario::usage =
"validateScenario[s] checks that the scenario s has all required Model keys and that \
the Model value is an Association. Returns s unchanged on success, or \
Failure[\"ScenarioValidation\", <|\"Message\" -> msg, \"MissingKeys\" -> {...}|>] \
on failure.";

completeScenario::usage =
"completeScenario[s] fills in derived fields: sets \"contentHash\" in the Identity \
block (SHA256 of the canonical Model string), and supplies default Benchmark values \
(\"Tier\" -> \"core\", \"Timeout\" -> 300) when missing. Returns a new scenario object.";

ScenarioData::usage =
"ScenarioData[s, key] returns the value associated with key in the scenario s, or \
Missing[\"KeyAbsent\", key] if absent. ScenarioData[s] returns the underlying \
Association.";

Begin["`Private`"];

(* Required keys that must appear inside the "Model" block. *)
$RequiredModelKeys = {
    "Vertices List",
    "Adjacency Matrix",
    "Entrance Vertices and Flows",
    "Exit Vertices and Terminal Costs",
    "Switching Costs"
};

(* Default benchmark settings applied by completeScenario. *)
$DefaultBenchmarkTier    = "core";
$DefaultBenchmarkTimeout = 300;

(* --- Type predicate --- *)

scenarioQ[scenario[_Association]] := True;
scenarioQ[_]                       := False;

(* --- Accessor --- *)

ScenarioData[scenario[assoc_Association]]           := assoc;
ScenarioData[scenario[assoc_Association], key_]     := Lookup[assoc, key, Missing["KeyAbsent", key]];

(* --- Validate --- *)

validateScenario[s_scenario] :=
    Module[{assoc, model, missing},
        assoc = ScenarioData[s];
        model = Lookup[assoc, "Model", Missing["KeyAbsent", "Model"]];
        Which[
            MissingQ[model],
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Model\" key is absent.",
                      "MissingKeys" -> {"Model"}|>],
            !AssociationQ[model],
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Model\" value must be an Association.",
                      "MissingKeys" -> {}|>],
            True,
                missing = Select[$RequiredModelKeys, !KeyExistsQ[model, #] &];
                If[missing === {},
                    s,
                    Failure["ScenarioValidation",
                        <|"Message" ->
                            "Missing required Model keys: " <> StringRiffle[missing, ", "],
                          "MissingKeys" -> missing|>]
                ]
        ]
    ];

(* Reject non-scenario inputs immediately. *)
validateScenario[x_] :=
    Failure["ScenarioValidation",
        <|"Message" -> "Input is not a scenario object.",
          "MissingKeys" -> {}|>];

(* --- Complete --- *)

completeScenario[s_scenario] :=
    Module[{assoc, identity, benchmark, model, hash, newAssoc},
        assoc     = ScenarioData[s];
        model     = Lookup[assoc, "Model", <||>];
        identity  = Lookup[assoc, "Identity", <||>];
        benchmark = Lookup[assoc, "Benchmark", <||>];

        (* Compute content hash from the canonical string of the Model. *)
        hash = Hash[ToString[model, InputForm], "SHA256", "HexString"];

        (* Computed hash always wins: placed on the right so it overrides any user-supplied value. *)
        identity  = Join[identity, <|"contentHash" -> hash|>];
        benchmark = Join[
            <|"Tier" -> $DefaultBenchmarkTier, "Timeout" -> $DefaultBenchmarkTimeout|>,
            benchmark
        ];

        newAssoc = Join[assoc, <|"Identity" -> identity, "Benchmark" -> benchmark|>];
        scenario[newAssoc]
    ];

completeScenario[x_] := x;   (* pass-through for non-scenario values *)

(* --- Constructor --- *)

makeScenario[rawAssoc_Association] :=
    Module[{wrapped, validated, completed},
        wrapped   = scenario[rawAssoc];
        validated = validateScenario[wrapped];
        If[FailureQ[validated],
            validated,
            completed = completeScenario[validated];
            completed
        ]
    ];

makeScenario[_] :=
    Failure["ScenarioValidation",
        <|"Message" -> "makeScenario requires an Association as input.",
          "MissingKeys" -> {}|>];

End[];
