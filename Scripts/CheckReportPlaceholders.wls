#!/usr/bin/env wolframscript
(* Fails when committed history/report files still contain placeholder prose. *)

targets = {
  "docs/history/DNF_PERFORMANCE_HISTORY.md",
  "docs/history/PARALLEL_PERFORMANCE_HISTORY.md",
  "History/DEVELOPMENT_HISTORY.md",
  "Results/compare_dnf.md",
  "Results/preprocessing_profile_report.md",
  "Results/reduction_strategies_profile.md"
};

placeholderPatterns = {
  "_TODO_",
  "_Describe the code or benchmark configuration change being measured in this run._",
  "_Summarize the main performance signal, regressions, or follow-up actions from this run._"
};

todoFiles = Select[
  targets,
  FileExistsQ[#] && StringContainsQ[Import[#, "Text"], Alternatives @@ placeholderPatterns] &
];

If[todoFiles === {},
  Print["No report/history placeholder tokens found."];
  Exit[0]
];

Print["Found placeholder prose in committed report/history files:"];
Scan[Print["  ", #] &, todoFiles];
Exit[1];
