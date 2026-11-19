(* Shared guard, sourced first by every script that takes a per-solve time
   limit (not a standalone CLI).

   wolframscript acts on a --timeout (or -timeout) argument itself, even when
   the argument follows the script name: it stops the whole run after that
   many seconds, prints a stray token and exits 0. A script that read the same
   flag as its own per-solve limit could therefore never report a timeout
   (wolframscript's clock starts first, so it always fires first), and any run
   longer than the limit in total ended without saving anything.

   The scripts' own flag is --solve-timeout. This guard turns the old spelling
   into a usage error, and runs before the package loads so that it prints
   ahead of wolframscript's cut-off. *)

If[MemberQ[$ScriptCommandLine, "--timeout" | "-timeout"],
  Print["--timeout is taken by wolframscript itself: it ends the whole run after that many seconds and exits 0."];
  Print["Use --solve-timeout <seconds> for this script's per-solve limit."];
  Exit[2]
];
