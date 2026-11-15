(* Wolfram Language package *)
(* primitives.wl — shared symbolic primitives and runtime flags for MFGraphs.

   This is the base layer. All other sub-packages depend on this one.

   Provides:
     - Symbolic heads: j, u. Both are Protected: primitives` sits ahead of
       Global` on $ContextPath, so without protection a notebook assignment
       such as u = 5 would silently rewrite every variable the package builds.
     - Runtime flag: $MFGraphsVerbose
     - Gated print helper: mfgPrint
*)

BeginPackage["primitives`"]

(* Lets a re-evaluation of this file reset the usage messages below. *)
Unprotect[j, u];

j::usage = "j[a, b] is the flow on edge {a,b} from a to b. j[r, i, w] is the transition flow at junction i from edge {r,i} onto edge {i,w}: summing over w gives j[r, i], and summing over r gives j[i, w].";
u::usage = "u[a, b] is the value of the value function at vertex b of edge {a,b}.";

$MFGraphsVerbose::usage =
"$MFGraphsVerbose controls whether progress and timing messages are printed. \
Default is False; set $MFGraphsVerbose = True to see them.";

mfgPrint::usage = "mfgPrint[args___] prints args only when $MFGraphsVerbose is True.";

Begin["`Private`"]

$MFGraphsVerbose = False;

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
mfgPrint[args___] := If[TrueQ[$MFGraphsVerbose], Print[args]];
(* :!CodeAnalysis::EndBlock:: *)

End[]

Protect[j, u];

EndPackage[]
