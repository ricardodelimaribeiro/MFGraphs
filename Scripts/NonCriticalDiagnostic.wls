#!/usr/bin/env wolframscript
(* NonCriticalDiagnostic.wls -- Convergence diagnostic for the iterative
   non-critical (alpha != 1) solver on the triangular network of the paper
   "An Exact Solver for Stationary Mean-Field Games on Networks"
   (Section 5.4, "Discussion: challenges and limitations").

   ----------------------------------------------------------------------------
   RECONSTRUCTION NOTE
   ----------------------------------------------------------------------------
   The exact solver in the MFGraphs package targets the critical case (alpha=1).
   The non-critical case is not part of the loaded package surface, so this
   diagnostic is SELF-CONTAINED: it implements the stationary non-critical
   fixed-point map for the triangular network directly from the model equations,
   independent of the package solver. It does not call CriticalCongestionSolver
   or any archived non-critical routine.

   This script was reconstructed to restore the artifact cited by the paper
   (the original Scripts/NonCriticalDiagnostic.wls was missing). Validate the
   numeric output locally before relying on the regenerated table; the last
   digits of the residuals are tolerance-dependent and may differ slightly from
   the committed figure.

   ----------------------------------------------------------------------------
   Model (triangular network, equal potential on all edges)
   ----------------------------------------------------------------------------
   Triangle v1,v2,v3: entry at v1 (inflow iota1), exit at v3 (cost phi3=0),
   zero switching costs. Two routes: 1->2->3 and 1->3.

   On each active edge carrying constant net flow j, the stationary MFG gives
       j^2 / (2 m(x)^(2-alpha)) + V(x) = g(m(x)),          (density)
       u_x = -j m(x)^(alpha-1),                            (transport)
   so the value drop across an edge is
       Delta u = j * R(j),   R(j) := Integrate[ m(x;j)^(alpha-1), {x,0,1} ].
   For alpha=1, R == 1 and Delta u = j (the critical, linear case).

   With equal potential on all edges R is the same function on each edge.
   Writing j12 = j23 = a and j13 = iota1 - a, route indifference (zero
   switching, both routes used) reads
       2 a R(a) = (iota1 - a) R(iota1 - a).
   The iterative solver freezes R at the current iterate and updates
       a_{k+1} = iota1 * R(iota1 - a_k) / ( 2 R(a_k) + R(iota1 - a_k) ),
   seeded from the critical solution a_0 = iota1/3. We record the flow-vector
   residual ||j_k - j_{k-1}||_inf = |a_k - a_{k-1}| at each step.

   Outputs:
     Results/noncritical_residuals.csv                (raw data)
     <paper>/figures/noncritical_residuals.tex        (LaTeX table, \input by paper)
     <paper>/figures/noncritical_residuals.pdf        (log-residual plot)

   Usage: wolframscript -file Scripts/NonCriticalDiagnostic.wls
*)

(* ----------------------------------------------------------------------------
   Parameters (Section 5.4)
   ---------------------------------------------------------------------------- *)
iota1   = 80;                                            (* inflow at v1 *)
Vfun    = Function[x, (1/2) Sin[2 Pi (x + 1/4)]^2];      (* potential V(x) *)
gfun    = Function[m, m];                                (* congestion g(m) = m *)
alphas  = {1.0, 1.2, 1.5};                               (* congestion exponents *)
maxIter = 8;                                             (* iteration cap *)
tol     = 10.^-12;                                       (* convergence tolerance *)

(* ----------------------------------------------------------------------------
   Density and edge resistance
   ---------------------------------------------------------------------------- *)
(* densityAt[j, x, alpha]: positive root m of  j^2/(2 m^(2-alpha)) + V(x) = g(m). *)
densityAt[j_?NumericQ, x_?NumericQ, alpha_?NumericQ] :=
  Module[{m, roots},
    If[Chop[j] == 0, Return[Missing["ZeroFlow"]]];
    roots = Quiet @ Check[
      m /. NSolve[{N[j]^2/(2 m^(2 - alpha)) + Vfun[x] == gfun[m], m > 0}, m, Reals],
      $Failed];
    If[roots === $Failed || roots === {} || !ListQ[roots],
      Missing["DensityNotSolved"],
      FirstCase[N @ roots, v_?NumericQ /; v > 0, Missing["DensityNotSolved"]]]
  ];

(* R[j, alpha] = Integrate[ m(x;j)^(alpha-1), {x,0,1} ]. For alpha==1, R==1. *)
edgeResistance[j_?NumericQ, alpha_?NumericQ] :=
  If[Abs[alpha - 1] < 10.^-14,
    1.,
    Quiet @ Check[
      NIntegrate[densityAt[j, x, alpha]^(alpha - 1), {x, 0, 1},
        Method -> "GlobalAdaptive", MaxRecursion -> 12],
      $Failed]
  ];

(* ----------------------------------------------------------------------------
   Fixed-point iteration on the flow split a (= j12 = j23), seeded from critical
   ---------------------------------------------------------------------------- *)
runDiagnostic[alpha_?NumericQ] :=
  Module[{a, aNew, Ra, Rb, residuals = {}, k, res},
    a = N[iota1]/3;                                      (* critical seed *)
    For[k = 1, k <= maxIter, k++,
      Ra = edgeResistance[a, alpha];
      Rb = edgeResistance[iota1 - a, alpha];
      If[Ra === $Failed || Rb === $Failed, AppendTo[residuals, Missing["IntegrationFailed"]]; Break[]];
      aNew = iota1 * Rb / (2 Ra + Rb);
      res = Abs[aNew - a];                               (* ||j_k - j_{k-1}||_inf *)
      AppendTo[residuals, res];
      a = aNew;
      If[res < tol, Break[]];
    ];
    <|"Alpha" -> alpha, "Flow" -> a, "Residuals" -> residuals|>
  ];

results = Association[(ToString[#["Alpha"]] -> #) & /@ (runDiagnostic /@ alphas)];

(* ----------------------------------------------------------------------------
   Assemble per-iteration table (rows = iteration index, cols = alpha)
   ---------------------------------------------------------------------------- *)
maxRows = Max[Length /@ (results[#]["Residuals"] & /@ ToString /@ alphas)];

cell[alpha_, row_] := Module[{rs = results[ToString[alpha]]["Residuals"]},
  If[row <= Length[rs], rs[[row]], Missing["Converged"]]];

(* ----------------------------------------------------------------------------
   Write raw CSV
   ---------------------------------------------------------------------------- *)
repoRoot = If[StringQ[$InputFileName] && $InputFileName =!= "",
  FileNameDrop[$InputFileName, -2],                      (* .../Scripts/this.wls -> repo root *)
  Directory[]];
resultsDir = FileNameJoin[{repoRoot, "Results"}];
If[!DirectoryQ[resultsDir], CreateDirectory[resultsDir]];

csvHeader = Prepend[("alpha=" <> ToString[#]) & /@ alphas, "iter"];
csvRows = Table[
  Prepend[
    Table[With[{c = cell[a, row]}, If[MissingQ[c], "", c]], {a, alphas}],
    row],
  {row, maxRows}];
Export[FileNameJoin[{resultsDir, "noncritical_residuals.csv"}],
  Prepend[csvRows, csvHeader]];
Print["Wrote Results/noncritical_residuals.csv"];

(* ----------------------------------------------------------------------------
   Write LaTeX table (matches the format the paper \input's)
   ---------------------------------------------------------------------------- *)
fmtCell[c_] := Which[
  MissingQ[c], "--",
  Chop[c, 10.^-30] == 0, "$0$",
  True, Module[{mant, ex},
    ex = Floor[Log10[Abs[c]]];
    mant = c/10.^ex;
    "$" <> ToString[NumberForm[mant, {3, 2}]] <> " \\times 10^{" <> ToString[ex] <> "}$"]];

texLines = Join[
  {"% Auto-generated by Scripts/NonCriticalDiagnostic.wls",
   "\\begin{tabular}{r" <> StringJoin[ConstantArray["c", Length[alphas]]] <> "}",
   "\\hline",
   "iter & " <> StringRiffle[("$\\alpha=" <> ToString[NumberForm[#, {2, 1}]] <> "$") & /@ alphas, " & "] <> " \\\\",
   "\\hline"},
  Table[
    ToString[row] <> " & " <>
      StringRiffle[fmtCell[cell[#, row]] & /@ alphas, " & "] <> " \\\\",
    {row, maxRows}],
  {"\\hline", "\\end{tabular}"}];

paperFigDir = "/Users/ribeirrd/Dropbox/Apps/Overleaf/Stationary solutions of a congested MFG on Networks/figures";
If[DirectoryQ[paperFigDir],
  Export[FileNameJoin[{paperFigDir, "noncritical_residuals.tex"}], StringRiffle[texLines, "\n"], "Text"];
  Print["Wrote figures/noncritical_residuals.tex"],
  Print["Paper figures/ not found; skipped .tex (looked in: ", paperFigDir, ")"]];

(* ----------------------------------------------------------------------------
   Write log-residual plot
   ---------------------------------------------------------------------------- *)
(* One {k, log10 residual} series per alpha, keeping only positive residuals
   (log of a zero residual is undefined; alpha=1 converges immediately). *)
plotSeries[alpha_] := Select[
  MapIndexed[{First[#2], #1} &, results[ToString[alpha]]["Residuals"]],
  NumericQ[Last[#]] && Last[#] > 0 &] /. {k_, r_} :> {k, Log10[r]};

activeAlphas = Select[alphas, plotSeries[#] =!= {} &];

If[DirectoryQ[paperFigDir],
  If[activeAlphas === {},
    Print["All residuals converged at step 1; no positive residuals to plot (.pdf left unchanged)."],
    Module[{plot},
      plot = ListLinePlot[
        plotSeries /@ activeAlphas,
        PlotMarkers -> Automatic,
        Frame -> True,
        FrameLabel -> {"iteration k", Row[{Subscript["log", 10], " ||j_k - j_{k-1}||\[Infinity]"}]},
        PlotLegends -> (("\[Alpha] = " <> ToString[#]) & /@ activeAlphas),
        ImageSize -> 360];
      Export[FileNameJoin[{paperFigDir, "noncritical_residuals.pdf"}], plot];
      Print["Wrote figures/noncritical_residuals.pdf"]]],
  Print["Paper figures/ not found; skipped .pdf"]];

Print["Flows (a = j12 = j23) per alpha: ",
  Association[(# -> results[ToString[#]]["Flow"]) & /@ alphas]];
Print["Done."];
