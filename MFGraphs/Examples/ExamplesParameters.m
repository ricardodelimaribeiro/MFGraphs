(* Wolfram Language package *)
(*TODO turn this into an association DONE *)



Parameters := Module[{alpha, beta, g, W, V, H},
 alpha = 1;
 beta = 2;
 g[m_] := If[beta == 0, Log[m], m^beta];
 W[x_, A_] := A Sin[2 Pi (x + 1/4)]^2;
 V[x_] := W[x, 1/2];
 H[x_, p_, m_] := p^2/(2 m^alpha) + V[x] - g[m];
 Association[{"alpha" -> alpha, "beta" -> beta, "g[m]" -> g[m], 
   "W[x,A]" -> W[x, A], "V[x]" -> V[x], "H[x,p,m]" -> H[x, p, m]}]]
   
(*alpha = 0;

beta = 4;

g[m_] :=
    If[ beta == 0,
        Log[m],
        m^beta
    ]
W[x_, A_] :=
    A Sin[2 Pi (x + 1/4)]^2;
V[x_] :=
    W[x, 1/2];
H[x_, p_, m_] :=
    p^2/(2 m^alpha) + V[x] - g[m];*)

