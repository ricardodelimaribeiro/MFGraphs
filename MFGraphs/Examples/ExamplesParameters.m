(* Wolfram Language package *)
(* this file looks like what we want.*)
Parameters :=
    Module[ {alpha, beta, g, W, V, H},
        alpha = 1;
        beta = 2;
        g[m_] :=
            If[ beta == 0,
                Log[m],
                m^beta
            ];
        W[x_, A_] :=
            A Sin[2 Pi (x + 1/4)]^2;
        V[x_] :=
            W[x, 1/2];
        H[x_, p_, m_] :=
            p^2/(2 m^alpha) + V[x] - g[m];
        Association[{"alpha" -> alpha, "beta" -> beta, "g[m]" -> g, 
          "W[x,A]" -> W, "V[x]" -> V, "H[x,p,m]" -> H}]
    ]