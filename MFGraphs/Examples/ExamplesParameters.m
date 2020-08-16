(* Wolfram Language package *)
(* this file looks like what we want.*)
   alpha = 0.5;
   beta = 1;
Parameters =
    Module[ {alpha, beta, g, W, V, H},
        g = If[ beta == 0, 
                Log,
                Function[{m},m^beta]
            ];
        W = Function[{x,A},
            A Sin[2 Pi (x + 1/4)]^2];
        V = Function[{x},     
        	W[x, 1/2]];
        H = Function[{x,p,m},
            p^2/(2 m^alpha) + V[x] - g[m]];
        Association[{"alpha" -> alpha, "beta" -> beta, "g[m]" -> g, 
          "W[x,A]" -> W, "V[x]" -> V, "H[x,p,m]" -> H}]
    ]
