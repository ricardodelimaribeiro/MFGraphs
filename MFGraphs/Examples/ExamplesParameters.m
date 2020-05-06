(* Wolfram Language package *)
(*TODO turn this into an association*)
alpha = 0;

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
    p^2/(2 m^alpha) + V[x] - g[m];

