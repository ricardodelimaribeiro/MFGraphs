(* Wolfram Language package *)
(* this file looks like what we want.*)
(* maybe the parameters association was not a good idea
Parameters =
    	Module[ {g, W, V, H},
        	g = If[ beta == 0, 
                Log,
                Function[{m},
                	m^beta]
            ];
        	W = Function[{x,a},
            	a Sin[2 Pi (x + 1/4)]^2];
        	V = Function[{x},     
        		W[x, A]];
        	H = Function[{x,p,m},
            	p^2/(2 m^alpha) + V[x] - g[m]];
        	Association[{
        		"alpha" -> alpha, 
        		"beta" -> beta, 
        		"g[m]" -> g, 
          		"W[x,A]" -> W, 
          		"V[x]" -> V, 
          		"H[x,p,m]" -> H
          	}]
    	]*)
g::usage =
""    	
g = If[ beta == 0,
		Log,
        Function[{m}, m^beta]
    ];

W::usage =
""
W = Function[{y, a}, a Sin[2 Pi (y + 1/4)]^2];

V::usage = 
""
V = Function[{x, edge}, W[x, A]];

H::usage = 
""
H = Function[{xi,p,m, edge}, p^2/(2 m^alpha) + V[xi, edge] - g[m]];
