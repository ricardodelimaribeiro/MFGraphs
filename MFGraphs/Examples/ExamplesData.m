(*2 edges entrance in the middle*)
    17 -> {
        (*VL=*) {1, 2, 3}, 
        (*AM*) {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}, 
        (*DataIn*) {{2, I1}}, 
        (*FinalCost*) {{1, U1}, {3,U2}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {3, 2, 1, S2}}},
    
     (*2 edges*)
    18 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){ {1, 2, 3, S1}, {3, 2, 1, S2}}},
        
    (*Y 1-in 2-out 4 vertices with switching*)
    19 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1},{4, I2}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {1, 2, 4, S2}, {3, 2, 1, S3}, {3, 2, 4, S4}, {4, 2, 1, S5}, {4, 2, 3, S6}}}
 (* 
    (*Y 1-in 2-out 3 vertices no switching*)
   18 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{2, U1}, {3, U2}},
        (*SwitchingCostsData=*){}}
      
    19 -> {{}, {}, {}, {}, {}},
    
 
    20 -> {
     	{1}, 
     	{{0}}, 
     	{{}}, 
     	{{1, U1}}, 
     	{}}
  
        
        *)      
];
    

DataG[n_] :=
    AssociationThread[{
    "Vertices List", 
    "Adjacency Matrix", 
    "Entrance Vertices and Currents", 
    "Exit Vertices and Terminal Costs", 
    "Switching Costs"}, 
    Test[
        n
    ]
    ](*To have a visual of the available examples, run: Table[DataToGraph[DataG[n]],{n,1,Length[Test]}]*)
                        
          


          

