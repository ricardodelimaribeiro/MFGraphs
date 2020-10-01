(* Wolfram Language package *)

Test = Association[{
    (*One vertex*)
    1 -> {(*VL=*){1},(*AM=*){{0}},(*DataIn=*){{1, 
        I1}},(*FinalCosts=*){{1, U1}},(*SwitchingCostsData=*){}},
    (*1 edge*)
    2 -> {(*VL=*){1, 
       2},(*AM=*){{0, 1}, {0, 0}},(*DataIn=*){{1, 
        I1}},(*FinalCosts=*){{2, U1}},(*SwitchingCostsData=*){}},
    (*Y 1-in 2-out 3 vertices no switching*)
    3 -> {(*VL=*){1, 2, 
       3},(*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},(*DataIn=*){{1, 
        I1}},(*FinalCosts=*){{2, U1}, {3, 
        U2}},(*SwitchingCostsData=*){}},
    (*Y 1-in 2-out 3 vertices with switching*)
    4 -> {(*VL=*){1, 2, 
       3},(*AM=*){{0, 0, 0}, {1, 0, 1}, {0, 0, 0}},(*DataIn=*){{2, 
        I1}},(*FinalCosts=*){{1, U1}, {3, 
        U2}},(*SwitchingCostsData=*){{1, 2, 3, S1}, {3, 2, 1, S2}}},
    (*Y 2-in 1-out 3 vertices no switching*)
    5 -> {(*VL=*){1, 2, 
       3},(*AM=*){{0, 0, 0}, {1, 0, 0}, {1, 0, 0}},(*DataIn=*){{2, 
        I1}, {3, I2}},(*FinalCosts=*){{1, 
        U1}},(*SwitchingCostsData=*){}},
    (*Y 2-in 1-out 3 vertices with switching*)
    6 -> {(*VL=*){1, 2, 
       3},(*AM=*){{0, 0, 0}, {1, 0, 0}, {1, 0, 0}},(*DataIn=*){{2, 
        I1}, {3, I2}},(*FinalCosts=*){{1, 
        U1}},(*SwitchingCostsData=*){{1, 2, 3, S1}, {3, 2, 1, S2}}},
    (*Y 1-in 2-out 4 vertices no switching*)
    7 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, I1}},(*FinalCosts=*){{3, U1}, {4, 
        U2}},(*SwitchingCostsData=*){}},
    (*Y 1-in 2-out 4 vertices with switching*)
    8 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, 2}},(*FinalCosts=*){{3, 5}, {4, 
        2}},(*SwitchingCostsData=*){{1, 2, 4, 4}, {1, 2, 3, 0}, {3, 2,
         1, 10}, {3, 2, 4, 10}, {4, 2, 1, 10}, {4, 2, 3, 10}}},
    (*Y 2-in 1-out 4 vertices no swithing*)
    9 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, I1}, {3, I2}},(*FinalCosts=*){{4, 
        U1}},(*SwitchingCostsData=*){}},
    (*Y 2-in 1-out 4 vertices with switching*)
    10 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, 1}, {3, 2}},(*FinalCosts=*){{4, 
        1}},(*SwitchingCostsData=*){{1, 2, 4, 1}, {1, 2, 3, 2}, {3, 2,
         1, 0}, {3, 2, 4, 1}, {4, 2, 1, 2}, {4, 2, 3, 0}}},
    (*Attraction Problem*)
    11 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, 0, 
        0}},(*DataIn=*){{1, I1}},(*FinalCosts=*){{4, 
        U2}}, (*SwitchingCostsData=*){{1, 2, 4, S1}, {1, 2, 3, 
        S2}, {3, 2, 1, S3}, {3, 2, 4, S4}, {4, 2, 1, S5}, {4, 2, 3, 
        S6}, {1, 3, 4, S7}, {4, 3, 1, S8}, {1, 3, 2, S9}, {2, 3, 1, 
        S10}, {3, 4, 2, S11}, {2, 4, 3, S12}, {2, 3, 4, S13}, {4, 3, 
        2, S14}, {3, 1, 2, S15}, {2, 1, 3, S16}}},
    (*Error in the switching costs*)
    12 -> {{1, 2, 
       3}, {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}, {{1, I1}}, {{3, 
        U1}}, {{1, 3, 2, S1}}},
    13 -> {{}, {}, {}, {}, {}},
    (*No entry currents*)
    14 -> {{1}, {{0}}, {{}}, {{1, U1}}, {}},
    (*Attraction without switching and no edge in the midle*)
    15 -> {{1, 2, 3, 
       4}, {{0, 1, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 
        0}}, {{1, I1}}, {{4, U1}}, {}},
    16 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, 0, 
        0}},(*DataIn=*){{1, I1}},(*FinalCosts=*){{4, 
        U2}}, (*SwitchingCostsData=*){}},
    (*Y 1-in 2-out 4 vertices no switching*)
    17 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, 2}},(*FinalCosts=*){{3, 2}, {4, 
        1}},(*SwitchingCostsData=*){}},
    (*Y 1-in 2-out 4 vertices no switching*)
    18 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, 2}},(*FinalCosts=*){{3, 5}, {4, 
        2}},(*SwitchingCostsData=*){}},
    (*Y 1-in 2-out 4 vertices no switching*)
    19 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, 2}},(*FinalCosts=*){{3, 2}, {4, 
        5}},(*SwitchingCostsData=*){}},
    (*Y 1-in 2-out 4 vertices no switching*)
    20 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, 2}},(*FinalCosts=*){{3, 2}, {4, 
        3}},(*SwitchingCostsData=*){}},
    (*Y 1-in 2-out 4 vertices no switching*)
    21 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, 1}},(*FinalCosts=*){{3, 2}, {4, 
        2.5}},(*SwitchingCostsData=*){}},
    22 -> {(*VL=*){1, 2, 
       3},(*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},(*DataIn=*){{1, 
        I1}},(*FinalCosts=*){{3, U1}},(*SwitchingCostsData=*){{1, 2, 
        3, S1}, {3, 2, 1, S2}}},
    23 -> {(*VL=*){1, 2, 
       3},(*AM=*){{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},(*DataIn=*){{1, 
        I1}},(*FinalCosts=*){{3, U1}},(*SwitchingCostsData=*){{1, 2, 
        3, S1}, {3, 2, 1, S2}}},
    (*Y 2-in 1-out 4 vertices no swithing*)
    24 -> {(*VL=*){1, 2, 3, 
       4},(*AM=*){{0, 1, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 
        0}},(*DataIn=*){{1, 1}, {3, 2}},(*FinalCosts=*){{4, 
        1}},(*SwitchingCostsData=*){}},
    
    (*2 edges*)
    25 -> {(*VL=*) {1, 2, 3}, 
    	(*AM*) {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}, 
    	(*DataIn*) {{2, 2}}, 
    	(*FinalCost*) {{1, 5}, {3,5}},
    	(*SwitchingCostsData=*){}},
    
    (*2 edges*)
    28 -> {(*VL=*) {1, 2, 3}, 
    	(*AM*) {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}, 
    	(*DataIn*) {{2, 2}}, 
    	(*FinalCost*) {{1, 4}, {3, 5}},
    	(*SwitchingCostsData=*){}},
    
    (*2 edges*)
    29 -> {(*VL=*) {1, 2, 3}, 
    	(*AM*) {{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}, 
    	(*DataIn*) {{2, 2}}, 
    	(*FinalCost*) {{1, 5}, {3, 9}},
    	(*SwitchingCostsData=*){}},
    
    (*1 edge*)
    26 -> {(*VL=*){1,2},(*AM=*){{0, 1}, {0, 0}},
    	(*DataIn=*){{1,80}},(*FinalCosts=*){{2, 15}},(*SwitchingCostsData=*){}},
    27 -> {(*VL=*){1,2,3,4,5,6,7,8,9,10},(*AM=*){
    	{0,1,0,0,0,0,0,0,0,0},
    	{0,0,1,0,0,0,0,0,0,0},
    	{0,0,0,1,0,0,0,0,0,0},
    	{0,0,0,0,1,0,0,0,0,0},
    	{0,0,0,0,0,1,0,0,0,0},
    	{0,0,0,0,0,0,1,0,0,0},
    	{0,0,0,0,0,0,0,1,0,0},
    	{0,0,0,0,0,0,0,0,1,0},
    	{0,0,0,0,0,0,0,0,0,1},
    	{0,0,0,0,0,0,0,0,0,0}},
    	(*DataIn=*){{1,80}},(*FinalCosts=*){{10, 15}},(*SwitchingCostsData=*){}}
    }];
    

DataG[n_] :=
    AssociationThread[{
    "Vertices List", 
    "Adjacency Matrix", 
    "Entrance Vertices and Currents", 
    "Exit Vertices and Terminal Costs", 
    "Switching Costs"}, 
    Test[n]
    ](*To have a visual of the available examples, run: Table[DataToGraph[DataG[n]],{n,1,Length[Test]}]*)
                        
          
