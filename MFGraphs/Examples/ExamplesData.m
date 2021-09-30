test = Association[
    (*One vertex*)
    1 -> {
        (*VL=*){1},
        (*AM=*){{0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{1, U1}},
        (*SwitchingCostsData=*){}},
        
    (*1 edge*)
    2 -> {
        (*VL=*){1, 2},
        (*AM=*){{0, 1}, {0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{2, U1}},
        (*SwitchingCostsData=*){}},
        
    (*2 edges*)
    3 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){}},
    
    (*3 edges*)
   4 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){}},
    
    (*4 edges*)    
    5 -> {
        (*VL=*){1, 2, 3, 4, 5},
        (*AM=*){{0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, 
        	    {0, 0, 0, 1, 0}, {0, 0, 0, 0 ,1}, {0, 0, 0, 0 ,0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{5, U1}},
        (*SwitchingCostsData=*){}},
   
   (*9 edges*)
    6 -> {
        (*VL=*){1,2,3,4,5,6,7,8,9,10},
        (*AM=*){
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
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{10, U1}},
        (*SwitchingCostsData=*){}},     
       
        
     (**********)
     
        
    (*Y 1-in 2-out 4 vertices no switching*)
    7 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}, {4, U2}},
        (*SwitchingCostsData=*){}},
        
    (*Y 1-in 2-out 4 vertices with switching*)
    8 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}, {4, U2}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {1, 2, 4, S2}, {3, 2, 1, S3}, {3, 2, 4, S4}, {4, 2, 1, S5}, {4, 2, 3, S6}}},
        
    (*Y 2-in 1-out 4 vertices no swithing*)
    9 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}, {3, I2}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){}},
        
    (*Y 2-in 1-out 4 vertices with switching*)
    10 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}, {3, I2}},
        (*FinalCosts=*){{4, U1}},
        (*SwitchingCostsData=*){{1, 2, 4, S1}, {1, 2, 3, S2}, {3, 2, 1, S3}, {3, 2, 4, S4}, {4, 2, 1, S5}, {4, 2, 3, S6}}},
   
   (********)
   
        
    (*Attraction Problem*)
    11 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{4, U2}}, 
        (*SwitchingCostsData=*){
            {1, 2, 4, S1}, 
            {1, 2, 3, S2}, 
            {3, 2, 1, S3},
            {3, 2, 4, S4}, 
            {4, 2, 1, S5}, 
            {4, 2, 3, S6}, 
            {1, 3, 4, S7}, 
            {4, 3, 1, S8}, 
            {1, 3, 2, S9}, 
            {2, 3, 1, S10}, 
            {3, 4, 2, S11}, 
            {2, 4, 3, S12}, 
            {2, 3, 4, S13}, 
            {4, 3, 2, S14}, 
            {3, 1, 2, S15}, 
            {2, 1, 3, S16}}},
   
   (* Attraction without switching*)
    12 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{4, U1}}, 
        (*SwitchingCostsData=*){}},
   
    (* Attraction without edge in the middle*)
    13 -> {
    	(*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 1}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{4, U1}}, 
        (*SwitchingCostsData=*){{1,2,4,S1},{1,3,4,S2},{4,2,1,S3},{4,3,1,S4}}}/. {I1 -> 1, U1 -> 0, S1 -> 0, S2 -> 0, S3 -> 0, S4 -> 0},
   
   (*********)
 
   (*triangle*)
    14 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {3, 2, 1, S2}}},
            
        
    (*Y 2-in 1-out 3 vertices with switching*)
    15 -> {
        (*VL=*){1, 2, 3},
        (*AM=*){{0, 0, 0}, {1, 0, 0}, {1, 0, 0}},
        (*DataIn=*){{2, I1}, {3, I2}},
        (*FinalCosts=*){{1, U1}},
        (*SwitchingCostsData=*){{2, 1, 3, S1}, {3, 1, 2,S2}}},
        
              
    (*Error in the switching costs*)
    16 -> {
        (*VL*){1, 2, 3}, 
        (*AM*){{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}, 
        (*DataIn*){{1, I1}}, 
        (*FinalCosts*){{3, U1}}, 
        (*SwiotchingCostsData*){{1, 3, 2, S1}}},
        

  
       
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
     
         20 -> {
        (*VL=*){1,2,3,4,5,6,7,8,9},
        (*AM=*){
        {0,0,1,0,0,0,0,0,0},
        {0,0,1,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0},
        {0,0,0,0,1,0,0,0,0},
        {0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,0,1,1,1},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{7, U1},{8, U2},{9, U3}},
        (*SwitchingCostsData=*){}}, 
        
         21 -> {
        (*VL=*){1,2,3,4,5,6,7,8,9,10,11,12},
        (*AM=*){
        {0,0,1,0,0,0,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0,0,0,0},
        {0,0,0,1,1,0,0,0,0,0,0,0},
        {0,0,0,0,0,1,0,0,0,0,0,0},
        {0,0,0,0,0,1,1,0,0,0,0,0},
        {0,0,0,0,0,0,0,1,0,0,0,0},
        {0,0,0,0,0,0,0,1,1,1,0,0},
        {0,0,0,0,0,0,0,0,1,0,1,0},
        {0,0,0,0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{10, U1},{11, U2},{12, U3}},
        (*SwitchingCostsData=*){}},     
        
        22 -> {
        (*VL=*){1,2,3,4,5,6,7},
        (*AM=*){
        {0,1,1,0,0,0,0},
        {0,0,0,1,0,0,0},
        {0,0,0,1,1,0,0},
        {0,0,0,0,0,1,0},
        {0,0,0,0,0,1,1},
        {0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{5, U1},{6, U2},{7, U3}},
        (*SwitchingCostsData=*){}},
       
       "Jamaratv9"->
       {
        (*VL=*){1,2,3,4,5,6,7,8,9},
        (*AM=*){
        {0,1,1,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0},
        {0,0,0,1,1,0,0,0,0},
        {0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,1,1,0,0},
        {0,0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0}
        },
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{7, U1},{8, U2},{9, U3}},
        (*SwitchingCostsData=*){}},
       
       
       
       
            23 -> {
        (*VL=*){1,2,3,4,5,6},
        (*AM=*){
        {0,1,1,0,0,0},
        {0,0,1,0,0,0},
        {0,0,0,1,1,0},
        {0,0,0,0,1,1},
        {0,0,0,0,0,1},
        {0,0,0,0,0,0}},
        (*DataIn=*){{1, I1},{2, I2}},
        (*FinalCosts=*){{4, U1},{5, U2},{6, U3}},
        (*SwitchingCostsData=*){}},
   
    (*Y 1-in 2-out 4 vertices with switching*)
    19 -> {
        (*VL=*){1, 2, 3, 4},
        (*AM=*){{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}},
        (*DataIn=*){{1, I1},{4, I2}},
        (*FinalCosts=*){{3, U1}},
        (*SwitchingCostsData=*){{1, 2, 3, S1}, {1, 2, 4, S2}, {3, 2, 1, S3}, {3, 2, 4, S4}, {4, 2, 1, S5}, {4, 2, 3, S6}}},
        27 -> {
        (*VL=*){1, 2},
        (*AM=*){{0, 1}, {1, 0}},
        (*DataIn=*){{1, I1}},
        (*FinalCosts=*){{2, U1}},
        (*SwitchingCostsData=*){}},
        
        "Braess split" -> {
        (*VL=*){1,2,3,4,5,6,7,8},
        (*AM=*){
        {0,1,1,0,0,0,0,0},
        {0,0,0,1,0,0,0,0},
        {0,0,0,0,1,0,0,0},
        {0,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->2},
        (*FinalCosts=*){{8, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,2,4,S1},{5,7,8,S2}}/.{S1->10,S2->10}},
        
        "Braess congest" -> {
        (*VL=*){1,2,3,4,5,6,7},
        (*AM=*){
        {0,1,1,0,0,0,0},
        {0,0,0,1,0,0,0},
        {0,0,0,1,0,0,0},
        {0,0,0,0,1,1,0},
        {0,0,0,0,0,0,1},
        {0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->2},
        (*FinalCosts=*){{7, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,2,4,S1},{4,6,7,S2}}/.{S1->10,S2->10}},
        
        "New Braess" -> {
        (*VL=*){1,2,3,4,5,6},
        (*AM=*){
        {0,1,0,1,0,0},
        {0,0,1,0,0,0},
        {0,0,0,0,0,1},
        {0,0,0,0,1,0},
        {0,0,0,0,0,1},
        {0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->4000},
        (*FinalCosts=*){{6, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,2,3,S1}, {3,2,1,S1},{4,5,6,S2}, {6,5,4,S2}}/.{S1->45,S2->45},
        (*a=*)Function[{j,edge},
        		Which[
        			edge === DirectedEdge[3,6] || edge === DirectedEdge[1,4],
        				j/100,
        			True, 0
        		]
        ]
        },
        
        "Big Braess split" -> {
        (*VL=*){1,2,3,4,5,6,7,8,9,10},
        (*AM=*){
        {0,1,1,0,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0,0},
        {0,0,0,0,0,1,0,0,0,0},
        {0,0,0,0,1,0,0,0,0,0},
        {0,0,0,0,0,0,1,0,0,0},
        {0,0,0,0,0,0,0,1,0,0},
        {0,0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->2},
        (*FinalCosts=*){{10, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,3,6,S1},{5,7,10,S2}}/.{S1->10,S2->10}},
        
        "Big Braess congest" -> {
        (*VL=*){1,2,3,4,5,6,7,8,9},
        (*AM=*){
        {0,1,1,0,0,0,0,0,0},
        {0,0,0,1,0,0,0,0,0},
        {0,0,0,0,1,0,0,0,0},
        {0,0,0,0,1,0,0,0,0},
        {0,0,0,0,0,1,1,0,0},
        {0,0,0,0,0,0,0,1,0},
        {0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,1},
        {0,0,0,0,0,0,0,0,0}},
        (*DataIn=*){{1, I1}/.I1->2},
        (*FinalCosts=*){{9, U1}/.U1->0},
        (*SwitchingCostsData=*){{1,3,5,S1},{5,7,9,S2}}/.{S1->10,S2->10}},
        
        "Grid1020" -> {
        	(*VL=*)VertexList[GridGraph[{10,20}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{10,20},DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{200, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        }(*,
        
        "Gridemen" -> {
        	(*VL=*)VertexList[GridGraph[{eme,ene}]],
        (*AM=*)AdjacencyMatrix[GridGraph[{eme,ene},DirectedEdges->True]],
        (*DataIn=*){{1, I1}/.I1->400},
        (*FinalCosts=*){{eme*ene, U1}/.U1->0},
        (*SwitchingCostsData=*){}
        }*)
        
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
Which[
	Length[test[[n]]] === 5,
    AssociationThread[{
    "Vertices List", 
    "Adjacency Matrix", 
    "Entrance Vertices and Currents", 
    "Exit Vertices and Terminal Costs", 
    "Switching Costs"}, 
    test[
        n
    ]
    ],
    Length[test[[n]]] === 6,
    AssociationThread[{
    "Vertices List", 
    "Adjacency Matrix", 
    "Entrance Vertices and Currents", 
    "Exit Vertices and Terminal Costs", 
    "Switching Costs","a"},
    test[
        n
    ]
    ],
    Length[test[[n]]] === 7,
    AssociationThread[{
    "Vertices List", 
    "Adjacency Matrix", 
    "Entrance Vertices and Currents", 
    "Exit Vertices and Terminal Costs", 
    "Switching Costs","a",
    "alpha"},
    test[
        n
    ]
    ]
    (*To have a visual of the available examples, run: Table[DataToGraph[DataG[n]],{n,1,Length[Test]}]*)
]                     
          



          

