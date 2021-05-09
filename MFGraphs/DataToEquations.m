(* Wolfram Language package *)

(*Assembles the stationary mfg equations using the notion of current*)
DataToEquations::usage = 
"DataToEquations[<|\"Vertices List\" -> {1, 2, 3}, \"Adjacency Matrix\" -> {{0, 0, 0}, {1, 0, 1}, {0, 0, 0}}, 
 \"Entrance Vertices and Currents\" -> {{2, I1}}, 
 \"Exit Vertices and Terminal Costs\" -> {{1, U1}, {3, U2}}, 
 \"Switching Costs\" -> {{1, 2, 3, S1}, {3, 2, 1, S2}}|>] returns the equations for the stationary mean-field game on the network, except for the coupling equations.";

Begin["`Private`"]

DataToEquations[Data_?AssociationQ] :=
    Module[ {showAll = True, BG, EntranceVertices, InwardVertices, InEdges, ExitVertices, 
      OutwardVertices, OutEdges, AuxiliaryGraph, FG, VL, EL, BEL, jargs, 
      js, jvars, FVL, AllTransitions, jts, jtvars, uargs, us, uvars, 
      EqPosCon, EqCurrentCompCon, EqTransitionCompCon, NoDeadEnds, 
      NoDeadStarts, EqBalanceSplittingCurrents, EqBalanceGatheringCurrents, EntryArgs, 
      EntryDataAssociation, EqEntryIn, NonZeroEntryCurrents, ExitCosts, EqExitValues, SwitchingCosts, OutRules, InRules, 
      EqSwitchingConditions, EqCompCon, EqValueAuxiliaryEdges, EqAllComp, EqAll, SignedCurrents, Nrhs, Nlhs, RulesEntryIn, MinimalTimeRhs,
      RulesExitValues, EqAllRules, EqAllCompRules, EqAllAll, EqAllAllRules, reduced1,reduced2,
      EqCriticalCase, BoundaryRules, criticalreduced1, criticalreduced2, EqGeneralCase, EqCcs, (*EqAllAllSimple,*) sol, system, rules ,TOL, time,nocasesystem,nocaserules(*, result*)},
      (*Set the tolerance for Chop*)
      TOL = 10^(-6);
      (*Checking consistency on the swithing costs*)
      ConsistentSwithingCosts[{a_, b_, c_, S_}] :=
      	Module[{sc = Data["Switching Costs"], or, de, bounds},
  			or = Cases[sc, {a, b, _, _}];
  			de = Cases[sc, {_, b, c, _}];
  			bounds = Outer[Plus, Last /@ or, Last /@ de] // Flatten;
  			And @@ (NonNegative[#-S] &) /@ bounds
  		];
  		EqCcs = ConsistentSwithingCosts /@ Data["Switching Costs"];
  		Print["DataToEquations: Triangle inequalities for switching costs: ", EqCcs,
  			"\nDataToEquations: Reduced: ",Reduce[EqCcs]];
  		If[And @@ EqCcs,
  			(*this happens when we have non-false condition.*)
  			Print["DataToEquations: The switching costs are compatible."],
  			Print["DataToEquations: The switching costs are ",Style["incompatible", Red],". \nStopping!"];
  			Return[]
  		];
      (*Begin Internal functions for DataToEquations: *)
        IncomingEdges[k_] :=
            {k, #1} & /@ IncidenceList[FG,k]; (*all edges "oriented" towards the vertex k*)
        OutgoingEdges[k_] :=
            OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]); (*all edges "oriented" away from the vertex k*)
        CurrentCompCon[a_ \[DirectedEdge] b_] :=
            jvars[{a, a \[DirectedEdge] b}] == 0 || jvars[{b, a \[DirectedEdge] b}] == 0;
        (*CurrentCompCon[edge_] := 
        	jvars[AtTail[edge]] || jvars[AtHead[edge]];*)
        CurrentSplitting[{c_, a_ \[DirectedEdge] b_}] :=
            Select[AllTransitions, (Take[#, 2] == {c, a \[DirectedEdge] b}) &];
        CurrentGathering[{c_, a_ \[DirectedEdge] b_}] :=
            Select[AllTransitions, (Part[#, {1, 3}] == OtherWay[{c, a \[DirectedEdge] b}]) &];
        TransitionCompCon[{v_, edge1_, edge2_}] :=
            jtvars[{v, edge1, edge2}] == 0 || jtvars[{v, edge2, edge1}] == 0;
        ExitValues[a_ \[DirectedEdge] b_] :=
            Total[uvars /@ {{b, DirectedEdge[a, b]}}] == ExitCosts[b];
        Transu[{v_, edge1_, edge2_}] :=
        	NonNegative[uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}]];(*u1<= u2+S(1,2) for agents going from edge 1 to 2 at some vertex*)
        Compu[{v_, edge1_, edge2_}] :=
            (jtvars[{v, edge1, edge2}] == 0) || uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}] == 0;
        a(*[SignedCurrents[edge], edge]:*)=
        	Lookup[Data, "a" , Function[{j, edge}, j]];(*Default corresponds to the classic critical congestion case*)
        (*End*)
       
        BG = AdjacencyGraph[Data["Vertices List"], Data["Adjacency Matrix"], VertexLabels -> "Name", DirectedEdges -> True, GraphLayout -> "SpringEmbedding"];
        EntranceVertices = First /@ Data["Entrance Vertices and Currents"];
        InwardVertices = AssociationThread[EntranceVertices, Unique["en"] & /@ EntranceVertices]; (*Defines auxiliary vertices for the entrance vertices*)
        InEdges = MapThread[DirectedEdge, {InwardVertices /@ EntranceVertices, EntranceVertices}];
        ExitVertices = First /@ Data["Exit Vertices and Terminal Costs"];
        OutwardVertices = AssociationThread[ExitVertices, Unique["ex"] & /@ ExitVertices];
        OutEdges = MapThread[DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
        AuxiliaryGraph = Graph[Join[InEdges, OutEdges], VertexLabels -> "Name", GraphLayout -> "SpringEmbedding"];
        FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
        VL = Data["Vertices List"];
        EL = EdgeList[FG];
        BEL = EdgeList[BG];
        jargs = Join[AtHead /@ EL, AtTail /@ EL];
        js = Unique["j"] & /@ jargs;
        jvars = AssociationThread[jargs, js];
        FVL = VertexList[FG];
        SignedCurrents =   AssociationThread[BEL, (jvars[AtHead[#]] - jvars[AtTail[#]] &) /@ BEL];
        AllTransitions = TransitionsAt[FG, #] & /@ FVL // Catenate (*at vertex from first edge to second edge*);
        jts = Unique["jt"] & /@ AllTransitions;
        jtvars = AssociationThread[AllTransitions, jts];
        uargs = Join[AtTail /@ EL, AtHead /@ EL];
        us = Unique["u"] & /@ uargs;
        uvars = AssociationThread[uargs, us];
        EqPosCon = And @@ (NonNegative /@ Join[jvars, jtvars]);
        EqCurrentCompCon = And @@ (CurrentCompCon /@ EL);
        EqTransitionCompCon = And @@ ((Sort /@ TransitionCompCon /@ AllTransitions) // Union);
        NoDeadEnds = IncomingEdges /@ VL // Flatten[#, 1] &;
        NoDeadStarts = OutgoingEdges /@ VL // Flatten[#, 1] &;
        EqBalanceSplittingCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentSplitting[#]]) & /@ NoDeadEnds);
        EqBalanceGatheringCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentGathering[#]]) & /@ NoDeadStarts);
        EqAll = EqPosCon && EqBalanceSplittingCurrents && EqBalanceGatheringCurrents;(*takes too long to Reduce, so don't!*)
        EntryArgs = AtHead /@ ((EdgeList[AuxiliaryGraph, _ \[DirectedEdge] #] & /@ (First /@ Data["Entrance Vertices and Currents"])) // Flatten[#, 1] &);
        EntryDataAssociation = AssociationThread[EntryArgs, Last /@ Data["Entrance Vertices and Currents"]];
        EqEntryIn = And @@ ((jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ InEdges));
        RulesEntryIn = ToRules[EqEntryIn];
        NonZeroEntryCurrents = And @@ (Positive[EntryDataAssociation[#]] & /@ (AtHead /@ InEdges)); (*useful for the general case...*)
        EqAll = EqAll && EqEntryIn;
        ExitCosts = AssociationThread[OutwardVertices /@ (First /@ Data["Exit Vertices and Terminal Costs"]), Last /@ Data["Exit Vertices and Terminal Costs"]];
        EqExitValues = And @@ (ExitValues /@ IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices]);
        RulesExitValues = ToRules[EqExitValues];
        EqAll = EqAll && EqExitValues;
        SwitchingCosts = AssociationThread[Join[AllTransitions, triple2path[Take[#, 3], FG] & /@ Data["Switching Costs"]], Join[Table[0, Length[AllTransitions]], Last[#] & /@ Data["Switching Costs"]]];(*Swithing cost is initialized with 0 AssociationThread associates the last association!*)
        (*Infinite switching costs here prevent the network from sucking agents from the exits.*)
        OutRules = Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, EL] // Flatten[#, 1] &);
        AssociateTo[SwitchingCosts, Association[OutRules]];
        InRules = Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // Flatten[#, 2] &);
        AssociateTo[SwitchingCosts, Association[InRules]];
        EqSwitchingConditions = Reduce[(And @@ Transu /@ AllTransitions), Reals];(*Reducing before joining with other equations*)
        EqAll = EqAll && EqSwitchingConditions && And @@ EqCcs;(*includes transition inequalities for the values and the triangle inequalities, too.*)
        EqCompCon = And @@ Compu /@ AllTransitions;
        EqValueAuxiliaryEdges = And @@ ((uvars[AtHead[#]] - uvars[AtTail[#]] == 0) & /@ EdgeList[AuxiliaryGraph]);
        EqAllComp = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        EqAll = EqAll && EqValueAuxiliaryEdges;
        (*Pre-processing: substitute rules in all equations and hand-sides...*)
        BoundaryRules = Join[RulesEntryIn, RulesExitValues];
        EqAllAll = EqAll && EqAllComp;
        EqAllAllRules = EqAllCompRules && EqAllRules;
        
        Print["DataToEquations: Finished assembling strucural equations. Reducing the structural system ... "];
               
        MinimalTimeRhs = Flatten[-a[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];

        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        Print["DataToEquations: Critical case ... "];
        EqCriticalCase = And @@ ((# == 0) & /@ Nlhs);
        
        EqMinimalTime = And @@ (MapThread[(#1 == #2) &, {Nlhs, MinimalTimeRhs}]);
        	
        {system, rules} = FixedPoint[EqEliminatorX, {(EqAllAll && EqCriticalCase), {}}];
        
        {nocasesystem, nocaserules} = FixedPoint[EqEliminatorX, {EqAllAll, {}}]; (*does this step help speed up fixedreducex1?*)
        
        {time,nocasesystem} = AbsoluteTiming @ NewReduce[nocasesystem]; (*does this step help speed up fixedreducex1?*)
        Print["DataToEquations: It took ", time, " seconds to reduce the general case with NewReduce.","\nThe system is :", nocasesystem];
 
        (*Print["DataToEquations: The system is:\n", system,
        	"\nand the rules are:\n", rules];*)
        	
        	(*There could be some alternatives and inequalities left, so we BooleanConvert and solve/replace throught the tree*)
        {time,system} = AbsoluteTiming @ NewReduce[system];
        Print["DataToEquations: It took ", time, " seconds to reduce the critical congestion case with NewReduce!","\nThe system is ", system];
        (*Print["DataToEquations: After NewReduce, the system is:\n", system,
        	"\nand the rules are:\n", rules];*)
        Which[system === True ||
        	Head[system] === Equal || 
        	(Head[system] === And && DeleteDuplicates[Head /@ system] === Equal),
        	sol = Solve[system, Reals];
        	If[Length[sol] == 1,
        		{system, rules} = {system, AssociateTo[ rules, First@ sol]} /. First@ sol,
        		Print[sol];
        	],
        	system === False,
        		Print["DataToEquations: ",Style["Incompatible system", Red],"."],
        	True,
        	Print["DataToEquations: The system does not have the original structure: ", system];
        	{system,rules} = FixedPoint[EqEliminatorX, {system, rules}
        		(*,SameTest -> (Length[#1[[2]]]===Length[#2[[2]]]&)*)];
        	system = Reduce[system, Reals];
        	Print["DataToEquations: ", Style["Multiple solutions", Red]," \n\t\t", system, "\n\tFinding one instance..."];
        	(*FindInstance!!!*)
        	Module[{solved = And @@ (Outer[Equal,Keys[rules],Values[rules]]//Diagonal) , onesol},
        		onesol = First@FindInstance[system && solved, Join[js, us, jts], Reals];	
        		(*TODO: should we choose a method to select a solution?*)
        		rules = Association[onesol]; 
        		system = system/.rules
        	];
        ];
        criticalreduced1 = {system, Simplify /@ rules};
        If[system === True,
        	Print["DataToEquations: Critical congestion solved."];,
        	Print["DataToEquations: There are ", Style["multiple solutions", Red]," for the given data."];
        ];
        Nrhs =  Flatten[IntM[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        EqGeneralCase = And @@ (MapThread[(#1 == #2) &, {Nlhs , Nrhs}]);
        
        
        (*Print["DataToEquations: Done."];*)
        If[showAll == True,
        	Association[{
        		"Data" -> Data,
          (*Graph structure*)
          "BG" -> BG, 
          "EntranceVertices" -> EntranceVertices, 
          "InwardVertices" -> InwardVertices, 
          "InEdges" -> InEdges, 
          "ExitVertices" -> ExitVertices, 
          "OutwardVertices" -> OutwardVertices, 
          "OutEdges" -> OutEdges, 
          "AuxiliaryGraph" -> AuxiliaryGraph, 
          "FG" -> FG, 
          "VL" -> VL, 
          "EL" -> EL, 
          "BEL" -> BEL, 
          "FVL" -> FVL, 
          "AllTransitions" -> AllTransitions, 
          "NoDeadEnds" -> NoDeadEnds, 
          "NoDeadStarts" -> NoDeadStarts, 
          (*variables*)
          "jargs" -> jargs, 
          "js" -> js, 
          "jvars" -> jvars, 
          "jts" -> jts, 
          "jtvars" -> jtvars, 
          "uargs" -> uargs, 
          "us" -> us, 
          "uvars" -> uvars,
          "SwitchingCosts" -> SwitchingCosts, 
          "OutRules" -> OutRules, 
          "InRules" -> InRules, 
          
          
          
          "EntryArgs" -> EntryArgs, 
          "EntryDataAssociation" -> EntryDataAssociation, 
          "jays" -> SignedCurrents, 
          "NonZeroEntryCurrents" -> NonZeroEntryCurrents, 
          "ExitCosts" -> ExitCosts, 
          
          
           
          (*equations*)
          (*complementarity*)
          "EqCurrentCompCon" -> EqCurrentCompCon, 
          "EqTransitionCompCon" -> EqTransitionCompCon, 
          "EqCompCon" -> EqCompCon, 
          "EqAllComp" -> EqAllComp, (*union of all complementarity conditions*)
          
          (*linear equations (and inequalities)*)
          "EqPosCon" -> EqPosCon, 
          "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
          "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
          "EqEntryIn" -> EqEntryIn, 
          "EqExitValues" -> EqExitValues, 
          "EqSwitchingConditions" -> EqSwitchingConditions, 
          "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
          "EqAll" -> EqAll, 
          (*"EqAllCompRules" -> EqAllCompRules,*)
          (*"EqAllRules" -> EqAllRules,*)
          "EqAllAll" -> EqAllAll,
          "BoundaryRules" -> BoundaryRules,
          (*"reduced1" -> reduced1,
          "reduced2" -> reduced2,*)
          (*"EqAllAllSimple" -> EqAllAllSimple,*)
          "RulesEntryIn"-> RulesEntryIn,
          "RulesExitValues" -> RulesExitValues,
          (*"EqAllAllRules" -> EqAllAllRules,*)
          "Nlhs" -> Nlhs,
          "nocasesystem" -> nocasesystem, 
          "nocaserules" -> nocaserules,
          "MinimalTimeRhs" -> MinimalTimeRhs,
          "EqCriticalCase" -> EqCriticalCase ,
          "criticalreduced1" -> criticalreduced1,
          (*"criticalreduced2" -> criticalreduced2,*)
          "Nrhs" -> Nrhs,
          "EqGeneralCase" -> EqGeneralCase,
          "TOL"->TOL
          }]
        	,
        Association[{
        	"Data" -> Data,
          (*Graph structure*)
          (*"BG" -> BG, 
          "EntranceVertices" -> EntranceVertices, 
          "InwardVertices" -> InwardVertices, 
          "InEdges" -> InEdges, 
          "ExitVertices" -> ExitVertices, 
          "OutwardVertices" -> OutwardVertices, 
          "OutEdges" -> OutEdges, 
          "AuxiliaryGraph" -> AuxiliaryGraph,*) 
          "FG" -> FG, (*
          "VL" -> VL, 
          "EL" -> EL, 
          "BEL" -> BEL, 
          "FVL" -> FVL, 
          "AllTransitions" -> AllTransitions, 
          "NoDeadEnds" -> NoDeadEnds, 
          "NoDeadStarts" -> NoDeadStarts, 
          (*variables*)
          "jargs" -> jargs, 
          "js" -> js, 
          "jvars" -> jvars, 
          "jts" -> jts, 
          "jtvars" -> jtvars, 
          "uargs" -> uargs, 
          "us" -> us, 
          "uvars" -> uvars,
          "SwitchingCosts" -> SwitchingCosts, 
          "OutRules" -> OutRules, 
          "InRules" -> InRules, 
          
          
          
          "EntryArgs" -> EntryArgs, 
          "EntryDataAssociation" -> EntryDataAssociation, 
          "jays" -> SignedCurrents, 
          "NonZeroEntryCurrents" -> NonZeroEntryCurrents, 
          "ExitCosts" -> ExitCosts, 
          
          
          *) 
          (*equations*)
          (*complementarity*)
(*          "EqCurrentCompCon" -> EqCurrentCompCon, 
          "EqTransitionCompCon" -> EqTransitionCompCon, 
          "EqCompCon" -> EqCompCon, 
          "EqAllComp" -> EqAllComp, (*union of all complementarity conditions*)*)
          
          (*linear equations (and inequalities)*)
          (*"EqPosCon" -> EqPosCon, 
          "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
          "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
          "EqEntryIn" -> EqEntryIn, 
          "EqExitValues" -> EqExitValues, 
          "EqSwitchingConditions" -> EqSwitchingConditions, 
          "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
          "EqAll" -> EqAll, *)
         (* "EqAllCompRules" -> EqAllCompRules,
          "EqAllRules" -> EqAllRules,*)
          "EqAllAll" -> EqAllAll,
          "BoundaryRules" -> BoundaryRules,
          "reduced1" -> reduced1,
          "reduced2" -> reduced2,
(*          "EqAllAllSimple" -> EqAllAllSimple,*)
(*          "RulesEntryIn"-> RulesEntryIn,
          "RulesExitValues" -> RulesExitValues,*)
(*          "EqAllAllRules" -> EqAllAllRules,*)
          (*"reduced system" -> reduced[[1]],
          "reducing rules" -> reduced[[2]],*)
          "Nlhs" -> Nlhs,
          "EqCriticalCase" -> EqCriticalCase ,
          "criticalreduced1" -> criticalreduced1,
          "criticalreduced2" -> criticalreduced2,
          "Nrhs" -> Nrhs,
          "EqGeneralCase" -> EqGeneralCase,
          "TOL" -> TOL
          }]
        ]
    ]

End[]