(* Wolfram Language package *)
(*TODO use module, block or with to localize the variables! DONE*)
AtHead[a_ \[DirectedEdge] b_] := {b, a \[DirectedEdge] b};
AtTail[a_ \[DirectedEdge] b_] := {a, a \[DirectedEdge] b};
TransitionsAt[G_, k_] := 
  Prepend[#, k] & /@ Permutations[IncidenceList[G, k], {2}];
OtherWay[{c_, a_ \[DirectedEdge] b_}] := {If[c === a, b, a], 
   a \[DirectedEdge] b} (*example: \
OtherWay[{a,a\[DirectedEdge]b}]={b,a\[DirectedEdge]b} *);
triple2path[{a_, b_, c_}, 
   G_] :=(*Takes 3 ordered vertices and gives the pair of edges \
conecting the first two and the last two. If this is not feasible, 
  error message with {a, b, c}!*)
  (
   EL = EdgeList[FG]; 
   If[SubsetQ[AdjacencyList[G, b], {a, c}], 
    If[MemberQ[EL, a \[DirectedEdge] b], 
     If[MemberQ[EL, b \[DirectedEdge] c], {b, a \[DirectedEdge] b, 
       b \[DirectedEdge] c}, {b, a \[DirectedEdge] b, 
       c \[DirectedEdge] b}], 
     If[MemberQ[EL, b \[DirectedEdge] a], 
      If[MemberQ[EL, b \[DirectedEdge] c], {b, b \[DirectedEdge] a, 
        b \[DirectedEdge] c}, {b, b \[DirectedEdge] a, 
        c \[DirectedEdge] b}]]], 
    StringJoin["\n There is no path from ", ToString[a], " to ", 
     ToString[b], " to ", ToString[c]]]
   );

   (*Reevaluate the importance of this function.*)
   DifferenceU[edge_] :=
 (
  H[x_, p_, m_, alpha_] := p^2/(2 m^alpha) + V[x] - g[m];
  HJ = (H[x, D[u[x], x], m[x], alpha] == 0);
  curr = ((j == -m[x]*D[H[x, p, m[x], alpha], p]) /. 
     p -> D[u[x], x]);
  ux = (Solve[curr, u'[x]] // Flatten);
  (*R=(HJ/.(curr/.Equal\[Rule]Rule)/.ux);*)
  R = (HJ /. ux);
  sol = RSolve[R, m[x], x] // Flatten;
  M = Select[Select[sol, Im[Values[# /. {x -> .4, j -> 1}]] == 0 &], 
    Re[Values[# /. {x -> .4, j -> 1}]] >= 
      0 &];(*Looks good!*)(*One way of selecting, 
  is there a better one?*)
  
  mfunc[edge] = 
   M /. {j -> jvars[AtHead[edge]] - jvars[AtTail[edge]]};
  fux[edge] = 
   ux /. M /. {j -> jvars[AtHead[edge]] - jvars[AtTail[edge]]};
  (*Integrate[(u'[x]/.fux),{x,0,1}]*)
  
  h[y_] := Values[fux[edge]] /. x -> y;
  ufunc[edge] = {u[y] -> 
     uvars[
       AtTail[edge]] - (jvars[AtHead[edge]] - 
         jvars[AtTail[edge]]) Integrate[
        m[x]^(alpha - 1) /. mfunc[edge], {x, 0, y}]};
  -(jvars[AtHead[edge]] - jvars[AtTail[edge]]) Integrate[
    m[x]^(alpha - 1) /. mfunc[edge], {x, 0, 1}]
  )

EE = Exists[#1, #2] &;
UU = Solve[#1, #2, Reals] &;
RR = Reduce[#1, #2, Reals] &;
