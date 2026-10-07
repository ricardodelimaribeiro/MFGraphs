(* Exact content comparison for artifacts, independent of mathematical validity.
   Wolfram 15.0.1 may compare Associations containing Graph objects unequally
   after WXF even when Normal[association] is identical. We compare EVERY key,
   value, context, and metadata field using an explicit tree representation.
   Graphs are represented by ordered vertices/edges and all explicit options;
   SparseArray is represented by dimensions and complete ArrayRules (including
   the default). Neither numeric tolerance nor a content hash defines equality.
   SHA256 in run manifests is only a subsequent file-corruption check. *)
ClearAll[artifactContent, artifactContentHeld, artifactHeldChildren, artifactContentReport];
artifactContent[expr_] := artifactContentHeld[HoldComplete[expr]];
artifactContentHeld[HoldComplete[a_Association]] := artifactAssociation[artifactContent[Normal[a]]];
artifactContentHeld[HoldComplete[g_Graph]] := artifactGraph[
    artifactContent[VertexList[g]], artifactContent[EdgeList[g]], artifactContent[Options[g]]];
artifactContentHeld[HoldComplete[a_SparseArray]] := artifactSparseArray[
    artifactContent[Dimensions[a]], artifactContent[ArrayRules[a]]];
(* Never strip a holding head with Apply[List, ...]: doing so evaluates held
   arithmetic, local Function parameters, and even side effects. Pattern-based
   decomposition hands each child directly to HoldComplete instead. Atoms also
   stay held, so a symbol's OwnValue cannot alter a function's stored body. *)
SetAttributes[artifactHeldChildren, HoldAllComplete];
artifactHeldChildren[] := {};
artifactHeldChildren[first_, rest___] :=
    Prepend[artifactHeldChildren[rest], artifactContentHeld[HoldComplete[first]]];
artifactContentHeld[HoldComplete[head_[args___]]] := artifactExpression[
    artifactContentHeld[HoldComplete[head]], artifactHeldChildren[args]];
artifactContentHeld[HoldComplete[atom_]] := HoldComplete[atom];

artifactContentReport[original_Association, reloaded_Association] := Module[{keys, checks},
    keys = Union[Keys[original], Keys[reloaded]];
    checks = AssociationMap[Function[key,
        KeyExistsQ[original, key] && KeyExistsQ[reloaded, key] &&
        SameQ[artifactContent[KeyTake[original, {key}]],
            artifactContent[KeyTake[reloaded, {key}]]]], keys];
    <|"Integrity" -> If[Keys[original] === Keys[reloaded] && AllTrue[Values[checks], TrueQ],
        "ExactContentMatch", "ContentMismatch"],
      "RawSameQ" -> SameQ[original, reloaded], "FieldChecks" -> checks,
      "Comparison" -> "Every association key/value; ordered graph vertices, edges and explicit options; sparse dimensions, entries and default; held expression trees and symbol contexts without evaluating held contents"|>
];
artifactContentReport[_, _] := <|"Integrity" -> "InvalidArtifactFormat"|>;
