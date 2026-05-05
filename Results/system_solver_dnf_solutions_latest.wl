(* Created with the Wolfram Language : www.wolfram.com *)
<|"Timestamp" -> "Mon 4 May 2026 16:33:05", "Commit" -> "0cbeec3", 
 "Solver" -> "dnf", "Cases" -> 
  <|"chain-2v" -> <|"Status" -> "OK", "Kind" -> "Rules", 
     "TransitionFlowStatus" -> "Unique", "TransitionFlowResidualCount" -> 0, 
     "Valid" -> True, "Solution" -> {u["auxExit2", 2] -> 0, j[1, 2] -> 10, 
       j[2, 1] -> 0, j[2, "auxExit2"] -> 10, j["auxEntry1", 1] -> 10, 
       j[1, 2, "auxExit2"] -> 10, j["auxEntry1", 1, 2] -> 10, u[1, 2] -> 0, 
       u[2, 1] -> 10, u["auxEntry1", 1] -> 10}|>, 
   "chain-3v-1exit" -> <|"Status" -> "OK", "Kind" -> "Rules", 
     "TransitionFlowStatus" -> "Unique", "TransitionFlowResidualCount" -> 0, 
     "Valid" -> True, "Solution" -> {u["auxExit3", 3] -> 0, u[3, 2] -> 10, 
       j[1, 2] -> 10, j[2, 1] -> 0, j[2, 3] -> 10, j[3, 2] -> 0, 
       j[3, "auxExit3"] -> 10, j["auxEntry1", 1] -> 10, j[1, 2, 3] -> 10, 
       j[2, 3, "auxExit3"] -> 10, j[3, 2, 1] -> 0, j["auxEntry1", 1, 2] -> 
        10, u[1, 2] -> 10, u[2, 1] -> 20, u[2, 3] -> 0, 
       u["auxEntry1", 1] -> 20}|>, "chain-3v-2exit" -> 
    <|"Status" -> "OK", "Kind" -> "Branched", "TransitionFlowStatus" -> 
      "Underdetermined", "TransitionFlowResidualCount" -> 3, "Valid" -> True, 
     "Solution" -> <|"Rules" -> {u["auxExit2", 2] -> 0, u["auxExit3", 3] -> 
          10, u[3, 2] -> u[1, 2], j[1, 2] -> 120, j[2, 1] -> 0, j[3, 2] -> 0, 
         j["auxEntry1", 1] -> 120, j[3, 2, 1] -> 0, j[3, 2, "auxExit2"] -> 0, 
         j["auxEntry1", 1, 2] -> 120}, "Equations" -> 
        (j[2, 3] == 0 && ((j[2, 3] == j[1, 2, 3] && j[2, 3] == 
             j[2, 3, "auxExit3"] && j[2, "auxExit2"] == 
             j[1, 2, "auxExit2"] && j[3, "auxExit3"] == 
             j[2, 3, "auxExit3"] && j[1, 2, "auxExit2"] == 0 && 
            j[1, 2, 3] + j[1, 2, "auxExit2"] == 120 && 120 + u[1, 2] == 
             u[2, 1] && u[2, 1] == u["auxEntry1", 1] && j[2, "auxExit2"] >= 
             0 && j[3, "auxExit3"] >= 0 && j[1, 2, 3] >= 0 && u[1, 2] <= 0 && 
            ((j[2, 3, "auxExit3"] == 0 && u[2, 3] <= 10) || (u[2, 3] == 10 && 
              j[2, 3, "auxExit3"] >= 0))) || (j[2, "auxExit2"] == 120 && 
            j[3, "auxExit3"] == 0 && j[1, 2, 3] == 0 && 
            j[1, 2, "auxExit2"] == 120 && j[2, 3, "auxExit3"] == 0 && 
            u[1, 2] == 0 && u[2, 1] == 120 && u["auxEntry1", 1] == 120 && 
            (u[2, 3] == 0 || u[2, 3] == 10)))) || (j[2, "auxExit2"] >= 0 && 
          j[3, "auxExit3"] >= 0 && j[1, 2, 3] >= 0 && j[2, 3] >= 0 && 
          j[2, 3, "auxExit3"] == 0 && j[1, 2, "auxExit2"] == 0 && 
          u[2, 1] == u["auxEntry1", 1] && u[1, 2] <= 0 && 
          j[1, 2, 3] + j[1, 2, "auxExit2"] == 120 && j[2, 3] == 
           j[2, 3, "auxExit3"] && j[2, 3] == j[1, 2, 3] && 
          j[2, "auxExit2"] == j[1, 2, "auxExit2"] && j[3, "auxExit3"] == 
           j[2, 3, "auxExit3"] && 120 + u[1, 2] == u[2, 1] && 
          u[2, 3] <= 10 && j[2, 3] + u[2, 3] == u[1, 2])|>|>, 
   "example-7" -> <|"Status" -> "OK", "Kind" -> "Branched", 
     "TransitionFlowStatus" -> "Underdetermined", 
     "TransitionFlowResidualCount" -> 8, "Valid" -> True, 
     "Solution" -> <|"Rules" -> {u["auxExit3", 3] -> 0, u["auxExit4", 4] -> 
          10, u[3, 2] -> u[1, 2], u[4, 2] -> u[1, 2], j[1, 2] -> 100, 
         j[2, 1] -> 0, j[3, 2] -> 0, j[4, 2] -> 0, j["auxEntry1", 1] -> 100, 
         j["auxEntry1", 1, 2] -> 100}, "Equations" -> 
        (j[3, 2, 4] == 0 && ((u[2, 4] == 10 && ((u[2, 3] == 0 && 
              ((j[4, 2, 1] == 0 && ((j[3, 2, 1] == 0 && ((j[4, 2, 3] == 0 && 
                    ((j[2, 3] == 0 && j[2, 4] == 100 && j[3, "auxExit3"] == 
                       0 && j[4, "auxExit4"] == 100 && j[1, 2, 3] == 0 && 
                      j[1, 2, 4] == 100 && j[2, 3, "auxExit3"] == 0 && 
                      j[2, 4, "auxExit4"] == 100 && u[1, 2] == 110 && 
                      u[2, 1] == 210 && u["auxEntry1", 1] == 210) || 
                     (j[2, 3] == 100 && j[2, 4] == 0 && j[3, "auxExit3"] == 
                       100 && j[4, "auxExit4"] == 0 && j[1, 2, 3] == 100 && 
                      j[1, 2, 4] == 0 && j[2, 3, "auxExit3"] == 100 && 
                      j[2, 4, "auxExit4"] == 0 && u[1, 2] == 100 && 
                      u[2, 1] == 200 && u["auxEntry1", 1] == 200) || 
                     (j[2, 3] == 55 && j[2, 4] == 45 && j[3, "auxExit3"] == 
                       55 && j[4, "auxExit4"] == 45 && j[1, 2, 3] == 55 && 
                      j[1, 2, 4] == 45 && j[2, 3, "auxExit3"] == 55 && 
                      j[2, 4, "auxExit4"] == 45 && u[1, 2] == 55 && 
                      u[2, 1] == 155 && u["auxEntry1", 1] == 155))) || 
                   (j[2, 3] == 0 && j[2, 4] == 0 && j[3, "auxExit3"] >= 0 && 
                    j[4, "auxExit4"] >= 0 && j[1, 2, 3] >= 0 && j[1, 2, 4] >= 
                     0 && j[1, 2, 3] + j[1, 2, 4] == 100 && j[2, 3] == 
                     j[2, 3, "auxExit3"] && j[3, "auxExit3"] == j[2, 3, 
                      "auxExit3"] && j[2, 3, "auxExit3"] >= 0 && j[2, 4] == 
                     j[2, 4, "auxExit4"] && j[4, "auxExit4"] == j[2, 4, 
                      "auxExit4"] && j[2, 4, "auxExit4"] >= 0 && j[2, 4] == 
                     j[1, 2, 4] + j[3, 2, 4] && j[3, 2, 1] + j[3, 2, 4] == 
                     0 && j[3, 2, 1] + j[4, 2, 1] == 0 && j[2, 3] == 
                     j[1, 2, 3] + j[4, 2, 3] && j[4, 2, 3] >= 0 && 
                    j[4, 2, 1] + j[4, 2, 3] == 0 && 100 + u[1, 2] == 
                     u[2, 1] && u[2, 1] == u["auxEntry1", 1]))) || 
                 (j[3, "auxExit3"] >= 0 && j[4, "auxExit4"] >= 0 && 
                  j[1, 2, 3] == 0 && j[1, 2, 4] >= 0 && j[1, 2, 3] + 
                    j[1, 2, 4] == 100 && j[2, 3] == j[2, 3, "auxExit3"] && 
                  j[3, "auxExit3"] == j[2, 3, "auxExit3"] && 
                  j[2, 3, "auxExit3"] >= 0 && j[2, 4] == j[2, 4, 
                    "auxExit4"] && j[4, "auxExit4"] == j[2, 4, "auxExit4"] && 
                  j[2, 4, "auxExit4"] >= 0 && j[3, 2, 1] >= 0 && j[2, 4] == 
                   j[1, 2, 4] + j[3, 2, 4] && j[3, 2, 1] + j[3, 2, 4] == 0 && 
                  j[3, 2, 1] + j[4, 2, 1] == 0 && j[2, 3] == j[1, 2, 3] + 
                    j[4, 2, 3] && j[4, 2, 3] >= 0 && j[4, 2, 1] + j[4, 2, 
                     3] == 0 && 100 + u[1, 2] == u[2, 1] && u[2, 1] == 
                   u["auxEntry1", 1] && ((j[2, 4] == 0 && (j[2, 3] == 0 || 
                     (j[2, 3] > 0 && j[2, 3] + u[2, 3] == u[1, 2]))) || 
                   (j[2, 3] >= 0 && j[2, 4] >= 0 && j[2, 3] + u[2, 3] == 
                     u[1, 2] && j[2, 4] + u[2, 4] == u[1, 2]))))) || (
                j[3, "auxExit3"] >= 0 && j[4, "auxExit4"] >= 0 && 
                j[1, 2, 4] == 0 && j[1, 2, 3] + j[1, 2, 4] == 100 && 
                j[2, 3] == j[2, 3, "auxExit3"] && j[3, "auxExit3"] == 
                 j[2, 3, "auxExit3"] && j[2, 3, "auxExit3"] >= 0 && 
                j[2, 4] == j[2, 4, "auxExit4"] && j[4, "auxExit4"] == 
                 j[2, 4, "auxExit4"] && j[2, 4, "auxExit4"] >= 0 && 
                j[2, 4] == j[1, 2, 4] + j[3, 2, 4] && j[3, 2, 1] + 
                  j[3, 2, 4] == 0 && j[4, 2, 1] >= 0 && j[3, 2, 1] + 
                  j[4, 2, 1] == 0 && j[2, 3] == j[1, 2, 3] + j[4, 2, 3] && 
                j[4, 2, 3] >= 0 && j[4, 2, 1] + j[4, 2, 3] == 0 && 
                100 + u[1, 2] == u[2, 1] && u[2, 1] == u["auxEntry1", 1] && 
                ((j[2, 3] == 0 && ((j[1, 2, 3] == 0 && j[3, 2, 1] >= 0) || 
                   (j[1, 2, 3] >= 0 && j[3, 2, 1] == 0)) && (j[2, 4] == 0 || 
                   (j[2, 4] > 0 && j[2, 4] + u[2, 4] == u[1, 2]))) || 
                 (j[2, 3] >= 0 && j[2, 3] + u[2, 3] == u[1, 2] && 
                  ((j[1, 2, 3] == 0 && j[3, 2, 1] >= 0 && (j[2, 4] == 0 || 
                     (j[2, 4] > 0 && j[2, 4] + u[2, 4] == u[1, 2]))) || 
                   (j[2, 4] >= 0 && j[1, 2, 3] >= 0 && j[3, 2, 1] == 0 && 
                    j[2, 4] + u[2, 4] == u[1, 2]))))))) || 
             (j[3, "auxExit3"] >= 0 && j[4, "auxExit4"] >= 0 && 
              j[1, 2, 3] + j[1, 2, 4] == 100 && j[2, 3] == j[2, 3, 
                "auxExit3"] && j[3, "auxExit3"] == j[2, 3, "auxExit3"] && 
              j[2, 3, "auxExit3"] == 0 && j[2, 4] == j[2, 4, "auxExit4"] && 
              j[4, "auxExit4"] == j[2, 4, "auxExit4"] && j[2, 4, 
                "auxExit4"] >= 0 && j[2, 4] == j[1, 2, 4] + j[3, 2, 4] && 
              j[3, 2, 1] + j[3, 2, 4] == 0 && j[3, 2, 1] + j[4, 2, 1] == 0 && 
              j[2, 3] == j[1, 2, 3] + j[4, 2, 3] && j[4, 2, 3] >= 0 && 
              j[4, 2, 1] + j[4, 2, 3] == 0 && 100 + u[1, 2] == u[2, 1] && 
              u[2, 3] <= 0 && u[2, 1] == u["auxEntry1", 1] && 
              ((j[2, 4] == 0 && ((j[1, 2, 4] == 0 && j[4, 2, 1] >= 0) || 
                 (j[1, 2, 4] >= 0 && j[4, 2, 1] == 0))) || (j[2, 4] >= 0 && 
                j[1, 2, 4] == 0 && j[4, 2, 1] >= 0 && j[2, 4] + u[2, 4] == 
                 u[1, 2])) && ((j[1, 2, 3] == 0 && j[3, 2, 1] >= 0) || (
                j[1, 2, 3] >= 0 && j[3, 2, 1] == 0)) && (j[2, 3] == 0 || (
                j[2, 3] > 0 && j[2, 3] + u[2, 3] == u[1, 2]))))) || 
           (j[3, "auxExit3"] >= 0 && j[4, "auxExit4"] >= 0 && 
            j[1, 2, 3] + j[1, 2, 4] == 100 && j[2, 3] == 
             j[2, 3, "auxExit3"] && j[3, "auxExit3"] == 
             j[2, 3, "auxExit3"] && j[2, 4] == j[2, 4, "auxExit4"] && 
            j[4, "auxExit4"] == j[2, 4, "auxExit4"] && j[2, 4, "auxExit4"] == 
             0 && j[2, 4] == j[1, 2, 4] + j[3, 2, 4] && 
            j[3, 2, 1] + j[3, 2, 4] == 0 && j[3, 2, 1] + j[4, 2, 1] == 0 && 
            j[2, 3] == j[1, 2, 3] + j[4, 2, 3] && j[4, 2, 3] >= 0 && 
            j[4, 2, 1] + j[4, 2, 3] == 0 && 100 + u[1, 2] == u[2, 1] && 
            u[2, 4] <= 10 && u[2, 1] == u["auxEntry1", 1] && 
            ((j[1, 2, 4] == 0 && j[4, 2, 1] >= 0) || (j[1, 2, 4] >= 0 && 
              j[4, 2, 1] == 0)) && (j[2, 4] == 0 || (j[2, 4] > 0 && 
              j[2, 4] + u[2, 4] == u[1, 2])) && ((j[2, 3] == 0 && 
              ((j[1, 2, 3] == 0 && j[3, 2, 1] >= 0) || (j[1, 2, 3] >= 0 && 
                j[3, 2, 1] == 0)) && ((j[2, 3, "auxExit3"] >= 0 && 
                u[2, 3] == 0) || (j[2, 3, "auxExit3"] == 0 && u[2, 3] <= 
                 0))) || (j[2, 3] + u[2, 3] == u[1, 2] && j[2, 3] >= 0 && 
              ((j[1, 2, 3] == 0 && j[3, 2, 1] >= 0 && 
                ((j[2, 3, "auxExit3"] >= 0 && u[2, 3] == 0) || 
                 (j[2, 3, "auxExit3"] == 0 && u[2, 3] <= 0))) || (
                j[1, 2, 3] >= 0 && j[2, 3, "auxExit3"] == 0 && j[3, 2, 1] == 
                 0 && u[2, 3] <= 0))))))) || (j[3, "auxExit3"] >= 0 && 
          j[4, "auxExit4"] >= 0 && j[1, 2, 3] + j[1, 2, 4] == 100 && 
          j[2, 3] == j[2, 3, "auxExit3"] && j[3, "auxExit3"] == 
           j[2, 3, "auxExit3"] && j[2, 4] == j[2, 4, "auxExit4"] && 
          j[4, "auxExit4"] == j[2, 4, "auxExit4"] && j[2, 4] == 
           j[1, 2, 4] + j[3, 2, 4] && j[3, 2, 4] >= 0 && 
          j[3, 2, 1] + j[3, 2, 4] == 0 && j[3, 2, 1] + j[4, 2, 1] == 0 && 
          j[2, 3] == j[1, 2, 3] + j[4, 2, 3] && j[4, 2, 3] == 0 && 
          j[4, 2, 1] + j[4, 2, 3] == 0 && 100 + u[1, 2] == u[2, 1] && 
          u[2, 1] == u["auxEntry1", 1] && ((j[2, 3] == 0 && 
            ((j[1, 2, 3] == 0 && j[3, 2, 1] >= 0) || (j[1, 2, 3] >= 0 && 
              j[3, 2, 1] == 0)) && ((j[2, 3, "auxExit3"] >= 0 && 
              u[2, 3] == 0) || (j[2, 3, "auxExit3"] == 0 && u[2, 3] <= 0)) && 
            ((j[2, 4] == 0 && ((j[1, 2, 4] == 0 && j[4, 2, 1] >= 0) || (
                j[1, 2, 4] >= 0 && j[4, 2, 1] == 0)) && 
              ((j[2, 4, "auxExit4"] >= 0 && u[2, 4] == 10) || (
                j[2, 4, "auxExit4"] == 0 && u[2, 4] <= 10))) || 
             (j[2, 4] >= 0 && j[2, 4] + u[2, 4] == u[1, 2] && 
              ((j[1, 2, 4] == 0 && j[4, 2, 1] >= 0 && 
                ((j[2, 4, "auxExit4"] >= 0 && u[2, 4] == 10) || 
                 (j[2, 4, "auxExit4"] == 0 && u[2, 4] <= 10))) || (
                j[1, 2, 4] >= 0 && j[2, 4, "auxExit4"] == 0 && j[4, 2, 1] == 
                 0 && u[2, 4] <= 10))))) || (j[2, 3] >= 0 && 
            j[2, 3] + u[2, 3] == u[1, 2] && ((j[2, 4] == 0 && 
              ((j[1, 2, 4] == 0 && j[4, 2, 1] >= 0) || (j[1, 2, 4] >= 0 && 
                j[4, 2, 1] == 0)) && ((j[2, 4, "auxExit4"] >= 0 && 
                u[2, 4] == 10) || (j[2, 4, "auxExit4"] == 0 && u[2, 4] <= 
                 10)) && ((j[1, 2, 3] == 0 && j[3, 2, 1] >= 0 && 
                ((j[2, 3, "auxExit3"] >= 0 && u[2, 3] == 0) || 
                 (j[2, 3, "auxExit3"] == 0 && u[2, 3] <= 0))) || (
                j[1, 2, 3] >= 0 && j[2, 3, "auxExit3"] == 0 && j[3, 2, 1] == 
                 0 && u[2, 3] <= 0))) || (j[2, 4] >= 0 && j[2, 4] + 
                u[2, 4] == u[1, 2] && ((j[1, 2, 3] == 0 && j[3, 2, 1] >= 0 && 
                ((j[1, 2, 4] == 0 && j[4, 2, 1] >= 0 && 
                  ((j[2, 3, "auxExit3"] >= 0 && u[2, 3] == 0) || 
                   (j[2, 3, "auxExit3"] == 0 && u[2, 3] <= 0)) && 
                  ((j[2, 4, "auxExit4"] >= 0 && u[2, 4] == 10) || 
                   (j[2, 4, "auxExit4"] == 0 && u[2, 4] <= 10))) || 
                 (j[4, 2, 1] == 0 && j[1, 2, 4] >= 0 && 
                  ((j[2, 4, "auxExit4"] == 0 && u[2, 4] <= 10 && 
                    ((j[2, 3, "auxExit3"] == 0 && u[2, 3] <= 0) || 
                     (u[2, 3] == 0 && j[2, 3, "auxExit3"] >= 0))) || 
                   (u[2, 3] == 0 && u[2, 4] == 10 && j[2, 3, "auxExit3"] >= 
                     0 && j[2, 4, "auxExit4"] >= 0))))) || (j[1, 2, 3] >= 
                 0 && j[3, 2, 1] == 0 && ((j[1, 2, 4] == 0 && j[4, 2, 1] >= 
                   0 && ((j[2, 4, "auxExit4"] >= 0 && u[2, 4] == 10 && 
                    ((j[2, 3, "auxExit3"] >= 0 && u[2, 3] == 0) || 
                     (j[2, 3, "auxExit3"] == 0 && u[2, 3] <= 0))) || 
                   (j[2, 3, "auxExit3"] == 0 && j[2, 4, "auxExit4"] == 0 && 
                    u[2, 3] <= 0 && u[2, 4] <= 10))) || (j[1, 2, 4] >= 0 && 
                  j[2, 3, "auxExit3"] == 0 && j[2, 4, "auxExit4"] == 0 && 
                  j[4, 2, 1] == 0 && u[2, 3] <= 0 && u[2, 4] <= 
                   10)))))))))|>|>, "chain-5v-1exit" -> 
    <|"Status" -> "OK", "Kind" -> "Rules", "TransitionFlowStatus" -> 
      "Unique", "TransitionFlowResidualCount" -> 0, "Valid" -> True, 
     "Solution" -> {u["auxExit5", 5] -> 0, u[3, 2] -> 30, u[4, 3] -> 20, 
       u[5, 4] -> 10, j[1, 2] -> 10, j[2, 1] -> 0, j[2, 3] -> 10, 
       j[3, 2] -> 0, j[3, 4] -> 10, j[4, 3] -> 0, j[4, 5] -> 10, 
       j[5, 4] -> 0, j[5, "auxExit5"] -> 10, j["auxEntry1", 1] -> 10, 
       j[1, 2, 3] -> 10, j[2, 3, 4] -> 10, j[3, 2, 1] -> 0, j[3, 4, 5] -> 10, 
       j[4, 3, 2] -> 0, j[4, 5, "auxExit5"] -> 10, j[5, 4, 3] -> 0, 
       j["auxEntry1", 1, 2] -> 10, u[1, 2] -> 30, u[2, 1] -> 40, 
       u[2, 3] -> 20, u[3, 4] -> 10, u[4, 5] -> 0, u["auxEntry1", 1] -> 
        40}|>, "grid-2x3" -> <|"Status" -> "TIMEOUT", "Kind" -> "-", 
     "TransitionFlowStatus" -> Missing["KeyAbsent", "TransitionFlowStatus"], 
     "TransitionFlowResidualCount" -> Missing["KeyAbsent", 
       "TransitionFlowResidualCount"], "Valid" -> Missing["TimedOut"], 
     "Solution" -> $TimedOut|>, "grid-3x2" -> <|"Status" -> "TIMEOUT", 
     "Kind" -> "-", "TransitionFlowStatus" -> Missing["KeyAbsent", 
       "TransitionFlowStatus"], "TransitionFlowResidualCount" -> 
      Missing["KeyAbsent", "TransitionFlowResidualCount"], 
     "Valid" -> Missing["TimedOut"], "Solution" -> $TimedOut|>, 
   "example-12" -> <|"Status" -> "TIMEOUT", "Kind" -> "-", 
     "TransitionFlowStatus" -> Missing["KeyAbsent", "TransitionFlowStatus"], 
     "TransitionFlowResidualCount" -> Missing["KeyAbsent", 
       "TransitionFlowResidualCount"], "Valid" -> Missing["TimedOut"], 
     "Solution" -> $TimedOut|>|>|>
