(* Created with the Wolfram Language : www.wolfram.com *)
<|"Timestamp" -> "Tue 5 May 2026 09:54:03", "Commit" -> "f480a1b", 
 "Cases" -> <|"chain-2v" -> <|"Status" -> "OK", "Kind" -> "Rules", 
     "TransitionFlowStatus" -> "Unique", "TransitionFlowResidualCount" -> 0, 
     "Valid" -> True, "Solution" -> <|"Rules" -> {u["auxExit2", 2] -> 0, 
         j[1, 2] -> 10, j[2, 1] -> 0, j[2, "auxExit2"] -> 10, 
         j["auxEntry1", 1] -> 10, j[1, 2, "auxExit2"] -> 10, 
         j["auxEntry1", 1, 2] -> 10, u[1, 2] -> 0, u[2, 1] -> 10, 
         u["auxEntry1", 1] -> 10}, "Equations" -> True|>|>, 
   "chain-3v-1exit" -> <|"Status" -> "OK", "Kind" -> "Rules", 
     "TransitionFlowStatus" -> "Unique", "TransitionFlowResidualCount" -> 0, 
     "Valid" -> True, "Solution" -> <|"Rules" -> {u["auxExit3", 3] -> 0, 
         u[3, 2] -> 10, j[1, 2] -> 10, j[2, 1] -> 0, j[2, 3] -> 10, 
         j[3, 2] -> 0, j[3, "auxExit3"] -> 10, j["auxEntry1", 1] -> 10, 
         j[1, 2, 3] -> 10, j[2, 3, "auxExit3"] -> 10, j[3, 2, 1] -> 0, 
         j["auxEntry1", 1, 2] -> 10, u[1, 2] -> 10, u[2, 1] -> 20, 
         u[2, 3] -> 0, u["auxEntry1", 1] -> 20}, "Equations" -> True|>|>, 
   "chain-3v-2exit" -> <|"Status" -> "OK", "Kind" -> "Parametric", 
     "TransitionFlowStatus" -> "Unique", "TransitionFlowResidualCount" -> 0, 
     "Valid" -> True, "Solution" -> <|"Rules" -> {u["auxExit2", 2] -> 0, 
         u["auxExit3", 3] -> 10, u[3, 2] -> 0, j[1, 2] -> 120, j[2, 1] -> 0, 
         j[2, 3] -> 0, j[2, "auxExit2"] -> 120, j[3, 2] -> 0, 
         j["auxEntry1", 1] -> 120, j[1, 2, 3] -> 0, j[1, 2, "auxExit2"] -> 
          120, j[2, 3, "auxExit3"] -> 0, j[3, 2, 1] -> 0, 
         j[3, 2, "auxExit2"] -> 0, j["auxEntry1", 1, 2] -> 120, u[1, 2] -> 0, 
         u["auxEntry1", 1] -> 120, j[3, "auxExit3"] -> 0, u[2, 1] -> 120}, 
       "Equations" -> u[2, 3] <= 10|>|>, "example-7" -> 
    <|"Status" -> "OK", "Kind" -> "Branched", "TransitionFlowStatus" -> 
      "Underdetermined", "TransitionFlowResidualCount" -> 1, "Valid" -> True, 
     "Solution" -> <|"Rules" -> {u["auxExit3", 3] -> 0, u["auxExit4", 4] -> 
          10, u[3, 2] -> -100 + u[2, 1], u[4, 2] -> -100 + u[2, 1], 
         j[1, 2] -> 100, j[2, 1] -> 0, j[2, 3] -> 100 - j[4, "auxExit4"], 
         j[2, 4] -> j[4, "auxExit4"], j[3, 2] -> 0, j[3, "auxExit3"] -> 
          100 - j[4, "auxExit4"], j[4, 2] -> 0, j["auxEntry1", 1] -> 100, 
         j[1, 2, 3] -> 100 - j[1, 2, 4], j[2, 3, "auxExit3"] -> 
          100 - j[4, "auxExit4"], j[2, 4, "auxExit4"] -> j[4, "auxExit4"], 
         j[3, 2, 1] -> -j[4, "auxExit4"] + j[1, 2, 4], 
         j[3, 2, 4] -> j[4, "auxExit4"] - j[1, 2, 4], j[4, 2, 1] -> 
          j[4, "auxExit4"] - j[1, 2, 4], j[4, 2, 3] -> -j[4, "auxExit4"] + 
           j[1, 2, 4], j["auxEntry1", 1, 2] -> 100, 
         u[1, 2] -> -100 + u[2, 1], u["auxEntry1", 1] -> u[2, 1]}, 
       "Equations" -> (j[4, "auxExit4"] == 0 && j[1, 2, 4] == 0 && 
          u[2, 1] == 200 && u[2, 3] == 0 && u[2, 4] <= 10) || 
         (j[4, "auxExit4"] == 45 && j[1, 2, 4] == 45 && u[2, 1] == 155 && 
          u[2, 3] == 0 && u[2, 4] == 10) || (j[4, "auxExit4"] == 100 && 
          j[1, 2, 4] == 100 && u[2, 1] == 210 && u[2, 3] <= 0 && 
          u[2, 4] == 10)|>|>, "chain-5v-1exit" -> 
    <|"Status" -> "OK", "Kind" -> "Rules", "TransitionFlowStatus" -> 
      "Unique", "TransitionFlowResidualCount" -> 0, "Valid" -> True, 
     "Solution" -> <|"Rules" -> {u["auxExit5", 5] -> 0, u[3, 2] -> 30, 
         u[4, 3] -> 20, u[5, 4] -> 10, j[1, 2] -> 10, j[2, 1] -> 0, 
         j[2, 3] -> 10, j[3, 2] -> 0, j[3, 4] -> 10, j[4, 3] -> 0, 
         j[4, 5] -> 10, j[5, 4] -> 0, j[5, "auxExit5"] -> 10, 
         j["auxEntry1", 1] -> 10, j[1, 2, 3] -> 10, j[2, 3, 4] -> 10, 
         j[3, 2, 1] -> 0, j[3, 4, 5] -> 10, j[4, 3, 2] -> 0, 
         j[4, 5, "auxExit5"] -> 10, j[5, 4, 3] -> 0, j["auxEntry1", 1, 2] -> 
          10, u[1, 2] -> 30, u[2, 1] -> 40, u[2, 3] -> 20, u[3, 4] -> 10, 
         u[4, 5] -> 0, u["auxEntry1", 1] -> 40}, "Equations" -> True|>|>, 
   "grid-2x3" -> <|"Status" -> "TIMEOUT", "Kind" -> "-", 
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
