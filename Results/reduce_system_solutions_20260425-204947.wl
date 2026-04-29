(* Created with the Wolfram Language : www.wolfram.com *)
<|"Timestamp" -> "Sat 25 Apr 2026 20:51:22", "Commit" -> "278f64f", 
 "Cases" -> <|"chain-2v" -> <|"Status" -> "OK", "Kind" -> "Rules", 
     "Valid" -> True, "Solution" -> {j[1, 2] -> 10, j[2, 1] -> 0, 
       j[2, "auxExit2"] -> 10, j["auxEntry1", 1] -> 10, 
       j[1, 2, "auxExit2"] -> 10, j["auxEntry1", 1, 2] -> 10, u[1, 2] -> 0, 
       u[2, 1] -> 10, u["auxEntry1", 1] -> 10, u["auxExit2", 2] -> 0}|>, 
   "chain-3v-1exit" -> <|"Status" -> "OK", "Kind" -> "Rules", 
     "Valid" -> True, "Solution" -> {j[1, 2] -> 10, j[2, 1] -> 0, 
       j[2, 3] -> 10, j[3, 2] -> 0, j[3, "auxExit3"] -> 10, 
       j["auxEntry1", 1] -> 10, j[1, 2, 3] -> 10, j[2, 3, "auxExit3"] -> 10, 
       j[3, 2, 1] -> 0, j["auxEntry1", 1, 2] -> 10, u[1, 2] -> 10, 
       u[2, 1] -> 20, u[2, 3] -> 0, u[3, 2] -> 10, u["auxEntry1", 1] -> 20, 
       u["auxExit3", 3] -> 0}|>, "chain-3v-2exit" -> 
    <|"Status" -> "OK", "Kind" -> "Underdetermined", "Valid" -> True, 
     "Solution" -> <|"Rules" -> {j[1, 2] -> 120, j[2, 1] -> 0}, 
       "Equations" -> (j[2, 3] == 0 && j[2, "auxExit2"] == 120 && 
          j[3, 2] == 0 && j[3, "auxExit3"] == 0 && j["auxEntry1", 1] == 
           120 && j[1, 2, 3] == 0 && j[1, 2, "auxExit2"] == 120 && 
          j[2, 3, "auxExit3"] == 0 && j[3, 2, 1] == 0 && 
          j[3, 2, "auxExit2"] == 0 && j["auxEntry1", 1, 2] == 120 && 
          u[1, 2] == 0 && u[2, 1] == 120 && ((u[2, 3] < 0 && 
            u[3, 2] == u[2, 3] && u["auxEntry1", 1] == 120 && 
            u["auxExit2", 2] == 0) || (u[2, 3] == 0 && u[3, 2] == 0 && 
            u["auxEntry1", 1] == 120 && u["auxExit2", 2] == 0) || 
           (u[2, 3] > 0 && u[3, 2] == u[2, 3] && u["auxEntry1", 1] == 120 && 
            u["auxExit2", 2] == 0))) || (Inequality[0, Less, j[2, 3], Less, 
           120] && j[2, "auxExit2"] == 120 - j[2, 3] && j[3, 2] == 0 && 
          j[3, "auxExit3"] == j[2, 3] && j["auxEntry1", 1] == 120 && 
          j[1, 2, 3] == 120 - j[2, "auxExit2"] && j[1, 2, "auxExit2"] == 
           120 - j[1, 2, 3] && j[2, 3, "auxExit3"] == j[2, 3] && 
          j[3, 2, 1] == 0 && j[3, 2, "auxExit2"] == 0 && 
          j["auxEntry1", 1, 2] == 120 && u[1, 2] == 0 && u[2, 1] == 120 && 
          u[2, 3] == -j[2, 3] && u[3, 2] == 0 && u["auxEntry1", 1] == 120 && 
          u["auxExit2", 2] == 0 && u["auxExit3", 3] == u[2, 3]) || 
         (j[2, 3] == 120 && j[2, "auxExit2"] == 0 && j[3, 2] == 0 && 
          j[3, "auxExit3"] == 120 && j["auxEntry1", 1] == 120 && 
          j[1, 2, 3] == 120 && j[1, 2, "auxExit2"] == 0 && 
          j[2, 3, "auxExit3"] == 120 && j[3, 2, 1] == 0 && 
          j[3, 2, "auxExit2"] == 0 && j["auxEntry1", 1, 2] == 120 && 
          u[1, 2] == 0 && u[2, 1] == 120 && u[2, 3] == -120 && 
          u[3, 2] == 0 && u["auxEntry1", 1] == 120 && u["auxExit3", 3] == 
           -120)|>|>, "example-7" -> <|"Status" -> "OK", 
     "Kind" -> "Underdetermined", "Valid" -> True, 
     "Solution" -> <|"Rules" -> {j[1, 2] -> 100, j[2, 1] -> 0}, 
       "Equations" -> (j[2, 3] == 0 && j[2, 4] == 100 && j[3, 2] == 0 && 
          j[3, "auxExit3"] == 0 && j[4, 2] == 0 && j[4, "auxExit4"] == 100 && 
          j["auxEntry1", 1] == 100 && j[1, 2, 3] == 0 && j[1, 2, 4] == 100 && 
          j[2, 3, "auxExit3"] == 0 && j[2, 4, "auxExit4"] == 100 && 
          j[3, 2, 1] == 0 && j[3, 2, 4] == 0 && j[4, 2, 1] == 0 && 
          j[4, 2, 3] == 0 && j["auxEntry1", 1, 2] == 100 && u[1, 2] == 110 && 
          u[2, 1] == 210 && u[2, 3] == 0 && u[2, 4] == 10 && u[3, 2] == 0 && 
          u[4, 2] == 110 && u["auxEntry1", 1] == 210 && u["auxExit4", 4] == 
           10) || (j[2, 3] == 55 && j[2, 4] == 45 && j[3, 2] == 0 && 
          j[3, "auxExit3"] == 55 && j[4, 2] == 0 && j[4, "auxExit4"] == 45 && 
          j["auxEntry1", 1] == 100 && j[1, 2, 3] == 55 && j[1, 2, 4] == 45 && 
          j[2, 3, "auxExit3"] == 55 && j[2, 4, "auxExit4"] == 45 && 
          j[3, 2, 1] == 0 && j[3, 2, 4] == 0 && j[4, 2, 1] == 0 && 
          j[4, 2, 3] == 0 && j["auxEntry1", 1, 2] == 100 && u[1, 2] == 55 && 
          u[2, 1] == 155 && u[2, 3] == 0 && u[2, 4] == 10 && u[3, 2] == 55 && 
          u[4, 2] == 55 && u["auxEntry1", 1] == 155 && u["auxExit3", 3] == 
           0 && u["auxExit4", 4] == 10) || (j[2, 3] == 100 && j[2, 4] == 0 && 
          j[3, 2] == 0 && j[3, "auxExit3"] == 100 && j[4, 2] == 0 && 
          j[4, "auxExit4"] == 0 && j["auxEntry1", 1] == 100 && 
          j[1, 2, 3] == 100 && j[1, 2, 4] == 0 && j[2, 3, "auxExit3"] == 
           100 && j[2, 4, "auxExit4"] == 0 && j[3, 2, 1] == 0 && 
          j[3, 2, 4] == 0 && j[4, 2, 1] == 0 && j[4, 2, 3] == 0 && 
          j["auxEntry1", 1, 2] == 100 && u[1, 2] == 100 && u[2, 1] == 200 && 
          u[2, 3] == 0 && u[2, 4] == 10 && u[3, 2] == 100 && u[4, 2] == 10 && 
          u["auxEntry1", 1] == 200 && u["auxExit3", 3] == 0)|>|>, 
   "chain-5v-1exit" -> <|"Status" -> "OK", "Kind" -> "Rules", 
     "Valid" -> True, "Solution" -> {j[1, 2] -> 10, j[2, 1] -> 0, 
       j[2, 3] -> 10, j[3, 2] -> 0, j[3, 4] -> 10, j[4, 3] -> 0, 
       j[4, 5] -> 10, j[5, 4] -> 0, j[5, "auxExit5"] -> 10, 
       j["auxEntry1", 1] -> 10, j[1, 2, 3] -> 10, j[2, 3, 4] -> 10, 
       j[3, 2, 1] -> 0, j[3, 4, 5] -> 10, j[4, 3, 2] -> 0, 
       j[4, 5, "auxExit5"] -> 10, j[5, 4, 3] -> 0, j["auxEntry1", 1, 2] -> 
        10, u[1, 2] -> 30, u[2, 1] -> 40, u[2, 3] -> 20, u[3, 2] -> 30, 
       u[3, 4] -> 10, u[4, 3] -> 20, u[4, 5] -> 0, u[5, 4] -> 10, 
       u["auxEntry1", 1] -> 40, u["auxExit5", 5] -> 0}|>, 
   "example-12" -> <|"Status" -> "TIMEOUT", "Kind" -> "-",
     "Valid" -> Missing["TimedOut"], "Solution" -> $TimedOut|>|>|>
