(* Created with the Wolfram Language : www.wolfram.com *)
<|"Timestamp" -> "Sat 25 Apr 2026 17:47:15", "Commit" -> "278f64f", 
 "Cases" -> <|"chain-2v" -> <|"Status" -> "OK", "Kind" -> "Underdetermined", 
     "Solution" -> <|"Rules" -> {j[1, 2] -> 10, j[2, 1] -> 0, 
         j[2, "auxExit2"] -> 10, j["auxEntry1", 1] -> 10, 
         j[1, 2, "auxExit2"] -> 10, j["auxEntry1", 1, 2] -> 10}, 
       "Equations" -> u[2, 1] == 10 + u[1, 2] && u["auxEntry1", 1] == 
          u[2, 1] && u["auxExit2", 2] == u[1, 2]|>|>, 
   "chain-3v-1exit" -> <|"Status" -> "OK", "Kind" -> "Underdetermined", 
     "Solution" -> <|"Rules" -> {j[1, 2] -> 10, j[2, 1] -> 0, j[2, 3] -> 10, 
         j[3, 2] -> 0, j[3, "auxExit3"] -> 10, j["auxEntry1", 1] -> 10, 
         j[1, 2, 3] -> 10, j[2, 3, "auxExit3"] -> 10, j[3, 2, 1] -> 0, 
         j["auxEntry1", 1, 2] -> 10}, "Equations" -> 
        u[2, 1] == 10 + u[1, 2] && u[2, 3] == -10 + u[1, 2] && 
         u[3, 2] == u[1, 2] && u["auxEntry1", 1] == u[2, 1] && 
         u["auxExit3", 3] == u[2, 3]|>|>, "chain-3v-2exit" -> 
    <|"Status" -> "OK", "Kind" -> "Underdetermined", 
     "Solution" -> <|"Rules" -> {j[1, 2] -> 120, j[2, 1] -> 0}, 
       "Equations" -> (j[2, 3] == 0 && j[2, "auxExit2"] == 120 && 
          j[3, 2] == 0 && j[3, "auxExit3"] == 0 && j["auxEntry1", 1] == 
           120 && j[1, 2, 3] == 0 && j[1, 2, "auxExit2"] == 120 && 
          j[2, 3, "auxExit3"] == 0 && j[3, 2, 1] == 0 && 
          j[3, 2, "auxExit2"] == 0 && j["auxEntry1", 1, 2] == 120 && 
          u[2, 1] == 120 + u[1, 2] && u[3, 2] == u[2, 3] && 
          u["auxEntry1", 1] == u[2, 1] && u["auxExit2", 2] == u[1, 2]) || 
         (Inequality[0, Less, j[2, 3], Less, 120] && j[2, "auxExit2"] == 
           120 - j[2, 3] && j[3, 2] == 0 && j[3, "auxExit3"] == j[2, 3] && 
          j["auxEntry1", 1] == 120 && j[1, 2, 3] == 120 - j[2, "auxExit2"] && 
          j[1, 2, "auxExit2"] == 120 - j[1, 2, 3] && j[2, 3, "auxExit3"] == 
           j[2, 3] && j[3, 2, 1] == 0 && j[3, 2, "auxExit2"] == 0 && 
          j["auxEntry1", 1, 2] == 120 && u[2, 1] == 120 + u[1, 2] && 
          u[2, 3] == -j[2, 3] + u[1, 2] && u[3, 2] == j[2, 3] + u[2, 3] && 
          u["auxEntry1", 1] == u[2, 1] && u["auxExit2", 2] == u[1, 2] && 
          u["auxExit3", 3] == u[2, 3]) || (j[2, 3] == 120 && 
          j[2, "auxExit2"] == 0 && j[3, 2] == 0 && j[3, "auxExit3"] == 120 && 
          j["auxEntry1", 1] == 120 && j[1, 2, 3] == 120 && 
          j[1, 2, "auxExit2"] == 0 && j[2, 3, "auxExit3"] == 120 && 
          j[3, 2, 1] == 0 && j[3, 2, "auxExit2"] == 0 && 
          j["auxEntry1", 1, 2] == 120 && u[2, 1] == 120 + u[1, 2] && 
          u[2, 3] == -120 + u[1, 2] && u[3, 2] == 120 + u[2, 3] && 
          u["auxEntry1", 1] == u[2, 1] && u["auxExit3", 3] == u[2, 3])|>|>, 
   "example-7" -> <|"Status" -> "OK", "Kind" -> "Underdetermined", 
     "Solution" -> <|"Rules" -> {j[1, 2] -> 100, j[2, 1] -> 0}, 
       "Equations" -> (j[2, 3] == 0 && j[2, 4] == 100 && j[3, 2] == 0 && 
          j[3, "auxExit3"] == 0 && j[4, 2] == 0 && j[4, "auxExit4"] == 100 && 
          j["auxEntry1", 1] == 100 && j[1, 2, 3] == 0 && j[1, 2, 4] == 100 && 
          j[2, 3, "auxExit3"] == 0 && j[2, 4, "auxExit4"] == 100 && 
          j[3, 2, 1] == 0 && j[3, 2, 4] == 0 && j[4, 2, 1] == 0 && 
          j[4, 2, 3] == 0 && j["auxEntry1", 1, 2] == 100 && 
          u[2, 1] == 100 + u[1, 2] && u[2, 4] == -100 + u[1, 2] && 
          u[3, 2] == u[2, 3] && u[4, 2] == 100 + u[2, 4] && 
          u["auxEntry1", 1] == u[2, 1] && u["auxExit4", 4] == u[2, 4]) || 
         (Inequality[0, Less, j[2, 3], Less, 100] && j[2, 4] == 
           100 - j[2, 3] && j[3, 2] == 0 && j[3, "auxExit3"] == j[2, 3] && 
          j[4, 2] == 0 && j[4, "auxExit4"] == j[2, 4] && j["auxEntry1", 1] == 
           100 && j[1, 2, 3] == 100 - j[2, 4] && j[1, 2, 4] == 
           100 - j[1, 2, 3] && j[2, 3, "auxExit3"] == j[2, 3] && 
          j[2, 4, "auxExit4"] == j[2, 4] && j[3, 2, 1] == 0 && 
          j[3, 2, 4] == 0 && j[4, 2, 1] == 0 && j[4, 2, 3] == 0 && 
          j["auxEntry1", 1, 2] == 100 && u[2, 1] == 100 + u[1, 2] && 
          u[2, 3] == -j[2, 3] + u[1, 2] && u[2, 4] == -j[2, 4] + u[1, 2] && 
          u[3, 2] == j[2, 3] + u[2, 3] && u[4, 2] == j[2, 4] + u[2, 4] && 
          u["auxEntry1", 1] == u[2, 1] && u["auxExit3", 3] == u[2, 3] && 
          u["auxExit4", 4] == u[2, 4]) || (j[2, 3] == 100 && j[2, 4] == 0 && 
          j[3, 2] == 0 && j[3, "auxExit3"] == 100 && j[4, 2] == 0 && 
          j[4, "auxExit4"] == 0 && j["auxEntry1", 1] == 100 && 
          j[1, 2, 3] == 100 && j[1, 2, 4] == 0 && j[2, 3, "auxExit3"] == 
           100 && j[2, 4, "auxExit4"] == 0 && j[3, 2, 1] == 0 && 
          j[3, 2, 4] == 0 && j[4, 2, 1] == 0 && j[4, 2, 3] == 0 && 
          j["auxEntry1", 1, 2] == 100 && u[2, 1] == 100 + u[1, 2] && 
          u[2, 3] == -100 + u[1, 2] && u[3, 2] == 100 + u[2, 3] && 
          u[4, 2] == u[2, 4] && u["auxEntry1", 1] == u[2, 1] && 
          u["auxExit3", 3] == u[2, 3])|>|>, "chain-5v-1exit" -> 
    <|"Status" -> "OK", "Kind" -> "Underdetermined", 
     "Solution" -> <|"Rules" -> {j[1, 2] -> 10, j[2, 1] -> 0, j[2, 3] -> 10, 
         j[3, 2] -> 0, j[3, 4] -> 10, j[4, 3] -> 0, j[4, 5] -> 10, 
         j[5, 4] -> 0, j[5, "auxExit5"] -> 10, j["auxEntry1", 1] -> 10, 
         j[1, 2, 3] -> 10, j[2, 3, 4] -> 10, j[3, 2, 1] -> 0, 
         j[3, 4, 5] -> 10, j[4, 3, 2] -> 0, j[4, 5, "auxExit5"] -> 10, 
         j[5, 4, 3] -> 0, j["auxEntry1", 1, 2] -> 10}, 
       "Equations" -> u[2, 1] == 10 + u[1, 2] && u[2, 3] == -10 + u[1, 2] && 
         u[3, 2] == 10 + u[2, 3] && u[3, 4] == -10 + u[2, 3] && 
         u[4, 3] == 10 + u[3, 4] && u[4, 5] == -10 + u[3, 4] && 
         u[5, 4] == 10 + u[4, 5] && u["auxEntry1", 1] == u[2, 1] && 
         u["auxExit5", 5] == u[4, 5]|>|>, "example-12" -> 
    <|"Status" -> "TIMEOUT", "Kind" -> "-",
     "Solution" -> $TimedOut|>|>|>
