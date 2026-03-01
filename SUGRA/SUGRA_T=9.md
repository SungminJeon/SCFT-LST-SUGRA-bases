| ii | Single LST + Single External | Two Externals | 2-3 External |
|---:|:---|:---|---|
| 1 | $`{\color{red}C_{1}}, 2, 1, \{{\color{red}C_{4}}, 2, \overset{{\color{red}C_{5}}}{1}, \{{\color{red}C_{6}}, 2, \overset{{\color{red}C_{7}}}{1}, 8\}\}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 2 | $`1, 3, 1, \{{\color{red}C_{5}}, 3, \overset{{\color{red}C_{6}}}{1}, \overset{{\color{red}C_{7}}}{\overset{1}{6}}\}, 1, 3, {\color{red}C_{9}}`$ |  | |
| 3 | $`{\color{red}C_{1}}, 2, 2, 1, \{{\color{red}C_{5}}, 2, \overset{{\color{red}C_{6}}}{1}, \overset{{\color{red}C_{7}}}{\overset{1}{8}}\}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 4 | $`2, \overset{{\color{red}C_{3}}}{\overset{1}{3}}, 2, 1, \underset{{\color{red}C_{7}}}{\underset{1}{\overset{{\color{red}C_{8}}}{\overset{1}{7}}}}, 1, {\color{red}C_{9}}`$ |  | |
| 5 | $`{\color{red}C_{1}}, 1, 2, 3, 2, 1, \underset{{\color{red}C_{7}}}{\underset{1}{\overset{{\color{red}C_{8}}}{\overset{1}{6}}}}, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{7}} -> 2}`$<br>$`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{8}} -> 2}`$<br>$`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 6 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, \underset{{\color{red}C_{5}}}{\underset{1}{\overset{{\color{red}C_{6}}}{\overset{1}{6}}}}, 1, \overset{{\color{red}C_{8}}}{3}, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 7 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, \{{\color{red}C_{5}}, 2, \overset{{\color{red}C_{6}}}{1}, 6\}, 1, \overset{{\color{red}C_{8}}}{3}, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 8 | $`1, 3, 1, \{{\color{red}C_{5}}, 3, \overset{{\color{red}C_{6}}}{1}, 6\}, 1, 2, 3, {\color{red}C_{9}}`$ |  | |
| 9 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 2, 1, \underset{{\color{red}C_{6}}}{\underset{1}{\overset{{\color{red}C_{7}}}{\overset{1}{7}}}}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 10 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 2, 1, \{{\color{red}C_{6}}, 2, \overset{{\color{red}C_{7}}}{1}, 7\}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 11 | $`1, 3, 2, 1, \{{\color{red}C_{6}}, 3, \overset{{\color{red}C_{7}}}{1}, 6\}, 1, 3, {\color{red}C_{9}}`$ |  | |
| 12 | $`{\color{red}C_{1}}, 1, 3, 2, 2, 1, \underset{{\color{red}C_{7}}}{\underset{1}{\overset{{\color{red}C_{8}}}{\overset{1}{7}}}}, 1, {\color{red}C_{9}}`$ |  | |
| 13 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{3}}}{\overset{1}{4}}, 1, 4, 1, \overset{{\color{red}C_{8}}}{\overset{1}{4}}, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{8}} -> 2}`$<br>$`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$<br>$`\overset{{\color{red}C_{3}} -> 2}{{\color{red}C_{8}} -> 2}`$<br>$`\overset{{\color{red}C_{3}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 14 | $`1, \overset{{\color{red}C_{3}}}{\overset{1}{4}}, 1, 4, 1, \overset{2}{3}, 2`$ |  | |
| 15 | $`{\color{red}C_{1}}, 1, \underset{{\color{red}C_{3}}}{\underset{1}{\overset{{\color{red}C_{4}}}{\overset{1}{5}}}}, 1, 3, 2, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 1, {External[3], 1, {External[4], 1, 5}}, 1, 3, 2, 2, 1, External[9]}, True]<br>mdCell[{External[3] -> 2, External[9] -> 3}, {External[3] -> 3, External[9] -> 2}, {External[1], 1, {External[3], 1, {External[4], 1, 5}}, 1, 3, 2, 2, 1, External[9]}, True]<br>mdCell[{External[4] -> 2, External[9] -> 3}, {External[4] -> 3, External[9] -> 2}, {External[1], 1, {External[3], 1, {External[4], 1, 5}}, 1, 3, 2, 2, 1, External[9]}, True] | |
| 16 | $`{\color{red}C_{1}}, 2, 2, 1, \underset{{\color{red}C_{5}}}{\underset{1}{\overset{{\color{red}C_{6}}}{\overset{1}{8}}}}, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 17 | $`{\color{red}C_{1}}, 2, 2, 1, \{{\color{red}C_{5}}, 2, \overset{{\color{red}C_{6}}}{1}, 8\}, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 18 | $`2, 2, 2, 1, \underset{{\color{red}C_{6}}}{\underset{1}{\overset{{\color{red}C_{7}}}{\overset{1}{8}}}}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 19 | $`2, 2, 2, 1, \{{\color{red}C_{6}}, 2, \overset{{\color{red}C_{7}}}{1}, 8\}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 20 | $`2, 2, 2, 2, 1, \underset{{\color{red}C_{7}}}{\underset{1}{\overset{{\color{red}C_{8}}}{\overset{1}{8}}}}, 1, {\color{red}C_{9}}`$ |  | |
| 21 | $`2, \overset{{\color{red}C_{3}}}{\overset{1}{3}}, 2, 1, \overset{{\color{red}C_{7}}}{\overset{1}{7}}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 22 | $`{\color{red}C_{1}}, 3, 1, 3, 2, 1, \underset{{\color{red}C_{7}}}{\underset{1}{\overset{{\color{red}C_{8}}}{\overset{1}{7}}}}, 1, {\color{red}C_{9}}`$ |  | |
| 23 | $`{\color{red}C_{1}}, 3, 1, \{{\color{red}C_{4}}, 2, \overset{{\color{red}C_{5}}}{1}, 6\}, 1, 3, 1, 4`$ |  | |
| 24 | $`{\color{red}C_{1}}, 3, 1, \{{\color{red}C_{4}}, 3, \overset{{\color{red}C_{5}}}{1}, 5\}, 1, 3, 2, 1`$ |  | |
| 25 | $`{\color{red}C_{1}}, 3, 1, \{{\color{red}C_{4}}, 3, \overset{{\color{red}C_{5}}}{1}, 6\}, 1, 3, 1, 3`$ |  | |
| 26 | $`{\color{red}C_{1}}, 3, 2, 1, \underset{{\color{red}C_{5}}}{\underset{1}{\overset{{\color{red}C_{6}}}{\overset{1}{7}}}}, 1, 2, 3, {\color{red}C_{9}}`$ |  | |
| 27 | $`{\color{red}C_{1}}, 3, 2, 1, \{{\color{red}C_{5}}, 2, \overset{{\color{red}C_{6}}}{1}, 7\}, 1, 2, 3, {\color{red}C_{9}}`$ |  | |
| 29 | $`{\color{red}C_{1}}, 1, 2, \overset{{\color{red}C_{3}}}{3}, 1, \overset{{\color{red}C_{6}}}{\overset{1}{5}}, 1, \overset{{\color{red}C_{8}}}{3}, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{6}} -> 2}`$<br>$`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{8}} -> 2}`$<br>mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 1, 2, {External[3], 3}, 1, {External[6], 1, 5}, 1, {External[8], 3}, 1, External[9]}, True] | |
| 30 | $`{\color{red}C_{1}}, 1, 2, 3, 2, 1, \overset{{\color{red}C_{7}}}{\overset{1}{6}}, 1, 2, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{7}} -> 2}`$<br>$`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 31 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 6, 1, 2, \overset{{\color{red}C_{8}}}{\overset{1}{3}}, 2`$ |  | |
| 32 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, \overset{{\color{red}C_{5}}}{\overset{1}{6}}, 1, 2, \overset{{\color{red}C_{8}}}{3}, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 33 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, \overset{{\color{red}C_{5}}}{\overset{1}{6}}, 1, 3, 1, 3, {\color{red}C_{9}}`$ |  | |
| 34 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 2, 1, \overset{{\color{red}C_{6}}}{\overset{1}{7}}, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 35 | $`{\color{red}C_{1}}, 1, 3, 2, 2, 1, \overset{{\color{red}C_{7}}}{\overset{1}{7}}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 36 | $`{\color{red}C_{1}}, 1, 4, 1, 4, 1, \overset{{\color{red}C_{7}}}{\overset{1}{4}}, 1, 4`$ |  | |
| 37 | $`{\color{red}C_{1}}, 1, 5, 1, 2, 3, 1, \overset{{\color{red}C_{8}}}{\overset{1}{5}}, 1, {\color{red}C_{9}}`$ |  | |
| 38 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{3}}}{\overset{1}{4}}, 1, 4, 1, \overset{{\color{red}C_{7}}}{3}, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 1, {External[3], 1, 4}, 1, 4, 1, {External[7], 3}, 2, 1, External[9]}, True]<br>mdCell[{External[3] -> 2, External[9] -> 3}, {External[3] -> 3, External[9] -> 2}, {External[1], 1, {External[3], 1, 4}, 1, 4, 1, {External[7], 3}, 2, 1, External[9]}, True] | |
| 39 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{3}}}{\overset{1}{4}}, 1, 4, 1, 4, 1, 2, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$<br>$`\overset{{\color{red}C_{3}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 40 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{3}}}{\overset{1}{5}}, 1, 3, 2, 2, 1, 5`$ |  | |
| 41 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{3}}}{\overset{1}{6}}, 1, 2, 3, 1, 4, 1, {\color{red}C_{9}}`$ |  | |
| 42 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{3}}}{\overset{1}{6}}, 1, 2, 3, 2, 1, 4`$ |  | |
| 43 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{3}}}{\overset{1}{6}}, 1, 3, 1, 4, 1, 3, {\color{red}C_{9}}`$ |  | |
| 44 | $`{\color{red}C_{1}}, 2, 1, 5, 1, 3, 1, \overset{{\color{red}C_{8}}}{\overset{1}{5}}, 1, {\color{red}C_{9}}`$ |  | |
| 45 | $`{\color{red}C_{1}}, 2, 1, \overset{{\color{red}C_{4}}}{\overset{1}{5}}, 1, 3, 2, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 2, 1, {External[4], 1, 5}, 1, 3, 2, 2, 1, External[9]}, True]<br>mdCell[{External[4] -> 2, External[9] -> 3}, {External[4] -> 3, External[9] -> 2}, {External[1], 2, 1, {External[4], 1, 5}, 1, 3, 2, 2, 1, External[9]}, True] | |
| 46 | $`{\color{red}C_{1}}, 2, 1, \overset{{\color{red}C_{4}}}{\overset{1}{6}}, 1, 3, 1, 4, 1, {\color{red}C_{9}}`$ |  | |
| 47 | $`2, 2, 2, 1, \overset{{\color{red}C_{6}}}{\overset{1}{8}}, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 48 | $`2, 2, 2, 2, 1, \overset{{\color{red}C_{7}}}{\overset{1}{8}}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 49 | $`2, 2, 2, 2, 2, 1, \overset{{\color{red}C_{8}}}{\overset{1}{8}}, 1, {\color{red}C_{9}}`$ |  | |
| 50 | $`2, 3, 1, 3, 2, 1, \overset{{\color{red}C_{8}}}{\overset{1}{7}}, 1, {\color{red}C_{9}}`$ |  | |
| 51 | $`2, 3, 1, \overset{{\color{red}C_{5}}}{\overset{1}{5}}, 1, 3, 1, 5`$ |  | |
| 52 | $`{\color{red}C_{1}}, 2, \overset{{\color{red}C_{3}}}{\overset{1}{3}}, 1, 5, 1, \overset{{\color{red}C_{7}}}{3}, 2, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 53 | $`2, \overset{{\color{red}C_{3}}}{\overset{1}{3}}, 2, 1, 7, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 56 | $`{\color{red}C_{1}}, 3, 1, 3, 2, 1, \overset{{\color{red}C_{7}}}{\overset{1}{7}}, 1, 2, {\color{red}C_{9}}`$ |  | |
| 57 | $`{\color{red}C_{1}}, 3, 1, 3, 2, 2, 1, \overset{{\color{red}C_{8}}}{\overset{1}{7}}, 1, {\color{red}C_{9}}`$ |  | |
| 58 | $`{\color{red}C_{1}}, 3, 1, \overset{{\color{red}C_{4}}}{\overset{1}{5}}, 1, 3, 1, 5, 1, {\color{red}C_{9}}`$ |  | |
| 59 | $`{\color{red}C_{1}}, 3, 1, \overset{{\color{red}C_{4}}}{\overset{1}{5}}, 1, 3, 2, 1, 5`$ |  | |
| 60 | $`{\color{red}C_{1}}, 3, 1, \overset{{\color{red}C_{4}}}{\overset{1}{6}}, 1, 2, 3, 1, 4`$ |  | |
| 61 | $`{\color{red}C_{1}}, 3, 2, 1, \overset{{\color{red}C_{5}}}{\overset{1}{6}}, 1, 3, 1, 4`$ |  | |
| 62 | $`3, 2, 2, 1, \overset{{\color{red}C_{6}}}{\overset{1}{7}}, 1, 2, 3, {\color{red}C_{9}}`$ |  | |
| 63 | $`{\color{red}C_{1}}, 1, 2, 2, 2, 2, 2, 2, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 6}, {External[1] -> 3, External[9] -> 5}, {External[1] -> 4, External[9] -> 4}, {External[1] -> 5, External[9] -> 3}, {External[1] -> 6, External[9] -> 2}, {External[1], 1, 2, 2, 2, 2, 2, 2, 2, 1, External[9]}, True] |$``$ |
| 64 | $`{\color{red}C_{1}}, 1, 2, 2, 3, 1, 3, 2, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 5}, {External[1] -> 3, External[9] -> 4}, {External[1] -> 4, External[9] -> 3}, {External[1] -> 5, External[9] -> 2}, {External[1], 1, 2, 2, 3, 1, 3, 2, 2, 1, External[9]}, True] |$`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 3}, \overset{{\color{red}C_{1}} -> 3}{{\color{red}C_{9}} -> 2}`$ |
| 65 | $`{\color{red}C_{1}}, 1, 2, \overset{{\color{red}C_{3}}}{3}, 1, 4, 1, \overset{{\color{red}C_{7}}}{3}, 2, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{7}} -> 2}`$<br>mdCell[{External[1] -> 2, External[9] -> 4}, {External[1] -> 3, External[9] -> 3}, {External[1] -> 4, External[9] -> 2}, {External[1], 1, 2, {External[3], 3}, 1, 4, 1, {External[7], 3}, 2, 1, External[9]}, True]<br>$`\overset{{\color{red}C_{3}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 66 | $`{\color{red}C_{1}}, 1, 2, \overset{{\color{red}C_{3}}}{3}, 1, 4, 1, 4, 1, 2, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 1, 2, {External[3], 3}, 1, 4, 1, 4, 1, 2, External[9]}, True] | |
| 67 | $`{\color{red}C_{1}}, 1, 2, \overset{{\color{red}C_{3}}}{3}, 1, 5, 1, 3, 1, 3, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 68 | $`{\color{red}C_{1}}, 1, 2, 3, 2, 1, 6, 1, 2, 2, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 69 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 4, 1, 3, 2, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 4}, {External[1] -> 3, External[9] -> 3}, {External[1] -> 4, External[9] -> 2}, {External[1], 1, {External[2], 3}, 1, 4, 1, 3, 2, 2, 1, External[9]}, True]<br>mdCell[{External[2] -> 2, External[9] -> 3}, {External[2] -> 3, External[9] -> 2}, {External[1], 1, {External[2], 3}, 1, 4, 1, 3, 2, 2, 1, External[9]}, True] | |
| 70 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 4, 1, 4, 1, \overset{{\color{red}C_{8}}}{3}, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{8}} -> 2}`$<br>mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 1, {External[2], 3}, 1, 4, 1, 4, 1, {External[8], 3}, 1, External[9]}, True]<br>$`\overset{{\color{red}C_{2}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 71 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 5, 1, 2, 3, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 1, {External[2], 3}, 1, 5, 1, 2, 3, 2, 1, External[9]}, True]<br>$`\overset{{\color{red}C_{2}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 72 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 5, 1, 3, 1, 4, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 73 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 5, 1, 3, 2, 1, 4`$ |  | |
| 74 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 6, 1, 2, 2, 3, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 75 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 6, 1, 2, 3, 1, 3, {\color{red}C_{9}}`$ |  | |
| 76 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 1, 6, 1, 3, 1, 3, 2, {\color{red}C_{9}}`$ |  | |
| 77 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 2, 1, 5, 1, \overset{{\color{red}C_{7}}}{3}, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 1, {External[2], 3}, 2, 1, 5, 1, {External[7], 3}, 2, 1, External[9]}, True]<br>$`\overset{{\color{red}C_{2}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 78 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 2, 1, 6, 1, 2, \overset{{\color{red}C_{8}}}{3}, 1, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 79 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 2, 1, 6, 1, 3, 1, 3, {\color{red}C_{9}}`$ |  | |
| 80 | $`{\color{red}C_{1}}, 1, \overset{{\color{red}C_{2}}}{3}, 2, 1, 7, 1, 2, 2, 2`$ |  | |
| 81 | $`{\color{red}C_{1}}, 1, 3, 2, 2, 1, 7, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 82 | $`{\color{red}C_{1}}, 1, 4, 1, 3, 1, 6, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 83 | $`{\color{red}C_{1}}, 1, 4, 1, 3, 2, 1, 6, 1, 2, {\color{red}C_{9}}`$ |  | |
| 84 | $`{\color{red}C_{1}}, 1, 4, 1, 3, 2, 2, 1, 6, 1, {\color{red}C_{9}}`$ |  | |
| 85 | $`{\color{red}C_{1}}, 1, 4, 1, 4, 1, 3, 1, 5, 1, {\color{red}C_{9}}`$ |  | |
| 86 | $`{\color{red}C_{1}}, 1, 4, 1, 4, 1, 3, 2, 1, 5`$ |  | |
| 87 | $`{\color{red}C_{1}}, 1, 5, 1, 2, 3, 2, 1, 5, 1, {\color{red}C_{9}}`$ |  | |
| 88 | $`{\color{red}C_{1}}, 1, 6, 1, 2, 3, 1, 4, 1, 3, {\color{red}C_{9}}`$ |  | |
| 89 | $`{\color{red}C_{1}}, 1, 6, 1, 3, 1, 4, 1, 3, 2, {\color{red}C_{9}}`$ |  | |
| 90 | $`{\color{red}C_{1}}, 2, 1, 4, 1, 4, 1, 4, 1, 2, {\color{red}C_{9}}`$ | $`\overset{{\color{red}C_{1}} -> 2}{{\color{red}C_{9}} -> 2}`$ | |
| 91 | $`{\color{red}C_{1}}, 2, 1, 5, 1, 3, 1, 5, 1, 2, {\color{red}C_{9}}`$ |  | |
| 92 | $`{\color{red}C_{1}}, 2, 1, 5, 1, 3, 2, 1, 5, 1, {\color{red}C_{9}}`$ |  | |
| 93 | $`{\color{red}C_{1}}, 2, 1, 5, 1, 3, 2, 2, 1, 5`$ |  | |
| 94 | $`{\color{red}C_{1}}, 2, 2, 1, 5, 1, 3, 2, 2, 1, {\color{red}C_{9}}`$ | mdCell[{External[1] -> 2, External[9] -> 3}, {External[1] -> 3, External[9] -> 2}, {External[1], 2, 2, 1, 5, 1, 3, 2, 2, 1, External[9]}, True] | |
| 96 | $`2, 2, 2, 2, 1, 8, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 97 | $`2, 2, 2, 2, 2, 1, 8, 1, 2, {\color{red}C_{9}}`$ |  | |
| 98 | $`2, 2, 2, 2, 2, 2, 1, 8, 1, {\color{red}C_{9}}`$ |  | |
| 100 | $`2, 2, 3, 1, 3, 2, 1, 7, 1, {\color{red}C_{9}}`$ |  | |
| 102 | $`2, 3, 1, 3, 2, 1, 7, 1, 2, {\color{red}C_{9}}`$ |  | |
| 103 | $`2, 3, 1, 3, 2, 2, 1, 7, 1, {\color{red}C_{9}}`$ |  | |
| 106 | $`{\color{red}C_{1}}, 3, 1, 3, 1, 6, 1, 3, 1, 3, {\color{red}C_{9}}`$ |  | |
| 107 | $`{\color{red}C_{1}}, 3, 1, 3, 2, 1, 7, 1, 2, 2, {\color{red}C_{9}}`$ |  | |
| 108 | $`{\color{red}C_{1}}, 3, 1, 3, 2, 2, 1, 7, 1, 2, {\color{red}C_{9}}`$ |  | |
| 109 | $`{\color{red}C_{1}}, 3, 1, 4, 1, 3, 1, 6, 1, 2, {\color{red}C_{9}}`$ |  | |
| 110 | $`{\color{red}C_{1}}, 3, 1, 4, 1, 4, 1, 3, 1, 5`$ |  | |
| 111 | $`{\color{red}C_{1}}, 3, 1, 5, 1, 2, 3, 1, 5, 1, {\color{red}C_{9}}`$ |  | |
| 112 | $`{\color{red}C_{1}}, 3, 1, 5, 1, 2, 3, 2, 1, 5`$ |  | |
| 113 | $`{\color{red}C_{1}}, 3, 1, 5, 1, 3, 1, 4, 1, 4`$ |  | |
| 114 | $`{\color{red}C_{1}}, 3, 2, 1, 5, 1, 3, 1, 5, 1, {\color{red}C_{9}}`$ |  | |
| 115 | $`{\color{red}C_{1}}, 3, 2, 1, 5, 1, 3, 2, 1, 5`$ |  | |
| 117 | $`4, 1, 2, 3, 2, 1, 6, 1, 2, {\color{red}C_{9}}`$ |  | |
| 120 | $`4, 1, 3, 2, 1, 6, 1, 2, 3, {\color{red}C_{9}}`$ |  | |
| 121 | $`4, 1, 3, 2, 2, 1, 6, 1, 3, {\color{red}C_{9}}`$ |  | |
| 122 | $`6, 1, 2, 2, 3, 1, 4, 1, 3, {\color{red}C_{9}}`$ |  | |
| 123 | $`6, 1, 2, 3, 1, 4, 1, 3, 2, {\color{red}C_{9}}`$ |  | |