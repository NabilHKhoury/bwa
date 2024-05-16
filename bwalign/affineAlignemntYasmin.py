import sys
from typing import List, Dict, Iterable, Tuple
import numpy as np

# Please do not remove package declarations because these are used by the autograder.

# Insert your AffineAlignment function here, along with any subroutines you need
def AffineAlignment(match_reward: int, mismatch_penalty: int,
                    gap_opening_penalty: int, gap_extension_penalty: int,
                    s: str, t: str) -> Tuple[int, str, str]:
    len_s, len_t = len(s), len(t)

    # Initialize matrices
    lower = np.full((len_s + 1, len_t + 1),float('-inf'))#was -inf
    middle = np.full((len_s + 1, len_t + 1),float('-inf'))
    upper = np.full((len_s + 1, len_t + 1),float('-inf'))
    middle[0, 0] = 0  # Starting point
    lower[0,0] = float('-inf')#-gap_opening_penalty
    upper[0,0] = float('-inf')#-gap_opening_penalty
    backtrack = []
    
    #initialize backtrack:
    for i in range(len(s)+1):
        row = []
        for j in range(len(t)+1):
            row.append("")
        backtrack.append(row)

    # Initialize the first row and column of the middle matrix to account for leading gaps
    # Initialize the first row of the upper matrix and the first column of the lower matrix
    for i in range(1, len_s + 1):
        lower[i, 0] = 0 - ((i - 1) * gap_extension_penalty) - gap_opening_penalty
        middle[i, 0] = 0 - ((i - 1) * gap_extension_penalty) - gap_opening_penalty
        upper[i, 0] = 0 - ((i - 1) * gap_extension_penalty) - 2*gap_opening_penalty #rem *2
        backtrack[i][0] = 'up'               
    # No need to initialize `middle[i, 0]` here since it's already done above

    for j in range(1, len_t + 1):
        upper[0, j] = 0 - (j - 1) * gap_extension_penalty - gap_opening_penalty
        lower[0, j] = 0 - ((j - 1) * gap_extension_penalty) - 2*gap_opening_penalty#remived x2
        middle[0, j] = 0 - (j - 1) * gap_extension_penalty - gap_opening_penalty
        backtrack[0][j] = 'left'
        
    # Fill in the matrices based on the recursive formulas
    for i in range(1, len_s + 1):
        for j in range(1, len_t + 1):
            char_s = s[i-1]
            char_t = t[j-1]

            # Calculate scores for each matrix
            lower[i, j] = max(lower[i-1, j] - gap_extension_penalty,
                          middle[i-1, j] - gap_opening_penalty)

            # Update upper matrix: Max of extending a gap from upper or opening a new gap from middle
            upper[i, j] = max(upper[i, j-1] - gap_extension_penalty,
                              middle[i, j-1] - gap_opening_penalty)

            # Update middle matrix: Max of diagonal move (match/mismatch), ending a gap from lower or upper
            middle[i, j] = max(middle[i-1, j-1] + (match_reward if s[i-1] == t[j-1] else -mismatch_penalty),
                               lower[i, j],
                               upper[i, j])
            
            #do Backtrack
            if(middle[i, j] == middle[i-1, j-1] - mismatch_penalty):
                backtrack[i][j] = "diag"
            elif(middle[i, j] == upper[i,j]):
                backtrack[i][j] = "left"
            elif(middle[i, j] == lower[i,j]):
                backtrack[i][j] = "up"
            elif(middle[i, j] == middle[i-1, j-1] + match_reward):
                backtrack[i][j] = "diag"
            
                

    # Traceback to find the optimal alignment
    alignment_s = ""
    alignment_t = ""
    i, j = len_s, len_t

    while i > 0 or j > 0:
        if backtrack[i][j] == 'diag' and i > 0 and j > 0 :
            alignment_s = s[i-1] + alignment_s
            alignment_t = t[j-1] + alignment_t
            i, j = i-1, j-1
        elif backtrack[i][j] == 'up' and i > 0:
            alignment_s = s[i-1] + alignment_s
            alignment_t = '-' + alignment_t
            i -= 1
        elif backtrack[i][j] == 'left' and j > 0:
            alignment_s = '-' + alignment_s
            alignment_t = t[j-1] + alignment_t
            j -= 1

    # Return the score of the best alignment and the alignments themselves
    return int(middle[len_s, len_t]), alignment_s, alignment_t
