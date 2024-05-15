from typing import List, Dict, Iterable, Tuple
import burrow_wheeler as bw
def linear_space_index(s1, s2, match_score, mismatch_penalty, indel_penalty):
    m, n = len(s1), len(s2)
    # Initialize the two rows for storing the scores.
    previous_row = [0] * (n + 1)
    current_row = [0] * (n + 1)
    
    # Variables to track the maximum score and its position.
    max_score = 0
    max_position = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i-1] == s2[j-1]:
                score = match_score
            else:
                score = mismatch_penalty
            
            # Calculate the scores for the three possibilities.
            diagonal = previous_row[j-1] + score
            up = previous_row[j] + indel_penalty
            left = current_row[j-1] + indel_penalty
            
            # Choose the best score.
            current_row[j] = max(0, diagonal, up, left)
            
            # Update the maximum score and its position.
            if current_row[j] > max_score:
                max_score = current_row[j]
                max_position = (i, j)

        # Prepare for the next iteration.
        previous_row, current_row = current_row, [0] * (n + 1)
    
    # Traceback to get the aligned sequences.
            
    # Convert the alignments to strings.
    return max_score, max_position
#this function is from cse 181 that i wrote.
def GobalAlignment(s, t, matchReward, mismatchPenalty, indelPenalty):

    m, n = len(s), len(t)
    dp = [[0 for i in range(n + 1)] for j in range(m + 1)]
    direction = [[None for i in range(n + 1)] for j in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = i * -indelPenalty
        direction[i][0] = 'up'
    for j in range(1, n + 1):
        dp[0][j] = j * -indelPenalty
        direction[0][j] = 'left'

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                across = dp[i - 1][j - 1] + matchReward
            else:
                across = dp[i - 1][j - 1] - mismatchPenalty
            down = dp[i - 1][j] - indelPenalty
            right = dp[i][j - 1] - indelPenalty
            dp[i][j] = max(across, down, right)
            
            if (dp[i][j] == right):
                direction[i][j] = 'left'
            elif (dp[i][j] == down):
                direction[i][j] = 'up'
            else:
                direction[i][j] = 'diag'
    sA, tA = traceback(direction, s, t)
    return dp[m][n], sA, tA
def traceback(direction, s, t):
    s_aligned, t_aligned = '', ''
    i, j = len(s), len(t)
    
    while (i > 0 or j > 0):
        if direction[i][j] == 'diag':
            s_aligned = s[i-1] + s_aligned
            t_aligned = t[j-1] + t_aligned
            i, j = i-1, j-1
        elif (direction[i][j] == 'left'):
            s_aligned = '-' + s_aligned
            t_aligned = t[j-1] + t_aligned
            j -= 1
        elif (direction[i][j] == 'up'):
            s_aligned = s[i-1] + s_aligned
            t_aligned = '-' + t_aligned
            i -= 1
            
    return s_aligned, t_aligned


# Insert your LocalAlignment function here, along with any subroutines you need
def LocalAlignment(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                    s: str, t: str) -> Tuple[int, str, str]:
        max_score1, end_max_position = linear_space_index(s, t, match_reward, -mismatch_penalty, -indel_penalty)  # Adjusted this line
        i_end, j_end = end_max_position
        rev_seq1 = s[: i_end][::-1]
        rev_seq2 = t[: j_end][::-1]
        max_score2, start_max_position = linear_space_index(rev_seq1, rev_seq2, match_reward, -mismatch_penalty, -indel_penalty)
        i_start, j_start = start_max_position
        i_start = i_end - i_start
        j_start = j_end - j_start
        alignment1 = s[i_start: i_end]
        alignment2 = t[j_start: j_end]
        return GobalAlignment(alignment1, alignment2, match_reward, mismatch_penalty, indel_penalty)

