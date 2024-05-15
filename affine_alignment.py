import sys

"""
Summary:

This is a dynamic programming algorithm that returns the highest scoring global alignment of
two strings using gap opening and gap extension penalties, therefore using the idea that
longer gaps in an alignment should be penalized less harshly than multiple shorter gaps.

The algorithm utilizes 3 matrices - the lower matrix, for which each [i][j] represents the
highest scoring alignment of s and t ending with a gap in s, the upper matrix, which is the
same as lower but with a gap in t at the end, and the middle matrix, for which [i][j] represents
the highest scoring alignment between s and t ending in s[i] and t[j] being aligned.

The recurrence relation is defined as:

lower[i][j] is the maximum of lower[i-1][j] - gap_extension penalty and middle[i-1][j] - gap_opening penalty.
    This represents the penalty for either opening a gap or extending it.
    
upper[i][j] is max(upper[i][j-1] - gap_extension_pen, middle[i][j-1] - gap_opening_pen)
    Same logic as lower.
    
middle[i][j] is max(lower[i][j], upper[i][j], middle[i-1][j-1] + value of aligning s[i], t[j])
    lower[i][j] and upper[i][j] are a factor here because it costs nothing to close a gap.
"""

# insert your AffineAlignment function here, along with any subroutines you need
def AffineAlignment(match_reward: int, mismatch_penalty: int,
                    gap_opening_penalty: int, gap_extension_penalty: int,
                    s: str, t: str) -> tuple[int, str, str]:
    """Generates the affine alignment of two strings s and t."""
    sys.setrecursionlimit(1500)
    lower = [[0] * (len(t) + 1) for _ in range(len(s) + 1)] # insertions matrix
    upper = [[0] * (len(t) + 1) for _ in range(len(s) + 1)] # deletions matrix
    middle = [[0] * (len(t) + 1) for _ in range(len(s) + 1)] # match/mismatch matrix
    b_lower = [[None] * (len(t)) for _ in range(len(s))] # backtrack insertions
    b_upper = [[None] * (len(t)) for _ in range(len(s))] # backtrack deletions
    b_middle = [[None] * (len(t)) for _ in range(len(s))] # backtrack ms/mms
    for i in range(0, len(s) + 1):
        for j in range(0, len(t) + 1):
            # set base cases for middle using gap opening and extension penalties, base for lower and upper is -infinity.
            if i == 0 and j == 0:
                lower[i][j], upper[i][j] = float('-inf'), float('-inf')
                continue
            elif i == 0:
                lower[i][j], middle[i][j] = float('-inf'), -gap_opening_penalty - ((j-1) * gap_extension_penalty)
                continue
            elif j == 0:
                upper[i][j], middle[i][j] = float('-inf'), -gap_opening_penalty - ((i-1) * gap_extension_penalty)
                continue
            
            # test for match
            match = -mismatch_penalty
            if s[i-1] == t[j-1]:
                match = match_reward
            
            lower[i][j] = max(lower[i-1][j] - gap_extension_penalty, middle[i-1][j] - gap_opening_penalty)
            upper[i][j] = max(upper[i][j-1] - gap_extension_penalty, middle[i][j-1] - gap_opening_penalty)
            middle[i][j] = max(lower[i][j], upper[i][j], middle[i-1][j-1] + match)
            if lower[i][j] == lower[i-1][j] - gap_extension_penalty:
                b_lower[i-1][j-1] = 'd'
            elif lower[i][j] == middle[i-1][j] - gap_opening_penalty:
                b_lower[i-1][j-1] = 'open'
            if upper[i][j] == upper[i][j-1] - gap_extension_penalty:
                b_upper[i-1][j-1] = 'r'
            elif upper[i][j] == middle[i][j-1] - gap_opening_penalty:
                b_upper[i-1][j-1] = 'open'
            if middle[i][j] == lower[i][j]:
                b_middle[i-1][j-1] = 'in_close'
            elif middle[i][j] == upper[i][j]:
                b_middle[i-1][j-1] = 'del_close'
            elif middle[i][j] == middle[i-1][j-1] + match:
                b_middle[i-1][j-1] = 'dr'
    s_align, t_align = backtrack(s, t, b_lower, b_upper, b_middle, len(s) - 1, len(t) - 1, 'middle')
    return max(lower[len(s)][len(t)], middle[len(s)][len(t)], upper[len(s)][len(t)]), s_align, t_align
                
def backtrack(s: str, t: str, b_lower: list[list[str]], b_upper: list[list[str]], b_middle: list[list[str]], i: int, j: int, LEVEL: str) -> tuple[str,str]:
    if i < 0 and j < 0: # if both i and j reach -1 at the same time, we can just return empty string.
        return '',''
    elif i < 0: # if i reaches -1 first, need to append the rest of t up to index j to t prime, and a number of '-' equal to the length of that string to s prime.
        return '-' * (j + 1), t[0:j+1]
    elif j < 0: # same as i < 0 condition but for j.
        return s[0:i+1], '-' * (i + 1)
    # my use of the LEVEL variable tells us which matrix this call is backtracking in.
    if LEVEL == 'middle':
        if b_middle[i][j] == 'in_close':
            s_prime, t_prime = backtrack(s, t, b_lower, b_upper, b_middle, i, j, 'lower')
            return s_prime, t_prime
        elif b_middle[i][j] == 'del_close':
            s_prime, t_prime = backtrack(s, t, b_lower, b_upper, b_middle, i, j, 'upper')
            return s_prime, t_prime
        elif b_middle[i][j] == 'dr':
            s_prime, t_prime = backtrack(s, t, b_lower, b_upper, b_middle, i-1, j-1, 'middle')
            return s_prime + s[i], t_prime + t[j]
    elif LEVEL == 'lower':
        if b_lower[i][j] == 'd':
            s_prime, t_prime = backtrack(s, t, b_lower, b_upper, b_middle, i-1, j, 'lower')
            return s_prime + s[i], t_prime + '-'
        elif b_lower[i][j] == 'open':
            s_prime, t_prime = backtrack(s, t, b_lower, b_upper, b_middle, i-1, j, 'middle')
            return s_prime + s[i], t_prime + '-'
    elif LEVEL == 'upper':
        if b_upper[i][j] == 'r':
            s_prime, t_prime = backtrack(s, t, b_lower, b_upper, b_middle, i, j-1, 'upper')
            return s_prime + '-', t_prime + t[j]
        elif b_upper[i][j] == 'open':
            s_prime, t_prime = backtrack(s, t, b_lower, b_upper, b_middle, i, j-1, 'middle')
            return s_prime + '-', t_prime + t[j]
        
def main():
    print(AffineAlignment(1,5,2,1, 'CCAT', 'GAT'))

if __name__ == '__main__':
    main()