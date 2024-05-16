import sys
from typing import List, Dict, Iterable, Tuple

# Please do not remove package declarations because these are used by the autograder.

def burrows_wheeler_transform(text: str) -> str:
    """
    Generate the Burrows-Wheeler Transform of the given text.
    """
    cyclic_rotations = []
    for i in range(0,len(text)):
        cyclic_rotations.append(text[i:] + text[:i])
    M = sorted(cyclic_rotations)
    bwt = ''
    for rotation in M:
        bwt = bwt + rotation[len(rotation)-1]
    return bwt

def suffix_array(text: str) -> List[int]:
    """
    Generate the suffix array for the given text.
    """
    s_array = []
    suffixes = []
    for i in range(0, len(text)):
        suffixes.append(text[i:])
    sorted_suffixes = sorted(suffixes)
    for suffix in sorted_suffixes:
        s_array.append(text.index(suffix))
    return s_array

def partial_suffix_array(text: str, k: int) -> dict[int, int]:
    """
    Generate a partial suffix array for the given text and interval K.
    """
    full_suffix_array = suffix_array(text)
    partial = dict()
    for idx in full_suffix_array:
        if full_suffix_array[idx] % k == 0:
            partial[idx] = full_suffix_array[idx]
    return partial

def compute_rank_arr(bwt: str) -> list[int]:
    """
    This function generates the rank of each position in the last column given by the bwt.
    The rank is the number of occurrences of whatever character is at that position, up to
    that position. This can be done in linear time by iterating through the bwt. The ranks
    will be returned in the form of a list of ranks, obviously indices will be in-built.
    """
    rank = []
    counts = dict()
    for char in bwt:
        if char not in counts: counts[char] = 0
        counts[char] += 1
        rank.append(counts[char])
    return rank

def compute_first_occurrences(bwt: str) -> dict[str, int]: # the mapping of c in C is the row at which the character c appears in the first column for the first time.
    """
    Generate a dict where each 'character' is mapped to the index in first column where
    these characters first appeared. In other words because the first column is in alphabetical
    order, we can count the ascii code of each character in the last column, then iterate
    from 0 to 255 to get the count of each ascii character in ascending lexicographic order.
    This is done in linear time.
    """
    C = dict()
    counts = [0 for _ in range(256)]
    for char in bwt:
        C[char] = 0
        counts[ord(char)] += 1
    curr_idx = 0
    for i in range(0,256):
        if counts[i] != 0:
            C[chr(i)] = curr_idx
        for _ in range(counts[i]):
            curr_idx += 1
    return C

def compute_checkpoint_arrs(bwt: str) -> dict[int, list[int]]:
    """
    Similar to ranks, but instead the list stored contains the rank of every character up to
    that index, if the index % C is 0. More memory efficient.
    """
    C = 5
    ranks = dict()
    rank = [0 for _ in range(256)]
    for i in range(0, len(bwt)):
        rank[ord(bwt[i])] += 1
        if i % C == 0:
            ranks[i] = rank[:]
    return ranks

def bw_better_match_pattern(bwt: str, pattern: str, first_occurrences: dict[str,int], ranks: list[list[int]]) -> tuple[int,int]:
    C = 5
    top = 0
    bot = len(bwt) - 1
    for i in range(len(pattern) - 1, -1, -1):
        symbol = pattern[i]
        if symbol not in first_occurrences: # in the case the symbol is not in text at all
            return (0,0)
        top_rank = compute_rank(bwt, top, ranks, symbol, C) # use checkpoint arrs to get the rank
        bot_rank = compute_rank(bwt, bot, ranks, symbol, C)
        marker = False
        if bwt[top] == symbol:
            marker = True
        top = first_occurrences[symbol] + top_rank
        if marker:
            top -= 1
        bot = first_occurrences[symbol] + bot_rank - 1
        if bot - top < 0:
            return (0,0)
    return (top, bot + 1)

def compute_idxes_from_top_bot(start: int, end: int, partial_s_array: dict[int, int], bwt: str, rank: list[int], occurrences: list[int]) -> list[int]:
    pattern_idxes = []
    for i in range(start, end):
        p = i
        plus_count = 0
        while True:
            predecessor = bwt[p]
            p = occurrences[predecessor] + rank[p] - 1
            plus_count += 1
            if p in iter(partial_s_array.keys()):
                break
        pattern_idxes.append((partial_s_array[p] + plus_count) % len(bwt))
    return pattern_idxes

def compute_rank(bwt: str, idx: int, ranks: dict[int, list[int]], symbol: str, C: int) -> int: # idx can be either top or bot
    idx_dist_from_chkpnt = idx % C
    idx_rank = ranks[idx - idx_dist_from_chkpnt][ord(symbol)]
    for j in range(idx - idx_dist_from_chkpnt + 1, idx + 1):
        if bwt[j] == symbol:
            idx_rank += 1
    return idx_rank

# Insert your multiple_pattern_matching function here, along with any subroutines you need
def multiple_pattern_matching(text: str, patterns: List[str]) -> Dict[str, List[int]]:
    """
    Find all starting positions in text where each string from patterns appears as a substring.
    """
    text = text + '$'
    pattern_indices = dict()
    bwt = burrows_wheeler_transform(text)
    first_occurrences = compute_first_occurrences(bwt)
    ranks = compute_checkpoint_arrs(bwt)
    rank = compute_rank_arr(bwt)
    K = 5
    partial_s_array = partial_suffix_array(text, K)
    for pattern in patterns:
        start, end = bw_better_match_pattern(bwt, pattern, first_occurrences, ranks)
        idxes = compute_idxes_from_top_bot(start, end, partial_s_array, bwt, rank, first_occurrences)
        pattern_indices[pattern] = idxes
    return pattern_indices

def main():
    text = 'AATCGGGTTCAATCGGGGT'
    patterns = ['ATCG', 'GGGT']
    print(multiple_pattern_matching(text, patterns))
    
    text = 'ATATATATAT'
    patterns = ['GT', 'AGCT', 'TAA', 'AAT', 'AATAT']
    print(multiple_pattern_matching(text, patterns))
    
    text = 'bananas'
    patterns = ['ana', 'as']
    print(multiple_pattern_matching(text, patterns))
    
    text = 'AAACAA'
    patterns = ['AA']
    print(multiple_pattern_matching(text, patterns))

if __name__ == '__main__':
    main()