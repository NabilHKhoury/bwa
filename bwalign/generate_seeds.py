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

def partial_suffix_array(text: str, k: int) -> dict[int, int]:
    """
    Generate a partial suffix array for the given text and interval K.
    """
    full_suffix_array = sorted(range(len(text)), key=lambda i: text[i:]) # taken from Nabil's bw code
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

def compute_rank(bwt: str, idx: int, ranks: dict[int, list[int]], symbol: str, C: int) -> int: # idx can be either top or bot
    idx_dist_from_chkpnt = idx % C
    idx_rank = ranks[idx - idx_dist_from_chkpnt][ord(symbol)]
    for j in range(idx - idx_dist_from_chkpnt + 1, idx + 1):
        if bwt[j] == symbol:
            idx_rank += 1
    return idx_rank

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

def generate_seeds(read: str, bwt: str, k: int, psa: dict[int, int],
                                                first_occurrences: dict[str, int],
                                                checkpoint_arrs: dict[int, list[int]],
                                                ranks: list[int]) -> list[list[int]]:
    """
    Takes a read and bwt created from reference genome, and generates a list of lists
    with each index being an index i in the read from 0 to len(read) - k + 1. The corresponding list at each
    index is a list of exact match indices of the kmer at read[i:i+k] located in the reference
    genome. Note that even if two kmers are identical, their indices in the read are not.
    """
    seed_idxes = []
    for i in range(0, len(read) - k + 1):
        kmer = read[i:i+k]
        start, end = bw_better_match_pattern(bwt, kmer, first_occurrences, checkpoint_arrs)
        idxes = compute_idxes_from_top_bot(start, end, psa, bwt, ranks, first_occurrences)
        seed_idxes.append(idxes)
    return seed_idxes

def main():
    """Example string finder preprocessing and usage"""
    ref = 'AATCGGGTTCAATCGGGGTAATCGGGTTCAATCGGGGT' # both ref and read will be passed in at another time
    read = 'TCGGGTTCAATCGG'
    # read2 = 'asdfasdfasdfas' # test chaff read
    
    K = 5 # partial suffix array constant
    text = ref + '$'
    bwt = burrows_wheeler_transform(text)
    psa = partial_suffix_array(text, K)
    first_occurrences = compute_first_occurrences(bwt)
    checkpoint_arrs = compute_checkpoint_arrs(bwt)
    ranks = compute_rank_arr(bwt)
    seeds = generate_seeds(read, bwt, 8, psa, first_occurrences, checkpoint_arrs, ranks)
    print(seeds)
    
if __name__ == '__main__':
    main()