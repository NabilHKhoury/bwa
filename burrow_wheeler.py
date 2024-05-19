from typing import List, Dict, Iterable, Tuple

class BWT:
    def __init__(self, bwt):
        self.bwt = bwt.upper()
        self.alphabet = sorted(set(self.bwt))  # all elements in the bwt
        self.lastColumn = [char for char in self.bwt]  
        self.firstColumn = sorted(self.lastColumn)  
        self.firstOccurence = self.FirstOccurence()
        self.count = self.Count()
        self.suffix_array = self.suffix_array(self.bwt)
    def FirstOccurence(self):
        return {char: self.firstColumn.index(char) for char in self.alphabet}

    def Count(self):
        count = {char: [0] * (len(self.bwt) + 1) for char in self.alphabet}
        #  increment counts over last column
        for i, char in enumerate(self.bwt):
            for acid in self.alphabet:
                count[acid][i + 1] = count[acid][i] + (1 if acid == char else 0)
        return count
    
    def BetterBWMatching(self, pattern):
        top = 0
        bottom = len(self.lastColumn) - 1

        while top <= bottom:
            if pattern:
                symbol = pattern[-1]
                pattern = pattern[: len(pattern) - 1]
                if symbol in self.alphabet and self.count[symbol][top] < self.count[symbol][bottom + 1]:
                    top = self.firstOccurence[symbol] + self.count[symbol][top]
                    bottom = self.firstOccurence[symbol] + self.count[symbol][bottom + 1] - 1
                else:
                    return 0
            else:
                return bottom - top + 1
        return 0
        
    def suffix_array(sequence):
        return sorted(range(len(sequence)), key=lambda i: sequence[i:])
    
    def find_seeds(bwt, sa, read, k):
        seeds = []
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            l, r = 0, len(bwt)
            for char in reversed(kmer):
                l = bwt.find(char, l, r)
                r = bwt.find(char, l, r) + bwt.count(char, l, r)
                if l == r:
                    break
            if l < r:
                seeds.append((kmer, sa[l:r]))
        return seeds
    
    def fm_index(sequence):
        bwt = bwt.bwt 
        sa = bwt.suffix_array()
        return bwt, sa

def better_bw_matching(bwt: str, patterns: List[str]) -> List[int]:
    bwMatching = []
    bw = BWT(bwt)
    for pattern in patterns:
        bwMatching.append(bw.BetterBWMatching(pattern))
    return bwMatching

