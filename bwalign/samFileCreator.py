
#first install using pip install pysam biopython

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Tuple, List, Dict
import numpy as np

#AffineAlignment function
def AffineAlignment(match_reward: int, mismatch_penalty: int,
                    gap_opening_penalty: int, gap_extension_penalty: int,
                    s: str, t: str) -> Tuple[int, str, str]:
    len_s, len_t = len(s), len(t)

    # Initialize matrices
    lower = np.full((len_s + 1, len_t + 1),float('-inf'))
    middle = np.full((len_s + 1, len_t + 1),float('-inf'))
    upper = np.full((len_s + 1, len_t + 1),float('-inf'))
    middle[0, 0] = 0  # Starting point
    lower[0,0] = float('-inf')
    upper[0,0] = float('-inf')
    backtrack = []

    #initialize backtrack:
    for i in range(len(s)+1):
        row = []
        for j in range(len(t)+1):
            row.append("")
        backtrack.append(row)

    # Initialize the first row and column of the middle matrix to account for leading gaps
    for i in range(1, len_s + 1):
        lower[i, 0] = 0 - ((i - 1) * gap_extension_penalty) - gap_opening_penalty
        middle[i, 0] = 0 - ((i - 1) * gap_extension_penalty) - gap_opening_penalty
        upper[i, 0] = 0 - ((i - 1) * gap_extension_penalty) - 2*gap_opening_penalty
        backtrack[i][0] = 'up'

    for j in range(1, len_t + 1):
        upper[0, j] = 0 - (j - 1) * gap_extension_penalty - gap_opening_penalty
        lower[0, j] = 0 - ((j - 1) * gap_extension_penalty) - 2*gap_opening_penalty
        middle[0, j] = 0 - (j - 1) * gap_extension_penalty - gap_opening_penalty
        backtrack[0][j] = 'left'

    # Fill in the matrices based on the recursive formulas
    for i in range(1, len_s + 1):
        for j in range(1, len_t + 1):
            char_s = s[i-1]
            char_t = t[j-1]

            lower[i, j] = max(lower[i-1, j] - gap_extension_penalty,
                          middle[i-1, j] - gap_opening_penalty)

            upper[i, j] = max(upper[i, j-1] - gap_extension_penalty,
                              middle[i, j-1] - gap_opening_penalty)

            middle[i, j] = max(middle[i-1, j-1] + (match_reward if s[i-1] == t[j-1] else -mismatch_penalty),
                               lower[i, j],
                               upper[i, j])

            if(middle[i, j] == middle[i-1, j-1] - mismatch_penalty):
                backtrack[i][j] = "diag"
            elif(middle[i, j] == upper[i,j]):
                backtrack[i][j] = "left"
            elif(middle[i, j] == lower[i,j]):
                backtrack[i][j] = "up"
            elif(middle[i, j] == middle[i-1, j-1] + match_reward):
                backtrack[i][j] = "diag"

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

    return int(middle[len_s, len_t]), alignment_s, alignment_t

#reference genome
reference = SeqIO.read("reference.fasta", "fasta")

#paired-end reads
reads1 = list(SeqIO.parse("reads_1.fastq", "fastq"))
reads2 = list(SeqIO.parse("reads_2.fastq", "fastq"))

#create the SAM file header
header = {
    'HD': {'VN': '1.0', 'SO': 'unsorted'},
    'SQ': [{'SN': reference.id, 'LN': len(reference)}]
}

#open the SAM file for writing
with pysam.AlignmentFile("output.sam", "w", header=header) as samfile:
    for read1, read2 in zip(reads1, reads2):
        #align read1
        score1, alignment1_s, alignment1_t = AffineAlignment(2, -1, -2, -1, str(read1.seq), str(reference.seq))
        pos1 = alignment1_t.find(alignment1_t.replace('-', '')) + 1  # Position in reference

        #align read2
        score2, alignment2_s, alignment2_t = AffineAlignment(2, -1, -2, -1, str(read2.seq), str(reference.seq))
        pos2 = alignment2_t.find(alignment2_t.replace('-', '')) + 1  # Position in reference

        #create SAM entry for read1
        a1 = pysam.AlignedSegment()
        a1.query_name = read1.id
        a1.query_sequence = str(read1.seq)
        a1.flag = 0
        a1.reference_id = 0
        a1.reference_start = pos1 - 1
        a1.mapping_quality = 60
        a1.cigarstring = f"{len(alignment1_s.replace('-', ''))}M"
        a1.next_reference_id = 0
        a1.next_reference_start = pos2 - 1
        a1.template_length = abs(pos2 - pos1)
        a1.query_qualities = pysam.qualitystring_to_array(read1.letter_annotations["phred_quality"])

        #create SAM entry for read2
        a2 = pysam.AlignedSegment()
        a2.query_name = read2.id
        a2.query_sequence = str(read2.seq)
        a2.flag = 16  # Reverse strand
        a2.reference_id = 0
        a2.reference_start = pos2 - 1
        a2.mapping_quality = 60
        a2.cigarstring = f"{len(alignment2_s.replace('-', ''))}M"
        a2.next_reference_id = 0
        a2.next_reference_start = pos1 - 1
        a2.template_length = -a1.template_length
        a2.query_qualities = pysam.qualitystring_to_array(read2.letter_annotations["phred_quality"])

        #write alignments to SAM file
        samfile.write(a1)
        samfile.write(a2)
