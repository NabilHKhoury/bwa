import pysam
from Bio import SeqIO
import numpy as np
from compute_max_seed import compute_max_seed
from generate_seeds import (burrows_wheeler_transform,
                            partial_suffix_array,
                            compute_first_occurrences,
                            compute_checkpoint_arrs,
                            compute_rank_arr,
                            generate_seeds)

def main():
    reference = SeqIO.read("reference.fasta", "fasta")
    reads = list(SeqIO.parse("reads.fastq", "fastq"))

#BWT and SA
    K = 5 
    ref_text = str(reference.seq) + '$'
    bwt = burrows_wheeler_transform(ref_text)
    psa = partial_suffix_array(ref_text, K)
    first_occurrences = compute_first_occurrences(bwt)
    checkpoint_arrs = compute_checkpoint_arrs(bwt)
    ranks = compute_rank_arr(bwt)

    #header for the SAM file (based on ref)
    header = {
        'HD': {'VN': '1.0', 'SO': 'unsorted'},
        'SQ': [{'SN': reference.id, 'LN': len(reference)}]
    }

    #write to sam FIle
    with pysam.AlignmentFile("output.sam", "w", header=header) as samfile:
        for read in reads:
            seed_idxes = generate_seeds(str(read.seq), bwt, 8, psa, first_occurrences, checkpoint_arrs, ranks)
            best_idx, score, alignment_s, alignment_t = compute_max_seed(str(reference.seq), str(read.seq), seed_idxes, 2, -1, -2, -1)

            #create a SAM entry for the aligned read
            a = pysam.AlignedSegment()
            a.query_name = read.id
            a.query_sequence = str(read.seq)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = best_idx
            a.mapping_quality = 60
            a.cigarstring = f"{len(alignment_s.replace('-', ''))}M"
            a.query_qualities = pysam.qualitystring_to_array(read.letter_annotations["phred_quality"])
            
            samfile.write(a)

if __name__ == '__main__':
    main()
