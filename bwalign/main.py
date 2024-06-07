from .utils import *
import pysam
import argparse
from tqdm import tqdm

"""
    File:          main.py
    Description:   Executes entire bwalign workflow

    Names:         Adrian Layer, Nabil Khoury, Yasmin A. Jaber
    Emails:        alayer@ucsd.edu, nkhoury@ucsd.edu, yjaber@ucsd.edu
    Project:       bwalign (final project for CSE 185 @ UC San Diego)
    Repository:    https://github.com/NabilHKhoury/bwalign
"""

def main():
    """
    Workflow can be broken down into three distinct parts:
    
    1. Parsing of fasta and fastq files, and storage in auxiliary data structures.
    2. Using Burrows-Wheeler/suffix array algorithms to seed reads, then using
       banded Needleman-Wunsch global alignment to extend those reads and determine
       best matches.
    3. Writing these best matches into SAM file format.
    """
    ### PART 1 OF WORKFLOW - PARSING FILES
    
    parser = argparse.ArgumentParser(description="Run BWA alignment with specified reference genome and FASTQ file.")
    parser.add_argument("reference_genome", type=str, help="Path to the reference genome file")
    parser.add_argument("fastq_file", type=str, help="Path to the FASTQ file")
    args = parser.parse_args()
    
    # Parse the FASTQ file
    reads = parse_fastq(args.fastq_file)
    _, first_read, _ = reads[0]
    read_length = len(first_read)
    
    # Parse the reference genome
    ref_id, reference = parse_reference_genome(args.reference_genome)
    
    ### PART 2 OF WORKFLOW - SEED AND EXTEND
    
    #BWT and SA
    K = 5
    ref_text = str(reference) + '$'
    bwt, psa = bwt_psa_out(ref_text, K)
    first_occurrences = compute_first_occurrences(bwt)
    checkpoint_arrs = compute_checkpoint_arrs(bwt)
    ranks = compute_rank_arr(bwt)
    
    #header for the SAM file (based on ref)
    header = {
        'HD': {'VN': '1.0', 'SO': 'unsorted'},
        'SQ': [{'SN': ref_id, 'LN': len(reference)}]
    }
    
    ### PART 3 OF WORKFLOW - WRITE TO SAM FILE
    
    #write to sam FILE
    with pysam.AlignmentFile("output.sam", "w", header=header) as samfile:
        for read_id, read_seq, qual_scores in tqdm(reads):
            
            seed_idxes = generate_seeds(str(read_seq), bwt, 19, psa, first_occurrences, checkpoint_arrs, ranks)
            best_idx, score, alignment_s, alignment_t = compute_max_seed(str(reference), str(read_seq), seed_idxes, 2, 2, 2, 10, read_length)
            
            #create a SAM entry for the aligned read
            a = pysam.AlignedSegment()
            a.query_name = read_id
            a.query_sequence = read_seq
            a.flag = 0
            a.reference_id = 0
            if score != float('-inf'):
                a.reference_start = best_idx
                a.mapping_quality = 60
                a.cigarstring = calculate_cigar(alignment_s, alignment_t)
            else: # in the case read could not be aligned
                a.reference_start = 0
                a.mapping_quality = 0
                a.cigarstring = '0M'
            a.query_qualities = qual_scores
            
            samfile.write(a)

if __name__ == "__main__":
    main()