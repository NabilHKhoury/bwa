from .utils import *
import pysam
import argparse

def main():
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
    
    #write to sam FILE
    with pysam.AlignmentFile("output.sam", "w", header=header) as samfile:
        for read_id, read_seq, qual_scores in reads:
            
            seed_idxes = generate_seeds(str(read_seq), bwt, 19, psa, first_occurrences, checkpoint_arrs, ranks)
            best_idx, score, alignment_s, alignment_t = compute_max_seed(str(reference), str(read_seq), seed_idxes, 2, 1, 2, 1, read_length)

            #create a SAM entry for the aligned read
            a = pysam.AlignedSegment()
            a.query_name = read_id
            a.query_sequence = read_seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = best_idx
            a.mapping_quality = calculate_mapping_quality(score, read_length)
            a.cigarstring = calculate_cigar(alignment_s, alignment_t)
            a.query_qualities = qual_scores
            
            samfile.write(a)

if __name__ == "__main__":
    main()