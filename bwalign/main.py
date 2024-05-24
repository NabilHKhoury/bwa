import utils
import argparse

def main(reference_genome_path, fastq_file_path):
    # Parse the FASTQ file
    sequences = utils.parse_fastq(fastq_file_path)
    
    # Print the first sequence
    print("Sequences from FASTQ:")
    for seq in sequences:
        print(seq)
    # Parse the reference genome
    reference_sequence = utils.parse_reference_genome(reference_genome_path)

    # Print the reference sequence ID and sequence
    print("Reference genome sequence:")
    print(reference_sequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BWA alignment with specified reference genome and FASTQ file.")
    parser.add_argument("reference_genome", type=str, help="Path to the reference genome file")
    parser.add_argument("fastq_file", type=str, help="Path to the FASTQ file")

    args = parser.parse_args()
    main(args.reference_genome, args.fastq_file)
