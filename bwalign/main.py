# This file is the main file for the project. It will be used to run the project and test the different algorithms.
import utils
import argparse

def main(reference_genome_path, fastq_file_path):
    # Read the input FASTQ file
    with open(fastq_file_path) as f:
        fastq = f.readlines()

    # Parse the FASTQ file
    sequences = utils.parse_fastq(fastq)
    
    # Print the first sequence
    print(sequences[0])

    # Read the reference genome file
    with open(reference_genome_path) as f:
        reference_genome = f.readlines()

    # Parse the reference genome
    reference_sequence = utils.parse_reference_genome(reference_genome)

    # Print the reference sequence ID and sequence
    print(reference_sequence)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BWA alignment with specified reference genome and FASTQ file.")
    parser.add_argument("reference_genome", type=str, help="Path to the reference genome file")
    parser.add_argument("fastq_file", type=str, help="Path to the FASTQ file")

    args = parser.parse_args()
    main(args.reference_genome, args.fastq_file)
