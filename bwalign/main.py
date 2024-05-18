# This file is the main file for the project. It will be used to run the project and test the different algorithms.
import utils
def main():
    # Read the input file
    with open(r"C:\Users\Nabil\Desktop\School\Spring2024\CSE185 - Bioinfo Lab\project\bwalign\data\test.fastq") as f:
        fastq = f.readlines()

    # Parse the FASTQ file
    sequences = utils.parse_fastq(fastq)
    
    # Print the first sequence
    print(sequences[0])


    # Read the reference genome file
    with open(r"C:\Users\Nabil\Desktop\School\Spring2024\CSE185 - Bioinfo Lab\project\bwalign\data\test_reference.fa") as f:
        reference_genome = f.readlines()

    # Parse the reference genome
    reference_sequence = utils.parse_reference_genome(reference_genome)

    # Print the reference sequence ID and sequence
    print(reference_sequence)

if __name__ == "__main__":
    main()
