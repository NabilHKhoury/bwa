from Bio import SeqIO
from typing import List, Tuple

def parse_fastq(fastq_path: str) -> List[Tuple[str, str, List[int]]]:
    """
    Parse a FASTQ file and return a list of tuples where each tuple contains the
    sequence ID, the sequence itself, and the quality scores.
    
    :param fastq_path: Path to the FASTQ file
    :return: A list of tuples where each tuple contains the sequence ID, the sequence itself, and the quality scores
    """
    fastq_list = []
    with open(fastq_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            fastq_list.append((record.id, str(record.seq), record.letter_annotations["phred_quality"]))
    return fastq_list

def parse_reference_genome(fasta_path: str) -> Tuple[str, str]:
    """
    Parse a FASTA file containing a reference genome and return a tuple containing
    the sequence ID and the sequence itself.

    :param fasta_path: Path to the FASTA file
    :return: A tuple containing the sequence ID and the sequence itself
    """
    with open(fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return record.id, str(record.seq)

def random_sequence(length: int) -> str:
    """
    Generate a random DNA sequence of the specified length.

    :param length: The length of the sequence
    :return: A random DNA sequence
    """
    import random
    return ''.join(random.choices('ACGT', k=length))

def random_quality_scores(length: int) -> List[int]:
    """
    Generate random quality scores for a sequence of the specified length.

    :param length: The length of the sequence
    :return: A list of random quality scores
    """
    import random
    return [random.randint(0, 40) for _ in range(length)]

def random_id(length: int) -> str:
    """
    Generate a random sequence ID of the specified length.

    :param length: The length of the sequence ID
    :return: A random sequence ID
    """
    import random
    return ''.join(random.choices('ACGT', k=length))

def generate_fasta_file(fasta_path: str, sequence_id: str, sequence: str):
    """
    Generate a FASTA file with the specified sequence ID and sequence.

    :param fasta_path: Path to the FASTA file to be generated
    :param sequence_id: The sequence ID
    :param sequence: The sequence
    """
    with open(fasta_path, "w") as handle:
        handle.write(f">{sequence_id}\n")
        handle.write(sequence)

def generate_fastq_file(fastq_path: str, sequence_id: str, sequence: str, quality_scores: List[int]):
    """
    Generate a FASTQ file with the specified sequence ID, sequence, and quality scores.

    :param fastq_path: Path to the FASTQ file to be generated
    :param sequence_id: The sequence ID
    :param sequence: The sequence
    :param quality_scores: The quality scores
    """
    with open(fastq_path, "w") as handle:
        handle.write(f"@{sequence_id}\n")
        handle.write(f"{sequence}\n")
        handle.write("+\n")
        handle.write("".join(chr(score + 33) for score in quality_scores))
def main():
    # Define parameters for random sequence and quality scores
    sequence_length = 100
    sequence_id_length = 10

    # Generate random reference genome
    ref_sequence_id = random_id(sequence_id_length)
    ref_sequence = random_sequence(sequence_length)
    generate_fasta_file(r"C:\Users\Nabil\Desktop\School\Spring2024\CSE185 - Bioinfo Lab\project\bwalign\data\random_reference.fasta", ref_sequence_id, ref_sequence)

    # Generate random FASTQ reads
    num_reads = 10
    for i in range(num_reads):
        read_id = random_id(sequence_id_length)
        read_sequence = random_sequence(sequence_length)
        read_quality = random_quality_scores(sequence_length)
        generate_fastq_file(f"random_read_{i+1}.fastq", read_id, read_sequence, read_quality)

if  __name__ == '__main__':
    main()

