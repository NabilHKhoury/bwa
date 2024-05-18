from typing import List, Tuple, Iterable

def parse_fastq(fastq: Iterable[str]) -> List[Tuple[str, str]]:
    """
    Parse a FASTQ file and return a list of tuples where each tuple contains the
    sequence ID, the sequence itself, and the quality scores.

    :param fastq: An iterable containing the lines of a FASTQ file
    :return: A list of tuples where each tuple contains the sequence ID, the sequence itself, and the quality scores

    Example:
    >>> fastq = ['@seq1', 'ATCG', '+', 'ABCD'] # seq1, ATCG, ABCD are the sequence ID, sequence, and quality scores respectively in this example
    >>> parse_fastq(fastq)
    [('seq1', 'ATCG', 'ABCD')]
    """
    fastq_list = []
    for i in range(0, len(fastq), 4):
        sequence_id = fastq[i][1:].strip()
        sequence = fastq[i+1].strip()
        quality_scores = fastq[i+3].strip()
        fastq_list.append((sequence_id, sequence, quality_scores))
    return fastq_list


def parse_reference_genome(fasta: Iterable[str]) -> Tuple[str, str]:
    """
    Parse a FASTA file containing a reference genome and return a tuple containing
    the sequence ID and the sequence itself.

    :param fasta: An iterable containing the lines of a FASTA file
    :return: A tuple containing the sequence ID and the sequence itself

    Example:
    >>> fasta = ['>ref_genome', 'ATCG'] # ref_genome and ATCG are the sequence ID and sequence respectively in this example
    >>> parse_reference_genome(fasta)
    ('ref_genome', 'ATCG')
    """
    sequence_id = fasta[0][1:].strip()
    sequence = fasta[1].strip()
    return sequence_id, sequence