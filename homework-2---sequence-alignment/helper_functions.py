from Bio.Align import substitution_matrices
from typing import Generator

def blosum62_score(el1, el2):
    if 'blosum62' not in globals():
        global blosum62
        blosum62 = substitution_matrices.load("BLOSUM62")
    if el1 == '-':
        el1 = '*'
    if el2 == '-':
        el2 = '*'
    return blosum62[el1][el2]


def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    """
    # Setup sequences
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    # Get lengths
    len1 = len(seq1)
    len2 = len(seq2)
    # Prepare scores and backtrace matrices
    score_mat = [[0 for i in range(len2)] for j in range(len1)]
    backtrace_mat = [[' ' for i in range(len2)] for j in range(len1)]
    # Starting position
    backtrace_mat[0][0] = 'S'
    # Setup first row and column with appropriate scores and directions for both matrices
    for i in range(1, len1):
        score_mat[i][0] = score_mat[i - 1][0] + scoring_function(seq1[i], '-')
        backtrace_mat[i][0] = 'V'
    for i in range(1, len2):
        score_mat[0][i] = score_mat[0][i - 1] + scoring_function('-', seq2[i])
        backtrace_mat[0][i] = 'H'
    # Calculate remaining fields for both matrices
    scores = {}
    for i in range(1, len1):
        for j in range(1, len2):
            scores['D'] = score_mat[i - 1][j - 1] + scoring_function(seq1[i], seq2[j])
            scores['H'] = score_mat[i][j - 1] + scoring_function('-', seq2[j])
            scores['V'] = score_mat[i - 1][j] + scoring_function(seq1[i], '-')
            score_mat[i][j] = max(scores.values())
            backtrace_mat[i][j] = max(scores, key=scores.get)
    # Backtrace and rebuild best aligned sequences
    i = len1 - 1
    j = len2 - 1
    aligned_seq1 = ""
    aligned_seq2 = ""
    direction = backtrace_mat[i][j]
    while direction != 'S':
        if direction == 'D':  # we came to this cell diagonally
            aligned_seq1 = seq1[i] + aligned_seq1
            aligned_seq2 = seq2[j] + aligned_seq2
            i -= 1
            j -= 1
        elif direction == 'H':  # we came to this cell horizontally
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j] + aligned_seq2
            j -= 1
        elif direction == 'V':  # we came to this cell vertically
            aligned_seq1 = seq1[i] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        direction = backtrace_mat[i][j]
    return aligned_seq1, aligned_seq2, score_mat[len1 - 1][len2 - 1]



def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    """
    # Setup sequences
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    # Get lengths
    len1 = len(seq1)
    len2 = len(seq2)
    # Prepare scores and backtrace matrices
    score_mat = [[0 for i in range(len2)] for j in range(len1)]
    backtrace_mat = [[' ' for i in range(len2)] for j in range(len1)]
    # Starting position
    backtrace_mat[0][0] = 'S'
    # Setup first row and column with appropriate directions for backtrace matrix
    for i in range(1, len1):
        backtrace_mat[i][0] = 'V'
    for i in range(1, len2):
        backtrace_mat[0][i] = 'H'
    # Calculate remaining fields for both matrices
    scores = {}
    max_score = 0
    max_loc = (0, 0)
    for i in range(1, len1):
        for j in range(1, len2):
            scores['D'] = score_mat[i - 1][j - 1] + scoring_function(seq1[i], seq2[j])
            scores['H'] = score_mat[i][j - 1] + scoring_function('-', seq2[j])
            scores['V'] = score_mat[i - 1][j] + scoring_function(seq1[i], '-')
            scores['S'] = 0
            score_mat[i][j] = max(scores.values())
            if score_mat[i][j] > max_score:
                max_score = score_mat[i][j]
                max_loc = (i, j)
            backtrace_mat[i][j] = max(scores, key=scores.get)
    # Backtrace and rebuild best aligned sequence
    i, j = max_loc
    aligned_seq1 = ""
    aligned_seq2 = ""
    direction = backtrace_mat[i][j]
    while direction != 'S':
        if direction == 'D':  # we came to this cell diagonally
            aligned_seq1 = seq1[i] + aligned_seq1
            aligned_seq2 = seq2[j] + aligned_seq2
            i -= 1
            j -= 1
        elif direction == 'H':  # we came to this cell horizontally
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j] + aligned_seq2
            j -= 1
        elif direction == 'V':  # we came to this cell vertically
            aligned_seq1 = seq1[i] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        direction = backtrace_mat[i][j]
    return aligned_seq1, aligned_seq2, max_score


def codons(seq: str) -> Generator[str, None, None]:
    """Walk along the string, three nucleotides at a time. Cut off excess."""
    for i in range(0, len(seq) - 2, 3):
        yield seq[i:i + 3]


def translate_to_protein(seq):
    """Translate a nucleotide sequence into a protein sequence.

    Parameters
    ----------
    seq: str

    Returns
    -------
    str
        The translated protein sequence.

    """
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W',
    }
    protein_sequence = ''
    for codon in codons(seq):
        protein_sequence += table[codon]
    return protein_sequence
