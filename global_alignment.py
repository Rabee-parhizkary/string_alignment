import numpy as np
from blosum import BLOSUM

# Define scoring matrix and penalties
blosum62 = BLOSUM(62)
gap_open_penalty = -5
gap_extend_penalty = 0

def get_blosum62_score(aa1, aa2):
    """Fetches the BLOSUM62 score for two amino acids."""
    return blosum62[aa1][aa2]

def global_alignment(seq1, seq2):
    len_x, len_y = len(seq1) + 1, len(seq2) + 1

    # Initialize scoring and traceback matrices
    score_matrix = np.zeros((len_x, len_y))
    gap_x_matrix = np.full((len_x, len_y), float('-inf'))  # Gaps in seq1
    gap_y_matrix = np.full((len_x, len_y), float('-inf'))  # Gaps in seq2
    traceback_matrix = np.zeros((len_x, len_y), dtype='object')

    # Initialize matrices with gap penalties
    for i in range(1, len_x):
        score_matrix[i][0] = gap_open_penalty + (i - 1) * gap_extend_penalty
        traceback_matrix[i][0] = 'up'
        gap_x_matrix[i][0] = gap_open_penalty + (i - 1) * gap_extend_penalty

    for j in range(1, len_y):
        score_matrix[0][j] = gap_open_penalty + (j - 1) * gap_extend_penalty
        traceback_matrix[0][j] = 'left'
        gap_y_matrix[0][j] = gap_open_penalty + (j - 1) * gap_extend_penalty

    # Fill the scoring matrices with affine gap penalties
    for i in range(1, len_x):
        for j in range(1, len_y):
            aa1, aa2 = seq1[i - 1], seq2[j - 1]

            # Calculate match/mismatch score
            match_score = get_blosum62_score(aa1, aa2)
            diagonal = score_matrix[i - 1][j - 1] + match_score

            # Calculate gap penalties
            gap_x_matrix[i][j] = max(score_matrix[i - 1][j] + gap_open_penalty, gap_x_matrix[i - 1][j] + gap_extend_penalty)
            gap_y_matrix[i][j] = max(score_matrix[i][j - 1] + gap_open_penalty, gap_y_matrix[i][j - 1] + gap_extend_penalty)

            # Update score matrix with the maximum score among diagonal and gaps
            max_score = max(diagonal, gap_x_matrix[i][j], gap_y_matrix[i][j])
            score_matrix[i][j] = max_score

            # Traceback direction
            if max_score == diagonal:
                traceback_matrix[i][j] = 'diagonal'
            elif max_score == gap_x_matrix[i][j]:
                traceback_matrix[i][j] = 'up'
            elif max_score == gap_y_matrix[i][j]:
                traceback_matrix[i][j] = 'left'

    # Traceback to build the aligned sequences
    aligned_seq1, aligned_seq2 = [], []
    i, j = len_x - 1, len_y - 1
    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 'diagonal':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'up':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback_matrix[i][j] == 'left':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse the sequences to get the correct alignment
    aligned_seq1 = ''.join(aligned_seq1[::-1])
    aligned_seq2 = ''.join(aligned_seq2[::-1])

    return aligned_seq1, aligned_seq2, score_matrix[-1][-1]  # Return alignment and final score

fasta1 = input()
seq1 = input()
fasta2 = input()
seq2 = input()

aligned_seq1, aligned_seq2, alignment_score = global_alignment(seq1, seq2)
print(int(alignment_score))
