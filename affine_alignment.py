import numpy as np
from blosum import BLOSUM
# Define scoring matrix and penalties
blosum62 = BLOSUM(62)
gap_open_penalty = -11
gap_extend_penalty = -1


def get_blosum62_score(aa1, aa2):
    """Fetches the BLOSUM62 score for two amino acids."""
    return blosum62[aa1][aa2]


def affine_alignment(seq1, seq2):
    len_x, len_y = len(seq1) + 1, len(seq2) + 1

    # Initialize matrices
    match_mismatch = np.full((len_x, len_y), float('-inf'))
    gap_x_matrix = np.full((len_x, len_y), float('-inf'))
    gap_y_matrix = np.full((len_x, len_y), float('-inf'))

    # Initialize traceback matrices
    traceback_M = np.full((len_x, len_y), None, dtype='object')
    traceback_Ix = np.full((len_x, len_y), None, dtype='object')
    traceback_Iy = np.full((len_x, len_y), None, dtype='object')

    # Initialization of the matrices
    match_mismatch[0][0] = 0

    for i in range(1, len_x):
        gap_x_matrix[i][0] = gap_open_penalty + gap_extend_penalty * (i - 1)
        match_mismatch[i][0] = gap_x_matrix[i][0]
        traceback_Ix[i][0] = ('Ix', i - 1, 0)

    for j in range(1, len_y):
        gap_y_matrix[0][j] = gap_open_penalty + gap_extend_penalty * (j - 1)
        match_mismatch[0][j] = gap_y_matrix[0][j]
        traceback_Iy[0][j] = ('Iy', 0, j - 1)

    # Fill in the matrices
    for i in range(1, len_x):
        for j in range(1, len_y):
            aa1, aa2 = seq1[i - 1], seq2[j - 1]
            score = get_blosum62_score(aa1, aa2)

            # Update match/mismatch matrix
            match_score = max(
                match_mismatch[i - 1][j - 1] + score,
                gap_x_matrix[i - 1][j - 1] + score,
                gap_y_matrix[i - 1][j - 1] + score
            )
            match_mismatch[i][j] = match_score

            if match_score == match_mismatch[i - 1][j - 1] + score:
                traceback_M[i][j] = ('M', i - 1, j - 1)
            elif match_score == gap_x_matrix[i - 1][j - 1] + score:
                traceback_M[i][j] = ('Ix', i - 1, j - 1)
            else:
                traceback_M[i][j] = ('Iy', i - 1, j - 1)

            # Update gap_x_matrix (vertical gap)
            gap_x_matrix[i][j] = max(
                match_mismatch[i - 1][j] + gap_open_penalty + gap_extend_penalty,
                gap_x_matrix[i - 1][j] + gap_extend_penalty
            )
            if gap_x_matrix[i][j] == match_mismatch[i - 1][j] + gap_open_penalty + gap_extend_penalty:
                traceback_Ix[i][j] = ('M', i - 1, j)
            else:
                traceback_Ix[i][j] = ('Ix', i - 1, j)

            # Update gap_y_matrix (horizontal gap)
            gap_y_matrix[i][j] = max(
                match_mismatch[i][j - 1] + gap_open_penalty + gap_extend_penalty,
                gap_y_matrix[i][j - 1] + gap_extend_penalty
            )
            if gap_y_matrix[i][j] == match_mismatch[i][j - 1] + gap_open_penalty + gap_extend_penalty:
                traceback_Iy[i][j] = ('M', i, j - 1)
            else:
                traceback_Iy[i][j] = ('Iy', i, j - 1)

    # Determine the maximum score
    end_scores = [
        match_mismatch[len_x - 1][len_y - 1],
        gap_x_matrix[len_x - 1][len_y - 1],
        gap_y_matrix[len_x - 1][len_y - 1]
    ]
    max_score = max(end_scores)
    if max_score == match_mismatch[len_x - 1][len_y - 1]:
        matrix = 'M'
    elif max_score == gap_x_matrix[len_x - 1][len_y - 1]:
        matrix = 'Ix'
    else:
        matrix = 'Iy'

    # Traceback
    aligned_seq1, aligned_seq2 = [], []
    i, j = len_x - 1, len_y - 1
    score = 0
    gap_y = False
    gap_x = False
    while i > 0 or j > 0:

        if matrix == 'M':
            prev_matrix, prev_i, prev_j = traceback_M[i][j]
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            score += get_blosum62_score(seq1[i-1],seq2[j-1])
            gap_y = False
            gap_x = False
        elif matrix == 'Ix':
            prev_matrix, prev_i, prev_j = traceback_Ix[i][j]
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            if(gap_x):
                score += gap_extend_penalty
            else:
                score += gap_open_penalty

            gap_x = True
            gap_y = False
        elif matrix == 'Iy':
            prev_matrix, prev_i, prev_j = traceback_Iy[i][j]
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            if(gap_y):
                score += gap_extend_penalty
            else:
                score += gap_open_penalty
            gap_x = False
            gap_y = True
        else:
            break

        i, j, matrix = prev_i, prev_j, prev_matrix

    # Reverse the sequences
    aligned_seq1 = ''.join(aligned_seq1[::-1])
    aligned_seq2 = ''.join(aligned_seq2[::-1])

    return aligned_seq1, aligned_seq2, int(score)


# Read input sequences in FASTA format
fasta1 = input()
seq1 = input().strip()
fasta2 = input()
seq2 = input().strip()

# Get the alignment result
aligned_seq1, aligned_seq2, alignment_score = affine_alignment(seq1, seq2)

# Print the output
print(alignment_score)
print(aligned_seq1)
print(aligned_seq2)
