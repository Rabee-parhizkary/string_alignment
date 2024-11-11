import numpy as np
def max_gap(seq1, seq2,g,m,d):
    gap = g
    len_x, len_y = len(seq1) + 1, len(seq2) + 1
    # Initialize scoring and traceback matrices
    score_matrix = np.zeros((len_x, len_y))
    traceback_matrix = np.zeros((len_x, len_y), dtype='object')

    # Initialize matrices with gap penalties
    for i in range(1, len_x):
        score_matrix[i][0] = i*gap
        traceback_matrix[i][0] = 'up'

    for j in range(1, len_y):
        score_matrix[0][j] = j*gap
        traceback_matrix[0][j] = 'left'
    # Fill the scoring matrices with affine gap penalties
    for i in range(1, len_x):
        for j in range(1, len_y):
            aa1, aa2 = seq1[i - 1], seq2[j - 1]

            # Calculate match/mismatch score
            match_score = (m if aa1==aa2 else d)
            diagonal = score_matrix[i - 1][j - 1] + match_score

            # Calculate gap penalties
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap

            # Update score matrix with the maximum score among diagonal and gaps
            max_score = max(diagonal, up,left)
            score_matrix[i][j] = max_score

            # Traceback direction
            # We prioritize creating gap
            if max_score == up:
                traceback_matrix[i][j] = 'up'
            elif max_score == left:
                traceback_matrix[i][j] = 'left'
            if max_score == diagonal:
                traceback_matrix[i][j] = 'diagonal'
    gap_count = 0
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
            gap_count +=1
        elif traceback_matrix[i][j] == 'left':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
            gap_count += 1

    # Reverse the sequences to get the correct alignment
    aligned_seq1 = ''.join(aligned_seq1[::-1])
    aligned_seq2 = ''.join(aligned_seq2[::-1])

    return gap_count # Return alignment and final score

fasta1 = input()
seq1 = input().strip()
fasta2 = input()
seq2 = input().strip()

#The number of maximum gaps would be the maximum number out of three bellow numbers: it depends on the gap penalty
#compared to the mismatch penalty
gap_count1 = max_gap(seq1, seq2,-10,2,-4)
gap_count2 = max_gap(seq1, seq2,-3,2,-3)
gap_count3 = max_gap(seq1, seq2,-4,2,-10)
gap_count = max(gap_count1,gap_count2,gap_count3)
print(int(gap_count))
