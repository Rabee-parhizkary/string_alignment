# Define scoring parameters
match_score = 1
mismatch_penalty = -2
gap_penalty = -2

def overlap_alignment(seq1, seq2):
    len_x, len_y = len(seq1), len(seq2)

    # Initialize scoring and traceback matrices
    score_matrix = [[0] * (len_y + 1) for _ in range(len_x + 1)]
    traceback_matrix = [[''] * (len_y + 1) for _ in range(len_x + 1)]

    # Fill the scoring matrix with overlap alignment scoring
    for i in range(1, len_x + 1):
        for j in range(1, len_y + 1):
            # Calculate match/mismatch score
            if seq1[i - 1] == seq2[j - 1]:
                match = match_score
            else:
                match = mismatch_penalty

            diagonal = score_matrix[i - 1][j - 1] + match
            up = score_matrix[i - 1][j] + gap_penalty  # Gap in seq2 (insert in seq1)
            left = score_matrix[i][j - 1] + gap_penalty  # Gap in seq1 (insert in seq2)

            # Determine the maximum score with priority rules
            max_score_cell = max(diagonal, left, up)
            if diagonal == max_score_cell:
                pointer = 'diagonal'
            elif left == max_score_cell:
                pointer = 'left'
            else:
                pointer = 'up'

            score_matrix[i][j] = max_score_cell
            traceback_matrix[i][j] = pointer

    # Find the maximum score in the last row and last column
    max_score = float('-inf')
    max_i, max_j = len_x, 0

    # Check the last row for maximum score
    for j in range(len_y + 1):
        if score_matrix[len_x][j] > max_score:
            max_score = score_matrix[len_x][j]
            max_i, max_j = len_x, j
        elif score_matrix[len_x][j] == max_score:
            # Choose the one with the longest alignment
            if j > max_j:
                max_i, max_j = len_x, j

    # Check the last column for maximum score
    for i in range(len_x + 1):
        if score_matrix[i][len_y] > max_score:
            max_score = score_matrix[i][len_y]
            max_i, max_j = i, len_y
        elif score_matrix[i][len_y] == max_score:
            if i > max_i:
                max_i, max_j = i, len_y

    # Traceback to get the aligned sequences
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_i, max_j
    while i > 0 and j > 0:
        pointer = traceback_matrix[i][j]
        if pointer == 'diagonal':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif pointer == 'left':
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        elif pointer == 'up':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            break

    # Reverse the aligned sequences to get the correct order
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    return int(max_score), aligned_seq1, aligned_seq2

# Input DNA sequences
seq1 = input().strip()
seq2 = input().strip()

# Calculate overlap alignment
alignment_score, aligned_seq1, aligned_seq2 = overlap_alignment(seq1, seq2)

# Output results
print(alignment_score)
print(aligned_seq1)
print(aligned_seq2)
