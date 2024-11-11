import numpy as np

# Define linear gap penalty
gap_penalty = 5

# PAM250 scoring matrix
def scoringMatrix():
    sMatrixTxt = '''
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  2 -2  0  0 -3  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3
C -2 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0
D  0 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4
E  0 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4
F -3 -4 -6 -5  9 -5 -2  1 -5  2  0 -3 -5 -5 -4 -3 -3 -1  0  7
G  1 -3  1  0 -5  5 -2 -3 -2 -4 -3  0  0 -1 -3  1  0 -1 -7 -5
H -1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0
I -1 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1
K -1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4
L -2 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1
M -1 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2
N  0 -4  2  1 -3  0  2 -2  1 -3 -2  2  0  1  0  1  0 -2 -4 -2
P  1 -3 -1 -1 -5  0  0 -2 -1 -3 -2  0  6  0  0  1  0 -1 -6 -5
Q  0 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4
R -2 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4
S  1  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3
T  1 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3
V  0 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2
W -6 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0
Y -3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10
'''
    sMatrixList = sMatrixTxt.strip().split('\n')
    aaList = sMatrixList[0].split()
    sMatrix = {aa: {} for aa in aaList}
    for line in sMatrixList[1:]:
        parts = line.split()
        for i, aa in enumerate(aaList):
            sMatrix[parts[0]][aa] = int(parts[i + 1])
    return sMatrix

PAM250 = scoringMatrix()

def local_alignment(s, t):
    m, n = len(s), len(t)
    H = np.zeros((m + 1, n + 1))
    max_score = 0
    max_pos = (0, 0)

    # Fill the scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = PAM250[s[i - 1]][t[j - 1]]
            score = max(
                0,
                H[i - 1][j - 1] + match,
                H[i - 1][j] - gap_penalty,
                H[i][j - 1] - gap_penalty
            )
            H[i][j] = score

            if score > max_score:
                max_score = score
                max_pos = (i, j)

    # Traceback
    i, j = max_pos
    max_i, max_j = i, j
    aligned_s = []
    aligned_t = []

    while i > 0 and j > 0 and H[i][j] > 0:
        if H[i][j] == H[i - 1][j - 1] + PAM250[s[i - 1]][t[j - 1]]:
            aligned_s.append(s[i - 1])
            aligned_t.append(t[j - 1])
            i -= 1
            j -= 1
        elif H[i][j] == H[i - 1][j] - gap_penalty:
            aligned_s.append(s[i - 1])
            aligned_t.append('-')
            i -= 1
        elif H[i][j] == H[i][j - 1] - gap_penalty:
            aligned_s.append('-')
            aligned_t.append(t[j - 1])
            j -= 1
        else:
            break  # Should not happen in Smith-Waterman

    # The starting positions of the substrings
    start_i, start_j = i, j

    # Extract the substrings r and u
    r = s[start_i:max_i]
    u = t[start_j:max_j]

    return r, u, int(max_score)
fasta1 = input()
seq1 = input().strip()
fasta2 = input()
seq2 = input().strip()

aligned_seq1, aligned_seq2, alignment_score = local_alignment(seq1, seq2)
print(alignment_score)
print(aligned_seq1)
print(aligned_seq2)
