import numpy as np

# INPUTS
QUERY = ""
probabilities = "/Users/jonbroad/Google Drive/School/McGill/Fall 2019/Comp 561/Project/Chr22probabilities.txt"

# PARAMS
w = 11
MATCH_SCORE = 1
M = np.array((4, 4))
delta = 10
alpha = 5
beta = 10

with open(probabilities, 'r') as f:
    D_prob = f.readline()
D_prob = D_prob.split(' ')
D_prob = D_prob[:-1]
D_prob = [float(i) for i in D_prob]


def singleBaseCompare(seq1, seq2, i, j):
    if seq1[i] == seq2[j]:
        return 2
    else:
        return -1

def scoreSeed(seed):
    index = 0
    seed_score = 0
    for i in range(index: index + w):
        seed_score = D_prob[i] * MATCH_SCORE

    return seed_score

