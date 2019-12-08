
# coding: utf-8

# In[1]:


from collections import defaultdict
import sqlite3
import numpy as np


# In[2]:


# INPUTS
query_file = "test_sequence.fa"  # random query selected as a slice from seq: seq[100:400] len=300
probability_file = "chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
fasta_file = "chr22.maf.ancestors.42000000.complete.boreo.fa.txt"

# PARAMS
w = 11
MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_PENALITY = -2
M = np.array((4, 4))
delta = 10      # ungapped extension
alpha = 0.8*w   # Seed stop
beta = 50       # Ungapped extension stop
epsilon = 20    # How far extra to look in the DB for gapped extension with NW


# In[3]:


# Open the probability file and store it is as a list of floats
with open(probability_file, 'r') as f:
    prob = f.readline()
prob = prob.split(' ')
prob = prob[:-1]
prob = [float(i) for i in prob]

# Open the fasta file and store it as a string
with open(fasta_file, 'r') as f:
    db = f.readline()
    if db[0] == '>':
        db = f.readline()

# Open the query file and store it as a string
with open(query_file) as f:
    query = f.readline()
    if query[0] == '>':
        query = f.readline()
    
# Test import
print(db[:10])
print(prob[:10])
print(query[:10])


# In[4]:


conn = sqlite3.connect('blast.db')
c = conn.cursor()


# In[5]:


def get_table_name(word_size):
    return "preprocess_wordsize_" + str(word_size)


# In[6]:


nucleotide_num = {'A':0,'C':1,'G':2,'T':3}

def get_word_encoding(word):
    inverted_word = word[::-1]
    n = len(word)
    
    index = 0
    
    for i in range(n):
        index += pow(4,i) * nucleotide_num[inverted_word[i]] 
    
    return index + 1


# In[7]:


def get_indexes_for_word(word):
    word_size = (len(word))
    new_table_name = get_table_name(word_size)
    
    encoding = get_word_encoding(word)
    
    s = "SELECT sequence_index FROM {} where word_encoding = ?".format(new_table_name)
    return c.execute(s,(encoding,))


# In[8]:


def singleBaseCompare(base1, base2):
    if base1 == base2:
        return MATCH_SCORE
    else:
        return MISMATCH_SCORE


# In[9]:


def scoreSeed(index):
    """ Given the starting index in the database of a seed
    returns the score of that seed based on the product of MATCH_SCORE and probabilities"""
    seed_score = 0
    for i in range(index, index + w + 1):
        seed_score += prob[i] * MATCH_SCORE
    return seed_score


# In[10]:


# Get all possible words from query and store them in the list of strings words
words = []
i = 0
while(i+w <= len(query)):
    words.append(query[i:i+w])
    i += 1

words[:3]


# In[ ]:


def nw(S,T,gap_penalty,match_score,mismatch_score, db_index):
    # call it nw because needleman-wunsch is kind of hard to spell
    
    len_S = len(S)
    len_T = len(T)
    
    dp = [[0 for j in range(len_S+1)] for i in range(len_T+1)]
    backpointers = [[None for i in range(len_S+1)] for j in range(len_T+1)]
    
    # init start scores
    
    for j in range(1,len_S + 1):
        dp[0][j] = j * gap_penalty
        backpointers[0][j] = [0,j-1]
    
    for i in range(1,len_T + 1):
        dp[i][0] = i * gap_penalty
        backpointers[i][0] = [i-1,0]
    
    
    # run needleman winsch
    
    for i in range(1,len_T+1):
        for j in range(1,len_S+1):
            
            # going diagonal, need to offset by 1 to get index 0
            cur_best_score = dp[i-1][j-1] + (prob[db_index] * (match_score if T[i-1] == S[j-1] else mismatch_score))
            cur_backpointer = [i-1,j-1]
            
            # going down, so along T
            tmp_score = dp[i-1][j] + (gap_penalty * prob[db_index])
            if tmp_score > cur_best_score:
                cur_best_score = tmp_score
                cur_backpointer = [i-1,j]
            
            # going right, so along S
            tmp_score = dp[i][j-1] + (gap_penalty * prob[db_index])
            if tmp_score > cur_best_score:
                cur_best_score = tmp_score
                cur_backpointer = [i,j-1]
            
            dp[i][j] = cur_best_score
            backpointers[i][j] = cur_backpointer
        db_index += 1
    
    # best score should be at [len_T][len_S]
    score = dp[len_T][len_S]

    # traceback
    back_ptr = [len_T, len_S]
    rev_sol_T = list()
    rev_sol_S = list()
    
    while not (back_ptr[0] == 0 and back_ptr[1] == 0):
        
        next_back_ptr = backpointers[back_ptr[0]][back_ptr[1]]
        
        # going diagonal
        if next_back_ptr[0] == back_ptr[0] - 1 and next_back_ptr[1] == back_ptr[1] - 1:
            rev_sol_T.append(T[back_ptr[0]-1])
            rev_sol_S.append(S[back_ptr[1]-1])
            
            
        
        # going up, so along T. Gap in S
        elif next_back_ptr[0] == back_ptr[0] - 1 and next_back_ptr[1] == back_ptr[1]:
            rev_sol_T.append(T[back_ptr[0]-1])
            rev_sol_S.append('-')
            
        
        # going down, so along S. Gap in T
        elif next_back_ptr[0] == back_ptr[0] and next_back_ptr[1] == back_ptr[1] - 1:
            rev_sol_T.append('-')
            rev_sol_S.append(S[back_ptr[1]-1])
            
        back_ptr = next_back_ptr
    
    sol_T = rev_sol_T[::-1]
    sol_S = rev_sol_S[::-1]
    
    T_final = "".join(sol_T)
    S_final = "".join(sol_S)
    
    return (S_final, T_final, score)


# In[11]:


def ungappedExtensionRight(query_index, db_index, seed_score):
    """Takes the index of the query and db at the end the seed and the seed_score
    outputs the indices of the ungapped extension and its score"""
    max_score = seed_score
    maxscoring_qi = query_index
    maxscoring_dbi = db_index
    score = seed_score
    query_index += 1
    db_index += 1
    
    # While loop that exits when the difference between max_score acheived and score is greater than delta
    while max_score - score < delta and query_index < len(query) and db_index < len(db):
        score += singleBaseCompare(query[query_index], db[db_index]) * prob[db_index]
        if score > max_score:
            max_score = score
            maxscoring_qi = query_index
            maxscoring_dbi = db_index
        query_index += 1
        db_index += 1
        
    
    return (maxscoring_qi, maxscoring_dbi, max_score)

def ungappedExtensionLeft(query_index, db_index, seed_score):
    """Takes the index of the query and db at the start the seed and the seed_score
    outputs the indices of the ungapped extension and its score"""
    max_score = seed_score
    maxscoring_qi = query_index
    maxscoring_dbi = db_index
    score = seed_score
    
    # While loop that exits when the difference between max_score acheived and score is greater than delta
    while max_score - score < delta and query_index > 0 and db_index > 0:
        query_index -= 1
        db_index -= 1
        score += singleBaseCompare(query[query_index], db[db_index]) * prob[db_index] 
        if score > max_score:
            max_score = score
            maxscoring_qi = query_index
            maxscoring_dbi = db_index
    
    return (maxscoring_qi, maxscoring_dbi, max_score)


# In[14]:


# Putting it all together
qi = 0
total_seeds = 0
stop_before_ungap = 0
stop_before_gap = 0
HSPs = 0

alignments = []
for word in words:
    cursor = get_indexes_for_word(word)
    pos_list = cursor.fetchall()
    
    for pos in pos_list:
        total_seeds += 1
        
        # Score Seed
        seed_score = scoreSeed(pos[0])
        if seed_score < alpha:
            stop_before_ungap += 1
            # print(word + " at pos " + str(pos[0]) + " stopped before \t ungapped Extension \t (seed score: "+str(seed_score)+" )")
            continue
        
        # ungapped Extension
        right = ungappedExtensionRight(qi+w, pos[0]+w, seed_score)
        left = ungappedExtensionLeft(qi, pos[0], seed_score)
        HSP_score = right[2] + left[2] + seed_score
        if HSP_score < beta:
            stop_before_gap += 1
            # print(word + " at pos " + str(pos[0]) + " stopped before \t gapped Extension \t (HSP score: "+str(HSP_score)+")")
            continue
        
        # print("!> \t"+word + " at pos " + str(pos[0]) + " stopped before \t gapped Extension \t (HSP score: "+str(HSP_score)+")")
        HSPs += 1

        # Needleman Winch
        # Right hand side
        # print("left: ",left, "\tright: ", right)
        right_nw = None
        left_nw = None
        db_start_aln = left[1]
        db_end_aln = right[1]
        if len(query) > right[0] + 1:
            right_ext = query[right[0] + 1:]
            if len(right_ext) > 0:
                if len(db) > right[1] + 1 + len(right_ext) + epsilon:
                    right_nw = nw(right_ext, db[right[1] + 1: right[1] + 1 + len(right_ext) + epsilon], GAP_PENALITY, MATCH_SCORE, MISMATCH_SCORE, right[1] + 1)
                    db_end_aln = right[1] + 1 + len(right_ext) + epsilon
                else:
                    right_nw = nw(right_ext, db[right[1] + 1:], GAP_PENALITY, MATCH_SCORE, MISMATCH_SCORE, right[1] + 1)
                    db_end_aln = len(db) - 1
        # Left hand side
        left_ext = query[:left[0]]
        if len(left_ext) > 0:
            if left[1] - len(left_ext) - epsilon >= 0:
                left_nw = nw(left_ext, db[left[1] - len(left_ext) - epsilon:left[1]], GAP_PENALITY, MATCH_SCORE, MISMATCH_SCORE, left[1] - len(left_ext) - epsilon)
                db_start_aln = left[1] - len(left_ext) - epsilon
            else:
                left_nw = nw(left_ext, db[:left[1]], GAP_PENALITY, MATCH_SCORE, MISMATCH_SCORE, 0)
                db_start_aln = 0

        # Print out alignment
        print("Alignment to database at index: ", db_start_aln, " - ", db_end_aln)
        if left_nw is not None and right_nw is not None:
            final_score = HSP_score + left_nw[2] + right_nw[2]
            query_aln = "Query: \t" + left_nw[0] + "HSP: \t" + query[left[0]: right[0] + 1] + "gap: \t" + right_nw[0]
            db_aln = "DB: \t" + left_nw[1]+ "HSP: \t"+ db[left[1]: right[1] + 1]+ "gap: \t" + right_nw[1]
        elif right_nw is not None:
            final_score = HSP_score + right_nw[2]
            query_aln = "Query: \t" + "HSP: \t" + query[left[0]: right[0] + 1] +  "gap: \t" + right_nw[0]
            db_aln = "DB: \t" + "HSP: \t" + db[left[1]: right[1] + 1] + "gap: \t" + right_nw[1]
        elif left_nw is not None:
            final_score = HSP_score + left_nw[2]
            query_aln = "Query: \t" + left_nw[0] + "HSP: \t" + query[left[0]: right[0] + 1]
            db_aln = "DB: \t" + left_nw[1] + "HSP: \t" + db[left[1]: right[1] + 1]
        else:
            final_score = HSP_score
            query_aln = "Query: \t" + "HSP: \t" + query[left[0]: right[0] + 1]
            db_aln = "DB: \t" + "HSP: \t" + db[left[1]: right[1] + 1]
        print("Score: \t", final_score)
        print('\n')

        # Append alignment to matrix
        alignments.append([db_start_aln, db_end_aln, query_aln, db_aln, final_score])

    qi += 1
    
    
print("total_seeds: ", total_seeds)  
print("stopped before ungapped extension: ", stop_before_ungap)
print("stopped before gapped extension: ", stop_before_gap)
print("HSPs: ", HSPs)


# Find the highest Scoring Alignment
def sort_by_score(entry):
    return entry[-1]

sorted_alignments = sorted(alignments, key=sort_by_score)
best_scoring_alignment = sorted_alignments[-1]
print("BEST SCORING ALIGNMENT: ")

for i in best_scoring_alignment:
    print(i)




# In[ ]:


conn.close()

