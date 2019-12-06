
# coding: utf-8

# In[1]:


from collections import defaultdict
import sqlite3
import numpy as np


# In[2]:


# INPUTS
query_file = "random_query.fa"  # random query selected as a slice from seq: seq[100:400] len=300
probability_file = "chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
fasta_file = "chr22.maf.ancestors.42000000.complete.boreo.fa.txt"

# PARAMS
w = 11
MATCH_SCORE = 1
MISMATCH_SCORE = -1
M = np.array((4, 4))
delta = 10
alpha = 0.8*w
beta = 100


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

# Open the query file and store it as a string
with open(query_file) as f:
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


def nw(S,T,gap_penalty,match_score,mismatch_score):
    # call it nw because needleman-wunsch is kind of hard to spell
    
    len_S = len(S)
    len_T = len(T)
    
    dp = [[0 for j in range(len_S+1)] for i in range(len_T+1)]
    backpointers = [[None for i in range(len_S+1)] for j in range(len_T+1)]
    
    # init start scores
    
    for j in range(1,len_S):
        dp[0][j] = j * gap_penalty
        backpointers[0][j] = [0,j-1]
    
    for i in range(1,len_T):
        dp[i][0] = i * gap_penalty
        backpointers[i][0] = [i-1,0]
    
    
    # run needleman winsch
    
    for i in range(1,len_T+1):
        for j in range(1,len_S+1):
            
            # going diagonal, need to offset by 1 to get index 0
            cur_best_score = dp[i-1][j-1] + (match_score if T[i-1] == S[j-1] else mismatch_score)
            cur_backpointer = [i-1,j-1]
            
            # going down, so along T
            tmp_score = dp[i-1][j] + gap_penalty
            if tmp_score > cur_best_score:
                cur_best_score = tmp_score
                cur_backpointer = [i-1,j]
            
            # going right, so along S
            tmp_score = dp[i][j-1] + gap_penalty
            if tmp_score > cur_best_score:
                cur_best_score = tmp_score
                cur_backpointer = [i,j-1]
            
            dp[i][j] = cur_best_score
            backpointers[i][j] = cur_backpointer
    
    # best score should be at [len_T][len_S]
    
    # traceback
    back_ptr = [len_T,len_S]
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
    
    print(T_final)
    print(S_final)


# In[11]:


def ungappedExtensionRight(query_index, db_index, seed_score):
    """Takes the index of the query and db at the end the seed and the seed_score
    outputs the indices of the ungapped extension and its score"""
    max_score = seed_score
    maxscoring_qi = 0
    maxscoring_dbi = 0
    score = seed_score
    query_index += 1
    db_index += 1
    
    # While loop that exits when the difference between max_score acheived and score is greater than delta
    while max_score - score < delta and query_index < len(query) and db_index < len(db):
        score += singleBaseCompare(query[query_index], db[db_index])
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
    maxscoring_qi = 0
    maxscoring_dbi = 0
    score = seed_score
    
    # While loop that exits when the difference between max_score acheived and score is greater than delta
    while max_score - score < delta and query_index >= 0 and db_index >= 0:
        query_index -= 1
        db_index -= 1
        score += singleBaseCompare(query[query_index], db[db_index])
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
        
        print("!> \t"+word + " at pos " + str(pos[0]) + " stopped before \t gapped Extension \t (HSP score: "+str(HSP_score)+")")
        HSPs += 1
            
    qi += 1
    
    
print("total_seeds: ", total_seeds)  
print("stopped before ungapped extension: ", stop_before_ungap)
print("stopped before gapped extension: ", stop_before_gap)
print("HSPs: ", HSPs)


# In[55]:


len(query)


# In[ ]:


conn.close()

