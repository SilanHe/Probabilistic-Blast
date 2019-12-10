from collections import defaultdict
import sqlite3
import numpy as np

class PBlast:

	# CLASS VARIABLES (STATIC)
	nucleotide_num = {'A':0,'C':1,'G':2,'T':3}

	# class constructors
	def __init__(self):

		self.query_file = "test_sequence.fa"  # random query selected as a slice from seq: seq[100:400] len=300
		self.probability_file = "chr22.maf.ancestors.42000000.complete.boreo.conf.txt"
		self.fasta_file = "chr22.maf.ancestors.42000000.complete.boreo.fa.txt"

		self.__init_params__()
		self.__init_files__()
		self.connect_db()

	def __init_params__(self):

		# DEFAULT PARAMS
		self.w = 11 # word length
		self.MATCH_SCORE = 1
		self.MISMATCH_SCORE = -1
		self.GAP_PENALITY = -1
		self.delta = 10      # ungapped extension max decrease
		self.alpha = 0.8 * self.w   # Seed stop
		self.beta = 50       # Ungapped extension stop
		self.epsilon = 20    # How far extra to look in the DB for gapped extension with NW

	# our loading function for the files
	def __init_files__(self):
		# Open the probability file and store it is as a list of floats
		with open(self.probability_file, 'r') as f:
			prob = f.readline()
		prob = prob.split(' ')
		prob = prob[:-1]
		prob = [float(i) for i in prob]

		self.prob = prob

		# Open the fasta file and store it as a string
		with open(self.fasta_file, 'r') as f:
			db = f.readline()
			if db[0] == '>':
				db = f.readline()

		self.db = db

		# Open the query file and store it as a string
		with open(self.query_file) as f:
			query = f.readline()
			if query[0] == '>':
				query = f.readline()

		self.query = query

	def set_nw_params(self, match_score,mismatch_score,gap_penalty):
		self.MATCH_SCORE = match_score
		self.MISMATCH_SCORE = mismatch_score
		self.GAP_PENALITY = gap_penalty

	def set_blast_params(self, max_decrease_stop = 10, seed_stop = 0.8, ungapped_stop = 50, gapped_stop = 20):
		self.delta = max_decrease_stop
		self.alpha = seed_stop * self.w
		self.beta = ungapped_stop
		self.epsilon = gapped_stop

	def set_word_length(self,length):
		self.w = length

	def connect_db(self):
		self.conn = sqlite3.connect('blast.db')
		self.c = self.conn.cursor()

	def exit(self):
		self.conn.close()

	# load a query file into the global instance variable
	def load_query_file(self,query_filename):
		# Open the query file and store it as a string
		with open(query_filename) as f:
			query = f.readline()
			if query[0] == '>':
				query = f.readline()

		self.query = query
		self.query_file = query_filename

	def get_table_name(self,word_size):
		return "preprocess_wordsize_" + str(word_size)

	def get_word_encoding(self,word):
		inverted_word = word[::-1]
		n = len(word)
		
		index = 0
		
		for i in range(n):
			index += pow(4,i) * self.nucleotide_num[inverted_word[i]] 
	
		return index + 1

	def get_indexes_for_word(self,word):
		word_size = (len(word))
		new_table_name = self.get_table_name(word_size)
		
		encoding = self.get_word_encoding(word)
		
		s = "SELECT sequence_index FROM {} where word_encoding = ?".format(new_table_name)
		return self.c.execute(s,(encoding,))

	def singleBaseCompare(self,base1, base2):
		if base1 == base2:
			return self.MATCH_SCORE
		else:
			return self.MISMATCH_SCORE

	def scoreSeed(self,index):
		""" Given the starting index in the database of a seed
		returns the score of that seed based on the product of self.MATCH_SCORE and probabilities"""
		seed_score = 0
		for i in range(index, index + self.w + 1):
			seed_score += self.prob[i] * self.MATCH_SCORE
		return seed_score

	def nw(self,S,T,gap_penalty,match_score,mismatch_score, db_index):
		# call it self.nw because needleman-wunsch is kind of hard to spell
		
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
				cur_best_score = dp[i-1][j-1] + (self.prob[db_index] * (match_score if T[i-1] == S[j-1] else mismatch_score))
				cur_backpointer = [i-1,j-1]
				
				# going down, so along T
				tmp_score = dp[i-1][j] + gap_penalty * self.prob[db_index]
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

	def ungappedExtensionRight(self,query_index, db_index, seed_score):
		"""Takes the index of the self.query and self.db at the end the seed and the seed_score
		outputs the indices of the ungapped extension and its score"""
		max_score = seed_score
		maxscoring_qi = query_index
		maxscoring_dbi = db_index
		score = seed_score
		query_index += 1
		db_index += 1
		
		# While loop that exits when the difference between max_score acheived and score is greater than self.delta
		while max_score - score < self.delta and query_index < len(self.query) and db_index < len(self.db):
			score += self.singleBaseCompare(self.query[query_index], self.db[db_index]) * self.prob[db_index]
			if score > max_score:
				max_score = score
				maxscoring_qi = query_index
				maxscoring_dbi = db_index
			query_index += 1
			db_index += 1
			
		
		return (maxscoring_qi, maxscoring_dbi, max_score)

	def ungappedExtensionLeft(self,query_index, db_index, seed_score):
		"""Takes the index of the self.query and self.db at the start the seed and the seed_score
		outputs the indices of the ungapped extension and its score"""
		max_score = seed_score
		maxscoring_qi = query_index
		maxscoring_dbi = db_index
		score = seed_score
		
		# While loop that exits when the difference between max_score acheived and score is greater than self.delta
		while max_score - score < self.delta and query_index > 0 and db_index > 0:
			query_index -= 1
			db_index -= 1
			score += self.singleBaseCompare(self.query[query_index], self.db[db_index]) * self.prob[db_index] 
			if score > max_score:
				max_score = score
				maxscoring_qi = query_index
				maxscoring_dbi = db_index
		
		return (maxscoring_qi, maxscoring_dbi, max_score)

	def probabilistic_blast(self):

		words = []
		i = 0
		while(i+self.w <= len(self.query)):
			words.append(self.query[i:i+self.w])
			i += 1

		# Putting it all together
		qi = 0
		total_seeds = 0
		stop_before_ungap = 0
		stop_before_gap = 0
		HSPs = 0

		alignments = []
		for word in words:
			cursor = self.get_indexes_for_word(word)
			pos_list = cursor.fetchall()
			
			for pos in pos_list:
				total_seeds += 1
				
				# Score Seed
				seed_score = self.scoreSeed(pos[0])
				if seed_score < self.alpha:
					stop_before_ungap += 1
					# print(word + " at pos " + str(pos[0]) + " stopped before \t ungapped Extension \t (seed score: "+str(seed_score)+" )")
					continue
				
				# ungapped Extension
				right = self.ungappedExtensionRight(qi+self.w, pos[0]+self.w, seed_score)
				left = self.ungappedExtensionLeft(qi, pos[0], seed_score)
				HSP_score = right[2] + left[2] + seed_score
				if HSP_score < self.beta:
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
				if len(self.query) > right[0] + 1:
					right_ext = self.query[right[0] + 1:]
					if len(right_ext) > 0:
						if len(self.db) > right[1] + 1 + len(right_ext) + self.epsilon:
							right_nw = self.nw(right_ext, self.db[right[1] + 1: right[1] + 1 + len(right_ext) + self.epsilon], self.GAP_PENALITY, self.MATCH_SCORE, self.MISMATCH_SCORE, right[1] + 1)
							db_end_aln = right[1] + 1 + len(right_ext) + self.epsilon
						else:
							right_nw = self.nw(right_ext, self.db[right[1] + 1:], self.GAP_PENALITY, self.MATCH_SCORE, self.MISMATCH_SCORE, right[1] + 1)
							db_end_aln = len(self.db) - 1
				# Left hand side
				left_ext = self.query[:left[0]]
				if len(left_ext) > 0:
					if left[1] - len(left_ext) - self.epsilon >= 0:
						left_nw = self.nw(left_ext, self.db[left[1] - len(left_ext) - self.epsilon:left[1]], self.GAP_PENALITY, self.MATCH_SCORE, self.MISMATCH_SCORE, left[1] - len(left_ext) - self.epsilon)
						db_start_aln = left[1] - len(left_ext) - self.epsilon
					else:
						left_nw = self.nw(left_ext, self.db[:left[1]], self.GAP_PENALITY, self.MATCH_SCORE, self.MISMATCH_SCORE, 0)
						db_start_aln = 0

				# Print out alignment
				#print("Alignment to database at index: ", db_start_aln, " - ", db_end_aln)
				if left_nw is not None and right_nw is not None:
					final_score = HSP_score + left_nw[2] + right_nw[2]
					query_aln = "\nquery: \t" + left_nw[0] + "\nHSP: \t" + self.query[left[0]: right[0] + 1] + "\ngap: \t" + right_nw[0]
					db_aln = "\ndb: \t" + left_nw[1]+ "\nHSP: \t"+ self.db[left[1]: right[1] + 1]+ "\ngap: \t" + right_nw[1]
				elif right_nw is not None:
					final_score = HSP_score + right_nw[2]
					query_aln = "\nquery: \t" + "\nHSP: \t" + self.query[left[0]: right[0] + 1] +  "\ngap: \t" + right_nw[0]
					db_aln = "db: \t" + "\nHSP: \t" + self.db[left[1]: right[1] + 1] + "\ngap: \t" + right_nw[1]
				elif left_nw is not None:
					final_score = HSP_score + left_nw[2]
					query_aln = "\nquery: \t" + left_nw[0] + "\nHSP: \t" + self.query[left[0]: right[0] + 1]
					db_aln = "\ndb: \t" + left_nw[1] + "\nHSP: \t" + self.db[left[1]: right[1] + 1]
				else:
					final_score = HSP_score
					query_aln = "\nquery: \t" + "\nHSP: \t" + self.query[left[0]: right[0] + 1]
					db_aln = "\ndb: \t" + "\nHSP: \t" + self.db[left[1]: right[1] + 1]
				#print("Score: \t", final_score)
				#print('\n')

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

		return sorted_alignments
