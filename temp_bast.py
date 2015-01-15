import random
from pprint import pprint

class blast:
	def __init__(self,genome,word_size=11):
		self.genome = genome
		self.word_size = word_size
		self.make_index()
	def make_index(self):
		ws = self.word_size
		index = {}
		for gene, seq in self.genome.iteritems():
			for i in range(len(seq)-ws):
				word = seq[i:i+ws]
				if not word in index:
					index.update({word:[]})
				index[word].append((gene,i))
		self.index = index
	
	def geneSize(self,gene):
		return len(self.genome[gene])

	def score(self,blastList,size):
		onlyMatches = blastList[:]
		while None in onlyMatches: onlyMatches.remove(None)
		score = 0
		for i in range(len(onlyMatches)-1):
			if onlyMatches[i]<onlyMatches[i+1]:    score+=1
			if onlyMatches[i+1]==onlyMatches[i]+1: score+=1
		return score/float(size)

			
			

	def blast(self,query):
		ws = self.word_size
		index = self.index
		similars={}
		#scan through the query 
		for i in range(len(query)-ws):
			word = query[i:i+ws]
			if word in index:
				for gene, location in index[word]:
					if not gene in similars:
						similars.update({gene:[None]*self.geneSize(gene)})
					similars[gene][location]=i
		results = {}
		for gene,blastlist in similars.iteritems():
			size = min([len(query),self.geneSize(gene)])
			score = self.score(blastlist,size)
			results.update({gene:score})
		return results
genome={
	'a':'CACCAGCCGACGATGCTTGATAGAGCAGCGAGAGACGTTGTCACGTGCCGCGCCAAGGCGACACCCATCTGCGATATTCGATTGGACTGG',
	'b':'TCGAGAGGCTCCCACTCAACATGGCGCCCTTCAAACAAACCTGGAGGATCACACAAACTGGGAGTCTCCTTTTCTTCCTTTCTTCACAGG',
	'c':'TCCAGGGATCGGTATCATTAAACGCGGCATTGGGTTTATGATTTCATCTGTGGGGGGTACCAGCTGGGTCTGCAGAAGTCTGACGACTTA',
	'd':'CCCTGTGCCGTATTCAGTACAGATGCCGGCTGCCTCAGGTGGGCCCGCGGACTCGGTCGGACCACAGAGCCGCCTTTGGGATCGATGCTG',
	'e':'CCCTGTGCCGTATTAGTACAGATGCCGGCTGCCTCAGGTGGGCCCGCGGTCTCGGTCGGACCACAGAGCCGCCTTTGGGATCGATGCTG',
	}
b = blast(genome)
q = 'CCCTGTGCCGTATTAGTACAGATGCCGGCTGCCTCAGGTGGGCCCGCGGTCAGAGACGTTGTCAC'
print b.blast(q)

	
