from pprint import pprint
modules  = ['os','subprocess','cPickle as pickle','glob','fa']#,'gtf']
for module in modules:
	print 'importing',module
	exec('import '+module)
'''
import os
import subprocess
import cPickle as pickle
import glob
import fa
import gtf
'''
def file_len(fname):
	'''
	just returns the number of lines of a file.
	'''
	p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
						  stderr=subprocess.PIPE)
	result, err = p.communicate()
	if p.returncode != 0:
		raise IOError(err)
	return int(result.strip().split()[0])

def genGOB(gmt=None):
	'''
	given a gmtfile, it creates a .gob (GOBlast) file
	a row of BGO looks like:
	[go term]	[all exons of all genes with that go concated by ';']
	'''
	if gmt == None:
		gmt = glob.glob('annotations/*.gmt')[0]
		if gmt == '':
			print 'NO gmt present'
			print 'download dataset from http://www.go2msig.org/cgi-bin/prebuilt.cgi'
			quit()
	f = open(gmt)
	num_lines = file_len(gmt)
	newFile = gmt[:-3]+'gob'
	open(newFile, 'w').close()
	newFile = open(newFile,'a')
	#print 1
	genesSeen={}
	n=0.0
	for i in range(num_lines):

		os.system('clear')
		percentage = str(round(100*n/num_lines,2))+'%'
		params = {'step':'GOB','percentage': percentage,'left':num_lines - n}
		remainder ='%(step)s_%(percentage)s_%(left)s_left' 
		remainder = remainder % params
		#print remainder
		os.system('rm *.remainder')
		os.system('touch '+remainder+'.remainder')
		line = f.readline().split()
		GO = line.pop(0)
		url = line.pop(0)
		genes = line
		exons = []
		print GO,';'
		for gene in genes:
			if not gene in genesSeen:
				try:
					exonNames = gtf.transcriptNames(gene)
					for exon in exonNames:
						#print exon
						c,s,e,strand = gtf.getTranscriptCoords(exon)
						exons.append( exon+':'+fa.seq_coords(c,s,e,strand))
					genesSeen.update({gene:exons})
				except:
					genesSeen.update({gene:[]})
					#print 'fail','\t',gene
			#print '\t',gene
			exons = genesSeen[gene]
		exons = ';'.join(exons)
		newLine = GO + '\t' + url + '\t' + exons + '\n'
		newFile.write(newLine)
		n+=1

	os.system('rm *.remainder')
	os.system('touch DONE.remainder')

def genGIB(gob=None, word_size=11, b=False):
	'''
	makes a GIB (gob index), from a gob file.
	'''
	if gob == None:
		gob = glob.glob('annotations/*.gob')[0]
		if gob == '':
			print 'NO GOB present, make one?'
			ans = raw_input('[yes/no]')
			if ans =='yes':
				genGOB(gmt = None)
				gob = glob.glob('annotations/*.gob')[0]
			else: quit()	
	num_lines = file_len(gob)
	gibfname = gob[:-4]+'_ws_'+str(word_size)+'.gib'
	gob = open(gob)
	metaGIB=[]
	n=0.0
	for l in range(num_lines):
		if n%1000 == 0: metaGIB.append({})
		os.system('clear')
		percentage = str(round(100*n/num_lines,2))+'%'
		params = {'step':'GIB','percentage': percentage,'left':num_lines - n}
		remainder ='%(step)s_%(percentage)s_%(left)s_left' 
		remainder = remainder % params
		print remainder
                os.system('rm *.remainder')
                os.system('touch '+remainder+'.remainder')
		line = gob.readline()
		go, url, seqs = line.split('\t')
		goid = url.split(':')[-1]
		seqs = seqs.split(';')
		for seqIndex in range(len(seqs)):
			if not  seqs[seqIndex].replace('\n',"") == '':
				exon, seq = seqs[seqIndex].replace('\n',"").split(':')
				words = len(seq)-word_size
				exonSize = len(seq)
				for i in range(words):
					word = seq[i:i+word_size]
					if not word in metaGIB[-1]: metaGIB[-1].update({word:[]})
					metaGIB[-1][word].append({'go':go,
							'goID':goid,
							'url':url,
							'line':l,
							'seqNumber':seqIndex,
							'location':i,
							'exon':exon,
							'exonSize':exonSize})
		n+=1
		if b:
			print n
			break
	print 'Merging GIBS into one GIB'
	gibD={}
	for d in metaGIB:
		gibD = dict(gibD,**d)		
	pickle.dump(gibD,open(gibfname,'wb'))
	os.system('rm *.remainder')

def _sepWords(seq, WordSize):
	wordsNumber = len(seq)-WordSize
	words = []
	for i in range(wordsNumber):
		word = seq[i:i+WordSize]
		words.append(word)
	return words

class glast:
	def __init__(self,GIB=None):
		if GIB == None:
			GIB = glob.glob('annotations/*.gib')[0]
		if GIB == '':
			print 'NO GIB present, make one?'
			ans = raw_input('[yes/no]')
			if ans =='yes':
				WS = int(raw_input('what word size? '))
				genGIB(gob=none,word_size=WS,b=False)
				GIB = glob.glob('annotations/*.gib')[0]
			else: quit()
		self.index = pickle.load(open(GIB))
		WS = GIB.split('.')[0].split('_')[-1]
		self.WS = int(WS)

	def glastSeq(self, seq, b = False):
		'''
		given a sequence it returns glasting results
		'''
		index = self.index
		words = _sepWords(seq,self.WS)
		matches = []
		n = 0
		for word in words:
			if word in index:
				
				matches.append()						
				
			
	def glastExon(self, exon , b=False):
		'''
		given a transcript name it returns the glasting results
		'''
		print 'glasting',exon
		chromosome,start,end,strand = gtf.getTranscriptCoords(exon)
		seq = fa.seq_coords(chromosome,start,end,strand)
		return self.glastSeq(seq, b = b)
		
	def glastGene(self, gene):
		'''
		given a gene short name, it returns the glasting
		results of all transcripts in a dictionary format
		'''
		exons = gtf.getExons(gene)
		results = dict.fromkeys(exons)
		for exon in exons:
			results[exons] = self.glastExon(exon)
		return results

if __name__ == '__main__':
	#genGOB('annotations/Mus_musculus_GSEA_GO_sets_all_symbols_September_2013.gmt')
	genGIB('annotations/Mus_musculus_GSEA_GO_sets_all_symbols_September_2013.gob',b=False) 
	#g =  glast()
	#print g.glastExon('Malat1-001',b=True)#[:3]
#
