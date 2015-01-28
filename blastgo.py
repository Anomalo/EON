from pprint import pprint
from threading import Thread
modules  = [	'os',
		'eta',
		'sys',
		'cPickle as pickle',
		'glob',
		'fa',
		]
for module in modules:
	#print 'importing',module
	exec('import '+module)


def file_len(fname):
	'''
	just returns the number of lines of a file.
	'''
	try:
		x = os.popen('wc -l '+fname).read().split()[0]
		return int(x)
	except:
		f = open(fname)
		c = f.read().count('\n')
		return c	


def genGOB(gmt=None):
	'''
	given a gmtfile, it creates a .gob (GOBlast) file
	a row of BGO looks like:
	[go term]	[all exons of all genes with that go concated by ';']
	'''
	import gtf
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
	ETA = eta.ETA(num_lines)
	n=0.0
	for i in range(num_lines):
		ETA.touch_status(prefix='GOB')
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

def genGIB(gob=None, word_size=11, b=False, checkpoint=1000,load=True):
	'''
	makes a GIB (gob index), from a gob file. if no gob is stated it will look for one in 
	annotations/ folder. if no gob present there it will promt to create one.
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
	print 'Gibbing'
	num_lines = file_len(gob)
	gibfname = gob[:-4]+'_ws_'+str(word_size)+'.gib'
	gob = open(gob)
	#GIB = {}

	if load: mode = 'c'
	else:mode  = 'n'
	#GIB = shelve.open(gibfname,mode,writeback=True)
	if load:
		try:
			GIB = pickle.load(open(gibfname,'rb'))
		except:
			load = False
	if not load: GIB={'position':1,'exon_sizes':{}}
	position=GIB['position']-1
	position = 0
	n=0.0
	ETA = eta.ETA(num_lines)
	def _save(target,content):
		pickle.dump(content,open(target,'wb'))
		#print 'saved'
	for l in range(position,num_lines):
		try: ETA.touch_status(prefix='GIB')
		except: pass
		line = gob.readline()
		go, url, seqs = line.split('\t')
		goid = url.split(':')[-1]
		seqs = seqs.split(';')
		for seqIndex in range(len(seqs)):
			if not  seqs[seqIndex].replace('\n',"") == '':
				exon, seq = seqs[seqIndex].replace('\n',"").split(':')
				exonSize = len(seq)

				if not exon in GIB['exon_sizes']:
					exon_sizes = GIB['exon_sizes']
					exon_sizes.update({exon:exonSize})
					GIB['exon_sizes']=exon_sizes
				words = len(seq)-word_size
				for i in range(words):
					word = seq[i:i+word_size]
					if not word in GIB: GIB.update({word:[]})
					if not word in GIB:print word in GIB
					oldList = GIB[word]
					newEntry = {	'go':go,
							#'goID':goid,
							#'url':url,
							#'line':l,
							#'seqNumber':seqIndex,
							'location':i,
							'exon':exon,
							#'exonSize':exonSize
							}
					if not newEntry in oldList:
						oldList.append(newEntry)
						GIB[word] = oldList
		n+=1
		if l%checkpoint == 0: # saves the dictionary every checkpoint lines
			GIB['position']=l
			#print l
			a=Thread(target=_save,args=(gibfname,GIB))
			a.start()
			a.join()
				#pickle.dump(GIB,open(gibfname,'wb'))
			#GIB.sync()
		if b:
			print n
			break
	#GIB.close()
	a=Thread(target=_save,args=(gibfname,GIB))
	a.start()
	a.join()

	#pickle.dump(GIB,open(gibfname,'wb'))

	
def _sepWords(seq, WordSize):
	wordsNumber = len(seq)-WordSize
	words = []
	for i in range(wordsNumber):
		word = seq[i:i+WordSize]
		words.append(word)
	return words

class glast:
	def __init__(self,GIB=None):
		'''
		if no gib specified will chose one from annotations/ folder.
		if no gib is present in that directory it will prompt to create one.	
		'''
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
		possibilities = {}
		for word in words:
			if word in index:
				wordMatches = index[word]
				for wordMatch in wordMatches:
					exon = wordMatch['exon']
					if not exon in possibilities:
						exonSize = wordMatch['exonSize']
						url = wordMatch['url']
						go = wordMatch['go']
						goID = wordMatch['goID']
						mockStrand = [None]*(exonSize-self.WS)
						data = {'url':url,'go':go,'goID':goID}
						possibilities.update({exon:{'strand':mockStrand, 'data':wordMatch}})
					location = wordMatch['location']
					possibilities[exon]['strand'][location]=n
					print possibilities
			if b:break
			n+=1
			
	def glastExon(self, exon , b=False):
		'''
		given a transcript name it returns the glasting results
		'''
		import gtf
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
	genGIB('annotations/Mus_musculus_GSEA_GO_sets_all_symbols_September_2013.gob',b=False,checkpoint = 500) 
	#g =  glast()
	#seq ='''ACCTCACTTGAGCCACGAGTGGGGTCAGGCATGTGGGTTTAAAGAGTTTTCCTTTGCAGAGCCTCATTTCATCCTTCATGGAGCTGCTCAGGACTTTGCATATAAGCGCTTGCCTCTGTCTTCTGTTCTGCTAGTGAGTGTGTGATGTGAGACCTTGCAGTGAGTTTGTTTTTCCTGGAATGTGGAGGGAGGGGGGGATGGGGCTTACTTGTTCTAGCTTTTTTTTTACAGACCACACAGAATGCAGGTGTCTTGACTTCAGGTCATGTCTGTTCTTTGGCAAGTAATATGTGCAGTACTGTTCCAATCTGCTGCTATTAGAATGCATTGTGACGCGACTGGAGTATGATTAAAGAAAGTTGTGTTTCCCCAAGTGTTTGGAGTAGTGGTTGTTGGAGGAAAAGCCATGAGTAACAGGCTGAGTGTTGAGGAAATGGCTCTCTGCAGCTTTAAGTAACCCGTGTTTGTGATTGGAGCCGAGTCCCTTTGCTGTGCTGCCTTAGGTAAATGTTTTTGTTCATTTCTGGTGAGGGGGGTTGGGAGCACTGAAGCCTTTAGTCTCTTCCAGATTCAACTTAAAATCTGACAAGAAATAAATCAGACAAGCAACATTCTTGAAGAAATTTTAACTGGCAAGTGGAAATGTTTTGAACAGTTCCGTGGTCTTTAGTGCATTATCTTTGTGTAGGTGTTCTCTCTCCCCTCCCTTGGTCTTAATTCTTACATGCAGGAACATTGACAACAGCAGACATCTATCTATTCAAGGGGCCAGAGAATCCAGACCCAGTAAGGAAAAATAGCCCATTTACTTTAAATCGATAAGTGAAGCAGACATGCCATTTTCAGTGTGGGGATTGGGAAGCCCTAGTTCTTTCAGATGTACTTCAGACTGTAGAAGGAGCTTCCAGTTGAATTGAAATTCACCAGTGGACAAAATGAGGACAACAGGTGAACGAGCCTTTTCTTGTTTAAGATTAGCTACTGGTAATCTAGTGTTGAATCCTCTCCAGCTTCATGCTGGAGCAGCTAGCATGTGATGTAATGTTGGCCTTGGGGTGGAGGGGTGAGGTGGGCGCTAAGCCTTTTTTTAAGATTTTTCAGGTACCCCTCACTAAAGGCACTGAAGGCTTAATGTAGGACAGCGGAGCCTTCCTGTGTGGCAAGAATCAAGCAAGCAGTATTGTATCGAGACCAAAGTGGTATCATGGTCGGTTTTGATTAGCAGTGGGGACTACCCTACCGTAACACCTTGTTGGAATTGAAGCATCCAAAGAAAATACTTGAGAGGCCCTGGGCTTGTTTTAACATCTGGAAAAAAGGCTGTTTTTATAGCAGCGGTTACCAGCCCAAACCTCAAGTTGTGCTTGCAGGGGAGGGAAAAGGGGGAAAGCGGGCAACCAGTTTCCCCAGCTTTTCCAGAATCCTGTTACAAGGTCTCCCCACAAGTGATTTCTCTGCCACATCGCCACCATGGGCCTTTGGCCTAATCACAGACCCTTCACCCCTCACCTTGATGCAGCCAGTAGCTGGATCCTTGAGGTCACGTTGCATATCGGTTTCAAGGTAACCATGGTGCCAAGGTCCTGTGGGTTGCACCAGAAAAGGCCATCAATTTTCCCCTTGCCTGTAATTTAACATTAAAACCATAGCTAAGATGTTTTATACATAGCACCTATGCAGAGTAAACAAACCAGTATGGGTATAGTATGTTTGATACCAGTGCTGGGTGGGAATGTAGGAAGTCGGATGAAAAGCAAGCCTTTGTAGGAAGTTGTTGGGGTGGGATTGCAAAAATTCTCTGCTAAGACTTTTTCAGGTGGACATAACAGACTTGGCCAAGCTAGCATCTTAGTGGAAGCAGATTCGTCAGTAGGGTTGTAAAGGTTTTTCTTTTCCTGAGAAAACAACCTTTTGTTTTCTCAGGTTTTGCTTTTTGGCCTTTCCCTAG'''
	#print g.glastSeq(seq,b=False)
	#print g.glastExon('Malat1-001',b=True)#[:3]
#
