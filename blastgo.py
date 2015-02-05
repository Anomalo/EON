from pprint import pprint
from threading import Thread
modules  = [	'os',
		'eta',
		'operator',
		'sys',
		'cPickle as pickle',
		'glob',
		'pprint',
		'fa',
		'csv',
		]
for module in modules:
	#print 'importing',module
	exec('import '+module)

def DtoTSV(d,tsvName,sep='\t',extras=None,head=None):
	'''
	it saves a dictionary into a csv file
	'''
	out = []
	if head != None:
		out =[sep.join(head)]
	for x,y in d.iteritems():
		out.append(sep.join([str(x),str(y)]))
		if extras != None:
			try:
				out[-1]=out[-1]+sep+str(extras[x])
			except: 
				pass
	tsv = '\n'.join(out)
	f = open(tsvName,'w')
	f.write(tsv)
	f.close()


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
	def __init__(self,GIB=None, GOB=None):
		'''
		if no gib specified will chose one from annotations/ folder.
		if no gib is present in that directory it will prompt to create one.	
		'''
		if GOB == None:
			GOB = glob.glob('annotations/*.gob')[0]
		self.GOB = GOB
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
		'''
	def _glastSeqIndex(self, seq, b = False):
		'''
		given a sequence it returns glasting results via GIB index
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
			
	def _glastSeqScan(self, seq, b = False, ws = 11,loops = 100):
		'''
		given a sequence it returns glasting results via scanning a gob file
		'''
		wordsList=_sepWords(seq,ws)
		GOB=open(self.GOB)
		virtual_seqs={}
		qseq=seq
		num_lines = file_len(self.GOB)
		if b: num_lines = loops
		ETA = eta.ETA(num_lines)
		ETA.print_status()
		exon_go={}
		while True:	
			ETA.print_status()
	
			GOBline = GOB.readline()
			if GOBline == '':break
			attr = GOBline.split()
			if len(attr)==2:attr.append('None:'+'n'*ws)
			go, url, seqs = attr
			print go
			#print go, url
			seqs = seqs.split(';')
			for seq in seqs:
				exon,seq = seq.split(':')
				if not exon in exon_go: exon_go.update({exon:[]})
				exon_go[exon].append(go)
				scanWords = _sepWords(seq,ws)
				n=0
				for word in scanWords:
					if word in wordsList:
						if not exon in virtual_seqs:
							virtual_exon = [-1]*len(scanWords)
							virtual_seqs.update({exon:virtual_exon})
						virtual_seqs[exon][n]=wordsList.index(word)
					n+=1				
			if b:
				loops-=1
				if loops == 0 :break


		grades = self.gradeMatchesD(virtual_seqs,qseq,ws)
		
		filter = [k for k, v in grades.items() if v>0.05]
		DtoTSV(grades,'blast.tsv',extras = exon_go, head = ['exon',
								    'score',
								    'gos'])
		##make a GOS (goScore) list {go:sum(score),...}
		GOS={}
		for exon, score in grades.iteritems():
			goIDS = exon_go[exon]
			for go in goIDS:
				if not go in GOS: GOS.update({go:0})
				GOS[go] += score
		GOS_sorted = sorted(GOS.items(), key=operator.itemgetter(1))[::-1]
		pprint.pprint(GOS_sorted)
	
		DtoTSV(GOS,'GOS.tsv')
		

	def gradeMatches(self,matches,original,ws):
		'''	
		given a list of matches it retuns a value of similarity
		'''
		origWords = len(original)-ws+1
		size = min([origWords,len(matches)])
		sequenciality = ws-1
		for i in range(1,len(matches)):
			if matches[i-1]+1==matches[i]:
				sequenciality+=1.0
		return sequenciality/size

	def gradeMatchesD(self, matchD, original, ws):
		'''
		given a dict of exon:mathces it grades all the mathces
		'''
		out = dict.fromkeys(matchD.keys())
		for exon in matchD:
			out[exon]=self.gradeMatches(matchD[exon],original,ws)
		return out


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
	def glastSeq(self, seq, b=False, scan = True, ws = 11,loops=100):
		if scan: return self._glastSeqScan(seq=seq, b=b,ws=ws,loops=loops)
		else: return self._glastSeqIndex(seq=seq, b=b)


if __name__ == '__main__':
	#genGOB('annotations/Mus_musculus_GSEA_GO_sets_all_symbols_September_2013.gmt')
	#genGIB('annotations/Mus_musculus_GSEA_GO_sets_all_symbols_September_2013.gob',b=False,checkpoint = 500) 
	g =  glast()
	seq ='''ACCATGGATCTCTCTGCCATCTACGAGGTGAGTACCTGTTAGACAGCATCCCGGGATCCCCGACGCACCAAACTTAGGCCC'''
	print g.glastSeq(seq)#,b=True,loops=100)
	#print g.glastExon('Malat1-001',b=True)#[:3]
#

