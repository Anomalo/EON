from pprint import pprint
from threading import Thread
modules  = [	'os',
		'plot',
		'operator',
		'sys',
		'cPickle as pickle',
		'glob',
		'pprint',
		'go',
		'fa',
		'gtf',
		'csv',
		]

for module in modules:
	#print 'importing',module
	exec('import '+module)

def DtoTSV(d,tsvName,sep='\t',extras=None,head=None,header=None):
	'''
	it saves a dictionary into a csv file
	'''
	out = []
	if not header == None:out.append(header)
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


def _sepWords(seq, WordSize):
	wordsNumber = len(seq)-WordSize
	words = []
	for i in range(wordsNumber):
		word = seq[i:i+WordSize]
		words.append(word)
	return words


class glast:
	def __init__(self,gmt,v=False ,output='results',ws=11,distortion = 3):
		'''
		if no gib specified will chose one from annotations/ folder.
		if no gib is present in that directory it will prompt to create one.
		'''
		self.output = output
		self.v = v
		self.WS = ws
		self.exonSeqOfGO={}
		self.gmt=gmt
		self.distortion=distortion
	def exonSeqs(self, genes,v=False):
		'''returns a a list of tuples of [(exonname,seq),...]'''
		exonsOut=[]
		for i,gene in enumerate(genes,start=1):
			exons = gtf.transcriptNames(gene+'-001')
			genePercent = round(100*(float(i)/len(genes)),2)
			for j, exon in enumerate(exons,start=1):
				exonPercent = round(100*(float(j)/len(exons)),2)
				if v:
					print genePercent,'%',exonPercent,'%',exon,' '*20
					sys.stdout.write("\033[F")
				c, s, e,strand = gtf.getTranscriptCoords(exon)
				seq = fa.seq_coords(c, s, e,strand)
				exonsOut.append((exon,seq))
		return exonsOut
	def glastSeq(self, 
			  seq,
			  exon,
			  ws = 11,
			  loops = 100, 
			  fname='GOS.tsv',
			  b=False,
			  header=None):
		'''
		given a sequence it returns glasting results via scanning a gob file
		'''
		exonQuery=exon
		wordsList=_sepWords(seq,ws)
		virtual_seqs={}
		qseq=seq
		v=self.v
		exon_go={}
		import go
		GO = go.GO(self.gmt)
		GOlabel_ID={}
		geneName = exon.split('-')[0]
		GOlabels = GO.GOgeneNames(geneName)
		goSizes={}
		for golabelIndex, golabel in enumerate(GOlabels):
			goSizes[golabel]=0
			if not golabel in self.exonSeqOfGO:
				genes = GO.GenesWithGO(golabel)	
				GOlabel_ID[golabel]=GO.GOlabel_ID(golabel)
				seqs = self.exonSeqs(genes,v=v)
				self.exonSeqOfGO[golabel]=seqs
			else: seqs = self.exonSeqOfGO[golabel]
			for seq in seqs:
				exon,seq = seq
				goSizes[golabel]+=len(seq)
				if exon == exonQuery: continue
				if v:
					os.system('clear')
					print 'Glasting', header,
					print 'for', golabel,
					print round(100*(golabelIndex/float(len(GOlabels))),1),'%'
					print exon,'(%(length)s)'%{'length':len(seq)} 
				if not exon in exon_go: exon_go.update({exon:[]})
				exon_go[exon].append(golabel)
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
		#DtoTSV(grades,'blast.tsv',extras = exon_go, head = ['exon',
		#						    'score',
		#						    'gos'])
		##make a GOS (goScore) list {go:sum(score),...}
		GOS={}
		for exon, score in grades.iteritems():
			goIDS = exon_go[exon]
			for go in goIDS:
				if not go in GOS: GOS.update({go:0})
				GOS[go] += score
			#normalizes by exon size
			for go,score in GOS.iteritems():
				GOS[go]=score/len(qseq)
		GOS_sorted = sorted(GOS.items(), key=operator.itemgetter(1))[::-1]
		if v:pprint.pprint(GOS_sorted)
		
		DtoTSV(GOS,fname,header=header,extras=GOlabel_ID)
		

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


	def glastExon(self, exon ,dir=None,header=''):
		'''
		given a transcript name it returns the glasting results
		'''
		import gtf
		v=self.v
		chromosome,start,end,strand = gtf.getTranscriptCoords(exon)
		seq = fa.seq_coords(chromosome,start,end,strand)
		header = '%(exon)s\n%(chromosome)s:%(start)s-%(end)s (%(strand)s)'%locals()
		if dir ==  None:dir = 'GOS_'+exon.split('-')[0]
		if not os.path.exists(dir): os.makedirs(dir)
		fname = '/'.join([dir,exon+'.tsv'])
		#fname = '%(dir)s/%(exon)s.tsv'%{dir:dir,exon:exon}
		gene = exon.split('-')[0].upper()
		GO = go.GO(self.gmt)
		guide = GO.GOgeneNames(gene)
		if v:print header
		return self.glastSeq(seq, exon, 
					fname=fname,header=header,ws=self.WS)
			
	def glastGene(self, gene,exonsToBlast=[]):
		'''
		given a gene short name, it returns the glasting
		results of all transcripts in a dictionary format
		'''
		v=self.v
		if not '-' in gene:gene+='-001'
		if v:print 'glasting', gene, 'exons:',','.join(exonsToBlast)
		output=self.output
		print 'importing gtf'
		import gtf
		exons = gtf.transcriptNames(gene)
		results = dict.fromkeys(exons)
		dir = output+'/'+gene
		if not os.path.exists(dir): os.makedirs(dir)		
		#if v:print exons
		for exon in exons:
			if len(exonsToBlast) !=0:
				if not exon.split(':')[-1] in exonsToBlast:
					if v:print 'skipping',exon
					continue
			if v:print exon
			
			results[exon] = self.glastExon(exon,dir=dir,
							 )
			
		plot.plotDIR(dir,distortion = self.distortion)
		return results
	


def main():
	pass
	
if __name__ == '__main__':
	main()

