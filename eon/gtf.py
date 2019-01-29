from glob import glob
import urllib
import cPickle as pickle
import os
import os.path

def getGTF(taxon='Mus musculus',release=77, dir='annotations/'):
	''' finds the GTF needed from ensembl and downloads it'''
	
	taxon = taxon.lower().replace(' ','_')
	url = 'ftp://ftp.ensembl.org/pub/release-'+str(release)+'/gtf/'+taxon+'/'
	response = urllib.urlopen(url).read()
	responselines = response.splitlines()
	for line in responselines:
		if '.gtf.gz' in line:
			fname = line.split()[-1]
			URL = url+fname

	os.system('clear')
	os.system(' '.join(['wget',
		URL,
		'-P',
		dir]))
	os.system(' '.join(['gunzip',
		dir+'/'+fname,
		'-d',
		dir]))

class gtf:
	"""
	A class used to represent a GTF file optimized to get fast queries

	...

	Attributes
	----------
	f : str
	a GTF file name
	"""



	def __init__(self, f = None):
		'''
		init object to get fast queries regarding gtf data
		f = gtf file
		'''

		if f == None:
			# if no GTF file given then see if there is only one in the annotation directory and uses it
			
			f = glob('annotations/*gtf')
			if len(f) != 1:
				raise ValueError('there is either no gtf file in annotations or there is more than one')
			f = f[0]
			
		if os.path.isfile(f+'.p'):
			# if used GTF file has a pickle version of self then use self instead of rereading it
			P = pickle.load( open( f+".p", "rb" ) )
			self.bigGTFlistdict = P['bigGTFlistdict']
			self.chrDict 	= P['chrDict']
			self.bigGTFdict	= P['bigGTFdict']
			self.names_transcripts = P['names_transcripts']
			self.gene_IDS	= P['gene_IDS']
			# end __init__ here
			return None
		
		f = open(f,'r')
		self.header = ['seqname',
					   'source',
					   'feature',
					   'start',
					   'end',
					   'score',
					   'strand',
					   'frame',
					   'attribute']
		# turns the gtf into a list
		bigGTFlistdict = []
		n =0
		while True:
			line = f.readline()
			n+=1
			if line == '': break
			if line[0]!='#':
				if 'exon' == line.split('\t')[2]:
					bigGTFlistdict.append(self._splitLine(line, self.header))
		self.bigGTFlistdict = bigGTFlistdict


		#  separates gtf list into a chromosome indexed dictionary
		chrDict = {}
		for line in bigGTFlistdict:
			chromosome = line['seqname']
			if not chromosome in chrDict:
				chrDict.update({chromosome:[]})
			chrDict[chromosome].append(line)
		self.chrDict = chrDict
		# creates the 'bigGTFdict' dictionary which will hold the preread GTF file
		self.bigGTFdict={}
		# creates the 'names_transcripts' dictionary which will hold the exons transcripts of each gene
		self.names_transcripts = {}
		# creates the 'gene_IDS' dictionary which will hold the genes ID of each gene
		self.gene_IDS={}
		for line in self.bigGTFlistdict:
			if not 'transcript_name' in line: continue
			if not line['transcript_name'] in self.bigGTFdict:
				exon_number = line['exon_number']
				exon_number = '0'*(3-len(exon_number))+exon_number
				transcript_name = line['transcript_name'].upper()+':'+exon_number
				self.bigGTFdict.update({transcript_name:line})
			name  = line['transcript_name'].upper()
			if not name in self.names_transcripts:
				self.names_transcripts.update({ name:[]})
			self.names_transcripts[name].append(transcript_name)
			ID = line['gene_id']
			if not ID in self.gene_IDS:
				self.gene_IDS[ID]=[]
			self.gene_IDS[ID].append(transcript_name)
		# saves the dictionaries into a pickle file for faster loading on the next run
		P = {}
		P['bigGTFlistdict']	= self.bigGTFlistdict
		P['chrDict']	= self.chrDict
		P['bigGTFdict']	= self.bigGTFdict
		P['names_transcripts']	= self.names_transcripts
		P['gene_IDS']	= self.gene_IDS
		pickle.dump( P, open( gtfFile+".p", "wb" ) )
		#orders and cleans self.names_transcripts
		#for i in self.names_transcripts:
		#	self.names_transcripts[i] =  sorted(set(self.names_transcripts[i]))
		
	def _splitLine(self,line,header):
		'''
		Given a tab delaminated string, and a header list, it returns a 
		dictionary with headers as keys
		'''
		line = line.replace('"','')
		line = line.split('\t')
		if len(line[0])<3:
			line[0] = 'chr'+line[0]
		for i in range(len(line)):
			try: line[i] = int(line[i])
			except: pass
		line = dict(zip(header, line))
		attr = line['attribute'].split('; ')
		for at in attr:
			attName= at.split(' ')[0]
			attVal = at.replace(attName+' ','')
			line.update({attName: attVal})
		del line['attribute']
		return line
	
	def readCONFIG(self, fname = 'config.txt'):
		'''
		return GTF file acording to config file, a string of the .gmt file name from 
		the config.txt file
		'''
		f = open(fname, 'r')
		exec(f.read())
		f.close()
		GTF = 'annotations/' + GTF
		return GTF

	def getGene(self, chromosome, start, end):
		'''
		given chromosome, start, and end of a gene it returns the 
		gene attributes
		'''
		chromosomeGenes = self.chrDict[chromosome]
		for line in chromosomeGenes:
			lstart = line['start']
			lend = line['end']
			if lstart<=start<=lend or lstart<=end<=lend:
				return line   

	def getExon(self, transcriptName):
		'''
		given a transcript name it returns the GFT data of that transcript in a 
		dictionaty.
		'''
		return self.bigGTFdict[transcriptName.upper()]

	def transcriptNames(self, gene):
		'''
		given a gene name it returns a list of transcript names (exon specific)
		'''
		if gene.upper() in self.names_transcripts:
			 return (self.names_transcripts[gene.upper()])
		else: return []

	def getGeneCoords(self,ID ,avoid_start=0,avoid_end=0):
		'''
		given a gene name it tries to return a list of all of its exons 
		coordinates (without the start and end coordinates)
		'''
		indices = set()
		for everyID in ID.split('+'):
			for exon in self.transcriptNamesFromID(everyID):
				if not exon==None:
					chr,start,end,strand = self.getTranscriptCoords(exon)
					indices.update(set(range(start,end)))

			indices = indices-set(range(avoid_start,avoid_end+1))
			indices = sorted(list(indices))
			coords = []
			if indices == []:return None
			start = indices[0]
			for n,i in enumerate(indices):
				if n == 0:continue
				if i != indices[n-1]+1:
					 coords.append((chr,start,i,strand))
					 start = i
			coords.append((chr,start,i+1,strand))
			return coords

	def getTranscriptCoords(self,transcript):
		'''
		Given a transcript name, it returns a tuple (chromosome, start, end, strand)
		'''
		data = self.getExon(transcript)
		return (data['seqname'].replace('chr',''),
				int(data['start']),
				int(data['end']),
		data['strand'])

	def getGeneList(self):
		'''
		returns a list of all the genes
		'''
		genes = []
		for i in self.bigGFTlistdict:
			genes.append(i['gene_name'])
		return sorted(set(genes))

	def transcriptNamesFromID(self,ID):
		'''
		returns a list of all transcripts from the ENSMUSG ID
		'''
		if ID in self.gene_IDS:
			return self.gene_IDS[ID]
		else : return []

def purge():
	'''redoes the GTF object'''
	GTF = gtf()
	pickle.dump(GTF,open('annotations/gtf.p','wb'))   

	

def getGene(chromosome, start, end):
	'''
	given chromosome, start, and end of a gene it returns the
	gene attributes
	'''
	return GTF.getGene(chromosome, start, end)

def getGeneCoords( gene,avoid_start=0,avoid_end=0):
	'''
	given a gene name it tries to return a list of all of its exons coordinates
	'''
	return GTF.getGeneCoords(gene,avoid_start,avoid_end)

def getExon(transcriptName):
	'''
	given a transcript name it returns the GFT data of that transcript in a 
	dictionaty.
	'''
	return GTF.getExon(transcriptName)

def transcriptNames(gene):
	'''
	given a gene name it returns a list of transcript names (exon specific)
	'''
	return GTF.transcriptNames(gene)

def getTranscriptCoords(transcript):
	'''
	given a transcript name, it returns a tuple (chromosome, start, end)
	'''
	return GTF.getTranscriptCoords(transcript)

def getGeneList():
	'''
	returns a list of all the genes
	'''
	return GTF.getGeneList()

if __name__ == '__main__':
	while True:
		# if run as main then emulate a REPL
		# usefull for debugging
		try:exec(raw_input('>>>'))
		except Exception as e: print e


