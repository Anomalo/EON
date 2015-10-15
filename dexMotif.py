##!/usr/bin/env python2.7
import os
from optparse import OptionParser
import csv
import glob
from eon import fa
def check_files(taxon='mus musculus',dir = 'annotations/'):
	'''
	checks, downloads and installs annotation files
	'''
	try: os.makedirs(dir)
	except OSError: pass

	#check gtf
	f = glob.glob(dir+'*.gtf')
	if f == []: 
		from eon import gtf
		gtf.getGTF(taxon = taxon, dir = dir)
		reload(gtf)
	
		
def main():
	description='''
	this program tool will glast(go blasting) exons of genes given. outputting how 
	a score of how common the motifs found on the exon where present in genes within 
	go annotation categories. The assumption is that motifs important for a specific 
	gene ontology category would be present in other genes within that GO category.
	'''.replace('\n','').replace('\t','')
	parser = OptionParser(description=description)
	parser.add_option('-f','--file',
				action='store',type='string',
				dest='dexFile',default='',
				help='dexseq output file to analyze')

	parser.add_option('-A','--annotation',
				action='store', type='string',
				dest='annotations',default='',
				help='downloads and generates anotation files of defined taxon -A "mus musculus"')

	parser.add_option('-p','--purge',
				action='store_true',
				dest='purge', default=False,
				help='delete all anotation files (usefull for changing the taxon of interest)')
	
	parser.add_option('-v','--verbose',
				action='store_true',
				dest ='v', default=False,
				help='shows you what am I thinking')
	'''
	parser.add_option('-o','--output',
				action='store', type='string',
				dest='output',default='results',
				help='directory where to save the results')
	'''
	parser.add_option('-n','--fastaline',
				action='store', type='int',
				dest='fastaline',default=80,
				help='line length for fasta files')

	parser.add_option('-P','--prosite',
				action='store', type='string',
				dest='prosite',default='annotations/prosite.dat',
				help='prosite.dat file')
	
	(options, args) = parser.parse_args()
	v = options.v
	#output = options.output
	annotations = options.annotations
	prosite = options.prosite
	dexseq = options.dexFile
	if options.purge:
		commandOptions ='-rf'
		if v:commandOptions+='v'
		os.system('rm '+commandOptions+' annotations/*')
	
	if annotations != '':
		check_files(taxon = annotations)
	
	tempFasta ,ids= dexSeqToFasta(dexseq,verbose = v,linelength=options.fastaline)
	if v:print tempFasta	
	prositeCMD = 'perl ps_scan/ps_scan.pl -d %(prosite)s %(tempFasta)s > %(tempFasta)s.prosite '
	if v:print prositeCMD % locals()
	os.system(prositeCMD % locals())
	
	
	prositeToDexseq(dexseq,
			'%(dexseq)s.withMotifs.csv'%locals(),
			verbose=True)
	
	#this part just cleans the temp files
	os.system('rm %(tempFasta)s %(tempFasta)s.prosite'%locals())
	if v: print 'completed'

def readProsite(prosite):
	f = open(prosite)
	prosite = f.read()
	chunks = prosite.split('>')
	proD={}
	for chunk in chunks:
		if chunk=='':continue
		lines = chunk.split('\n')
		description = lines.pop(0)
		count = len(lines)
		name = description.split(' :')[0]
		motif = description.split()[3]
		if not name in proD: proD[name]=''
		proD[name]+='%(motif)s/%(count)s '%locals()
	f.close()
	return proD
	

def prositeToDexseq(dexseq,dexseqOut,sep=',',verbose=False):
	'''
	reads a dexseq file and its prosite output and saves the results in dexseqOut
	'''
	f = open(dexseq)
	csvfile = csv.reader(f,delimiter=sep)
	n=1
	newCSV= []
	prosite = '%(dexseq)s.tmp.fasta.prosite'%locals()
	proD = readProsite(prosite)
	for row in csvfile:
		if n==1:
			header = row
			n+=1
			newCSV.append(','.join(header+['prosite_motifs']))
			continue
		rowD = dict(zip(header,row))	
		line = row
		n+=1
		seqname	= rowD['genomicData.seqnames']
		strand  = rowD['genomicData.strand']
		start   = rowD['genomicData.start']
		end     = rowD['genomicData.end']

		ID = '%(seqname)s:%(start)s-%(end)s_%(strand)s' % locals()
		line.append(proD[ID])
		newCSV.append(','.join(line))

	f.close()
	newCSV='\n'.join(newCSV)
	f = open(dexseqOut,'w')
	f.write(newCSV)
	f.close()
	
		
def dexSeqToFasta(dexseq,sep=',',verbose = False,linelength=80):
	'''
	reads a dexseq output and produces a fasta file based on the coordinates of sequences altered
	returns the filename of the fasta produced and a list of ids
	'''
	f = open(dexseq)
	csvfile = csv.reader(f,delimiter=sep)
	n=0
	fasta = []
	ids = []
	for row in csvfile:
		if n==0:
			header = row
			n+=1
			continue
		rowD = dict(zip(header,row))	
		n+=1
		seqname	= rowD['genomicData.seqnames']
		strand  = rowD['genomicData.strand']
		start   = rowD['genomicData.start']
		end     = rowD['genomicData.end']
		#gene    = rowD['gene']
		seq     = fa.seq_coords(seqname,start,end,strand)
		choppedSeq = ''
		for i in range(0,len(seq),linelength):
			choppedSeq+=seq[i:i+linelength]+'\n'
		ID = '%(seqname)s:%(start)s-%(end)s_%(strand)s' % locals()
		newFasta= '>%(ID)s\n%(choppedSeq)s' % locals()
		fasta.append(newFasta)
		ids.append(ID)
	f.close()
	fasta = ''.join(fasta)
	fastaFname = '%(dexseq)s.tmp.fasta'%locals()
	f = open(fastaFname,'w')
	f.write(fasta)
	f.close()
	return fastaFname,ids

if __name__ == '__main__': 
	#check_files()
	main()
'''
        while True:
                try:exec(raw_input('>>>'))
                except Exception as e: print e

'''
