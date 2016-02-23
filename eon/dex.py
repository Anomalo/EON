##!/usr/bin/env python2.7
import os
import os.path
from optparse import OptionParser
import csv
import glob
import fa
import gtf
import sys
def err(*args):
	sys.stderr.write(' '.join(map(str,args))+'\n')
class dex:
	def __init__(self,dexseq,verbose=False,sep=',',taxon='Mus_musculus'
			,version='GRCm38',temp=False,annotationDir='annotations'):
		self.dexseq=dexseq
		self.prosite = 'ps_scan/prosite.dat'
		self.verbose=verbose
		self.sep=sep
		self.temps = temp
		self.annotationDir = annotationDir
		self.ps_scan = '/'.join(os.path.realpath(__file__).split('/')[:-2])+'/ps_scan/'
		#err(self.ps_scan)
		#return None
		gtf_file = glob.glob(annotationDir+'/*.gtf')
		if gtf_file ==[]:raise Exception('no gtf found in '+annotationDir)
		if len(gtf_file) >1: raise Exception('more than one gtf in %(annotationDir)s '%locals())
		self.gtf_file = gtf_file[0]
		fa.set_taxon(taxon,version,annotationDir = annotationDir)
#		if not os.path.isfile(self.prosite):
#			os.system('wget -O annotations/prosite.dat ftp://ftp.expasy.org/databases/prosite/prosite.dat')
	

	def addMotifs(self,min_score=2):
		dexseq = self.dexseq
		prosite = self.prosite
		verbose = self.verbose
		sep = self.sep
		ps_scan = self.ps_scan
		tempFasta ,ids= self.dexSeqToFasta()
		if verbose: err(tempFasta)
		prositeCMD = 'perl %(ps_scan)sps_scan.pl --pfscan %(ps_scan)spfscan -d %(ps_scan)sprosite.dat %(tempFasta)s > %(tempFasta)s.prosite '
		if verbose:err(prositeCMD % locals())
		os.system(prositeCMD % locals())
		if verbose: err('ps_scan done, reading results')
		
		self.prositeToDexseq(min_score=min_score)

		#this part just cleans the temp files
		if not self.temps: os.system('rm %(tempFasta)s %(tempFasta)s.prosite'%locals())
		if verbose: err( 'completed')
	
	def readPrositeOut(self,min_score=2):
		'''
		reads the prosite output and returns it in a dictionary format
		'''
		dexseq = self.dexseq
		prositeOutput = '%(dexseq)s.tmp.fasta.prosite'%locals()
		f = open(prositeOutput)
		prosite = f.read()
		chunks = prosite.split('>')
		proD={}
		for chunk in chunks:
			if chunk=='':continue
			lines = chunk.split('\n')
			description = lines.pop(0)
			count = 0
			for line in lines:
				if line =='':continue
				if len( line.split())>3:
					start,space,end = line.split()[:3]
					count += int(end)-int(start)
			#count = len(lines)
			name , ground, length,motif = description.split(':')
			motif = motif.split()[1]
			if not name in proD: proD[name]={}
			if not motif in proD[name]:
				proD[name][motif]={}
			proD[name][motif][ground]=map(int,(0,length))
			proD[name][motif][ground][0]+=round(count,2)
		f.close()
		
		for exon,motifs in proD.iteritems():
			out = []
			for motif,counts in motifs.iteritems():
				if 'foreground' in counts:
					f_counts, f_length=counts['foreground']
					if 'background' in counts:	
						b_counts, b_length=counts['background']
						points = (float(f_counts)/float(f_length))/(float(b_counts)/float(b_length))
						if points < min_score:continue
						points = '%(points).2f'%locals()
					else:
						points = 'N%(f_counts)s'%locals()
					out.append('%(motif)s:%(points)s'%locals())
			out = ' ; '.join(out)
			proD[exon]=out
		return proD
		
	
	def prositeToDexseq(self,min_score=2):
		'''
		reads a dexseq file and its prosite output and saves the results in dexseqOut
		'''
		dexseq = self.dexseq
		verbose = self.verbose
		sep = self.sep
		#dexseqOut = '%(dexseq)s.withMotifs.csv'%locals()

		f = open(dexseq)
		csvfile = csv.reader(f,delimiter=sep)
		n=1
		#newCSV= []
		prosite = '%(dexseq)s.tmp.fasta.prosite'%locals()
		proD = self.readPrositeOut(min_score=min_score)
		for row in csvfile:
			if n==1:
				header = row
				n+=1
				print ','.join(['prosite_motifs']+header)
				continue
			rowD = dict(zip(header,row))	
			line = row
			n+=1
			try:
				seqname	= rowD['genomicData.seqnames']
				strand  = rowD['genomicData.strand']
				start   = rowD['genomicData.start']
				end     = rowD['genomicData.end']
			except:continue
			ID = '%(seqname)s_%(start)s-%(end)s_%(strand)s' % locals()
			if ID in proD:line = [proD[ID]]+line
			else: line = ['-']+line
			print ','.join(line)
		f.close()
		#newCSV='\n'.join(newCSV)
		#f = open(dexseqOut,'w')
		#f.write(newCSV)
		#f.close()
		
			
	def dexSeqToFasta(self,linelength=80):
		'''
		reads a dexseq output and produces a fasta file based on the coordinates of sequences altered
		returns the filename of the fasta produced and a list of ids
		'''
		dexseq = self.dexseq
		verbose = self.verbose
		sep = self.sep
		GTF = gtf.gtf(self.gtf_file)
		if verbose:
			f = open(dexseq)
			numlines= float(len(f.readlines()))*2
			f.close()
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
#			try:
			seqname	= rowD['genomicData.seqnames']
			strand  = rowD['genomicData.strand']
			start   = rowD['genomicData.start']
			end     = rowD['genomicData.end']
			id      = rowD['groupID']
#			except:continue
			#gene    = rowD['gene']
			start = int(start)
			end = int(end)
			'''
			for mode in ['foreground','background','upstream','downstream']:
				if verbose:
					done = 100*((2*n-1)/numlines)
					err('%(done).2f%%\tretrving sequence for %(id)s %(mode)s'%locals())
				if mode == 'foreground':
					seq = fa.seq_coords(seqname,start,end,strand)
				elif mode == 'background':
					seq = ''.join(fa.seqs_coords(GTF.getGeneCoords(id,
										avoid_start=start,
										avoid_end  =end)))
				elif mode == 'upstream':
					continue
					if strand == '+':
						start,end = int(end),int(end)+200
					else:
						start,end = int(start)-200,int(start)
					seq = fa.seq_coords(seqname,start,end,strand)
				elif mode == 'downstream':
					continue
					if strand == '+':
						start,end = int(start)-200,int(start)
					else:
						start,end = int(end),int(end)+200
					seq = fa.seq_coords(seqname,start,end,strand)

				elid mode == 'downstream':
					continue
				else: continue  #this part avoids the script from running with downstream or upstream
						
						#in order for upstream and downstream to run a suitable background
						#is needed.
						#
						#the readPrositeOut() needs to be updated to handle this modes
				length=len(seq)
				seqs = sliceSeq(seq)
				for seq in seqs:			
					ID = '%(seqname)s_%(start)s-%(end)s_%(strand)s:%(mode)s:%(length)s' % locals()
					newFasta= '>%(ID)s\n%(seq)s' % locals()
					fasta.append(newFasta)
	
	
			'''
			if verbose:
				done = 100*((2*n-1)/numlines)
				err('%(done).2f%%\tretrving sequence for %(id)s foreground'%locals())
			seq = fa.seq_coords(seqname,start,end,strand)
			print seq
			length=len(seq)
			seqs = sliceSeq(seq)
			for seq in seqs:			
				ID = '%(seqname)s_%(start)s-%(end)s_%(strand)s:foreground:%(length)s' % locals()
				newFasta= '>%(ID)s\n%(seq)s' % locals()
				fasta.append(newFasta)


			if verbose:
				done = 100*((2*n)/numlines)
				err('%(done).2f%%\tretrving sequence for %(id)s background'%locals())
			seq = ''.join(fa.seqs_coords(GTF.getGeneCoords(id,
									avoid_start=start,
									avoid_end  =end)))
			length=len(seq)
			seqs = sliceSeq(seq)
			for seq in seqs:			
				ID = '%(seqname)s_%(start)s-%(end)s_%(strand)s:background:%(length)s' % locals()
				newFasta= '>%(ID)s\n%(seq)s' % locals()
				fasta.append(newFasta)




			ids.append(ID)
		f.close()
		fasta = ''.join(fasta)
		fastaFname = '%(dexseq)s.tmp.fasta'%locals()
		f = open(fastaFname,'w')
		f.write(fasta)
		f.close()
		return fastaFname,ids
def sliceSeq(seq,max_length=40000,linelength=80):
	'''
	returns a list of sequences, each sequence no longer than max_length
	and with each line no longer than linelength
	'''
	seqs = []
	for i in range(0,len(seq),max_length):
		subseq=seq[i:i+max_length]
		choppedSeq=''
		for j in range(0,len(subseq),linelength):
			choppedSeq+=subseq[j:j+linelength]+'\n'
		seqs.append(choppedSeq)
	return seqs
if __name__ == '__main__': 
	#check_files()
	main()
'''
        while True:
                try:exec(raw_input('>>>'))
                except Exception as e: print e

'''
