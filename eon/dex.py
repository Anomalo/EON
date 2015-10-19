##!/usr/bin/env python2.7
import os
from optparse import OptionParser
import csv
import glob
import fa
#import gtf
class dex:
	def __init__(self,dexseq,prosite,verbose=False,sep=','):
		self.dexseq=dexseq
		self.prosite=prosite
		self.verbose=verbose
		self.sep=sep

		

	def addMotifs(self):
		dexseq = self.dexseq
		prosite = self.prosite
		verbose = self.verbose
		sep = self.sep

		tempFasta ,ids= self.dexSeqToFasta()
		if verbose:print tempFasta
		prositeCMD = 'perl ps_scan/ps_scan.pl --pfscan ps_scan/pfscan -d %(prosite)s %(tempFasta)s > %(tempFasta)s.prosite '
		if verbose:print prositeCMD % locals()
		os.system(prositeCMD % locals())
		if verbose:print 'ps_scan done, reading results'
		
		self.prositeToDexseq()
		
		#this part just cleans the temp files
		os.system('rm %(tempFasta)s %(tempFasta)s.prosite'%locals())
		if verbose: print 'completed'
	
	def readPrositeOut(self):
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
			count = len(lines)
			name , ground, length,motif = description.split(':')
			motif = motif.split()[1]
			if not name in proD: proD[name]={}
			if not motif in proD[name]:proD[name][motif]={}
			proD[name][motif][ground]=map(int,(count,length))
		f.close()
		
		for exon,motifs in proD.iteritems():
			out = []
			for motif,counts in motifs.iteritems():
				if 'foreground' in counts:
					f_counts, f_length=counts['foreground']
					if 'background' in counts:	
						b_counts, b_length=counts['background']
						points = (float(f_counts)/float(f_length))/(float(b_counts)/float(b_length))
						points = '%(points).2f'%locals()
					else:
						points = 'N%(f_counts)s'%locals()
					out.append('%(motif)s:%(points)s'%locals())
			out = ' ; '.join(out)
			proD[exon]=out
		return proD
		
	
	def prositeToDexseq(self):
		'''
		reads a dexseq file and its prosite output and saves the results in dexseqOut
		'''
		dexseq = self.dexseq
		verbose = self.verbose
		sep = self.sep
		dexseqOut = '%(dexseq)s.withMotifs.csv'%locals()

		f = open(dexseq)
		csvfile = csv.reader(f,delimiter=sep)
		n=1
		newCSV= []
		prosite = '%(dexseq)s.tmp.fasta.prosite'%locals()
		proD = self.readPrositeOut()
		for row in csvfile:
			if n==1:
				header = row
				n+=1
				newCSV.append(','.join(['prosite_motifs']+header))
				continue
			rowD = dict(zip(header,row))	
			line = row
			n+=1
			seqname	= rowD['genomicData.seqnames']
			strand  = rowD['genomicData.strand']
			start   = rowD['genomicData.start']
			end     = rowD['genomicData.end']

			ID = '%(seqname)s_%(start)s-%(end)s_%(strand)s' % locals()
			if ID in proD:line = [proD[ID]]+line
			else: line = ['-']+line
			newCSV.append(','.join(line))
		f.close()
		newCSV='\n'.join(newCSV)
		f = open(dexseqOut,'w')
		f.write(newCSV)
		f.close()
		
			
	def dexSeqToFasta(self,linelength=80):
		'''
		reads a dexseq output and produces a fasta file based on the coordinates of sequences altered
		returns the filename of the fasta produced and a list of ids
		'''
		dexseq = self.dexseq
		verbose = self.verbose
		sep = self.sep
	
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
			id      = rowD['groupID']
			#gene    = rowD['gene']
			start = int(start)
			end = int(end)
			if verbose:
				print 'retrving sequence for ',id, 'foreground'
			seq = fa.seq_coords(seqname,start,end,strand)
			length=len(seq)
			seqs = sliceSeq(seq)
			for seq in seqs:			
				ID = '%(seqname)s_%(start)s-%(end)s_%(strand)s:foreground:%(length)s' % locals()
				newFasta= '>%(ID)s\n%(seq)s' % locals()
				fasta.append(newFasta)


			if verbose:
				print 'retrving sequence for ',id,'background'
			seq = ''.join(fa.seqs_coords(gtf.getGeneCoords(id,
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
	returns a list of sequences, each sequence no longer than mac_length
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
