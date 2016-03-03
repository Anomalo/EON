##!/usr/bin/env python2.7
import os
import os.path
from optparse import OptionParser
import csv
import glob
import fa
import gtf
import math
import sys
import time
import subprocess
import atexit

def processes():
	pl = subprocess.Popen('ps aux'.split(), stdout=subprocess.PIPE).communicate()[0]
	return map(lambda x: int(x.split()[1]),pl.split('\n')[1:-1])
def command(cmd):
	pl = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
	return pl.pid

def err(*args):
	sys.stderr.write(' '.join(map(str,args))+'\n')


class maltese:
	def __init__(self,dexseq,verbose=False,sep=',',taxon='Mus_musculus'
		,version='GRCm38',temp=False,annotationDir='annotations',
		output='',inputFormat="0,8,9,10,12,7",outputMode='w',PvalFilter = 0.01):
		self.dexseq=dexseq
		self.prosite = 'ps_scan/prosite.dat'
		self.verbose=verbose
		self.sep=sep
		self.temps = temp
		self.annotationDir = annotationDir
		self.ps_scan = '/'.join(os.path.realpath(__file__).split('/')[:-2])+'/ps_scan/'
		self.dexseqOut=output
		self.inputFormat = inputFormat
		self.outputMode = outputMode
		self.pvalFilter = PvalFilter

		if self.dexseqOut=='':
			self.dexseqOut = '%(dexseq)s.withMotifs.csv'%locals()

		#err(self.ps_scan)
		#return None
		gtf_file = glob.glob(annotationDir+'/*.gtf')
		if gtf_file ==[]:raise Exception('no gtf found in '+annotationDir)
		if len(gtf_file) >1: raise Exception('more than one gtf in %(annotationDir)s '%locals())
		self.gtf_file = gtf_file[0]
		fa.set_taxon(taxon,version,annotationDir = annotationDir)

	def addMotifs(self,skipProsite=False):
		#adds motifs to a dexseq file stated in __init__
		dexseq = self.dexseq
		prosite = self.prosite
		verbose = self.verbose
		sep = self.sep
		ps_scan = self.ps_scan
		if skipProsite==False:
			tempFasta ,ids= self.dexSeqToFasta()
			if verbose: 
				err(tempFasta)
				f = open('%(tempFasta)s'%locals()).read().split('>')
				f = set(map(lambda x: ':'.join(x.split(':')[:3]),f))
				fastaCount= len(f)

			prositeCMD = 'perl %(ps_scan)sps_scan.pl --pfscan %(ps_scan)spfscan -d %(ps_scan)sprosite.dat %(tempFasta)s > %(tempFasta)s.prosite'
			os.system(prositeCMD % locals())
			'''
			PID = command(prositeCMD % locals())
			self.PID = PID
			atexit.register(self.premature_exit)
			if verbose:
				print (prositeCMD % locals())
				print 'starting prosite, PID',PID

				print ''
			time.sleep(5)

			while PID in processes(): #while prosite is running

				if verbose: 
					f = open('%(tempFasta)s.prosite'%locals()).read().split('>')
					F = set(map(lambda x: ':'.join(x.split(':')[:3]),f))
					prositeCount= len(F)
					print  "\033[F",'checked',prositeCount, 'of',fastaCount,'sequences', 100*prositeCount/fastaCount,'% ',
					last = (map(lambda x: ':'.join(x.split(':')[:3]),f))[-1]
					size = last.split('_')[1]
					start,end = size.split('-')
					size = abs(int(start)-int(end))
					print last, 'size = ',size
				time.sleep(1)
			'''
			if verbose: err('ps_scan done, reading results')
		self.prositeToDexseq()

		#this part just cleans the temp files
		if not self.temps: os.system('rm %(tempFasta)s %(tempFasta)s.prosite'%locals())
		if verbose: err( 'completed')
	
	def readPrositeOut(self):
		'''
		reads the prosite output and returns it in a dictionary format
		{exon:"motif:logfold2change:exonMotifCounts:exonLength:backgroundMotifCounts:backgroundLength"}
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
			#count = 0
			virtual =[0]*10**6
			for line in lines:
				if line =='':continue
				if len( line.split())>3:
					start,space,end = line.split()[:3]
					#count += int(end)-int(start)
					for i in range(int(start),int(end)):
						virtual[i]=1
			count = sum(virtual)
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
						points = math.log((float(f_counts)/float(f_length))/(float(b_counts)/float(b_length)),2)
						#if points <= abs(min_score):continue
						points = '%(points).2f'%locals()
						out.append('%(motif)s:%(points)s:%(f_counts)s:%(f_length)s:%(b_counts)s:%(b_length)s'%locals())

					else:
						points = f_counts/ f_length
						points = 'N%(points)s'%locals()
						out.append('%(motif)s:%(points)s:%(f_counts)s:%(f_length)s:0:-'%locals())
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
		inputFormat=self.inputFormat

		IDi,CHR,START,END,STRAND,PVAL,CHANGE = map(int,inputFormat.replace('-','-1').split(','))

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
				newCSV.append(self.sep.join(['prosite_motifs']+header))
				continue
			rowD = dict(zip(header,row))	
			line = row
			n+=1
			try:
				seqname	= row[CHR]#D['genomicData.seqnames']
				strand  = row[STRAND]#D['genomicData.strand']
				start   = row[START]#D['genomicData.start']
				end     = row[END]#D['genomicData.end']
				pVal	= row[PVAL]

			except:continue
			if float(pVal)>self.pvalFilter:continue
			ID = '%(seqname)s_%(start)s-%(end)s_%(strand)s' % locals()
			if ID in proD:line = [proD[ID]]+line
			else: line = ['-']+line
			newCSV.append(self.sep.join(line))
		f.close()
		newCSV='\n'.join(newCSV)
		f = open(self.dexseqOut,self.outputMode)
		f.write(newCSV)
		f.close()
		if self.verbose: err('saved to '+self.dexseqOut)
		
	def dexSeqToFasta(self,linelength=80):
		'''
		reads a dexseq output and produces a fasta file based on the coordinates of sequences altered
		returns the filename of the fasta produced and a list of ids
		'''
		inputFormat=self.inputFormat
		dexseq = self.dexseq
		verbose = self.verbose
		sep = self.sep
		GTF = gtf.gtf(self.gtf_file)
		IDi,CHR,START,END,STRAND,PVAL,CHANGE = map(int,inputFormat.split(','))
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
			#rowD = dict(zip(header,row))	
			n+=1
			seqname	= row[CHR]#D['genomicData.seqnames']
			strand  = row[STRAND]#D['genomicData.strand']
			start   = row[START]#D['genomicData.start']
			end     = row[END]#D['genomicData.end']
			id      = row[IDi]#D['groupID']
			pVal	= row[PVAL]
			#except:continue
			#gene    = rowD['gene']
			start = int(start)
			end = int(end)
			if float(pVal)>self.pvalFilter:continue
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
				err('\033[F%(done).2f%%\tretrving sequence for %(id)s foreground'%locals())
			seq = fa.seq_coords(seqname,start,end,strand)
			length=len(seq)
			seqs = sliceSeq(seq)
			for seq in seqs:			
				ID = '%(seqname)s_%(start)s-%(end)s_%(strand)s:foreground:%(length)s' % locals()
				newFasta= '>%(ID)s\n%(seq)s' % locals()
				fasta.append(newFasta)


			if verbose:
				done = 100*((2*n)/numlines)
				err('\033[F%(done).2f%%\tretrving sequence for %(id)s background'%locals())

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

	def premature_exit(self):
		#kills ps_scan in case of keyboard interrupt
		print 'closing ', self.PID
		os.system('kill '+str(self.PID))

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
