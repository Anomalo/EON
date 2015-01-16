from pprint import pprint
import os
import subprocess
import cPickle as pickle

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
def genGOB(gmt):
	'''
	given a gmtfile, it creates a .gob (GOBlast) file
	a row of BGO looks like:
	[go term]	[all exons of all genes with that go concated by ';']
	'''
	print 'importing gtf'
	import gtf
	print 'importing fa'
	import fa

	f = open(gmt)
	num_lines = file_len(gmt)
	newFile = gmt[:-3]+'gob'
	open(newFile, 'w').close()
	newFile = open(newFile,'a')
	print 1
	genesSeen={}
	n=0.0
	for i in range(num_lines):

		os.system('clear')
		percentage = str(round(100*n/num_lines,2))+'%'
		params = {'step':'GOB','percentage': percentage,'left':num_lines - n}
		remainder ='%(step)s_%(percentage)s_%(left)s_left' 
		remainder = remainder % params
		print remainder
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
						exons.append( fa.seq_coords(c,s,e,strand))
					genesSeen.update({gene:exons})
				except:
					genesSeen.update({gene:[]})
					print 'fail','\t',gene
			print '\t',gene
			exons = genesSeen[gene]
		exons = ';'.join(exons)
		newLine = GO + '\t' + url + '\t' + exons + '\n'
		newFile.write(newLine)
		n+=1

	os.system('rm *.remainder')
	os.system('touch DONE.remainder')

def genGIB(gob, word_size=11):
	'''
	makes a GIB (gob index), from a gob file.
	'''

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
			seq = seqs[seqIndex]
			words = len(seq)-word_size
			for i in range(words):
				word = seq[i:i+word_size]
				if not word in metaGIB[-1]: metaGIB[-1].update({word:[]})
				metaGIB[-1][word].append({'go':go,
						'goID':goid,
						'url':url,
						'line':l,
						'seqNumber':seqIndex,
						'location':i})
		n+=1
	print 'Merging GIBS into one GIB'
	gibD={}
	for d in metaGIB:
		gibD = dict(gibD,**d)		
	pickle.dump(gibD,open(gibfname,'wb'))
	os.system('rm *.remainder')
if __name__ == '__main__':
	#genGOB('annotations/Mus_musculus_GSEA_GO_sets_all_symbols_September_2013.gmt')
	genGIB('annotations/Mus_musculus_GSEA_GO_sets_all_symbols_September_2013.gob') 

