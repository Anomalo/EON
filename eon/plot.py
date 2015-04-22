import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab

import glob

def plotDIR(dir,distortion=3):
	''' pieplots all the TSVs in that directory '''
	files = glob.glob(dir+'/*.tsv')
	for f in files:
		pieplot(f,distortion=distortion)
		#colplot(f)
		
def wrap(txt, linelength=15):
	'''turns a line of text into a wraped block'''
	words = txt.split()
	out = ''
	length = 0
	for word in words:
		if length > linelength:
			out +='\n'
			length = 0
		out += ' '+word
		length += len(word)+1
	return out

def color(txt):
	'''given a string returns a random number consistently'''
	h = str(hash(txt)**4)
	r = float(h[-3:-1])/100
	g = float(h[-5:-3])/100
	b = float(h[-7:-5])/100
	return r,g,b
def readTSV(tsv,wraplabel=10):

	f = open(tsv,'r')
	exon = f.read()
	f.close()
	lines = exon.splitlines()
	head,data = lines[:2], lines[2:]
	head ='\n'.join(head)
	vals , labels, colors = [],[],[]
	for line in data:
		label,val = line.split('\t')[:2]
		val = float(val)
		if val > 0:
			c = color(label)
			colors.append(c)
			labels.append(wrap(label.replace('_',' ')+' '+str(round(val,2)),wraplabel))
			vals.append(val)
	return vals,labels,colors
	
def pieplot(tsv,show=False,distortion=3,wraplabel=10):
	'''pieplots atsv file with a header'''
	f = open(tsv,'r')
	exon = f.read()
	f.close()
	lines = exon.splitlines()
	head,data = lines[:2], lines[2:]
	head ='\n'.join(head)
	vals , labels, colors = [],[],[]
	for line in data:
		label,val = line.split('\t')[:2]
		val = float(val)
		if val > 0:
			c = color(label)
			colors.append(c)
			labels.append(wrap(label.replace('_',' ')+' '+str(round(val,2)),wraplabel))
			vals.append(val**distortion)
	pylab.close('all')
	#pylab.figure(figsize=(13,13))
	ax = pylab.axes(aspect=1)
	#f, ax = pylab.subplots(figsize(10,10))
	pylab.pie(vals, labels=labels,colors=colors,)
	ax.set_title(head)
	fname = tsv[:-4]+'_pie.png'
	pylab.savefig(fname)
	if show:pylab.show()
	
def colplot(tsv):
	vals,labels,colors = readTSV(tsv,wraplabel=1000)
	pylab.close('all')
	pylab.figure(figsize=[.7+len(vals),10])
	ax = pylab.axes()
	index = np.arange(len(vals))
	barwidth = 1
	ax.bar(index,
		vals,
		width=barwidth,
		)
	pylab.xticks(index+barwidth,labels,rotation=90)
	fname = tsv[:-4]+'_bar.png'
	pylab.savefig(fname,dpi=200)

def main():
	import sys
	dir = sys.argv[-2]
	distortion = sys.argv[-1]
	plotDIR(dir)

if __name__=="__main__": main()
