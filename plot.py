import matplotlib.pyplot as plt
import glob

def plotDIR(dir):
	''' pieplots all the TSVs in that directory '''
	files = glob.glob(dir+'/*.tsv')
	for f in files:
		pieplot(f)
		
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
	
def pieplot(tsv):
	'''pieplots atsv file with a header'''
	f = open(tsv,'r')
	exon = f.read()
	f.close()
	lines = exon.splitlines()
	head,data = lines[:2], lines[2:]
	head ='\n'.join(head)
	vals , labels, colors = [],[],[]
	for line in data:
		label,val = line.split('\t')
		val = float(val)
		if val > 0:
			vals.append(val)
			c = color(label)
			colors.append(c)
			labels.append(wrap(label.replace('_',' '),10))
	plt.close('all')
	f, ax = plt.subplots(figsize(10,10))
	ax.pie(vals, labels=labels,colors=colors,autopct='%1,1f%%')
	ax.set_title(head)
	fname = tsv[:-3]+'png'
	plt.savefig(fname)

	
