import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import scipy.stats as st
import matplotlib.cm as cmx
import matplotlib.colors as colors
from mpl_toolkits.axes_grid.inset_locator import inset_axes
rcParams.update({'figure.autolayout': True})
from pprint import pprint


def bionmial_test(sample,background):
	mean = np.mean(background)
	std = np.std(background)
	z_score = (sample - mean)/std
	p_values = st.norm.sf(abs(z_score))
	return p_values

def summaryDex(fname,fnameOut='test',sep=',',FORMAT="0,8,9,10,12,7,-",plotFormat='pdf'):
	'''
	given the maltesers motif of dexseq output, it will generate a boxplots summarizing the motifs changes
	'''
	IDi,GENENAME,CHR,START,END,STRAND,PVAL,CHANGE = map(lambda x: int(x)+6,FORMAT.replace('-','-1').split(','))
	if GENENAME ==-1: GENENAME=IDi
	f = open(fname).read().split('\n')
	plotFormat = '.'+plotFormat.replace('.','')
	CHANGE-=1
	header = f.pop(0)
	header = header.replace('prosite_motifs',
							sep.join([	'motif',
									'logFold2',
									'motifExonCount',
									'exonLen',
									'motifGeneCount',
									'geneLen']))
	newLines = [header]
	for line in f:
		malteseData = line.split(sep)[0]
		dexseqData = line.replace(malteseData,'')
		#this loops gets a list of all changes to form a binomial distribution
		changes = []
		#for motif in malteseData.split(' ; '):
		#	if motif=='-':continue
		#	changes.append(float(motif.split(':')[1]))
		for motif in malteseData.split(' ; '):
			if motif=='-':
				newLines.append(sep.join(['-']*7)+'%(dexseqData)s'%locals())
				continue
			if motif=='':continue
			motif, logFold2,motifExonCount,exonLen,motifGeneCount,geneLen = motif.split(':')
			#print motif, logFold2,motifExonCount,exonLen,motifGeneCount,geneLen
			#pValue = bionmial_test(float(logFold2),changes)
			newLine = '%(motif)s%(sep)s%(logFold2)s%(sep)s%(motifExonCount)s%(sep)s%(exonLen)s%(sep)s%(motifGeneCount)s%(sep)s%(geneLen)s%(dexseqData)s'%locals()
			newLines.append(newLine)
	open(fnameOut+'.csv','w').write('\n'.join(newLines))
	##############################################################################################
	#makes plot###################################################################################
	##############################################################################################
	#joins all the motifs together
	#pprint(	map(lambda x: (x.split(sep),len(x.split(sep))),
	#		open(fnameOut).read().split('\n')[1:]))
	f = map(lambda x: x.split(sep)[:7]+[x.split(sep)[CHANGE]]+[x.split(sep)[GENENAME]],
			open(fnameOut+'.csv').read().split('\n')[1:])
	data = {}
	colChanges = {}
	#for i in enumerate(map(lambda x: x.split(sep),
	#		open(fnameOut).read().split('\n'))[0]):print i
	
	for motif, change, n,n,n,n, n,colChange,exon in f:
		#print motif, change,colChange,exon
		if change =='-':continue
		if 'N' in change:
			#print motif,change, exon
			continue
		change = float(change)
		if not motif in data:
			data[motif]=[]
			colChanges[motif]=[]
		data[motif].append(change)
		colChanges[motif].append(colChange)

	data = (sorted(data.items(),
			key=lambda x: (sum(x[1])/len(x[1]))))
	values = map(lambda x: x[1],data)
	labels = map(lambda x: x[0],data)
	fig = plt.figure(figsize = (10,0.5*len(data)))
	ax = fig.add_subplot(111)
	if values != []:
		bp = ax.boxplot(values, 0, '', 0,patch_artist=True)
		vminmax = 0.5
		cNorm  = colors.Normalize(vmin=-vminmax, vmax=vminmax)
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='bwr')
		for i in range(len(labels)): #this adds all the points
			y = data[i][1]
			colChange = colChanges[data[i][0]]
			changes = data[i][1]
			x = np.random.normal(1+i, 0.04, size=len(y))
			for X,Y,c in zip(x,changes,colChange):
				#c = '#e7298a'
				plt.plot(Y,X,'r.',marker='o', color=scalarMap.to_rgba(float(c)))#, alpha=0.5)

		## change outline color, fill color and linewidth of the boxes
		for box in bp['boxes']:
	    	# change outline color
			box.set( color='#7570b3', linewidth=2)
	    		# change fill color
			box.set( facecolor = '#1b9e77' )

		## change color and linewidth of the whiskers
		for whisker in bp['whiskers']:
			whisker.set(color='#7570b3', linewidth=2)

		## change color and linewidth of the caps
		for cap in bp['caps']:
			cap.set(color='#7570b3', linewidth=2)

		## change color and linewidth of the medians
		for median in bp['medians']:
			median.set(color='#b2df8a', linewidth=2)

		## change the style of fliers and their fill
		for flier in bp['fliers']:
			flier.set(marker='o', color='#e7298a', alpha=0.5)

		ax.set_xlabel('logFold2')
		ax.grid(True)
		LABELS = []
		for a,b in zip(labels, values):
			LABELS.append(a+' ('+str(len(b))+')')
		plt.yticks(range(1,1+len(LABELS)),
			  	LABELS)
		fig.savefig(fnameOut+'_motifs'+plotFormat)
		print fnameOut+'_motifs'+plotFormat
	
	
	##############################################################################################
	#makes plot###################################################################################
	##############################################################################################
	#joins all the motifs together

	f = map(lambda x: x.split(sep)[:7]+[x.split(sep)[CHANGE]]+[x.split(sep)[GENENAME]],
			open(fnameOut+'.csv').read().split('\n')[1:])
	data = {}
	colChanges = {}
	#for i in enumerate(map(lambda x: x.split(sep),
	#		open(fnameOut).read().split('\n'))[0]):print i
	#print CHANGE
	for motif, change, n,n,n,n, n,colChange,exon in f:
		if change =='-':continue
		if 'N' in change:
			#print motif,exon, change
			continue
		change = float(change)
		if not exon in data:
			data[exon]=[]
			colChanges[exon]=[]
		data[exon].append(change)
		colChanges[exon].append(colChange)

	data = (sorted(data.items(),
			key=lambda x: (sum(x[1])/len(x[1]))))
	values = map(lambda x: x[1],data)
	labels = map(lambda x: x[0],data)
	fig = plt.figure(figsize = (10,0.5*len(data)))
	ax = fig.add_subplot(111)
	if values != []:
		bp = ax.boxplot(values, 0, '', 0,patch_artist=True)
		for i in range(len(labels)): #this adds all the points
			y = data[i][1]
			colChange = colChanges[data[i][0]]
			changes = data[i][1]
			x = np.random.normal(1+i, 0.04, size=len(y))
			for X,Y,c in zip(x,changes,colChange):
				#c = '#e7298a'
				plt.plot(Y,X,'r.',marker='o', color=scalarMap.to_rgba(float(c)))#, alpha=0.5)

		## change outline color, fill color and linewidth of the boxes
		for box in bp['boxes']:
	    	# change outline color
			box.set( color='#7570b3', linewidth=2)
			# change fill color
			box.set( facecolor = '#1b9e77' )

		## change color and linewidth of the whiskers
		for whisker in bp['whiskers']:
			whisker.set(color='#7570b3', linewidth=2)

		## change color and linewidth of the caps
		for cap in bp['caps']:
			cap.set(color='#7570b3', linewidth=2)

		## change color and linewidth of the medians
		for median in bp['medians']:
			median.set(color='#b2df8a', linewidth=2)

		## change the style of fliers and their fill
		for flier in bp['fliers']:
			flier.set(marker='o', color='#e7298a', alpha=0.5)

		ax.set_xlabel('logFold2')
		ax.grid(True)
		LABELS = []
		for a,b in zip(labels, values):
			LABELS.append(a+' ('+str(len(b))+')')
		plt.yticks(range(1,1+len(LABELS)),
			  	LABELS)

		gradient = np.linspace(0, 1, 256)
		gradientlabels = map(lambda x: x/float(max(gradient))*vminmax*2-vminmax,gradient)
		gradient = np.vstack((gradientlabels, gradientlabels)).transpose()

		#legend = inset_axes(ax,
    	#                width=0.5, # width = 30% of parent_bbox
    	#                height="30%", # height : 1 inch
    	#                loc=4)
		#legend.imshow(gradient, aspect='auto', cmap='bwr')
		'''
		legend = mpl.colorbar.ColorbarBase(ax, cmap='bwr',
                                	# to use 'extend', you must
                                	# specify two extra boundaries:
                                	#boundaries=[0] + bounds + [13],
                                	#extend='both',
                                	#ticks=bounds,  # optional
                                	spacing='proportional',
                                	orientation='horizontal')
    	'''
		#legend.get_yaxis().set_visible(False)
		#legend.get_xaxis().set_visible(False)
		#legend.set_xlabel(' '.join([-vminmax,0,1.5vminmax]))
		#legend.set_yticklabels()
		#legend.set_yaxis([-vminmax,0,vminmax])
		fig.savefig(fnameOut+'_exons'+plotFormat)
		print fnameOut+'_exons'+plotFormat


	f = map(lambda x: x.split(sep)[:7]+[x.split(sep)[CHANGE]]+[x.split(sep)[GENENAME]],
			open(fnameOut+'.csv').read().split('\n')[1:])
	motifs = sorted(list(set([x[ 0] for x in f])))
	exons  = sorted(list(set([x[-1] for x in f])))
	table = [['']+list(motifs)]
	for exon in exons:
		table+=[[exon]]
		for motif in motifs:
			X = (filter(lambda x: exon in x and motif in x,f))
			if len(X)==0:	table[-1].append('0')
			else:		table[-1].append(str(X[0][1]))
	table = '\n'.join(map(lambda x: sep.join(x),table))
	f=open(fnameOut+'_motifGene.csv','w')
	print fnameOut+'_motifGene.csv'
	f.write(table)
	f.close()

	try:
		import pandas
		import seaborn as sns
		data = pandas.DataFrame.from_csv(fnameOut+'_motifGene.csv',sep=sep).applymap(lambda x: float(str(x).replace('N','')))
		f = sns.clustermap(data,
				figsize=(
					0.3	* len(data.columns),
					0.2	* len(data.index)
					))
		plt.setp(f.ax_heatmap.get_yticklabels(), rotation=0)
		f.savefig(fnameOut+"_motifGene"+plotFormat)
		print fnameOut+"_motifGene"+plotFormat
	except:
		print "pandas or seaborn python modules not available"
		print "these modules are needed for clustermap plotting only"

if __name__ == '__main__':
	summaryDex('smallDexseq.csv.withMotifs.csv')
    # while True:
    #     try:exec(raw_input('>>>'))
    #     except Exception as e: print e

