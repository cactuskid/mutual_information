#functions for dealing with prot families involved in the interaction graph
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import shlex
import colour

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import time

import csv
import math

import gc

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.SeqUtils import seq1

import glob

import mutualInfo
from multiprocessing import Pool
import pickleutils
import similarityMatrix

####multidimensional MI analysis based on sequence and physical props#####
####you can incrementally feed in more prots to build a network#####

#return a dictionary of props for each AA and gap
def init_props(propfile):
	handle = open(propfile, 'r')
	my_data = csv.reader(handle, delimiter=',')
	propdict = {}
	AAnum = 0
	for i, line in enumerate(my_data):
		if i == 0 :
			header = line
			for prop in header[2:]:
				propdict[prop] = {} 
			#propdict['AA'] = {}
		else :
			for index,entry in enumerate(line):
				if index == 1:
					AA = entry
					AAnum +=1
					##propdict['AA'][AA] = AAnum 
				if index > 1:
					propdict[header[index]][AA] = float(entry)
	for prop in header[2:]:
		propdict[prop]['-'] = 0
	#propdict['AA']['-'] =0
	return propdict

def column_to_mat(args):
	column,propdict = args
	Mat = np.zeros((len(column),len(propdict)))
	for i, letter in enumerate(column):
		for j,prop in enumerate(propdict):
			try:
				Mat[i,j] = propdict[prop][letter]
			except:
				Mat[i,j] = 0

	return Mat

def line_to_sig(args):
	seq, gaussian, propdict = args
	lines = []
	for i in range(len(propdict)):
		lines.append(np.zeros((len(seq),1)))
	#build vector of prot
	for i,prop in enumerate(propdict)
		for j,letter in enumerate(seq):
			lines[i][j] = propdict[prop][letter]
	#use mirror conidition convolution
	if gaussian != None:
		for line in lines:
			line = 

	return lines
	

def chop_lines(prot , lines , propdict):
	mats = []
	for i in range(0,maxlen):
		Mats.append([])
		newmat = np.zeros(len(allLines.keys()), len(propdict.keys()))			
		for code in allLines:
			for signal in enumerate(allLines[code]):
				


def aln_to_Mats(prot,propdict):
	#returns numpy mat of physical props for each sequence in an alignment to be passed to mutualinfo analysis

	#X is a matrix corresponding to one column of the alignment
	#It contains the amino acids' physical props
	# X.shape = nsequences * nphysicalprops
	#each column is a matrix in Mats
	columns = []
	Mats =[]
	#store the uniprot refs of each prot in the alignment
	refs =[]
	params = []

	#generate raw and filtered vectors for all prots in a fasta
	record_iterator = AlignIO.read(prot, "fasta")
	maxlen = record_iterator.get_alignment_length()
	print maxlen
	print 'columns in alignment'

	for i in range(0,maxlen):
		columns.append([])
		Mats.append([])

	for record in record_iterator:
		refs.append(record.id)
		#generate an array for each column
		for i,letter in enumerate(record.seq):
			columns[i].append(letter)


	for column in columns:
		params.append((column,propdict))
	
	pool = Pool()
	results = pool.map_async(column_to_mat, params, 25 ).get()
	for i,mat in enumerate(results):
		Mats[i] = mat


	pool.close()
	return refs,Mats

def makegaussian(a,l):
	#l is the length of the vector
	#a is alpha...
	gaussian = np.zeros((1,l))
	for index in range(l):
		x = index - l/2
		gaussian[0,index] = math.sqrt(a/math.pi) * math.exp(-a*x*x)
	return gaussian


def smooth_mats(mats):
	# aplly a convolution of gaussian vec and each mat representation of the columns
	for mat in mats:

	pass


def contact_map(args):
	subfasta1_aln, subfasta2_aln, propdict = args
	#perform mutual information analysis on two aligned fastas to generate a contact map.
	#uses protasvec and mutual information analysis based on physical properties.
	
	refs1,mts_1 = aln_to_Mats(subfasta1_aln,propdict)
	refs2,mts_2 = aln_to_Mats(subfasta2_aln,propdict)
	
	MI=np.zeros((len(mts_1),len(mts_2)))
	
	#prepare for multiporcessing run
	pool = Pool()

	selfMI = False
	if subfasta1_aln == subfasta2_aln:
		selfMI = True

	##calculate entropy for each column in 1 
	rundata = []
	for i,mat1 in enumerate(mts_1):
		rundata.append((mat1,3))
	results1 =  pool.map_async(mutualInfo.entropy,rundata, 25).get()
	
	rundata = []
	if selfMI == False:
		for i,mat2 in enumerate(mts_2):
			rundata.append((mat2,3))
		results2 = pool.map_async(mutualInfo.entropy,rundata, 25).get()
	else :
		results2 = results1

	rundata = []
	IJ = []

	print 'starting pairwise calculation'
	#divide this beast into chunks otherwise you'll run out of mem
	results3=[]
	rundata = []
	for i,mat1 in enumerate(mts_1):
		for j, mat2 in enumerate(mts_2):
			if i <= j:
				rundata.append((np.hstack((mat1,mat2)),3) )
				IJ.append((i,j))
				if len(rundata) > len(mts_1):
					results3 +=  pool.map_async(mutualInfo.entropy,rundata, 25).get()
					rundata = []
	#leftovers
	results3 +=  pool.map_async(mutualInfo.entropy,rundata, 25).get()
	avgMI = 0
	for k,pair in enumerate(IJ):
		i,j = pair
		MI[i,j]= results1[i] + results2[j] - results3[k] 
		avgMI += MI[i,j]
	#apply correction...Average Product Corection (APC), proposed by Dunn et al. (2008), that defines an APC, and subtracts this value from MI:
	#avgMI /= len(IJ)

	#for k, pair in enumerate(IJ):
	#	MI[i,j] -=  results1[i] * results2[j] / avgMI
	#symmetric MI matrix
	MI = MI.T + MI
	print 'MI calculation finished'
	return MI

def StartContactMap(fasta1,fasta2, MIname = None):
	start = time.time()
	propdict = init_props('physicalporpTable.csv')
	MI =contact_map((fasta1,fasta2,propdict))
	stop = time.time()
	print stop-start
	print 'seconds elapsed'
	if MIname != None:
		pickleutils.save_obj(MI, MIname)
	return MI
	#show global columnwise MI for the whole alignment

def showcontactmap(MIname):
	MI = pickleutils.load_obj(MIname)
	#similarityMatrix.plot_coo_matrix(MI)
	plt.imshow(MI)
	plt.colorbar()
	plt.show()

def addSecToGraph(G,seq):
	count = 0
	for i,resi in enumerate(seq.seq):
		if resi != '-':
			G.add_node(seq.id + '_'  + resi + str(count) , ID = seq.id , Number = count)
			#add connections for the peptide bonds
			if count > 0:
				G.add_edge(lastnode,seq.id + '_'  + resi + str(count), type = 'peptide')
			lastnode = seq.id + '_' + resi + str(count)
			count +=1

def gap_index(seq):
	index = []
	for i,resi in enumerate(seq.seq):
		if '-' not in resi:
			index.append(i)
	return index


def fetch_seq(subfasta1,representative1=None):
	record_iterator = AlignIO.read(subfasta1, "fasta")
	for record in record_iterator:
		if representative1 != None:
			if representative1 in record.id:
				print record.id
				seq1 = record
				print 'found ' + representative1
				break
		else:
			representative1 = record.id
			seq1 = record
	return representative1,seq1

def contact_graph( MI , cutoff ,  subfasta1,subfasta2 , representative1 = None, representative2=None , G = None):
	#verify that rep1 and rep2 are in the same organism before sending them to this function
	#if rep1 or rep2 is already in the graph then add the new prot to the already existing graph
	#sort the Mutual info scores and extract the top ones
	#take the two representative sequences from one organism. If this is a contact map for just one prot use the same code
	
	representative1,seq1 = fetch_seq(subfasta1,representative1)
	representative2,seq2 = fetch_seq(subfasta2,representative2)

	#fill an already existing graph with more prots
	if G != None:
		if seq1.id not in G.graph['prots']:	
			G.graph['prots'].append(seq1.id)
			addSeqToGraph(G,seq1)
			gap_index1 = gap_index(seq1)
		if seq2.id not in G.graph['prots']:	
			G.graph['prots'].append(seq2.id)
			addSeqToGraph(G,seq2)
			gap_index2 = gap_index(seq2)

	else:
		if seq1.id != seq2.id:
			G = nx.Graph(prots =[seq1.id,seq2.id] )
			addSecToGraph(G,seq1)
			addSecToGraph(G,seq2)
			gap_index1 = gap_index(seq1)
			gap_index2 = gap_index(seq2)

		else:
			G = nx.Graph(prots = [seq1.id] )
			addSecToGraph(G,seq1)
			gap_index1 = gap_index(seq1)
			gap_index2 = gap_index(seq1)

	
	#grab th rows and columns corresponding to the representative sequences
	MIreduced = MI[gap_index1,:]
	MIreduced = MIreduced[:,gap_index2]

	#grab the Ncutoff best interactions
	originalshape = MIreduced.shape
	flatMI = MIreduced.reshape((1,MIreduced.shape[0]*MIreduced.shape[1]))
	#descending order
	order = np.fliplr(np.argsort(flatMI))
	order = order.reshape(originalshape)
	MIreduced = MIreduced.reshape(originalshape)
	indices = np.argwhere(order<cutoff)

	print indices
	names = G.nodes()
	#check each pair in the MImap
	for number,index in enumerate(indices):
		i,j = index
		if math.fabs(MIreduced[i,j]) > 0:
		 	nodeName1 = seq1.id + '_' + seq1.seq[gap_index1[i]] + str(i)
			nodeName2 = seq2.id + '_' + seq2.seq[gap_index2[j]] + str(j)
			if nodeName1 in names and nodeName2 in names:
				G.add_edge(nodeName1, nodeName2)
				G[nodeName1][nodeName2]['MI'] = MIreduced[i,j]
				G[nodeName1][nodeName2]['type'] = 'MI'	


	return G


def contact_visualize(G):

	maxline = 4

	#ranges for graph
	linethickness = np.linspace(0,maxline, 100)

	red = colour.Color("red")
	blue = colour.Color("blue")
	if  len(G.graph['prots']) > 1:
		nodecolors = list(red.range_to(blue, len(G.graph['prots'])))
	else:
		nodecolors = [red]

	edgecolors = list(blue.range_to(red, 100 ))

	#update position dictionary
	rMin = 10
	rIncr = 1

	nprots = len(G.graph['prots'])
	#pos will be radial
	pos = {}
	graph_colors = []

	for node in G.nodes(data= True):
		i = G.graph['prots'].index(node[1]['ID'])
		x = (rMin + rIncr * node[1]['Number'] ) * math.sin(i * 2*math.pi/nprots) 
		y = (rMin + rIncr * node[1]['Number'] ) * math.cos( i* 2* math.pi/nprots)
		pos[node[0]] = (x,y)
		graph_colors.append(nodecolors[i].rgb)

	edgemin = 100000000
	edgemax = -1000000000

	for edge in G.edges(data = True):
		if edge[2]['type'] == 'MI':
			if math.fabs(edge[2]['MI']) < edgemin:
				edgemin = math.fabs(edge[2]['MI'])
				
			if math.fabs(edge[2]['MI']) > edgemax:
				edgemax = math.fabs(edge[2]['MI'])


	f, ax = plt.subplots()
	nodeLines = nx.draw_networkx_nodes(G, pos, nodelist=None, node_size=1, node_color='r', 
	node_shape='o', alpha=1.0, cmap=None, vmin=None, vmax=None, ax=None, 
	linewidths=None, label=None)

	ax.add_collection(nodeLines)

	edgeLines=[]
	for edge in G.edges(data=True):
		if edge[2]['type'] == 'peptide':
			#prot backbone is a thick line
			thickness = maxline
			colorindex = 0
			alpha = 1
		else:
			#thickness and color proportional to MI
			#maybe add an ln() to this...
			thickness = maxline * math.fabs(edge[2]['MI']) / edgemax
			alpha = 1 * math.fabs(edge[2]['MI']) / edgemax
			colorindex = int( math.floor( 99 * (math.fabs(edge[2]['MI']) / edgemax) ) )

		#set position to nodes
		posA = pos[edge[0]]
		posB = pos[edge[1]]
		#1/2 dist between the two nodes...
		rad= .5*((posA[0]-posB[0])**2 + (posA[1]-posB[1])**2 )**.5  * .005
		
		color = edgecolors[colorindex]
		arrow = patches.FancyArrowPatch(posA=posA, posB=posB, path=None, arrowstyle='simple', arrow_transmuter=None
		, connector=None, patchA=None, patchB=None, shrinkA=2.0, shrinkB=2.0,
		mutation_scale=1.0, mutation_aspect=None, dpi_cor=1.0,
		lw = thickness , color = color.rgb , alpha = alpha  )

		arrow.set_connectionstyle(connectionstyle='arc3' , rad = rad)
		ax.add_patch(arrow)
		edgeLines.append(arrow)
	
	##plot the arrows and dots...###

	plt.show()

	nx.draw_circular(G)
	plt.show()

	return ax,pos
	#linear cartoons of each prot to visualize contact info

	#networkx: each node is an amino acid of . each 

#MI = StartContactMap('cip4_reduced.fasta','cip4_reduced.fasta', MIname = 'cip4_selfMI')



#make an MI matrix for each interactor with cip4
MI = pickleutils.load_obj('cip4_selfMI')
aln1 = 'cip4_reduced.fasta'
cutoff = 500
showcontactmap('cip4_selfMI')

G= contact_graph(MI,cutoff,aln1,aln1)
contact_visualize(G)


alignments = glob.glob(infolder + '*aln.fasta')
for aln2 in aligments:
	if aln2 != aln1:
		MI = StartContactMap(aln1,aln2) 
		pickleutils.save_obj(MI , aln1 + aln2 )

for i,aln2 in enumerate(alignments):
	MI = pickleutils.load_obj(aln1+aln2)
	if i == 0:
		#init graph
		G= contact_graph(MI,cutoff,aln1,aln2)
	else:
		#add to graph
		G = contact_graph(MI,cutoff,aln1,aln2,G)
pickleutils.save_obj(G, 'cip4_graph')

print 'done!'