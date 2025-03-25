
from grakel import GraphKernel
from grakel import Graph
import grakel
import networkx as nx
from rdkit import Chem
import random
import time
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go






def get_grakel_graph_from_networkx(graph):
	edges = nx.get_edge_attributes(graph, "bond_type")
	node_labels = nx.get_node_attributes(graph, "atom_type")
	

	G = Graph(edges,edge_labels=edges, node_labels = node_labels)


	return G





def get_networkx_graph_from_smiles(smiles):
	
	mol = Chem.MolFromSmiles(smiles)

	

	edges = {}
	node_attributes = {}

	for bond in mol.GetBonds():
		src = bond.GetBeginAtomIdx()
		tar = bond.GetEndAtomIdx()
		bond_type = bond.GetBondType()
		
		node_attributes[src] = bond.GetBeginAtom().GetSymbol()
		node_attributes[tar] = bond.GetEndAtom().GetSymbol()
		edges[(src, tar)] = str(bond_type)


	graph = nx.Graph()

	for key in node_attributes.keys():
		graph.add_node(key, atom_type = node_attributes[key])

	for key in edges.keys():
		graph.add_edge(key[0], key[1], bond_type = edges[key])



	return graph




def greedy_prefix_ordering(graphs, kernel):

	grakel_graphs = [get_grakel_graph_from_networkx(get_networkx_graph_from_smiles( i)) for i in graphs]

	grakel_matrix = kernel.fit_transform(grakel_graphs)
	

	min_coord = (0,0)
	min = float("inf")

	for i in range(len(graphs)):
		for j in range(len(graphs)):
			if i<=j:
				continue
			if grakel_matrix[i,j]<min:
				min = grakel_matrix[i,j]
				min_coord = (i,j)

	

	min_i = float("inf")
	min_j = float("inf")

	for k in range(len(graphs)):
		if k==min_coord[0] or k==min_coord[1]:
			continue
		if grakel_matrix[min_coord[0],k]<min_i:
			min_i = grakel_matrix[min_coord[0],k]
		if grakel_matrix[min_coord[1],k]<min_j:
			min_j = grakel_matrix[min_coord[0],k]

	ordering = []

	if min_i<min_j:
		ordering = [min_coord[0], min_coord[1]]
	else:
		ordering = [min_coord[1], min_coord[0]]

	while len(ordering)<len(graphs):
		min = float("inf")
		min_coord = 0
		for i in range(len(graphs)):
			if i in ordering:
				continue

			sum = 0
			for j in ordering:
				#sum += grakel_matrix[j,i]
				if sum<grakel_matrix[j,i]:
					sum = grakel_matrix[j,i]
			if sum<min:
				min = sum
				min_coord = i
		ordering.append(min_coord)

	
	return ordering





def read_chembl():
	iter = 0
	f = open("chembl_22_clean_1576904_sorted_std_final.smi")
	lines = f.readlines() 
	lines = lines[1:]
	#random.seed(10)
	#random.shuffle(lines)
	smiles_list = []
	for line in lines[1:]:
		smiles = line.split("	")[0]
		smiles_list.append(smiles)

	f.close()
	return smiles_list




def create_sized_instances(kernel):

	smiles_list = read_chembl()

	random.seed(13)

	random.shuffle(smiles_list)

	num_graphs = 40

	done = 0

	sizes = [i for i in range(10, 150)]

	size_lists = {}

	for i in sizes:
		size_lists[i] = []

	i = 0
	while done<len(sizes):
		print(done)
		i += 1
		smi = smiles_list[i]
		mol = Chem.MolFromSmiles(smi)
		if not mol:
			continue
		if mol.GetNumBonds() in sizes:
			if len(size_lists[mol.GetNumBonds()])==num_graphs:
				continue
			elif len(size_lists[mol.GetNumBonds()])==num_graphs-1:
				done += 1
				size_lists[mol.GetNumBonds()].append(smi)
			else:
				size_lists[mol.GetNumBonds()].append(smi)

	

	instances = []

	for i in range(len(sizes)):
		instances.append(size_lists[i+10])

	file = open("instances_wl3.txt", "w")

	for smiles_strings in instances:

		order = greedy_prefix_ordering(smiles_strings, kernel)
		#order = range(len(smiles_strings))
		for i in order:

			file.write(smiles_strings[i] + "\n")

		file.write("#\n")






def create_instances(kernel):


	


	file = open("instances_wl3.txt", "w")

	smiles_list = read_chembl()

	good_smiles = []

	for smi in smiles_list[1200000:]:#1200000
		try:
			mol = Chem.MolFromSmiles(smi)
			print(mol.GetNumAtoms())
			if mol.GetNumAtoms() >= 35 and mol.GetNumAtoms()<=35:#35
				Chem.Kekulize(mol)
				good_smiles.append(Chem.MolToSmiles(mol, kekuleSmiles = True))

		except:
			continue
		print(len(good_smiles))
		if len(good_smiles)>7500:
			break

	num_graphs = 5

	instances = []

	random.seed(19)#random seed 18

	for i in range(5000):#5000
		instances.append(random.sample(good_smiles, num_graphs))

	start_time = time.time()

	for inst in instances:

		order = greedy_prefix_ordering(inst, kernel)

		smiles_strings = []

		for i in order:
			smiles_strings.append(inst[i])

		for smi in smiles_strings:

			file.write(smi + "\n")

		file.write("#\n")

	print(f"TIME: {time.time() - start_time}")




def create_csv_file():

	file = open("instances_before_ordering.csv", "w")

	smiles_list = read_chembl()

	good_smiles = []

	for smi in smiles_list[1200000:]:#1200000
		try:
			mol = Chem.MolFromSmiles(smi)
			print(mol.GetNumAtoms())
			if mol.GetNumAtoms() >= 35 and mol.GetNumAtoms()<=35:#35
				Chem.Kekulize(mol)
				good_smiles.append(Chem.MolToSmiles(mol, kekuleSmiles = True))

		except:
			continue
		print(len(good_smiles))
		if len(good_smiles)>7500:
			break

	num_graphs = 5

	instances = []

	random.seed(19)#random seed 18

	for i in range(5000):#5000
		instances.append(random.sample(good_smiles, num_graphs))

	start_time = time.time()

	for inst in instances:

		order = greedy_prefix_ordering(inst, kernel)

		smiles_strings = []

		for i in order:
			smiles_strings.append(inst[i])

		for smi in smiles_strings:

			file.write(smi + ",")

		file.write("\n")


def order_csv_file(kernel, file, output_file):

	f = open(file)

	lines = f.readlines()

	instances = []

	for line in lines:
		spl = line.split(",")
		instance = [i for i in spl if len(i)>1]
		instances.append(instance)

	file = open(output_file, "w")
	for inst in instances:

		order = greedy_prefix_ordering(inst, kernel)

		smiles_strings = []

		for i in order:
			smiles_strings.append(inst[i])

		for smi in smiles_strings:

			file.write(smi + ",")

		file.write("\n")



def visualize_runtime():
	x = ["VH/normal", "VH/pruned", "minmax/normal", "minmax/pruned", "WL/normal", "WL/pruned", "NSPD/normal", "NSPD/pruned"]
	ordering_time = [15.534, 15.534, 374.811, 374.811, 26.1809, 26.1809, 163.7517 , 163.7517]
	time = [7227.42 + 15.534, 2050.8 + 15.534, 1663.469 + 374.811, 899.779 + 374.811, 8718.75 +26.1809, 2501.13 +26.1809, 6700.92 + 163.7517, 1907.18 + 163.7517]


	data = {'Heuristics': ["VH/normal","VH/normal", "VH/pruned", "VH/pruned", "minmax/normal", "minmax/normal", "minmax/pruned","minmax/pruned", "WL/normal", "WL/normal", "WL/pruned", "WL/pruned", "NSPD/normal", "NSPD/normal", "NSPD/pruned", "NSPD/pruned"],
        'Legend': ["MCS time","Ordering time", "MCS time","Ordering time", "MCS time","Ordering time", "MCS time","Ordering time", "MCS time","Ordering time", "MCS time","Ordering time", "MCS time","Ordering time", "MCS time", "Ordering time"],
        'Time in seconds': [7227.42, 15.534, 2050.8, 15.534, 1663.469, 374.811, 899.779, 374.811, 8718.75, 26.1809, 2501.13, 26.1809, 6700.92, 163.7517, 1907.18 , 163.7517]}

	fig = px.bar(data, x="Heuristics", y="Time in seconds", color="Legend", title="Runtime analysis")
	fig.update_layout(
	    title=dict(text="Runtime Analysis", font=dict(size=50)),  # Title size
	    xaxis=dict(title=dict(text="Heuristics", font=dict(size=36)), tickfont=dict(size=28)),  # X-axis label size
	    yaxis=dict(title=dict(text="Time in seconds", font=dict(size=36)), tickfont=dict(size=28)),  # Y-axis label size
	    legend=dict(font=dict(size=28)),  # Legend text size
	    
	)
	fig.show()


def visualize_runtime_sizes():
	file = open("time_sizes.txt")

	lines = file.readlines()

	sizes = [i for i in range(10, 150)]

	times = []

	for line in lines[1:]:
		time = float(line.split(",")[0])
		times.append(time)


	plt.bar(sizes, times)
	plt.show()

"""
def visualize_runtime_hard():
	
	x = ["minmax/normal", "minmax/pruned"]
	ordering_time = [547.546,547.546]
	time = [?,3077.714]


	data = {'Heuristics': ["minmax/normal", "minmax/normal", "minmax/pruned","minmax/pruned"],
        'Legend': ["MCS time","Ordering time", "MCS time","Ordering time"],
        'Time in seconds': [7227.42, 15.534, 2050.8, 15.534, 1663.469, 374.811, 899.779, 374.811, 8718.75, 26.1809, 2501.13, 26.1809, 6700.92, 163.7517, 1907.18 , 163.7517]}

	fig = px.bar(data, x="Heuristics", y="Time in seconds", color="Legend", title="Runtime analysis")
	fig.update_layout(
	    title=dict(text="Runtime Analysis", font=dict(size=50)),  # Title size
	    xaxis=dict(title=dict(text="Heuristics", font=dict(size=36)), tickfont=dict(size=28)),  # X-axis label size
	    yaxis=dict(title=dict(text="Time in seconds", font=dict(size=36)), tickfont=dict(size=28)),  # Y-axis label size
	    legend=dict(font=dict(size=28)),  # Legend text size
	    
	)
	fig.show()

	
"""
import numpy as np
def get_q1_q3(data):
    q1 = np.percentile(data, 25)  # First quartile (Q1)
    q3 = np.percentile(data, 75)  # Third quartile (Q3)
    return q1, q3


def box_plots():
	# Example data
	data1 = [float(i) for i in open("times_vh_normal.txt").readlines()]
	data2 = [float(i) for i in open("times_vh_pruned.txt").readlines()]
	data3 = [float(i) for i in open("times_wl_normal.txt").readlines()]
	data4 = [float(i) for i in open("times_wl_pruned.txt").readlines()]
	data5 = [float(i) for i in open("times_nspd_normal.txt").readlines()]
	data6 = [float(i) for i in open("times_nspd_pruned.txt").readlines()]
	data7 = [float(i) for i in open("times_minmax_normal.txt").readlines()]
	data8 = [float(i) for i in open("times_minmax_pruned.txt").readlines()]

	data = [data1, data2, data3, data4, data5, data6, data7, data8]
	#fig, ax = plt.subplots()

	#bp = ax.boxplot(data)

	#plt.boxplot(data, whis=0, showfliers=False, showmeans=True, meanline=True)
	#plt.show()


	# Create figure and axis
	# Create figure and axis
	fig, ax = plt.subplots(figsize=(10, 6))

	# Create a list of 8 distinct colors
	colors = [
	    'lightblue', 'lightgreen', 'lightcoral', 'lightyellow', 
	    'lightpink', 'lightgray', 'lightskyblue', 'lightseagreen'
	]

	# Create a box plot with custom settings
	bp = ax.boxplot(data, whis=0, showfliers=False, showmeans=True, meanline=True, patch_artist=True)

	# Apply custom colors to each box
	for i, box in enumerate(bp['boxes']):
	    box.set_facecolor(colors[i])  # Set each box with a different color

	# Customize whiskers (keep at Q1/Q3)
	for whisker in bp['whiskers']:
	    whisker.set(color='gray', linewidth=1.5)

	# Customize medians (set line style)
	for median in bp['medians']:
	    median.set(color='black', linewidth=2)

	# Customize mean lines (bolder and different color)
	for mean in bp['means']:
	    mean.set(color='black', linestyle='dotted', linewidth=6)  # Bolder and black

	    # Extend the mean line
	    # Get the position of the mean line
	    y = mean.get_ydata()
	    x = mean.get_xdata()
	    
	    # Extend the mean line beyond the box plot range
	    mean.set_xdata([x[0] - 0.2, x[-1] + 0.2])  # Extend both sides of the line
	    mean.set_ydata([y[0], y[0]])  # Keep the mean constant across the line


	# Customize caps (hide if you want to remove them)
	for cap in bp['caps']:
	    cap.set_visible(False)

	# Set y-axis to logarithmic scale
	plt.ylabel("Time in seconds", fontsize=50, labelpad=40)

	# Customize y-axis tick labels to display floats with a specific format
	#ax.yaxis.set_major_formatter(ScalarFormatter())
	#ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.2f}'))  # Format as floats with 2 decimal places

	# Set title and labels
	ax.set_title("Runtime analysis", fontsize=60)
	
	plt.xlabel("Heuristics", fontsize=50, labelpad=40)
	ax.set_xticklabels(["VH/not pruned", "VH/pruned", "WL/not pruned", "WL/pruned", "NSPD/not pruned", "NSPD/pruned", "minmax/not pruned", "minmax/pruned"], fontsize=24)

	# Add gridlines for clarity
	ax.grid(True, axis='y', linestyle='--', alpha=0.7)

	ax.tick_params(axis='y', labelsize=40)
	ax.tick_params(axis='x', labelsize=40)
	plt.rcParams["figure.figsize"] = (20,3)

	plt.xticks(rotation=15)

	# Show the plot
	plt.show()


	"""
	q1_1, q3_1 = get_q1_q3(data1)
	q1_2, q3_2 = get_q1_q3(data2)
	q1_3, q3_3 = get_q1_q3(data3)
	q1_4, q3_4 = get_q1_q3(data4)
	q1_5, q3_5 = get_q1_q3(data5)
	q1_6, q3_6 = get_q1_q3(data6)
	q1_7, q3_7 = get_q1_q3(data7)
	q1_8, q3_8 = get_q1_q3(data8)
	

	fig = go.Figure()

	fig.add_trace(go.Box(y=data1, name="Set 1", boxmean=True, boxpoints=False, lowerfence=[q1_1], upperfence=[q3_1]))
	fig.add_trace(go.Box(y=data2, name="Set 2", boxmean=True, boxpoints=False, q1=[q1_2], q3=[q3_2], lowerfence=[q1_2], upperfence=[q3_2]))
	fig.add_trace(go.Box(y=data3, name="Set 3", boxmean=True, boxpoints=False, q1=[q1_3], q3=[q3_3], lowerfence=[q1_3], upperfence=[q3_3]))
	fig.add_trace(go.Box(y=data4, name="Set 4", boxmean=True, boxpoints=False, q1=[q1_4], q3=[q3_4], lowerfence=[q1_4], upperfence=[q3_4]))
	fig.add_trace(go.Box(y=data5, name="Set 5", boxmean=True, boxpoints=False, q1=[q1_5], q3=[q3_5], lowerfence=[q1_5], upperfence=[q3_5]))
	fig.add_trace(go.Box(y=data6, name="Set 6", boxmean=True, boxpoints=False, q1=[q1_6], q3=[q3_6], lowerfence=[q1_6], upperfence=[q3_6]))
	fig.add_trace(go.Box(y=data7, name="Set 7", boxmean=True, boxpoints=False, q1=[q1_7], q3=[q3_7], lowerfence=[q1_7], upperfence=[q3_7]))
	fig.add_trace(go.Box(y=data8, name="Set 8", boxmean=True, boxpoints=False, q1=[q1_8], q3=[q3_8], lowerfence=[q1_8], upperfence=[q3_8]))

	#fig.update_layout(yaxis_type="log", title="Box Plots with Logarithmic Scale")

	# Show the plot
	fig.show()
	"""


import sys
if __name__ == '__main__':
	kernel1 = grakel.VertexHistogram(n_jobs=None, normalize=True, verbose=False, sparse=False)
	kernel5 = grakel.NeighborhoodSubgraphPairwiseDistance(n_jobs=None, normalize=False, verbose=False, r=2, d=3)
	kernel8 = grakel.WeisfeilerLehmanOptimalAssignment(n_jobs=None, verbose=False, normalize=True, n_iter=4, sparse=False)

	

	#create_instances(kernel1)
	#create_sized_instances(kernel5)
	#visualize_runtime()
	#box_plots()
	#create_csv_file()	
	#order_csv_file(kernel1, "instances_before_ordering.csv", "ordered_instances.csv")

	args = sys.argv

	input = args[1]
	output = args[2]
	kernel_arg = int(args[3])

	if kernel_arg==1:
		kernel = grakel.VertexHistogram(n_jobs=None, normalize=True, verbose=False, sparse=False)
	elif kernel_arg==2:
		kernel = grakel.NeighborhoodSubgraphPairwiseDistance(n_jobs=None, normalize=False, verbose=False, r=2, d=3)
	elif kernel_arg==3:
		kernel = grakel.WeisfeilerLehmanOptimalAssignment(n_jobs=None, verbose=False, normalize=True, n_iter=4, sparse=False)
	else:
		print("Argument 3 specifying kernel type must be 1,2 or 3.")

	order_csv_file(kernel1, input, output)

