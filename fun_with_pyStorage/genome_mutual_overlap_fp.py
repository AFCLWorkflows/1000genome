import logging
from random import sample
import io
import tarfile
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import random
from Inspector import *

from storage.pyStorage import pyStorage

def lambda_handler(event, context):
	
	logging.info('Python HTTP trigger function processed a request.')
	inspector = Inspector()
	inspector.inspectAll()
	
	inspector.addAttribute("function_name", "genome_mutual_overlap_nv")
	inspector.addTimeStamp("Time Stamp at start")
	
	pyStorage.create_cred_file(event['aws_access_key_id'], event['aws_secret_key'], event['aws_session_token'], event['gcp_client_email'], event['gcp_private_key'], event['gcp_project_id'])
	
	
	bucket_name = event["output_buckets"][4]
	
	failure_prob = event['failure_prob']
	
	siftfile = event["sifted"]
	
	key_columnsfile = event['key_columnsfile']
	
	merged_result = event['merged_result']
	
	if random.random() < (failure_prob/100):
		raise Exception('failure')
	
	
	c = 22 
	
	POP_bucket = event['POP']
	POP = pyStorage.retrieve_file_name(POP_bucket)
	POP = POP.replace("input/","")
	
	n_runs=1000
	
	
	font = {'family':'serif', 'size':14   }
	plt.rc('font', **font)
	
	rd = ReadData()
	res = Results()
	wr = WriteData()
	pd = PlotData()
	
	half_indpairsfile = './output_no_sift/individual_half_pairs_overlap_chr' + str(c) + '_sNO-SIFT_' + POP + '.txt'
	total_indpairsfile = './output_no_sift/total_individual_pairs_overlap_chr' + str(c) + '_sNO-SIFT_' + POP + '.txt'
	genepairsfile = './output_no_sift/gene_pairs_count_chr' + str(c) + '_sNO-SIFT_' + POP + '.txt'
	random_indpairsfile = './output_no_sift/100_individual_overlap_chr' + str(c) + '_sNO-SIFT_' + POP + '.txt'
	
	colormap = './plots_no_sift/colormap_distribution_c' + str(c) + '_sNO-SIFT_' + POP + '.png'
	half_overlap = './plots_no_sift/half_distribution_c' + str(c) + '_sNO-SIFT_' + POP + '.png'
	total_overlap = './plots_no_sift/total_distribution_c' + str(c) + '_sNO-SIFT_' + POP + '.png'
	random_overlap = './plots_no_sift/100_distribution_c' + str(c) + '_sNO-SIFT_' + POP + '.png'
	
	total_mutations_filename = './output_no_sift/total_mutations_individual' + str(c) + '_sNO-SIFT' + '_' + POP + '.txt'
	random_mutations_filename = './output_no_sift/random_mutations_individual' + str(c) + '_sNO-SIFT' + '_' + POP 
	
	mutation_index_array_file = './output_no_sift/mutation_index_array' + str(c) + '_sNO-SIFT_' + POP + '.txt'
	
	map_variations_file = './output_no_sift/map_variations' + str(c) + '_sNO-SIFT_' + POP + '.txt'
	
	
	print("POP: " + POP)
	print("POP_bucket: " + POP_bucket)
	ids = rd.read_names(POP, POP_bucket, key_columnsfile, inspector)
	
	rs_numbers, map_variations = rd.read_rs_numbers(siftfile, inspector)
	mutation_index_array, total_mutations, total_mutations_list = rd.read_individuals(ids, c, rs_numbers, merged_result, inspector)
	
	
	tmp_tar = "/tmp/archive.tar.gz"
	tar = tarfile.open(tmp_tar,mode = "w:gz")
	
	wr.write_total_indiv(total_mutations_filename, total_mutations,tar)
	wr.write_total_indiv(map_variations_file, map_variations, tar)
	
	#cross-correlations mutations overlapping
	half_pairs_overlap = res.half_pair_individuals(mutation_index_array)
	total_pairs_overlap, simetric_overlap = res.total_pair_individuals(mutation_index_array)
	random_pairs_overlap = res.pair_individuals(mutation_index_array, n_runs)
	
	wr.write_mutation_index_array(mutation_index_array_file, mutation_index_array,tar)
	wr.write_pair_individuals(half_indpairsfile, half_pairs_overlap, tar)
	wr.write_pair_individuals(total_indpairsfile, total_pairs_overlap, tar)
	wr.write_pair_individuals(random_indpairsfile, random_pairs_overlap, tar)
	
	pd.individual_overlap(POP, c, n_runs, half_pairs_overlap, half_overlap, tar)
	pd.individual_overlap(POP, c, n_runs, simetric_overlap, total_overlap, tar)
	pd.individual_overlap(POP, c, n_runs, random_pairs_overlap, random_overlap, tar)
	pd.total_colormap_overlap(POP, total_pairs_overlap, colormap, tar)
	
	#list of frecuency of mutations in 26 individuals
	random_mutations_list=res.group_indivuals(total_mutations_list, n_runs)
	wr.write_random_mutations_list(n_runs, random_mutations_filename, random_mutations_list,tar)
	
	# gen overlapping
	gene_pair_list = res.gene_pairs(mutation_index_array)
	wr.write_total_indiv(genepairsfile, gene_pair_list,tar)
	
	tar.close()
	result_key= bucket_name+ 'output/chr%s-%s.tar.gz' % (c, POP)
	
	inspector.addTimeStamp("Time Stamp before Upload")
	pyStorage.copy(tmp_tar, result_key)
	
	inspector.addTimeStamp("Time Stamp after Upload")
	
	inspector.addTimeStamp("Time Stamp at end")
	inspector.addAttribute("output", result_key)
	inspector.addAttribute("aws_access_key_id", event['aws_access_key_id'])
	inspector.addAttribute("aws_secret_key", event['aws_secret_key'])
	inspector.addAttribute("aws_session_token", event['aws_session_token'])
	inspector.addAttribute("gcp_client_email", event['gcp_client_email'])
	inspector.addAttribute("gcp_private_key", event['gcp_private_key'])
	inspector.addAttribute("gcp_project_id", event['gcp_project_id'])
	
	inspector.inspectAllDeltas()
	return inspector.finish() 
	
	
class ReadData :
#reads POP + columns makes set
	def read_names(self, POP, POP_bucket, key_columnsfile, inspector) :
		pop_tmp = "/tmp/"+ POP
		
		inspector.addTimeStamp("Time Stamp before Download 1")
		pyStorage.copy(POP_bucket, pop_tmp)
		inspector.addTimeStamp("Time Stamp after Download 1")
		
		object = open(pop_tmp, "r").read()
		bytes(object, "utf-8")
		text = object.split()
		all_ids = text[0:]
		
		columns_tmp = "/tmp/columns.txt"
		
		inspector.addTimeStamp("Time Stamp before Download 2")
		pyStorage.copy(key_columnsfile, columns_tmp)
		inspector.addTimeStamp("Time Stamp after Download 2")
		
		object = open(columns_tmp, "r").read()
		bytes(object, "utf-8")
		genome_ids = object.split()
		
		ids = list(set(all_ids) & set(genome_ids))
		return ids

#reads siftfile
	def read_rs_numbers(self, siftfile, inspector) :
		rs_numbers = []
		map_variations = {}
		siftfile_tmp = "/" + pyStorage.retrieve_file_name(siftfile)
		
		inspector.addTimeStamp("Time Stamp before Download 3")
		pyStorage.copy(siftfile, siftfile_tmp)
		inspector.addTimeStamp("Time Stamp after Download 3")
		
		object = open(siftfile_tmp, "r").read()
		bytes(object, "utf-8")
		sift_file = object.split("\n")
		for item in sift_file:
			item = item.split()
			if len(item) > 2:
				rs_numbers.append(item[1])
				map_variations[item[1]] = item[2]
		print(rs_numbers)
		return rs_numbers, map_variations
	
#reads different cromes
	def read_individuals(self, ids, c, rs_numbers, merged_result, inspector) :
		mutation_index_array = []
		total_mutations={}  
		total_mutations_list =[]
		
		tmp_tar = "/tmp/file-%s.tar.gz" % (c)
		
		inspector.addTimeStamp("Time Stamp before Download 4")
		pyStorage.copy(merged_result, tmp_tar)
		inspector.addTimeStamp("Time Stamp after Download 4")
		
		object = open(tmp_tar, "rb").read()
		bytes = io.BytesIO(object)
		
		inputtar = tarfile.open(mode = "r:gz", fileobj = bytes)

		for member in inputtar.getmembers():
	
			f = inputtar.extractfile(member)
			if f==None:
				continue
			for name in ids:
				if name in member.name:
					break;
	
			text = f.read().decode("utf-8").split()
			
			sifted_mutations = list(set(rs_numbers).intersection(text))
			mutation_index_array.append(sifted_mutations)
			total_mutations[name]= len(sifted_mutations)
			total_mutations_list.append(len(sifted_mutations))
		return mutation_index_array, total_mutations, total_mutations_list

class Results :

	def group_indivuals(self, total_mutations_list, n_runs) :
		n_group = 26
		random_mutations_list= []
		for run in range(n_runs):
			random_mutations_list.append(sample(total_mutations_list, n_group))
		return random_mutations_list

	def pair_individuals(self, mutation_index_array, n_runs) :
		n_p = len(mutation_index_array)
		n_pairs = int(round(n_p/2))
		list_p = np.linspace(0, n_p - 1, n_p).astype(int)
		pairs_overlap = np.zeros((n_runs, n_pairs))
		for run in range(n_runs) :
			randomized_list = sample(list(list_p) , n_p)
			for pq in range(n_pairs) :
				array1 = mutation_index_array[randomized_list[2*pq]]
				array2 = mutation_index_array[randomized_list[2*pq]]
				pair_array = set(array1) & set(array2)
				pairs_overlap[run][pq] = len(pair_array)

		return pairs_overlap

	def total_pair_individuals (self, mutation_index_array) :
		n_p = len(mutation_index_array)
		total_pairs_overlap = np.zeros((n_p, n_p))
		simetric_overlap = np.zeros((n_p, n_p))
		for run in range(n_p):
						array1 = mutation_index_array[run]
						start = run +1
						for pq in range(start, n_p) :
								array2 = mutation_index_array[pq]
								pairs_array = set(array1) & set(array2)
								total_pairs_overlap[run][pq]=len(pairs_array)
								simetric_overlap[run][pq] = len(pairs_array)
								simetric_overlap[pq][run]= len(pairs_array)

		return total_pairs_overlap , simetric_overlap

	def half_pair_individuals(self, mutation_index_array) :
		n_p = len(mutation_index_array)
		n_pairs = int(round(n_p/2))
		pairs_overlap = np.zeros((n_pairs, n_pairs))
		for run in range(n_pairs):
			array1 = mutation_index_array[run]
			index =0
			for pq in range(n_pairs+1, n_p):
				array2 = mutation_index_array[pq]
				pairs_array = set(array1) & set(array2)
				pairs_overlap[run][index]=len(pairs_array)

		return pairs_overlap

	def gene_pairs(self, mutation_index_array):
		n_p = len(mutation_index_array)
		gene_pair_list = {}
		for pp in range(n_p) :  
			pairs = itertools.combinations(mutation_index_array[pp], 2)
			for pair in pairs :
				key = str(pair)
				if key not in gene_pair_list : gene_pair_list[key] = 1
				else : gene_pair_list[key] += 1

		return gene_pair_list

class PlotData :		

	def individual_overlap(self, POP, c, n_runs, pairs_overlap, outputFile, tar):
		print('plotting cross matched number of individuals:%s '% len(pairs_overlap))
		pairs_overlap = np.array(pairs_overlap)	 

		min_p = np.min(pairs_overlap)
		max_p = np.max(pairs_overlap)
		nbins = int(max_p) + 1
		n_runs = len(pairs_overlap)

		nbins = int(np.max(pairs_overlap))
		bin_centres = np.linspace(0, nbins, nbins)
		bin_edges = np.linspace(-0.5, nbins + 0.5, nbins + 1)

		fig = plt.figure(frameon=False, figsize=(10, 9))
		ax = fig.add_subplot(111)
		hists = []
		max_h = 0
		for run in range(n_runs) :
			h, edges = np.histogram(pairs_overlap[run], bins = bin_edges)
			ax.plot(bin_centres, h, alpha = 0.5)
			if len(h) > 0:
				max_h = max(max_h, max(h))

		plt.xlabel('Number of overlapping gene mutations', fontsize = 24)
		plt.ylabel(r'frequency', fontsize = 28)
		text1 = 'population ' + POP + '\n chromosome ' + str(c) + '\n SIFT < NO-SIFT \n' + str(n_runs) + ' runs'
		plt.text(.95, .95, text1, fontsize = 24, verticalalignment='top', horizontalalignment='right', transform = ax.transAxes)
		string = io.BytesIO()
		plt.savefig(string)
		plt.close()
		info = tarfile.TarInfo(name=outputFile)
		info.size=len(string.getbuffer())
		string.seek(0)
		tar.addfile(tarinfo=info, fileobj=string)

	def total_colormap_overlap(self, POP, total_pairs_overlap, outputFile, tar):
		print('plotting colormap number of individuals: %s' % len(total_pairs_overlap))
		fig = plt.figure()
		cmap = colors.ListedColormap(['blue','black','red', 'green', 'pink'])
		img = plt.imshow(total_pairs_overlap,interpolation='nearest', cmap = cmap, origin='lower')
		plt.colorbar(img,cmap=cmap)
		string = io.BytesIO()
		plt.savefig(string)  
		plt.close()
		info = tarfile.TarInfo(name=outputFile)
		info.size=len(string.getbuffer())
		string.seek(0)
		tar.addfile(tarinfo=info, fileobj=string)


class WriteData :
	def saveText(self, filename, text, tar):
		string = io.BytesIO(text.encode())
		info = tarfile.TarInfo(name=filename)
		info.size=len(string.getbuffer())
		tar.addfile(tarinfo=info, fileobj=string)

	def write_pair_individuals(self, indpairsfile, pairs_overlap, tar) : 
		string = io.BytesIO()
		np.savetxt(string, pairs_overlap, fmt = '%i')
		info = tarfile.TarInfo(name=indpairsfile)
		info.size=len(string.getbuffer())
		string.seek(0)
		tar.addfile(tarinfo=info, fileobj=string)

	def write_total_indiv(self, total_mutations_filename, total_mutations, tar) :
		s=""
		for key, count in total_mutations.items() :
			s=s+key + '\t' + str(count) + '\n'
		self.saveText(total_mutations_filename, s, tar)
	
	def write_random_mutations_list(self, n_runs, random_mutations_filename, random_mutations_list,tar) :
		for run in range(n_runs):
			filename= random_mutations_filename +'_run_' + str(run) + '.txt'
			text= '\n'.join([str(item) for item in random_mutations_list[run]])
			self.saveText(filename, text, tar)
	
	def write_mutation_index_array(self, mutation_index_array_file, mutation_index_array, tar):
		text = '\n'.join(str(x) for x in mutation_index_array)
		self.saveText(mutation_index_array_file, text, tar)