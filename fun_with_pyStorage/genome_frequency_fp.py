import logging
import numpy as np
import io
import tarfile
from random import sample
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
from collections import Counter
import random
from storage.pyStorage import pyStorage
from Inspector import *

def lambda_handler(event, context):

	logging.info('Python HTTP trigger function processed a request.')
	inspector = Inspector()
	inspector.inspectAll()
	
	inspector.addTimeStamp("Time Stamp at start")

	bucket_name = event["output_buckets"][5]
	pyStorage.create_cred_file(event['aws_access_key_id'], event['aws_secret_key'], event['aws_session_token'], event['gcp_client_email'], event['gcp_private_key'], event['gcp_project_id'])

	#Expected failure probability in percent
	failure_prob = event['failure_prob']

	#introduce failure probability
	if random.random() < (failure_prob/100):
		raise Exception('failure')

	key_columnsfile = event['key_columnsfile']

	siftfile = event["sifted"]
	
	merged_result = event["merged_result"]
	
	c = 22
	POP_bucket = event['POP']
	POP = pyStorage.retrieve_file_name(POP_bucket)
	POP = POP.replace("input/","")
	#c = event['chromNr']


	n_runs = 1000
	n_indiv= 52


	font = {'family':'serif',
		'size':14   }
	plt.rc('font', **font)

	rd = ReadData()
	res = Results()
	wr = WriteData()
	pd = PlotData()
	
	histogram_overlapfile = './output_no_sift/Histogram_mutation_overlap_chr' + str(c) + '_sNO-SIFT_' + POP + '_'
	mutation_overlapfile = './output_no_sift/Mutation_overlap_chr' + str(c) + '_sNO-SIFT_' + POP + '_'
	mutation_index_array_file = './output_no_sift/mutation_index_array' + str(c) + '_sNO-SIFT_' + POP + '.txt'
	histogram_overlap_plot = './plots_no_sift/Frequency_mutations' + str(c) + '_sNO-SIFT_' + POP 
	map_variations_file = './output_no_sift/map_variations' + str(c) + '_sNO-SIFT_' + POP + '.txt'
	randomindiv_file = './output_no_sift/random_indiv' + str(c) + '_sNO-SIFT_' + POP + '_'


	ids = rd.read_names(POP, POP_bucket, key_columnsfile, inspector)
	n_pairs = len(ids)/2
	

	rs_numbers, map_variations = rd.read_rs_numbers(siftfile, inspector)
	mutation_index_array = rd.read_individuals(ids, rs_numbers, c, merged_result, inspector)

	# gen final output
	tmp_tar = "/tmp/%s-%s-archive.tar.gz"% (c, POP)
	tar = tarfile.open(tmp_tar, mode = "w:gz")

	wr.write_map_variations(map_variations_file, map_variations, tar)
	wr.write_mutation_index_array(mutation_index_array_file, mutation_index_array, tar)
	
	mutation_overlap, random_indiv= res.overlap_ind(ids, n_runs, n_indiv, mutation_index_array)

	histogram_overlap= res.histogram_overlap(n_runs, mutation_overlap)
	wr.write_mutation_overlap(n_runs, mutation_overlapfile, mutation_overlap, tar)
	wr.write_histogram_overlap(n_runs, n_indiv, histogram_overlapfile, histogram_overlap, tar)
	wr.write_random_indiv(n_runs, randomindiv_file, random_indiv, tar)
	
	pd.plot_histogram_overlap(POP, n_runs, histogram_overlap, histogram_overlap_plot, tar)

	tar.close()
	result_key='/output/chr%s-%s-freq.tar.gz' % (c, POP)
	
	inspector.addTimeStamp("Time Stamp before Upload")
	pyStorage.copy(tmp_tar, bucket_name+result_key)
	inspector.addTimeStamp("Time Stamp after Upload")
	
	inspector.addTimeStamp("Time Stamp at end")
	
	inspector.addAttribute("output", result_key)
	inspector.inspectAllDeltas()
	
	return inspector.finish() 


	

class ReadData :
	def read_names(self, POP,POP_bucket, key_columnsfile, inspector) :
		
		pop_tmp = "/tmp/"+POP
		
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


	def read_rs_numbers(self, siftfile, inspector) :
		rs_numbers = []
		variations = {}
		map_variations = {}
		all_variations = []
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
	
	def read_individuals(self, ids, rs_numbers, c , merged_result, inspector) :
		mutation_index_array = []
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
			text = f.read().decode("utf-8").split()
			sifted_mutations = list(set(rs_numbers).intersection(text))
			mutation_index_array.append(sifted_mutations)
		return mutation_index_array


class Results :

	def overlap_ind(self, ids, n_runs, n_indiv, mutation_index_array):
		n_p = len(mutation_index_array)
		list_p = np.linspace(0, n_p - 1, n_p).astype(int)
		mutation_overlap = []
		random_indiv = []
		for run in range(n_runs) :
			randomized_list = sample(list(list_p), n_p)
			result = Counter()
			r_ids=[]
			for pq in range(n_indiv):
				if 2*pq >= len(randomized_list):
					break
				if randomized_list[2*pq] >= len(ids):
					break
				b_multiset = collections.Counter(mutation_index_array[randomized_list[2*pq]])
				r_ids.append(ids[randomized_list[2*pq]])
				result = result + b_multiset
			random_indiv.append(r_ids)
			mutation_overlap.append(result)
		return mutation_overlap, random_indiv
	
	def histogram_overlap(self, n_runs, mutation_overlap):
		histogram_overlap= []
		for run in range(n_runs):
			final_counts = [count for item, count in mutation_overlap[run].items()]
			histogram_overlap.append(collections.Counter(final_counts))
		return histogram_overlap			

class PlotData :		

	def plot_histogram_overlap(self, POP, n_runs, histogram_overlap, outputFile, tar):
		for run in range(n_runs):
			output = outputFile + str(run) + '.png'
			final_counts = [count for item, count in histogram_overlap[run].items()]
			N = len( final_counts )
			x = range( N )
			width = 1/1.5
			bar1=plt.bar( x, final_counts, width, color="grey" )
			plt.ylabel( 'Mutations' )
			plt.xlabel('Individuals')
			plt.xticks( np.arange( 1,N+1 ) )
			string = io.BytesIO()
			plt.savefig(string)
			plt.close()
			info = tarfile.TarInfo(name=output)
			info.size=len(string.getbuffer())
			string.seek(0)
			tar.addfile(tarinfo=info, fileobj=string)
	

class WriteData :
	def saveText(self, filename, text, tar):
		string = io.BytesIO(text.encode())
		info = tarfile.TarInfo(name=filename)
		info.size=len(string.getbuffer())
		tar.addfile(tarinfo=info, fileobj=string)

	def write_histogram_overlap(self, n_runs, n_indiv, histogram_overlapfile, histogram_overlap, tar) :	
		for run in range(n_runs):
			overlapfile = histogram_overlapfile + str(run) + '.txt'
			s='Number Individuals - Number Mutations  \n'
			for i in range(1,n_indiv+1):
				if i in histogram_overlap[run]:
					s=s+str(i) + '-' + str(histogram_overlap[run][i]) + '\n'
				else:
					s=s+str(i) + '-' + str(0) + '\n'
			self.saveText(overlapfile, s, tar)
				
	def write_mutation_overlap(self, n_runs, mutation_overlapfile, mutation_overlap, tar) :	
		for run in range(n_runs):
			overlapfile = mutation_overlapfile + str(run) + '.txt'
			s="Mutation Index- Number Overlapings \n"
			for key, count in mutation_overlap[run].items() :
				s=s+key + '-' + str(count) + '\n'
			self.saveText(overlapfile, s, tar)
	
	def write_random_indiv(self, n_runs, randomindiv_file, random_indiv, tar) :
		for run in range(n_runs):
			randomfile = randomindiv_file + str(run) + '.txt'
			text= '\n'.join([str(item) for item in random_indiv[run]])
			self.saveText(randomfile, text, tar)
	
	def write_mutation_index_array(self, mutation_index_array_file, mutation_index_array, tar):
		text = '\n'.join(str(x) for x in mutation_index_array)
		self.saveText(mutation_index_array_file, text, tar)
	
	def write_map_variations(self, map_variations_file, map_variations, tar) :
		s=""
		for key, count in map_variations.items() :
			s=s+key + '\t' + str(count) + '\n'
		self.saveText(map_variations_file, s, tar)


