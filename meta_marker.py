#Translates DNA sequences in all 6 reading frames, ignoring start / stop codons.

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import sys
import tempfile
import subprocess
from collections import Counter
from scipy.spatial import distance
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import time

def translate_6frames(input_file, min_size):
	input_handle = open(input_file, "rU")	
	f = tempfile.NamedTemporaryFile(delete=False)
	for record in SeqIO.parse(input_handle, "fasta") :
		if len(record.seq) >= min_size:
			#Frame 1
			original = record.seq
			f.write(">" + str(record.id) + "_1\n")
			f.write(str(record.seq.translate()).replace("*","") + "\n")
			#Frame 2
			f.write(">" + str(record.id) + "_2\n")
			record.seq = Seq(str(record.seq)[1:])
			f.write(str(record.seq.translate()).replace("*","") + "\n")
			#Frame 3
			f.write(">" + str(record.id) + "_3\n")
			record.seq = Seq(str(record.seq)[1:])
			f.write(str(record.seq.translate()).replace("*","") + "\n")

			record.seq = original.reverse_complement()

			#Frame -1
			f.write(">" + str(record.id) + "_-1\n")
			f.write(str(record.seq.translate()).replace("*","") + "\n")
			#Frame -2
			record.seq = Seq(str(record.seq)[1:])
			f.write(">" + str(record.id) + "_-2\n")
			f.write(str(record.seq.translate()).replace("*","") + "\n")
			#Frame -3
			record.seq = Seq(str(record.seq)[1:])
			f.write(">" + str(record.id) + "_-3\n")
			f.write(str(record.seq.translate()).replace("*","") + "\n")
	return f

def calc_tetra(seq_record):

	tetramers = {}
	for a in ['A', 'C', 'G', 'T']:
		for b in ['A', 'C', 'G', 'T']:
			for c in ['A', 'C', 'G', 'T']:
				for d in ['A', 'C', 'G', 'T']:
					tetramers[a+b+c+d] = 0

	start = 0
	end = 4	
	for i in range(0,len(str(seq_record.seq))):
		if len(str(seq_record.seq[start:end])) == 4:
			try:
				tetramers[str(seq_record.seq[start:end])] += 1
			except:
				pass	
		start += 1
		end += 1
	
	#Normalize 
	total = sum(tetramers.values())
	for k in tetramers.keys():
		tetramers[k] = float(tetramers[k]) / float(total)

	return tetramers

if __name__ == '__main__':

	start = time.time()

	input_file = sys.argv[1]
	#Keep only contigs bigger than X bp
	min_size = 5000

	#Translate contigs in all 6 frames
	print "Translating DNA..."
	tempfile = translate_6frames(input_file, min_size)

	#Search for markers using hmm file
	print "Searching for marker proteins..."
	output = subprocess.Popen(["hmmsearch", "-A", tempfile.name + ".aa", 'ribosomal.hmm', tempfile.name], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)	
	stdout, stderr = output.communicate()

	#Write out markers
	print "Outputting marker proteins..."
	fname = '.'.join(sys.argv[0].split(".")[:-1]) + ".markers"
	f = open(fname, 'a+')
	for al in AlignIO.parse(open(tempfile.name + ".aa"), "stockholm"):
		for seq in al:
			f.write(">" + '/'.join(str(seq.id).split("/")[:-1]) + "\n")
			f.write(str(seq.seq).replace('-', '') + "\n")
	f.close()

	#BLAST markers against blast db
	print "Blasting marker proteins against reference DB..."
	output = subprocess.check_output(['blastp', '-query', fname, '-db', 'markers', '-outfmt', '6 qseqid stitle pident evalue', '-max_target_seqs', '1'])
	lines = str(output).splitlines()


	print "Outputting marker protein BLAST results..."
	fname = '.'.join(sys.argv[0].split(".")[:-1]) + ".blast"
	f = open(fname, 'a+')
	for line in lines:
		f.write(line + "\n")
	f.close()

	#Calculate most common species for a contig.
	contigs = {}
	for line in lines:
		contig = '_'.join(line.split("\t")[0].split("_")[:-1])
		try:
			genus = line.split("\t")[1].split("[")[1].split("]")[0].replace("Candidatus","").lstrip().strip().split()[0]
		except:
			genus = ''
		if contig not in contigs.keys():
			contigs[contig] = [genus]
		else:
			contigs[contig].append(genus)

	names = {}
	for contig in contigs:
		count = Counter(contigs[contig])
		names[contig] = str(count.most_common(1)[0][0])

	#Calculate tetranucleotide frequencies for all marked contigs
	print "Calculating tetranucleotide frequencies for marker contigs..."
	tetramers = {}
	sizes = {}
	input_handle = open(input_file, "rU")
	for record in SeqIO.parse(input_handle, "fasta") :
		if record.id in names and len(record.seq) >= min_size:
			tetramers[record.id] = calc_tetra(record)
			sizes[record.id] = len(record.seq)

	#Run PCA on that
	print "Using Scikitlearn to run PCA on 4mer distances..."

	tetramer_array = []
	for t in sorted(tetramers.keys()):
		temp = []
		for tet in sorted(tetramers[t].keys()):
			temp.append(tetramers[t][tet])
		tetramer_array.append(temp)

	tetramers_np = np.array(tetramer_array)
	pca = PCA(n_components=2)
	fit = pca.fit(tetramers_np).transform(tetramers_np)

	#Create graph with matplotlib

	#normalize sizes on scale from 20 - 200
	for contig in sizes:
		sizes[contig] = float(sizes[contig]) / float(max(sizes.values())) * 250 + 35
	sizes_list = []
	for contig in sorted(sizes.keys()):
		sizes_list.append(sizes[contig])


	#color contig points by blast assignments
	print "Plotting PCA graph of marker contigs"
	color_assigned = {}
	color_list = []
	color_brewer = ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a', '#b15928', '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6', '#ffff99']
	i = 0
	#assign colors
	most_common_names = Counter(names.values()).most_common()
	for name in most_common_names:
		if len(name) > 1:
			n = name[0]
			if n not in color_assigned:
				if i < len(color_brewer)-1:
					color_assigned[n] = color_brewer[i]
					i += 1
				else:
					color_assigned[n] = color_brewer[i]
		else:
			color_assigned[n] = '#ffff99'
	#create color list
	for name in sorted(names.keys()):
		color_list.append(color_assigned[names[name]])
	colors = np.random.rand(len(fit[:,0]))
	plt.scatter(fit[:,0], fit[:,1], s=sizes_list, c=color_list, alpha=0.9)

	#Add labels
	positions = {}
	i = 0 
	avg_size = sum(sizes.values()) / len(sizes.values())
	for name in sorted(names.keys()):
		if len(names.keys()) > 5 and sizes[name] > avg_size:
			positions[names[name]] = [fit[i,0], fit[i,1]]
		i += 1

	for name in positions.keys():
		x = positions[name][0]
		y = positions[name][1]
		
		plt.annotate(
        name, 
        xy = (x, y), xytext = (-20, 20),
        textcoords = 'offset points', ha = 'right', va = 'bottom',
        fontsize = 10,
        bbox = dict(boxstyle = 'round,pad=0.3', fc = 'white', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

	end = time.time()
	print "Run time: "
	print(end - start)

	plt.show()
