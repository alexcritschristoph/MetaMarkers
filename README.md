# MetaMarkers
Super fast microbial taxonomic profiling and visualization of metagenome assemblies.

This software was created for initial exploratory analysis of the prokaryotic taxa in a metagenomic assembly. It translates DNA naively in all 6 reading frames, searches using HMMER for 11 highly conserved marker ribosomal proteins, finds the closest known references to those proteins, and plots a tetranucleotide frequency PCA for the labeled contigs containing marker proteins. The results should allow you to quickly identify key species and contigs for genome binning of microbial genomes.

Most importantly, run time on assemblies is generally less than 5 minutes, even on a personal computer.

*Python Dependencies: Python 2.7+, Biopython, Scikit-learn, Matplotlib.*
*Other Dependencies: Blast+ (blastn), HMMER3. These should be in your path.*

Installing all of the dependencies on Ubuntu is super easy, and can be done in one line: `sudo apt-get install hmmer ncbi-blast+ python-scikits-learn python-matplotlib python-biopython`.

###Usage

```
python meta_markers.py metagenome_assembly.fa
```

###Examples
