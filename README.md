# MetaMarkers
Super fast microbial taxonomic profiling and visualization of metagenome assemblies.

This software was created for initial exploratory analysis of the prokaryotic taxa in a metagenomic assembly. It translates DNA naively in all 6 reading frames, searches using HMMER for 11 highly conserved marker ribosomal proteins, finds the closest known references to those proteins, and plots a tetranucleotide frequency PCA for the labeled contigs containing marker proteins. The results should allow you to quickly identify key species and contigs for genome binning of microbial genomes.

Most importantly, run time on full assemblies is generally less than **5 minutes** on a personal computer. 

*Python Dependencies: Python 2.7+, Biopython, Scikit-learn, Matplotlib.*
*Other Dependencies: Blast+ (blastn), HMMER3. These should be in your path.*

Installing all of the dependencies on Ubuntu is super easy, and can be done in one line: `sudo apt-get install hmmer ncbi-blast+ python-scikits-learn python-matplotlib python-biopython`.

###Usage

```
python meta_markers.py metagenome_assembly.fa
```

###Examples

(Runtimes were calculated on laptop with 4 GB RAM, 7,200 RPM HD, Intel core i3 with nproc=2.)

**Halite Metagenome**
*runtime: 4.15 minutes*
![halite](http://i.imgur.com/ZByehLT.png)

**Skin metagenome assembly**
*runtime: 1.05 minutes*
![skin](http://i.imgur.com/Xmo2HuC.png)

**Gut metagenome assembly**
*runtime: 4.4 minutes*
![gut](http://i.imgur.com/wPwH73j.png)

**Calcite endolithic Metagenome**
*runtime: 4.15 minutes*
![calcite](http://i.imgur.com/31nGcK4.png)

**Throat metagenome assembly**
*runtime: 3.45 minutes*
![throat](http://i.imgur.com/fpCq1Ad.png)
