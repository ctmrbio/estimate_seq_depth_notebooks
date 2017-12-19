##################################################
Estimation of minimum metagenomic sequencing depth
##################################################
This report represents a draft summary of the outcomes of estimating the
required sequencing depth for different types of shotgun metagenome samples.

:Authors: Fredrik Boulund, <fredrik.boulund@ki.se>, Luisa Hugerth, <luisa.hugerth@ki.se>
:Date: 2017-12-19


Background
==========
The microbiome composition and microbial load of different sample types can
vary significantly.  For example, intestinal biopsies have been estimated to
contain about 1-3% microbial DNA, making it difficult (and potentially
expensive) to achieve decent resolution of microbial DNA using brute force
shotgun methods.

This work was a small pilot to provide good guesstimates to use as starting
points for assessing the costs associated with shotgun sequencing of different
types of microbiome samples.


Methods
=======
The steps to assess the required sequencing depth are described in detail below. 
All code related to this endeavour is located in our `Github repository`_.

.. _Github repository: https://github.com/ctmrbio/estimate_seq_depth

Simulated sample type compositions
----------------------------------
The microbial composition of four different samples types was simulated. The
basis of the metagenome simulation was a handmade microbial community,
represented by whole genome sequences. To produce a somewhat realistic
situation, some whole genome sequences were present in multiple copies, to
simulate a relatively larger presence of those specific species. 

**PLACEHOLDER, more text needed**


Simulated metagenome sequencing data
------------------------------------
Metagenome sequencing data was simulated using `BBMap`_'s ``reformat.sh``,
using ``samplereadstarget=N`` with ``N={10000000,1000000,100000,10000}``, to
produce simulated paired end metagenome sequences, with somewhat realistic
error profiles. Each "reference metagenome" was randomly sampled three times to
produce three technical replicates of each sample type. The code for the
metagenome simulation is available in the Nextflow workflow file
``simulate_metagenomes.nf``, in our `Github repository`_.

.. _BBMap: http://seqanswers.com/forums/showthread.php?t=41057


Taxonomic profile
-----------------
The taxonomic profiles of the triplicate simulated metagenomes were assessed
using `Kaiju`_, `MetaPhlAn2`_, and `Centrifuge`_. The methods use different
approaches to the taxonomic profiling problem. Kaiju uses a modified backwards
search using Burrows-Wheeler transform to find maximum exact matches,
MetaPhlAn2 uses Bowtie2 to align reads to a set of marker gene sequences, and
Centrifuge uses a specialized Burrows-Wheeler transform and Ferragina-Manzini
index that compresses redundant and non-species specific information in
complete reference genomes to produce an efficient index. 
Taxonomic profiles are visualized using `Krona`_. The taxonomic
profiling procedure is described in detail in the Nextflow workflow
``taxonomic_profile.nf`` available in our `Github repository`_.

.. _Kaiju: http://kaiju.binf.ku.dk/
.. _MetaPhlAn2: https://bitbucket.org/biobakery/metaphlan2
.. _Centrifuge: https://ccb.jhu.edu/software/centrifuge/manual.shtml
.. _Krona: https://github.com/marbl/Krona/wiki


Functional profile
------------------
The metagenomes' functional potential was assessed using `TIGRFAMs`_, which is
a collection of curated multiple sequence alignments in the form of hidden
Markov models (HMMs). The HMMs accurately represent conserved protein families.
We used `HMMER3`_'s ``hmmsearch`` to search our metagenome sequence replicates.
As the comparisons are made on protein level, the reads were first translated
into all six reading frames using `BBmap`_'s ``translate6frames.sh`` with
default settings.  Filtering and counting the TIGRFAM HMM matches was done
using a custom script, ``count_tigrfam_annotations.py``. Only matches with
scores above the trusted cutoffs were used in downstream profiling analysis.
The full code for functional profiling is available in ``annotate_reads.nf``,
``annotate_contigs.nf``, and ``annotate_reference_contigs.nf`` in our 
`Github repository`_.

.. _TIGRFAMs: http://www.jcvi.org/cgi-bin/tigrfams/index.cgi
.. _HMMER3: http://hmmer.org/download.html

The functional profiles of each sample was normalized in two ways: 

1. By sample read count
2. By HMM length

This way, the samples can be compared to the normalized counts of observed
TIGRFAMs in the reference (meta)genomes.


Results
=======
Overall, we produced about 36 GiB of gzipped FASTQ data, across three
replicates for each of the four samples types::

    Sample types: Biopsies, Faeces, Saliva, Vagina
    Sequencing depths: 10M, 1M, 100k, 10k
    Replicates: 3

This produced 48 "samples" in total.

Taxonomic profile
-----------------


Functional profile
------------------
Functional profiling data is hard to compare accurately.

.. figure:: saliva_Mainrole_diffs.png
    :figwidth: 75%
    :alt: Average TIGRFAM mainrole differences for saliva samples.

    Average TIGRFAM mainrole differences for saliva samples.

.. figure:: saliva_Subrole_diffs.png
    :figwidth: 75%
    :alt: Average TIGRFAM subrole differences for saliva samples.

    Average TIGRFAM subrole differences for saliva samples.

.. figure:: faeces_correlations.png
    :figwidth: 75%
    :alt: Correlation matrix for faeces samples

    Correlation matrix for faeces samples


.. figure:: biopsy_boxplots.png
    :figwidth: 75%
    :alt: Boxplots of biopsy samples.

    Boxplots of biopsy sample type.

Such text. Many figure. Much wow.



Discussion
==========
Based on Krona plots, it seems a fairly good representation of the original
community is achieved even at fairly low sequencing depths.
Kaiju ...
MetaPhlAn2 ...
Centrifuge ...

The functional profiles based on TIGRFAM annotation of reads seems to indicate
that the functional profile reaches decent detection coverage (>75%) somewhere
after 1M reads. It also shows some indications of overprediction at the 10M seq
depths, based on the detection coverage being above that of the reference
sequences. 

Performance-wise, taxonomic profiling is fairly light-weight and our
experiments were all run a fairly modest Linux server: 2x10 core Intel Xeon
E5-2630v4 CPUs @ 2.20 Ghz, with 64 GB RAM. Functional profiling, however, is
much more demanding. It just barely completed in over two weeks when run on the
lightweight Linux server. 




Conclusions
===========

For taxonomic profiling, shotgun sequencing appears to provide good results
already at sequencing depths around 100k reads. We expect sensitivity to
increase with increasing read depth, so if detailed resolution is required for
low abundance species, higher is generally better. 

For functional profiling, it is evident that higher sequencing depth leads to a
better reproduction of the actual functional profile. However, increasing read
depth also increases the likelihood of overpredicting the presence of TIGRFAMs.
