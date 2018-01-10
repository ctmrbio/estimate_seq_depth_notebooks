##################################################
Estimation of minimum metagenomic sequencing depth
##################################################
This report represents a draft summary of the outcomes of estimating the
required sequencing depth for different types of shotgun metagenome samples.

:Authors: Fredrik Boulund, <fredrik.boulund@ki.se>, Luisa Hugerth, <luisa.warchavchik.hugerth@ki.se>
:Date: 2017-12-19


Background
==========
The microbiome composition and microbial load of different sample types can
vary significantly.  For example, intestinal biopsies have been estimated to
contain about 1-3% microbial DNA, making it difficult (and potentially
expensive) to achieve decent resolution of microbial DNA using brute force
shotgun methods.

This work was a small pilot to provide descent guesstimates to use as starting
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
microbial composition profile of each sample type was designed based on
previous 16S results of similar samples. 10-50 samples of each type were summed
up to create a "metasample" of species, genera and higher order clades
characteristic of each sample type, associated to their total counts. Genomes
were then extracted from the RefSeq database in these proportions, while
adjusting the total number of genomes included per sample to a realistic value.
This means that rarely occuring genomes were included or not at random, which
reflects real ecological stochastic filtering. In the case of the vaginal swab,
a single genome from *Candida albicans* was added, corresponding to 10% of
cells in the sample.

.. TODO: more text needed here?


Simulated metagenome sequencing data
------------------------------------
First, we simulated large 50 million read metagenomes based on the "metasamples"
we produced earlier using `BBMap`_'s ``randomreads.sh`` with settings
``paired=true``, ``length=125``. Note that we did not use the ``metagenome``
feature of ``randomreads.sh`` due to our metasample genome files already
contain the expected proportion of all genomes already contain the expected
proportions of all genomes. The error profile of the reads was automatically
created by ``randomreads.sh`` and typic Illumina substitution errors were also
automatically added via the default settings of ``randomreads.sh``. 

Then, smaller metagenome sequencing samples were created using `BBMap`_'s
``reformat.sh``, using ``samplereadstarget=N`` with
``N={10000000,1000000,100000,10000}``, to produce simulated paired-end
metagenome sequences at several sequencing depths. Each "reference
metagenome" was randomly sampled three times to produce three technical
replicates of each sample type. The code for the metagenome simulation is
available in the Nextflow workflow file ``simulate_metagenomes.nf``, in our
`Github repository`_.

.. _BBMap: http://seqanswers.com/forums/showthread.php?t=41057


Taxonomic profile
-----------------
The taxonomic profiles of the triplicate simulated metagenomes were assessed
using `Centrifuge`_, `Kaiju`_, and `MetaPhlAn2`_. The methods use different
approaches to the taxonomic profiling problem.  Centrifuge uses a specialized
Burrows-Wheeler transform and Ferragina-Manzini index that compresses redundant
and non-species specific information in complete reference genomes to produce
an efficient index. Kaiju uses a modified backwards search using
Burrows-Wheeler transform to find maximum exact matches, and MetaPhlAn2 uses
Bowtie2 to align reads to a set of marker gene sequences.  Taxonomic profiles
are visualized using `Krona`_. The taxonomic profiling procedure is described
in detail in the Nextflow workflow ``taxonomic_profile.nf`` available in our
`Github repository`_.

.. _Kaiju: http://kaiju.binf.ku.dk/
.. _MetaPhlAn2: https://bitbucket.org/biobakery/metaphlan2
.. _Centrifuge: https://ccb.jhu.edu/software/centrifuge/manual.shtml
.. _Krona: https://github.com/marbl/Krona/wiki


.. TODO: More details about the comparison procedures needed. 


Functional profile
------------------
The metagenomes' functional potential was assessed using `TIGRFAMs`_, which is
a collection of curated multiple sequence alignments in the form of hidden
Markov models (HMMs). The HMMs accurately represent conserved protein families.
We used `HMMER3`_'s ``hmmsearch`` to search our metagenome sequence replicates.
This was done in two different ways: 
1) by directly translating each read into six reading frames, and searching for
matches to the TIGRFAM HMMs using the translated reads,
2) by first assembling each metagenomic sample individually, then making ORF
predictions using MetaGeneMark, and then running ``hmmsearch`` on the predicted
ORFs.
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
.. _MetaGeneMark: http://exon.gatech.edu/meta_gmhmmp.cgi

The functional profiles of each sample was normalized in two ways: 

1. By sample read count
2. By HMM length

This way, the samples can be compared to the normalized counts of observed
TIGRFAMs in the reference (meta)genomes.

.. TODO: More details about the comparison procedures needed. 


Results
=======
The basic statistics of the primary metasamples is summarized in the general
stats table below.

.. _general stats table:
.. table:: General statistics of the produced "samples". Statistics computed by summarizing FastQC output with MultiQC.
    :widths: auto
    :align: center

    +-------------+--------+--------+--------+
    |Sample Name  | % Dups | % GC   | M Seqs |
    +=============+========+========+========+
    |biopsy_50M_1 |   8.7% |  47%   |  50.0  |
    +-------------+--------+--------+--------+
    |biopsy_50M_2 |   8.7% |  47%   |  50.0  |
    +-------------+--------+--------+--------+
    |faeces_50M_1 |  12.7% |  48%   |  50.0  |
    +-------------+--------+--------+--------+
    |faeces_50M_2 |  12.8% |  48%   |  50.0  |
    +-------------+--------+--------+--------+
    |saliva_50M_1 |  10.6% |  43%   |  50.0  |
    +-------------+--------+--------+--------+
    |saliva_50M_2 |  10.7% |  43%   |  50.0  |
    +-------------+--------+--------+--------+
    |vag_50M_1    |  18.4% |  38%   |  50.0  |
    +-------------+--------+--------+--------+
    |vag_50M_2    |  18.5% |  38%   |  50.0  |
    +-------------+--------+--------+--------+


Overall, to create the smaller samples from the primary metasamples we produced
about 36 GiB of gzipped FASTQ data, across three replicates for each of the
four samples types::

    Sample types: Biopsies, Faeces, Saliva, Vagina
    Sequencing depths: 10M, 1M, 100k, 10k
    Replicates: 3

This produced 48 simulated shotgun metagenome samples in total.


Taxonomic profile
-----------------

The taxonomic profiles can be assess in multiple ways. First, we present
taxonomic composition profiles using Krona plots.

.. TODO: Here I want to put a few figures showing Krona plots.
   You need to help me make them look good Luisa :).

The proportion of unclassified reads shows a similar pattern between the
methods. 

.. figure:: taxonomy_plots/Centrifuge_unclassified_reads.png
    :figwidth: 50%
    :alt: Centrifuge unclassified reads.

    Centrifuge performs well on all sample types except vagina (possibly due to the presence of yeast).

.. figure:: taxonomy_plots/Kaiju_unclassified_reads.png
    :figwidth: 50%
    :alt: Kaiju unclassified reads.

    Kaiju performs well on all samples types except vagina (possibly due to the presence of yeast).

.. figure:: taxonomy_plots/MetaPhlAn2_unclassified_reads.png
    :figwidth: 50%
    :alt: MetaPhlAn2 unclassified reads.

    MetaPhlAn2 can only classify reads that align to the relatively
    short reference genes it uses, which is why the everage proportion
    of unclassified reads is much higher than Centrifuge or Kaiju.

A PCA on the abundance profiles show separation between the sample types at all
sequencing depths for Centrifuge and Kaiju. MetaPhlAn2 that appears to produce
profiles that are less clearly separated, especially for faeces, saliva, and
biopsy samples.

.. figure:: taxonomy_plots/Centrifuge_PCA_species.png
    :figwidth: 50%
    :alt: PCA of Centrifuge species composition 

    The taxonomic profiles are cleanly separated with Centrifuge, except for
    biopsy and saliva sample types that are located a bit closer together in
    the PCA plot.

.. figure:: taxonomy_plots/Kaiju_PCA_species.png
    :figwidth: 50%
    :alt: PCA of Kaiju species composition 

    The taxonomic profiles are cleanly separated with Kaiju.

.. figure:: taxonomy_plots/MetaPhlAn2_PCA_species.png
    :figwidth: 50%
    :alt: PCA of MetaPhlAn2 species composition 

    The taxonomic profiles are not cleanly separated across all sequencing
    depths, based on the MetaPhlAn2 output.

Pearson and Spearman correlations between samples for each method are available
as supplementary material.

We also compared the taxonomic profiles of Centrifuge with the true composition
profiles.

.. figure:: taxonomy_plots/Centrifuge_PCA_true.png
    :figwidth: 50%
    :alt: PCA of Centrifuge with true compositions.

    The taxonomic profiles of Centrifuge differ from the true compositions
    in all cases but for biopsies. The large discrepancy between the compositions 
    of vaginal samples is probably due to the large proportion of unclassified reads
    (likely yeast).


Functional profile
------------------
We begin with results from identifying TIGRFAMs directly from read data, and
then continue with results from mapping reads back to ORFs identified in
assembled contigs.

Identifying TIGRFAMs directly from read data
............................................
.. figure:: tigrfam_reads/saliva_Mainrole_diffs.png
    :figwidth: 75%
    :alt: Average TIGRFAM mainrole differences for saliva samples.

    Average TIGRFAM mainrole differences for saliva samples.

.. figure:: tigrfam_reads/saliva_Subrole_diffs.png
    :figwidth: 75%
    :alt: Average TIGRFAM subrole differences for saliva samples.

    Average TIGRFAM subrole differences for saliva samples.

.. figure:: tigrfam_reads/faeces_correlations.png
    :figwidth: 50%
    :alt: Correlation matrix for faeces samples

    Correlation matrix for faeces samples.

.. figure:: tigrfam_reads/biopsy_boxplots.png
    :figwidth: 75%
    :alt: Boxplots of biopsy samples.

    Boxplots of biopsy sample differences.


Assembled and reads mapped back to ORFs
.......................................
The proportion of reads that can be mapped back to assembled contigs increase
with increased sequencing depth. 

.. figure:: tigrfam_assembled_orfs/Proportion_reads_mapped_to_ORFs.png
    :figwidth: 75%
    :alt: Proportion of reads mapped to ORFs.
    
    Proportion of reads mapped to ORFs.

The low proportion of reads mapped to assembled ORFs gives low correlations between
predicted and expected TIGRFAM profiles.

.. figure:: tigrfam_assembled_orfs/faeces_correlations.png
    :figwidth: 50%
    :alt: Correlation matrix for faeces samples of reads mapped to assembled ORFs

    Correlation matrix for faeces samples of reads mapped to assembled ORFs.

The average error in estimated profile compared to true functional profile goes down with increasing
sequencing depth regardless of method, but the average difference and standard deviation 
looks best for high depth samples that are assembled. For read 

.. figure:: tigrfam_assembled_orfs/comparison_raw_reads_assembled_orfs.png
    :figwidth: 75%
    :alt: Comparison of average errors 

    Comparison of average absolute difference in functional profiles produced
    by annotating raw reads versus annotating ORFs from assemblies. 

Required sequencing depths for different sample types
-----------------------------------------------------
Here we summarize the minimum required sequencing depths for different sample
types, under the assumption that the expected bacterial content represents the
actual proportion of reads from the target community in the final sequencing
data. Our guesstimates are listed in the table below. 

.. table:: Estimated minimum required sequencing depths for different sample types.
    :widths: auto
    :align: center

    +-------------+-------------------------+----------------------------+----------------------------+
    | Sample type | Expected bacterial load | Minimum seq depth Taxonomy | Minimum seq depth Function |
    +=============+=========================+============================+============================+
    | Biopsy      |                    1-3% |                    50-100M |                         1G |
    +-------------+-------------------------+----------------------------+----------------------------+
    | Faecal      |                    >90% |                    500K-1M |                        10M |
    +-------------+-------------------------+----------------------------+----------------------------+
    | Saliva      |                    >90% |                    500K-1M |                        10M |
    +-------------+-------------------------+----------------------------+----------------------------+
    | Vaginal     |                  50-90% |                    500K-2M |                     10-20M |
    +-------------+-------------------------+----------------------------+----------------------------+

Discussion
==========
Based on Krona plots, it seems a decent representation of the original
community is achieved even at fairly low sequencing depths.
Kaiju ...
MetaPhlAn2 ...
Centrifuge ...

The functional profiles based on TIGRFAM annotation of reads seems to indicate
that when mapping individual translated reads, the functional profile reaches
decent detection coverage (>75%) somewhere after 1M reads. It also shows some
indications of overprediction at the 10M sequencing depths, based on the
observation that the average detection coverage being slightly higher that of
the reference sequences for the 10M sequencing depth samples. Our original
hypothesis regarding mapping individual translated reads versus mapping open
reading frames predicted from assembled metagenomes was that assembly and ORF
prediction would yield better sensitivity for samples with lower sequencing
depth (at least down to a lower limit). Our tests indicate that mapping reads
back to assembled predicted ORFs generally provide lower absolute sum of
differences to the normalized reference contig TIGRFAM counts than annotating
reads directly. When annotating reads directly it appears the gain in increased
sequencing depth diminishes as the average absolute sum of differences to 
the reference counts appears asymptotically converge to around 0.3. In the case 
of assembling contigs, calling ORFs, and mapping reads back, we see no
indication of this difference leveling out (but it might continue to improve
beyond sequencing depths of 10M). The quality of the functional profiles
of metagenomes when using an assembly-based method appears quite dependent on
having sequenced deep enough so that the assembly produces decent contigs to
map the reads back to. 

Performance-wise, taxonomic profiling is fairly light-weight and our
experiments were all run a fairly modest Linux server: 2x10 core Intel Xeon
E5-2630v4 CPUs @ 2.20 Ghz, with 64 GB RAM. Kaiju runs without issues with the
officially prepared ProGenomes database, even giving room for running multiple
processes in parallel. MetaPhlAn2 is very lightweight, both with regards to
memory and CPU usage. Centrifuge requires a lot of memory to prepare a custom
database, so we used the database of archea, bacteria, and viruses that the
authors' provide. 
Functional profiling, however, is much more demanding. Running TIGRFAM
annotation using ``hmmsearch`` directly on the annotated reads is a
computationally heavy task, to such an extent that some samples actually failed
to run to completion on our Linux server. It just barely completed in over two
weeks when run on the aforementioned server. To assemble reads, perform ORF
prediction, and annotate predicted ORFs with ``hmmsearch`` proved
computationally feasible. The time to assemble all samples using ``MegaHIT``
was less than 6 hours total on our server, using 40 cores and less than 64 GB
of memory.



Conclusions
===========

For taxonomic profiling, shotgun sequencing appears to provide good results
already at sequencing depths around 100k reads. We expect sensitivity to
increase with increasing read depth, so if detailed resolution is required for
study of particular low abundance species, higher sequencing depth is generally
better. We would hazard to conclude that a sequencing depth of 1M reads
from the target community is adequate. This means that for sample types such as
biopsies prepared with a standard protocol with an average bacterial content of
around 1-3%, a very high sequencing depth of about 50-100M is required, unless
the bacterial fraction of the sample can be separated from the biopsy.

For functional profiling, it is evident from our data that higher sequencing
depth leads to a better reproduction of the actual functional profile. However,
increasing read depth also increases the likelihood of overpredicting the
presence of TIGRFAMs, which is an important thing to consider if making an
analysis that depends on the relative abundances of identified TIGRFAMs.  All
in all, if you have a low coverage metagenome sample, our results indicate that
it is better to try to annotate raw reads directly, rather than trying to
assemble and identify ORFs. 

