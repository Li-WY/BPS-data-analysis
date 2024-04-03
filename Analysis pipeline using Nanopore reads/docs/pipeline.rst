.. _pipeline:

BPS Analysis Pipeline
#############################################################


First off, we're happy to support your use of BPS.
Please make an "Issue" in the `git repo's issues section`_.
Feel free to post an issue just asking for more info, or with
a question or complaint or suggestion.

.. _git repo's issues section: https://gitlab.com/darachm/bps-dev/-/issues

This BPS pipeline *is a bit* complicated - the complexity is exposed 
because this pipeline is designed for flexibility.
This means you can use the pipeline for diverse applications,
including those unanticipated by the authors, but it does also mean that
the interface for users is less refined right now.

This means this pipeline is in flux. We have plans to greatly simplify the
pipeline now that we understand how it is likely to be used. Stay tuned
for a Genbank-file-based auto-configuration, more lightweight Docker containers,
performance tuning, and QC reports!

You may find *the fastest way to test if this works for you* is to 
consider the :ref:`Examples <examples>` and just change
what is relevant for your application.
But if you're keen on digging into :ref:`the details <details>`, 
both those as well as :ref:`an overview of the pipeline <bioinf-overview>`
are available on this page.

.. contents::
    :depth: 2
    :local:
    :backlinks: top

Experimental details you need to collect before analysis
=======================================================================

To use this tool, you'll need to configure it. Here's what you'll need:

#.  Where are your **FASTQ files**? This is very important.
#.  What is the right **medaka configuration** for these runs?
    Consult `their documentation`_ to select `which model`_ to use.
    You need to specify the model name (but not including the ``_model.tar.gz``
    suffix please!).
#.  Which **barcode corresponds to which sample** (original tube of combined 
    plasmids) in which sequencing pool? Usually you'll have just one barcode 
    per sample in one pool, but this 
    can be useful if you have samples run in multiple pools.
#.  What plasmids are in this run that you want to analyze? 
    More specifically, what is **a sequence that identifies
    each plasmid** from the others (or just identifies the only plasmid) ?
    Best to pick about a 100bp stretch.
#.  For each plasmid, what is **the sequence around the positioning barcodes**?
    Here, looking for about 100bp in the upstream and downstream 1kb or so
    (for the coarse chop), and then about 10-12bp immediately 
    upstream and downstream (for the fine-scale trimming).
#.  What are the **known positioning barcodes**?
#.  What kind of assembly do you want to do? **MSA or flye?** 
    Can you assemble them all 
    (with flye) and cut out the region of interest after? Do you want the whole
    plasmid? If it's very short (<1kb), you need to use the either the 
    MSA method or assemble first and then cut it out.
#.  What kind of post-assembly processing do you want to do?
    Do you want it **trimmed or rotated** to start with similar sequence?
#.  Optional: what are your target sequences? You don't need to compare the 
    assemblies to your targets, but this pipeline will **compare each assembly 
    to a set of designated targets**, that is *if* you provide it.

.. _their documentation: https://github.com/nanoporetech/medaka#models

.. _which model: https://github.com/nanoporetech/medaka/tree/master/medaka/data


Using the pipeline
=====================

Installation / Dependencies
------------------------------------------------------------

This pipeline depends on having:

* ``nextflow`` - `install nextflow from here`_
* ``singularity`` - `install singularity from here`_

.. _install nextflow from here: https://www.nextflow.io/index.html#GetStarted
.. _install singularity from here: https://docs.sylabs.io/guides/3.10/admin-guide/installation.html#installation-on-windows-or-mac

You will want to clone the git repo and then move into it, then modify the
:ref:`config file <configfile>`, then run it like so:

.. code-block::

    git clone http://gitlab.com/darachm/bps ./
    cd bps 
    # edit the config
    nextflow run main.nf ...[rest of command]....

If you detect any other [#deps]_ needs/dependencies, please `let us know`_!

.. _let us know: https://gitlab.com/darachm/bps-dev/-/issues

.. [#deps] If you want to use the ``Makefile``, you'll need ``make``. 
    To build the docs you'll need ``python3`` and ``tox``, or use ``sphinx``
    yourself (see the ``.gitlab-ci.yml``.


First write the config file
-----------------------------------------------------------------------

.. _configfile:

We're going to edit a text file in the `YAML`_ format.
Copy the file ``example-config.yaml`` to a new file -
you might call it ``config.yaml``.

.. _YAML: https://yaml.org/

You need four sections:

#. ``plasmids:`` defines which plasmids and how to recognize each
#. ``pools:`` defines the library pools sequenced and how to demultiplex
#. ``runs:`` defines the sequencing data and where to find it
#. ``experiments:`` - defines how to analyze/assemble it

``plasmids:``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's start with ``plasmids``. The below sets ``plasmids`` to have two
plasmids, ``plasmid1`` and ``plasmid2``, each with their corresponding
signatures with which to recognize them in the sequencing data.

.. code-block::

    plasmids:

      plasmid1:
        signature: 'tctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccCAGGCGGGCTCACCTCCGTGtggGCGGCCATggcgcgcc'

      plasmid2:
        signature: 'gggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtggGCGGCCATCGAGGGCTAGAATTACCTACCGGCCTCCACCATGCCTGCGggcgcgcc'

That's done. 

``pools:``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Then we'll setup ``pools`` to contain one pool that contains these 
plasmids and has a particular demultiplexing scheme.
Below, barcode10 corresponds to sample1:

.. code-block::

    pools:

      pool1:
        plasmids:
          - plasmid1
          - plasmid2
        demux: |
          barcode10	sample1
          barcode11	sample2
          barcode12	sample3
          barcode13	sample4
          barcode14	sample5
          barcode15	sample6
          barcode16	sample7
          barcode17	sample8
          barcode18	sample9
          barcode19	sample10

**Note that there's a TAB** character in between barcode and sample! 
**Note** that the ``|`` character after ``demux:`` permits us to have it 
formatted like this, so keep that ``|``.

``runs:``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Then we can define the ``runs``, so this defines what FASTQ files are of each
pool, with appropriate settings on ``medaka-model``. The model is determined
by what flowcell and setting you used to sequence it.
The ``pool:`` is for which of the sequencing ``pools:`` you put on there. 


.. code-block::

    runs:

      run1:
        fastq: '/path/to/the/folder/where/fastq/lives/the.fastq.gz'
        medaka-model: 'r941_min_sup_g507'
        pool: 'pool1'


``experiment:``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, ``experiment`` - the complicated one.
It's complicated because it defines the analysis that are actually done...
which takes several steps. For an overview, puruse the 
:ref:`Overview <bioinf-overview>`.

There are two example ``experiment`` blocks below, to contrast some different
analyses done on the same data.

The configuration block below details an analysis that:

#.  Considers the data for ``pool1``
#.  Filters out anything matching ``weird-stuff.fasta``
#.  Considers plasmids ``plasmid1`` and ``plasmid2``
#.  Uses the sequences in the ``samlami`` [#samlami]_ 
    section to trim these reads on the left and right
#.  Uses ``itermae-known-codes.yaml`` as a config for extracting the barcode from
    the coarsely chopped sequence using ``itermae`` [#itermae]_ 
#.  Again uses ``samlami`` [#samlami]_ to coarsely trim out the payload 
#.  Uses a MSA-based approach to generate and polish a consensus payload sequence
#.  Then cuts out the "payload" sequence using ``itermae-payload.yaml`` 
    and compares then to the target's FASTA file [#shortie]_


.. code-block::

    experiments: 

      extract-and-polish-small-payload:
        pools:
          - pool1
        filter-out: 'weird-stuff.fasta'
        plasmids:
          - plasmid1
        extract-barcode:
          plasmid1:
            samlami: 
              - arg: '--cut-on-left'
                ref: 'ctcaagcaaggttttcagtataatgttacatgcgtacacgcgtctgtacagaaaaaaaagaaaaatttgaaatataaataacgttcttaatactaacata'
              - arg: '--cut-on-right'
                ref: 'ttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
            itermae: 'itermae-known-codes.yaml'
            knowncodes: 'known-barcodes.fasta' 
        pre-assembly:
          samlami: 
            - arg: '--cut-on-left'
              ref: 'cgagGGTACCttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
            - arg: '--cut-on-right'
              ref: 'cacatatacctgccgttcactattatttagtgaaatgagatattatgatattttctgaattgtgattaaaaaggcaactttatgcccatgcaacagaaac'
        assembly-method: 'msa'
        post-assembly:
          cutoutdonor:
            itermae: 'itermae-payload.yaml'
            target-fasta: 'donor-barcodes.fasta' 
            target-size: 'short'

To contrast, the below config block specifies something very similar,
except that instead:

- Uses ``flye`` to assemble de novo from each well (after sub-clustering 
    [#subcluster]_ !)
- Uses ``samlami`` and ``itermae`` to cut out the payload and compare it to
    the donor barcodes using the ``short`` comparator
- Also, orients the entire plasmid sequences to start in the same spot
    (modulus) and compares it to the parent backbone sequence
    (using the default ``long`` comparator)

``experiments`` can have multiple ``post-assembly``'s to perform, so you can
assemble once and cut out particular pieces for each position's assembly.

.. code-block::

    experiments: 

      assemble-polish-then-extract::
        pools:
          - pool1
        filter-out: 'weird-stuff.fasta'
        plasmids:
          - plasmid1
          - plasmid2
        extract-barcode:
          plasmid1:
            samlami: 
              - arg: '--cut-on-left'
                ref: 'ctcaagcaaggttttcagtataatgttacatgcgtacacgcgtctgtacagaaaaaaaagaaaaatttgaaatataaataacgttcttaatactaacata'
              - arg: '--cut-on-right'
                ref: 'ttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
            itermae: 'itermae-known-codes.yaml'
            knowncodes: 'known-barcodes.fasta' 
        assembly-method: 'flye'
        post-assembly:
          rec-bc: 
            samlami: 
              - arg: '--cut-on-left'
                ref: 'cgagGGTACCttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
              - arg: '--cut-on-right'
                ref: 'cacatatacctgccgttcactattatttagtgaaatgagatattatgatattttctgaattgtgattaaaaaggcaactttatgcccatgcaacagaaac'
            itermae: 'itermae-payload.yaml'
            target-fasta: 'donor-barcodes.fasta' 
            target-size: 'short'
          modulus: 
            samlami:
             - arg: '--modulus'
               ref: 'AGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCA'
            target-fasta: 'entire-plasmid-backbone.fasta'

Then you go about running it (next section).

.. [#mangled] This is just there because our basecalled mangles the header 
   to an underscore-delimited barcode in the header. Just use plain guppy.
   Ignore and the default is guppy.

.. [#itermae] This is a tool for applying fuzzy regular expressions, iteratively,
   to extract out sequence patterns on DNA sequencing datasets. 
   `More details available on the itermae gitlab.`_

.. _More details available on the itermae gitlab.: https://gitlab.com/darachm/itermae/

.. [#samlami] This is a tool internal to this repository that does something
   very simple: it takes a SAM file and returns the soft-clipped part.
   So if you align a particular sequence and ``--cut-on-left`` you remove
   that sequence and everything to left, vice versa on right, and 
   ``--modulus`` just "rotates" the sequence to start with the matching (no
   sequence cut off here!). Very simple, very handy way to use bwa/mimnimap2
   to do trimming for you! **(thx lh3)**

.. [#shortie] Here using the ``short`` comparator to use ``bwa mem`` instead of 
   ``minimap2`` - because it seems to work better for short seqs. And that's
   probably in the documentation of it. Anyways, ``long`` is default, using
   ``minimap2`` to compare to targets.

.. [#subcluster] Plasmids often are recovered as dimers from cells that
   are propogating the plasmid. A few reads of dimers can scaffold the monomers
   to assemble as the dimer. Multiple different species may be present in
   a well and we want to see that. Therefore, we fit a mixture-model of
   at least two peaks and one long-tail of fragmented reads (from the main
   peak), and attempt to assemble each of these. This is the read-length 
   clustering portion of the pipeline. Future work aims to clarify the
   interpretation of this clustering, but for now we use it to correctly
   assemble each cluster of read-lengths and give them to the user to
   decide how to interpet any artifactual or biological heterogeneity.


Running it
-----------------------------------------------------------------------

Check out the Makefile, or adapt the below codeblock to launch the pipeline.

#.  This one starts from ``main.nf`` (as it should!) and reading ``nextflow``
    configuration from ``nextflow.config`` (which makes sense). 
    You should edit this config file to reflect your computational resources,
    such as available CPUs and RAM, and to use the right
    `nextflow executor`_.
    You can use the named docker containers, and future versions will be 
    trimmed down. On the `todo list`_.
#.  Next is some arguments that I always stuff into my ``nextflow`` calls.
    I like these ones.
#.  Next is selecting the ``common`` and ``blocky`` profiles from the config 
    file. This `optionality in nextflow`_ is useful, but you will have to 
    configure for your system.
#.  The option ``slurm_queue`` specifies what queue to submit to. 
    Consult your sysadmin or HPC documentation for details, or omit.
#.  Finally, the ``params-file`` is the path to the YAML config file 
    detailed at the start.

.. _nextflow executor: https://www.nextflow.io/docs/latest/executor.html
.. _todo list: https://gitlab.com/darachm/bps-dev/-/issues
.. _optionality in nextflow: https://www.nextflow.io/docs/latest/config.html#config-profiles

.. code-block::

    nextflow run main.nf \
        -resume -ansi-log false -with-dag reports/dag.html \
        -profile common,blocky \
        --slurm_queue high \
        -params-file your-config.yaml

Alternatively you could use the gitlab repo to direct it running that, but
you would likely want to specify a custom ``nextflow.config``
and YAML config locally first before launching the pipeline it!

What to expect - the results
-----------------------------------------------------------------------

When run, the pipeline will take a few seconds to load, then print 
the details of the configuration back. If you have issues, please scroll up
to check these - if some of the values are wrong then please check your
config and/or open a GitLab issue.

You should see jobs being submitted. If re-running, you'll see the jobs that
are loaded from cache specified as such. 

Outputs should be linked in the ``output`` folder.
Reports (run time, resource usage, etc) in the ``report`` folder.
``work`` contains all the runs, so when it says something like ``a6/c42569``
then there's a folder like that inside of ``work`` that is the working folder
for that step. Check it out.

At the end, you should get a file in ``output`` folder that's a ``.tsv``
file with a name derived from your config file. Each row is a 
position/cluster/assembly, each column is information about that.
Here are the columns:

* experiment_1

    Name of the experiment

* medaka_2

    Which ``medaka`` configuration was used for this call - different configs
    are not yet merged.

* plasmid_3

    Which plasmid this is for - multiple plasmids can be separated based on
    the "signature" sequence you configure.

* sample_4

    Which sample - defined in the ``pools:`` section

* position_5

    Which position in the sample. This is just the FASTA ID of the barcode that
    it matches.

* positioning_code_6

    Copying back what that barcode is.

* coverage_7

    How many raw reads were used to generate this assembly.

* well_purity_8

    Of those raw reads, how many had >90% identity with the consensus sequence.

* assembly_method_9

    What method is used for assembly - [msa] or [flye].

* assembly_detail_10

    A detail output by each assembly method, for [msa] it's the votes for the
    most controversial base, and for [flye] it's the cluster/contig.

* post_assembly_11

    What post assembly output this is, this allows you to do one assembly and
    then both output the full length product but also sections of it that have
    been cut out of the full length assembly. All together.

* length_sequence_12

    The length of the assembly.

* cluster_13

    For [flye], which cluster this is (based on read lengths).

* cluster_n_reads_14

    For [flye], how many reads went into this read-length-cluster.

* cluster_median_length_15

    For [flye], what the median read-length is for this read-length-cluster.

* cluster_cv_length_16

    For [flye], what the CV (mean/sd) read-length is for this 
    read-length-cluster.

* matches_target_17

    What target is the best match for the assembly. 

* query_errors_18

    When we align to that best match, how many errors are made? This is the
    sum of the change lengths, insertion lengths, deletion legnths, and 
    soft-clipping at the end.

* query_changes_19

    How many changes, relative to reference? This is from the MD tag.

* query_change_length_20

    How long are all changes, relative to reference? This is from the MD tag.

* query_insertions_21

    How many insertions, relative to reference? This is from the CIGAR tag.

* query_insertion_length_22

    How long are all insertions, relative to reference? This is from the CIGAR tag.

* query_deletions_23

    How many deletions, relative to reference? This is from the CIGAR tag.

* query_deletion_length_24

    How long are all deletions, relative to reference? This is from the CIGAR tag.

* query_clipped_length_25

    How long are all soft-clipping off the ends, when aligning to reference
    best match? This is from the CIGAR tag.

* sequence_26

    Finally, the sequence of the assembly - whew. 

Also note there's a few plots output. And in the ``output/`` directory there
should be more subdirectories with more useful files:

* assemblies

    FASTA files of the assemblies for each position.

* aligns

    Alignment files of these assemblies back against the targets.

* cluster_plots

    Plots showing how each position's set of reads was clustered by read-length
    into read-length clusters. 

* readlens

    Files that have all the readlengths for each cluster, for later analysis.


.. _bioinf-overview:

Overview of the bioinformatics approach and concepts
=================================================================

This overview is intended to give bioinformaticians some technical insight, 
without having to dig into the code. 

This is a ``nextflow`` pipeline, executed by directing the ``nextflow`` 
`executable`_ to run the workflow as specified the file ``main.nf``. 
This imports
the various other ``*.nf`` files (in this repo) to add sub-pipeline modules, 
and relevant configuration profiles are defined in the ``nextflow.config``.

You *must* edit this config file so that it fits the resources and context you
are executing it in! This means if it's on an HPC, if using different Docker
containers, the available memory/cpus, etc. 
Documentation for this :ref:`executors <executor>` section is linked.

Nextflow then runs the steps of each process 
in different directories using ``singularity``
executing from copies of ``docker`` images [#docker]_ . 

There are several additional scripts that are used, these are 
written in python and R and are in the ``scripts/`` folder in this repo.

So, this should just need ``nextflow`` and ``singularity`` installed to
launch and run.

.. _executable: nextflow.io/index.html#GetStarted

.. [#docker] All of the benefits of developing in ``docker``, but without having
   to run it as root :sunglasses: .


Diagram
------------------------------------------------------------

Here's a flowchart of what is done.

.. mermaid::

    flowchart TB

        classDef optional fill:#fff,stroke:#666,stroke-width:2px,font-size:12pt;
        classDef step fill:#fff,stroke:#000,stroke-width:2px,font-size:12pt;
        classDef file fill:#ccf,stroke:#000,stroke-width:2px,font-size:12pt;

        input[input FASTQ files]:::file

        filterz([>500bp and non-contaminant reads]):::optional
        input -->|``chopper`` filters for length and<br>``minimap2`` against contaminants| filterz 

        demux([assigned sample label based on each pool and barcode]):::step
        filterz -->|shell commands| demux

        plasmids([separated reads for each plasmid]):::step
        demux -->|align to signatures using ``minimap2``,<br>split the SAM file into plasmid files| plasmids

        chopped([coarsely extracted position barcode]):::step
        plasmids -->|chop away backbone coarsely with<br>``minimap2`` and samlami.py| chopped

        barcode([each read's exact position barcode]):::step
        chopped -->|finely extract barcode using ``itermae``| barcode

        clustered([barcodes clustered and assigned to position]):::step
        barcode -->|cluster barcodes with ``starcode`` and <br>assign to position with barcode2well.py| clustered

        together([reads, and a table of each read assigned to a position]):::step
        plasmids & clustered --> together 

        separated([reads separated for each position]):::step
        together -->|pool all reads for sample from across runs and<br>separate reads per position with awk| separated

        preass([optionally trimmed before assembly]):::optional
        separated -->|optional trimming using<br>samlami.py and/or ``itermae``| preass

        collected([reads per position]):::step
        separated & preass --> collected

        cluster_length([reads clustered by length within each position]):::step
        collected -->|separate read-length clusters using mixture-model<br>in length-cluster.py| cluster_length

        flye([de novo assemblies]):::step
        cluster_length -->|de novo assembly with ``flye`` using ``trycycler``| flye

        flye_polished([polished de novo assembly]):::step
        flye -->|polish assembly with ``medaka``| flye_polished

        msa([multiple-sequence alignment consensus]):::step
        collected -->|aligned and merged with kalign3 and msafasta2consensus.py| msa
        msa_polished([polished consensus MSA]):::step
        msa -->|polish consensus with racon+``medaka``| msa_polished

        assembled([candidate assemblies]):::step
        flye_polished --> assembled
        msa_polished --> assembled



        postproc([processed/trimmed assembly/consensus]):::optional
        assembled -->|optional processing/re-orientation<br>with samlami and/or ``itermae``| postproc

        aln2ref([assembly/consensus aligned to reference target]):::step
        postproc -->|``minimap2``/``bwa`` align<br>to reference target| aln2ref
        assembled -->|``minimap2``/``bwa`` align<br>to reference target| aln2ref

        purity([assess position purity by alignment of<br>input reads to the consensus/assembly]):::step
        assembled -->|align raw reads<br>to each position's<br>assembly/consensus| purity

        analyze([analysis]):::step
        postproc & purity & aln2ref & assembled --> analyze

        results[result tsv files with positions called<br>as pure and/or matching target seq]:::file
        analyze -->|R script| results



Tools used
----------------------------------------------------------------------

The pipeline makes use of at least these tools

* Ubuntu GNU/Linux and various shell utilities
* GNU parallel
* awk
* chopper
* bwa
* minimap2
* itermae
* starcode
* kalign3
* flye
* trycycler
* samtools
* R
* racon
* medaka
* python3, with packages:
    * BioPython
    * matplotlib
    * numpy
    * pandas
    * pysam
    * rapidfuzz
    * scipy

and some custom scripts:

* samlami.py

    This tool simply returns the sequence from a SAM record, but removes or
    rearranges it based on the arguments. ``--cut-on-left`` removes the
    match and everything upstream of it. ``--cut-on-right`` is the opposite.
    ``--modulus`` re-arranges the sequence to start with the match (removing
    nothing, ie rotating it).

* pairaln2ref.py 

    Align a query to a ref using pairwise alignments in BioPython.
    Used to assess how pure the position is for a particular assembly/consensus.

* barcode2well.py

    Used to map extracted barcodes to a reference set to find the closest
    match.

* length-cluster.py

    For a particular FASTQ file, uses a mixture model to separate the 
    reads into at least two clusters (however many minimizes the AIC,
    above two clusters not including the tail of shorter subreads).
    Done via MLE with scipy.

* msafasta2consensus.py

    Reads a multiple sequence alignment and votes on the consensus base for
    each position based on the base's quality annotation. 
    Used to generate a draft consensus.

* puritysam2tsv.py

    A quick function to parse SAM files to assess the alignment statistics.
    Necessary because the SAM files were crashing my R script.

.. _executor:

Execution context
----------------------------------------------------------------------


You need to edit the ``nextflow.config`` to change the ``executor =``
to reflect how you're executing it. 
``nextflow`` `supports`_ a variety of executors, such as local (default),
SLURM, or even rented computers through companies like Amazon and Google.

.. _supports: https://www.nextflow.io/docs/latest/executor.html

The specific containers run are detailed in the ``nextflow.config``, 
and are specified in each ``process`` block using the ``label`` directive
and that corresponds to specific containers specified in the
``nextflow.config`` file.

.. _concepts:

Glossary of terms, concepts
-----------------------------------------------------------


For the purpose of this pipeline, here's my understanding and usage of some 
concepts, as that may help clarify how to configure and use the pipeline.

barcode

    A particular variable sequence surrounded by a particular fixed sequence.
    This variable sequence is relatively short, allowing it to stably and
    cheaply represent more complex genetic variation or operations.

plasmid

    Circular DNA sequence. The type of DNA we are (presently) 
    expected to be sequencing.

    For the purposes of this pipeline, we distinguish different plasmids using
    "signatures" of DNA sequence that are somewhat unique. This is using
    alignment, so about 100bp seems to work but it depends on your design.

position

    A well or position where there is a colony that was subject to this 
    barcode multiplexing. In essence, the smallest wetbench unit of 
    multiplexing, and the unit of our particular concern here.

    Example: well A1 on plate 2 in sample 14 is a particular position.

sample

    The particular mixture of plasmids extracted from a combined collection of 
    cells in various positions.

    Example: the plasmids extracted from all of plate 2 scraped 
    into a single tube.

pool

    A particular pool of DNA from samples that are multiplexed together, 
    such that another sample barcode has been introduced 
    (for instance during library preparation) that can distinguish the different
    samples. 

    Example: the samples from plate 2 and plate 3 are in a pool7 where they
    are represented by barcodes 2 and 3 respectively.

    Different mixtures of different samples are different pools, so if you
    add more multiplexed samples to a pool then it is a new pool.

run

    A batch of reads (here, pretty much exclusively Nanopore reads) that
    were generated from one run of a flowcell. This has one pool of prepared
    and (optionally) multiplexed sequences on it.

    Example: pool7 was run on a v9.4.1 flowcell to generate a FASTQ file,
    then later pool7 was run on a v10.4.1 flow to generate a different 
    FASTQ file.

experiment

    A set of analyses done, with the aim of separating plasmid long-reads by
    a particular barcode then to use those reads to assemble a consensus 
    sequence and determine
    how "pure" and "correct" that assembly is.

    Basically, the analysis choices to recognize, partition, and 
    assemble/extract the region of interest from the plasmid for each 
    position.

pure

    A colony is "pure" if it has mostly one genotype of concern. 
    No colony is truly pure, but here the term is used for what is likely
    to be operationally useful (ie there is no significant evidence that
    multiple genotypes are present in the colony).

correct

    A colony is "correct" if it has one of the desired target sequences, ie
    there aren't mismatches, insertions, or deletions in the region of interest.



.. _details:

Details
==============


In this page, I attempt to document all the options and details.
This is a reference page.

An quick primer on the YAML format
---------------------------------------------------------------------

The YAML format is a way of specifying data objects in plain text.
It's like JSON, but readable. Here's `the docs`_, 
here's `an easier to use guide`_ 
and below is a quick primer on the format.
You can `play around with it here`_ on that interactive demo parser.

.. _the docs: https://yaml.org/spec/1.2.2/
.. _an easier to use guide: https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html
.. _play around with it here: https://yaml-online-parser.appspot.com/

Indentation is important, and determines nesting (an indented object is 
contained within a previous less-indented object, roughly).
A colon ':' usually notates that it's a key-value / dictionary / map -type of
object.

Consider this YAML object:

.. code-block::

    a_key: 'a_value'
    b_key: 10
    c_key: 
      deeper_key: 'frank'

This is a dictionary with keys of 'a_key', 'b_key', and 'c_key'.
The value stored under 'a_key' is the string 'a_value'. 
The value stored under 'b_key' is the number 10.
The value stored under 'c_key' is another dictionary, this stored dictionary 
has a key 'deeper_key' that corresponds to a value of the string 'frank'.

An ordered list is notated by a hyphen '-'. These lists can be inside of
dictionaries or lists, and can contain dictionaries or lists. Or strings.
Or whatever.

Consider this YAML object:

.. code-block::

    here_be_a_list:
        - a string
        - a key:
            a value
    here-be-a-list:
        - 
            - these
            - are
            - items in a list inside the list
        - b key:
            - list in here
            - too

This is a dictionary with two keys. 
The value stored under 'here_be_a_list' is a list of two items.
The first item is the value of 'a string'.
The second item is a dictionary with a single key of 'a key', and
that corresponds to the value of the string 'a value'.

The value stored under 'here-be-a-list' is a list of two items, as well.
The first item is a list, this sub-list contains three values of
'these', 'are', and 'items in a list inside the list'.
The second item is a dictionary, where the key 'b key' corresponds to the
list of strings 'list in here' and 'too'.

That should be enough to understand the nomenclature of the YAML config
file. 

Config file details
---------------------------------------------------------------------

There are four parts. 
These must be the top-level keys. 
Other extra keys may be present, but will be ignored.

Below are the four top-level keys. 

plasmids:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each key in the plasmids dictionary is the name of plasmid.
You need to keep the name of the plasmid consistent across
other uses of it, so maybe keep it short and simple.

For each plasmid, here are the mandatory fields:

* signature:

        This is a DNA sequence used to identify the plasmid. This works by
        aligning reads using minimap2 to a reference of possible signatures
        (as specified in the ``pools:`` section). Maybe go with 100 bases.

Example:

.. code-block::

    plasmids:

      plasmid1:
        signature: 'tctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccCAGGCGGGCTCACCTCCGTGtggGCGGCCATggcgcgcc'

      plasmid2:
        signature: 'gggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtggGCGGCCATCGAGGGCTAGAATTACCTACCGGCCTCCACCATGCCTGCGggcgcgcc'


pools:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each key in the pools dictionary is the name of a pool.
You need to keep the name of the pool consistent as you refer to it in
experiments.

For each pool, here are the mandatory fields:

* plasmids:

        This is a list of which plasmids (as named above) are in the pool.
        This is important so that it will use minimap2 to try and use the
        signature sequence to separate out each plasmid.

* demux:

        This is the multiplexing/de-multiplexing information. Note the ``|``
        character that permits us to write a tab-delimited file in the config
        field. There must be a tab between the barcode and sample fields.
        One line per each. The barcode field (first) corresponds to the 
        barcode designator as output by ``guppy`` demultiplexer. 

.. code-block::

    pools:

      pool1:
        plasmids:
          - plasmid1
          - plasmid2
        demux: |
          barcode10	sample1
          barcode11	sample2
          barcode12	sample3
          barcode13	sample4
          barcode14	sample5
          barcode15	sample6
          barcode16	sample7
          barcode17	sample8
          barcode18	sample9
          barcode19	sample10


runs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each key in the runs dictionary is the name of a run.
The name of the run doesn't matter, because you specify which pool is in the
run. 

For each run, here are the mandatory fields:

* fastq:

        The path of a FASTQ (or FASTQZ) file, a list of these, or a `glob`_
        that expands to some of these files.

* pool:

        Which pool was run on this run. A different mix of barcodes is a 
        different pool, so you can't mix pools together - they become a new
        pool.

* medaka-model:

        The medaka model appropriate for polishing using these raw reads.

.. _glob: https://www.gnu.org/software/bash/manual/bash.html#Filename-Expansion

Optional fields:

* flowcell:

        Which flowcell is used. Not important anymore.

* header-parse:

        By default, the tool expects to find the FASTQ comment section similar
        to how ``guppy`` outputs it. However, if you've mangled it, perhaps
        have deleted the comment and put the barcode name appended to the ID
        using a "_", then setting this to "modified" will parse that correctly.

.. code-block::

    runs:

      run1:
        fastq: '/path/to/the/folder/where/fastq/lives/the.fastq.gz'
        pool: 'pool1'
        medaka-model: 'r941_min_sup_g507'


experiments:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The experiment names are unimportant, but they will be used to separate out
files through the pipeline and name the result files.

For each experiment, here are the mandatory fields:

* pools:

        What pools should be considered for the experiment?
        These names should match the names of the pools in previous
        sections.
        Specific samples can be included or excluded in the optional options
        below.

* plasmids:

        What plasmids are considered in this experiment?
        These names should match the names of the plasmids in previous
        sections.

* extract-barcode:

        Instructions for how to extract the positioning barcode.
        This is keyed for each plasmid name.
        Each plasmid-keyed entry needs to specify:

            * samlami: 

                A set of ``SAMlami`` slicer instructions. This consists of
                a list with fields of 'arg', 'ref', and optionally 
                'head' or 'tail'.
                'arg' specifies the mode, so '--cut-on-left' to delete the
                'ref' sequence and everything to the left, and '--cut-on-right'
                to delete the 'ref' sequence and everything to the right.
                '--modulus' rotates the sequence to start with the 'ref'
                sequence and does not delete it.
                'ref' is the sequence used in the comparison, about 100 bases
                works okay.
                'head' or 'tail' is an option to take the head or tail of
                'ref'. Totally optional. Not recommended.

            * itermae: 

                An ``itermae`` configuration file that extracts the barcode.
                Documenting that is beyond the scope of this document, and
                we hope to make it simpler in the future.

            * knowncodes: 

                A FASTA file of known barcodes, where IDs are the positions
                that each barcode represents.

* assembly-method:

        Two options, the string 'flye' or the string 'msa'.
        The 'flye' option uses ``trycycler`` and ``flye`` to do *de novo*
        complete assembly. The 'msa' option uses ``kalign`` to generate a 
        multiple sequence alignment, uses a custom script to merge them into
        a draft consensus, and uses ``racon`` and ``medaka`` to polish this
        consensus using the raw reads.

* post-assembly:

        Analyses to do, post-assembly. 
        These are keyed with the analysis name.
        Then, they can have *optionally* a ``SAMlami`` stage with multiple
        operations (such as cut on left or right) or a single operation
        (such as '--modulus' to rotate it to start with similar sequence.
        They can have an *optional* step of an ``itermae`` config to fine
        trim using that tool.
        If none of these are specified, then it just passes the assembly
        through to the next step (target comparison).

        Then we have the *option* of specifying a single or list of target
        fasta sequences with 'target-fasta'. If present, the product after
        applying ``SAMlami`` and ``itermae`` steps will be aligned to the
        target. By default, this is done with ``minimap2``, but if 
        'target-size' is set to 'short' then we use ``bwa mem``.

Optional:

* include-samples:

        A list of samples to include. If this is present, then 
        only these samples in the named pools are considered.

* exclude-samples:

        A list of samples to exclude. If this is present, all other samples
        in the pool are considered.

* pre-assembly:

        Optional pre-assembly processing. These are steps done before the
        assembly. A common one is to use SAMlami slicer to coarsely trim
        the reads down to a target region before doing the MSA assembly.
        Similar syntax as the 'post-assembly' section.

.. code-block::

    experiments: 

      extract-and-polish-small-payload:
        pools:
          - pool1
        filter-out: 'weird-stuff.fasta'
        plasmids:
          - plasmid1
        extract-barcode:
          plasmid1:
            samlami: 
              - arg: '--cut-on-left'
                ref: 'ctcaagcaaggttttcagtataatgttacatgcgtacacgcgtctgtacagaaaaaaaagaaaaatttgaaatataaataacgttcttaatactaacata'
              - arg: '--cut-on-right'
                ref: 'ttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
            itermae: 'itermae-known-codes.yaml'
            knowncodes: 'known-barcodes.fasta' 
        pre-assembly:
          samlami: 
            - arg: '--cut-on-left'
              ref: 'cgagGGTACCttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
            - arg: '--cut-on-right'
              ref: 'cacatatacctgccgttcactattatttagtgaaatgagatattatgatattttctgaattgtgattaaaaaggcaactttatgcccatgcaacagaaac'
        assembly-method: 'msa'
        post-assembly:
          cutoutdonor:
            itermae: 'itermae-payload.yaml'
            target-fasta: 'donor-barcodes.fasta' 
            target-size: 'short'

    experiments: 

      assemble-polish-then-extract::
        pools:
          - pool1
        filter-out: 'weird-stuff.fasta'
        plasmids:
          - plasmid1
          - plasmid2
        extract-barcode:
          plasmid1:
            samlami: 
              - arg: '--cut-on-left'
                ref: 'ctcaagcaaggttttcagtataatgttacatgcgtacacgcgtctgtacagaaaaaaaagaaaaatttgaaatataaataacgttcttaatactaacata'
              - arg: '--cut-on-right'
                ref: 'ttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
            itermae: 'itermae-known-codes.yaml'
            knowncodes: 'known-barcodes.fasta' 
        assembly-method: 'flye'
        post-assembly:
          rec-bc: 
            samlami: 
              - arg: '--cut-on-left'
                ref: 'cgagGGTACCttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
              - arg: '--cut-on-right'
                ref: 'cacatatacctgccgttcactattatttagtgaaatgagatattatgatattttctgaattgtgattaaaaaggcaactttatgcccatgcaacagaaac'
            itermae: 'itermae-payload.yaml'
            target-fasta: 'donor-barcodes.fasta' 
            target-size: 'short'
          modulus: 
            samlami:
             - arg: '--modulus'
               ref: 'AGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCA'
            target-fasta: 'entire-plasmid-backbone.fasta'






Input files
---------------------------------------------------------------------

This expects to read in FASTQ files as input, so after basecalling.
This also expects those FASTQ files to be basecalled by ``guppy`` as
published by Oxford Nanopore. You'll have to register with them to download
and set it up.

By default the pipeline expects that the sample-demultiplexing barcode is
present in the header as output by ``guppy``. 
If the header instead consists of the unique read ID, then and underscore,
then the barcode (such as ``@44803c81-5018-46f3-9fa1-b069f68a37fa_barcode84``),
then you need to put ``header-parse: 'modified'`` in the run's block in the
config file.





You don't need to do this, but if you want to specify plasmids that you don't 
intend to analyze, you can specify those to "subtract" them from being 
recognized as other plasmids.
This is all done by alignment of the "signature" sequence, so the best
match wins...



.. _examples:

Examples
================================

These are examples of analyses, specifically the analyses as run for the paper.
This is a big section!

First, we have the overall configuration, then the specific configs.

Common file - ``nextflow.config``
----------------------------------------------------

::

    profiles {
        common{ 
            singularity {
                enabled = true
                autoMounts = true
                runOptions = '--no-home'    // essential to prevent mounting
                                                    // the local HOME and thus
                                                    // pip doing stupid shit!
                cacheDir = '/home/l/.singularity/'
            }
            cache = 'lenient'
            trace.enabled       = true
            report.enabled      = true
            timeline.enabled    = true
            dag.enabled         = true
            trace.overwrite     = true
            report.overwrite    = true
            timeline.overwrite  = true
            dag.overwrite       = true
            trace.file          = "reports/nextflow_pipeline_trace.txt"
            report.file         = "reports/nextflow_pipeline_report.html"
            timeline.file       = "reports/nextflow_pipeline_timeline.html"
            dag.file            = "reports/nextflow_pipeline_dag.png"
        }
        blocky { // This is our computer. It is blocky.
            executor {
                name = 'slurm'
                cpus = 36
                memory = '60GB'
            }
            process {
                executor = 'slurm'
                clusterOptions = '--propagate=ALL'
                withLabel: 'all_core' { cpus = 36 }
                withLabel: 'half_core' { cpus = 18 }
                withLabel: 'quarter_core' { cpus = 6 }
                withLabel: 'one_core' { cpus = 1 }
                withLabel: 'all_mem' { memory = '58G' }
                withLabel: 'half_mem' { memory = '28G' }
                withLabel: 'smol_mem' { memory = '2G' }
                withLabel: 'bioinfmunger' { container = 'docker://darachm/bioinf:bioinf-sam-bedtools-emboss-ncbi-ucsc-genometools-htslib'}
                withLabel: 'lh3aligners'  { container = 'docker://darachm/lh3-aligners:minimap2-bwa-bwamem2'}
                withLabel: 'itermae' { container = 'docker://darachm/itermae:plus' }
                withLabel: 'starcode' { container = 'docker://darachm/starcode:latest' }
                withLabel: 'medaka'  { container = 'docker://darachm/nanopore:medaka-hack' }
                withLabel: 'chopper'  { container = 'docker://darachm/nanopore:chopper' }
                withLabel: 'assemble'  { container = 'docker://darachm/flye:flye-miniasm-mash-muscle-r-pkg-mm-try-canu' }
                withLabel: 'racon'  { container = 'docker://darachm/nanopore:racon' }
                withLabel: 'kalign' { container = 'docker://darachm/kalign:bioinf-kalign' }
                withLabel: 'r' { container = 'docker://darachm/rr:r-4.3.1-tidy-db-viz-mod-bio' }
                withLabel: 'jbrowse' { container = 'docker://darachm/jbrowse:serve' }
                withLabel: 'plannotate' { container = 'docker://darachm/plannotate:latest' }
            }
            mail {
                smtp.host = 'localhost'
                smtp.port = 25
                smtp.user = 'darachm'
            }
        }
    }




For each, we first have the config file, then how it was run, then the names
of the outputs.

1100 random oligos analysis
--------------------------------------------

Command run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    nextflow main.nf -c nextflow.config \
        -resume -ansi-log false -with-dag reports/dag.html \
        -profile common,blocky \
        --slurm_queue high \
        -params-file config-1100.yaml

BPS config yaml
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    pools:

      bc_pool_I:
        filter-out: 'nuconfig/filter-out-r6k-psc101-937.fasta'
        plasmids:
          - 438-937
          - 439-937
        demux: |
          barcode81	bc_bc_3
          barcode82	bc_bc_4
          barcode83	bc_bc_5
          barcode84	bc_bc_6
          barcode85	bc_bc_7
          barcode86	bc_bc_8
          barcode87	bc_bc_9
          barcode88	bc_bc_10

    plasmids:

      438-937:
        signature: 'ATGCATATGGGTTACCTGTACACGTACGTTCGAAGCCGGCGCCCCTAGGACTAGTACGCGTccaCACGGAGGTGAGCCCGCCTGCCCGGGGCTAGCGTCG'

      439-937:
        signature: 'ATGCATATGGGTTACCTGTACACGTACGTTCGAAGCCGGCGCCCCTAGGACTAGTACGCGTCCCGGGGCTAGCGTCGACTCTAGAGGATCGATCCTTTTT'

      1064-937-lasso:
        signature: 'GGGCCACTAGGGACAGGATTGGGCGGCCATCGAGGGCTAGAATTACCTACCGGCCTCCACCATGCCTGCG'

    runs:

    #  20221111_1544_MN21180_FAT59932_c3bac6f0:
    #    fastq: '/archive/l/seq_data/bps/20221111_1544_MN21180_FAT59932_c3bac6f0/basecalling/output/demuxed.fastq.gz'
    #    flowcell: 'FLO-MIN112'
    #    pool: '428_lasso_bc_pool_I'
    #    medaka-model: 'r104_e81_sup_g610'
    #    header-parse: 'modified'
    #
    #  20221113_1507_MN24328_FAT59932_3e555f90:
    #    fastq: '/archive/l/seq_data/bps/20221113_1507_MN24328_FAT59932_3e555f90/basecalling/output/demuxed.fastq.gz'
    #    flowcell: 'FLO-MIN112'
    #    pool: '428_lasso_bc_pool_I'
    #    medaka-model: 'r104_e81_sup_g610'
    #    header-parse: 'modified'

      20230718_plasmidsaurus:
        fastq: '/archive/l/seq_data/bps/20230718_plasmidsaurus/*gz'
        flowcell: 'FLO-PRO114M'
        pool: 'bc_pool_I'
        medaka-model: 'r1041_e82_400bps_sup_v4.2.0'
        header-parse: 'guppy'

    #####
    #####
    ##### TODO conceptually rewrite so there's a top level ontology category of
    ##### 'assembly', so all the pre-assembly and assembly steps can be pooled
    ##### across 'experiments', which maybe should just be 'analyses' then
    #####
    #####

    samlami-rec-bc: &samlami-rec-bc
      - arg: '--cut-on-left'
        ref: 'ctcaagcaaggttttcagtataatgttacatgcgtacacgcgtctgtacagaaaaaaaagaaaaatttgaaatataaataacgttcttaatactaacata'
        tail: 100
      - arg: '--cut-on-right'
        ref: 'ttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
        head: 100

    samlami-internal: &samlami-internal
      - arg: '--cut-on-left'
        ref: 'cgagGGTACCttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
        tail: 100
      - arg: '--cut-on-right'
        ref: 'cacatatacctgccgttcactattatttagtgaaatgagatattatgatattttctgaattgtgattaaaaaggcaactttatgcccatgcaacagaaac'
        head: 100


    experiments: 

      bc-recipient-msa-internal:
        pools: 
          - 'bc_pool_I'
        plasmids:
          - 438-937
          - 439-937
        extract-barcode:
          438-937:
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            knowncodes: 'nuconfig/recipient-barcodes.fasta' 
          439-937:
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            knowncodes: 'nuconfig/recipient-barcodes.fasta' 
        pre-assembly:
          samlami: *samlami-internal
        assembly-method: msa
        post-assembly:
          cutoutdonor:
            itermae: 'nuconfig/itermae_BC_donor.yaml'
            target-fasta: 
              - 'nuconfig/pSL438_donor_bc_corrected_v3_revcomp.fasta'
              - 'nuconfig/pSL439_donor_bc_corrected_v3_revcomp.fasta'
            target-size: 'short'

      bc-recipient-flye-backbone:
        pools: 
          - 'bc_pool_I'
        plasmids:
        - 438-937
        - 439-937
        extract-barcode:
          438-937:
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            knowncodes: 'nuconfig/recipient-barcodes.fasta' 
          439-937:
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            knowncodes: 'nuconfig/recipient-barcodes.fasta' 
        pre-assembly:
        assembly-method: flye
        post-assembly:
          modulus: 
            samlami:
             - arg: '--modulus'
               ref: 'AGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCA'
               head: 100
            target-fasta: 
              - 'nuconfig/psl439_bc_x_psl937_bc-changed-to-reflect-empirical-sequencing.fasta'
              - 'nuconfig/psl438_bc_x_psl937_bc-updated-to-reflect-empirical-sequencing.fasta'
          donor-bc: 
            samlami: *samlami-internal
            itermae: 'nuconfig/itermae_BC_donor.yaml'
            target-fasta: 
              - 'nuconfig/pSL438_donor_bc_corrected_v3_revcomp.fasta'
              - 'nuconfig/pSL439_donor_bc_corrected_v3_revcomp.fasta'
            target-size: 'short'

    ## below is using internal to separate
      bc-internal-msa-recipient:
        pools: 
          - 'bc_pool_I'
        plasmids:
        - 438-937
        - 439-937
        extract-barcode:
          438-937:
            samlami: *samlami-internal
            itermae: 'nuconfig/itermae_BC_donor.yaml'
            knowncodes: 
              - 'nuconfig/pSL438_donor_bc_corrected_v3_revcomp.fasta'
          439-937:
            samlami: *samlami-internal
            itermae: 'nuconfig/itermae_BC_donor.yaml'
            knowncodes: 
              - 'nuconfig/pSL439_donor_bc_corrected_v3_revcomp.fasta'
        pre-assembly:
          samlami: *samlami-rec-bc
        assembly-method: msa
        post-assembly:
          rec-bc: 
            itermae: 'nuconfig/itermae_BC_937.yaml'
            target-fasta: 'nuconfig/recipient-barcodes.fasta' 
            target-size: 'short'

    # below is using internal to separate
      bc-internal-flye-backbone-438:
        pools: 
          - 'bc_pool_I'
        plasmids:
        - 438-937
        extract-barcode:
          438-937:
            samlami: *samlami-internal
            itermae: 'nuconfig/itermae_BC_donor.yaml'
            knowncodes: 
              - 'nuconfig/pSL438_donor_bc_corrected_v3_revcomp.fasta'
        pre-assembly:
        assembly-method: flye
        post-assembly:
          rec-bc: 
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            target-fasta: 'nuconfig/recipient-barcodes.fasta' 
            target-size: 'short'
          modulus: 
            samlami:
             - arg: '--modulus'
               ref: 'AGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCA'
               head: 100
            target-fasta: 
              - 'nuconfig/psl438_bc_x_psl937_bc-updated-to-reflect-empirical-sequencing.fasta'

      bc-internal-flye-backbone-439:
        pools: 
          - 'bc_pool_I'
        plasmids:
        - 439-937
        extract-barcode:
          439-937:
            samlami: *samlami-internal
            itermae: 'nuconfig/itermae_BC_donor.yaml'
            knowncodes: 
              - 'nuconfig/pSL439_donor_bc_corrected_v3_revcomp.fasta'
        pre-assembly:
        assembly-method: flye
        post-assembly:
          rec-bc: 
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            target-fasta: 'nuconfig/recipient-barcodes.fasta' 
            target-size: 'short'
          modulus: 
            samlami:
             - arg: '--modulus'
               ref: 'AGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCA'
               head: 100
            target-fasta: 
              - 'nuconfig/psl439_bc_x_psl937_bc-changed-to-reflect-empirical-sequencing.fasta'



LASSO analysis
--------------------------------------------

Command run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    nextflow main.nf -c nextflow.config \
        -resume -ansi-log false -with-dag reports/dag.html \
        -profile common,blocky \
        --slurm_queue high \
        -params-file config-lasso.yaml

BPS config yaml
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    pools:

      455_lasso_pool_I:
        plasmids:
          - 1064-937-lasso
        demux: |
          barcode73	455_lasso_s1
          barcode74	455_lasso_s2
          barcode75	455_lasso_s3

      455_lasso_bc_pool_I:
        plasmids:
          - 438-937
          - 439-937
          - 1064-937-lasso
        demux: |
          barcode76	bc_bc_1
          barcode77	bc_bc_2
          barcode73	455_lasso_s1
          barcode74	455_lasso_s2
          barcode75	455_lasso_s3
          barcode78	455_lasso_s4
          barcode79	455_lasso_s5
          barcode80	455_lasso_s6
          barcode81	455_lasso_s7
          barcode82	455_lasso_s8
          barcode83	455_lasso_s9
          barcode84	455_lasso_s10
          barcode85	455_lasso_s11
          barcode86	455_lasso_s12
          barcode87	455_lasso_s13
          barcode88	455_lasso_s14

      455_lasso_pool_II:
        plasmids:
        - 1064-937-lasso
        demux: |
          barcode1	455_lasso_s1
          barcode2	455_lasso_s2
          barcode3	455_lasso_s3
          barcode4	455_lasso_s4
          barcode5	455_lasso_s5
          barcode6	455_lasso_s6
          barcode7	455_lasso_s7
          barcode8	455_lasso_s8
          barcode9	455_lasso_s9
          barcode10	455_lasso_s10
          barcode11	455_lasso_s11
          barcode12	455_lasso_s12
          barcode13	455_lasso_s13
          barcode14	455_lasso_s14

      455_lasso_pool_III:
        plasmids:
        - 1064-937-lasso
        demux: |
          barcode17	455_lasso_s15
          barcode18	455_lasso_s16
          barcode19	455_lasso_s17
          barcode20	455_lasso_s18
          barcode21	455_lasso_s19
          barcode22	455_lasso_s20
          barcode23	455_lasso_s21
          barcode24	455_lasso_s22
          barcode25	455_lasso_s23
          barcode26	455_lasso_s24
          barcode27	455_lasso_s25
          barcode28	455_lasso_s26
          barcode29	455_lasso_s27
          barcode30	455_lasso_s28
          barcode31	455_lasso_s29
          barcode32	455_lasso_s30

    plasmids:

      438-937:
        signature: 'TAGTACGCGTccaCACGGAGGTGAGCCCGCCTGCCCGGGGCTAGCGTCG'

      439-937:
    # go diff??
        signature: 'TAGTACGCGTCCCGGGGCTAGCGTCGACTCTAGAGGATCGATCCTTTTT'

      1064-937-lasso:
        signature: 'CATCGAGGGCTAGAATTACCTACCGGCCTCCACCATGCCTGCGggcgcgccaccgctaagctcaaggtcacaaaagcAgACGACGGCCAGTgtcgacATG'

    runs:

      20221111_1544_MN21180_FAT59932_c3bac6f0:
        fastq: '/archive/l/seq_data/bps/20221111_1544_MN21180_FAT59932_c3bac6f0/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: '455_lasso_bc_pool_I'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20221113_1507_MN24328_FAT59932_3e555f90:
        fastq: '/archive/l/seq_data/bps/20221113_1507_MN24328_FAT59932_3e555f90/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: '455_lasso_bc_pool_I'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20221012_1434_MN21180_FAT60123_2f1f2cff:
        fastq: '/archive/l/seq_data/bps/20221012_1434_MN21180_FAT60123_2f1f2cff/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: '455_lasso_pool_I'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20230124_1146_MN24357_FAV20471_aa3aca5f:
        fastq: '/archive/l/seq_data/bps/20230124_1146_MN24357_FAV20471_aa3aca5f/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN114'
        pool: '455_lasso_pool_II'
        medaka-model: 'r1041_e82_400bps_sup_v4.2.0'
        header-parse: 'modified'

      20230329_1348_MN24357_FAV33302_7c8b2367:
        fastq: '/archive/l/seq_data/bps/20230329_1348_MN24357_FAV33302_7c8b2367/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN114'
        pool: '455_lasso_pool_III'
        medaka-model: 'r1041_e82_400bps_sup_v4.2.0'
        header-parse: 'modified'

      20230330_1541_MN24357_FAV33302_a9025c19:
        fastq: '/archive/l/seq_data/bps/20230330_1541_MN24357_FAV33302_a9025c19/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN114'
        pool: '455_lasso_pool_III'
        medaka-model: 'r1041_e82_400bps_sup_v4.2.0'
        header-parse: 'modified'


    #####
    #####
    ##### TODO conceptually rewrite so there's a top level ontology category of
    ##### 'assembly', so all the pre-assembly and assembly steps can be pooled
    ##### across 'experiments', which maybe should just be 'analyses' then
    #####
    #####

    samlami-rec-bc: &samlami-rec-bc
      - arg: '--cut-on-left'
        ref: 'ctcaagcaaggttttcagtataatgttacatgcgtacacgcgtctgtacagaaaaaaaagaaaaatttgaaatataaataacgttcttaatactaacata'
        tail: 100
      - arg: '--cut-on-right'
        ref: 'ttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
        head: 100

    samlami-internal: &samlami-internal
      - arg: '--cut-on-left'
        ref: 'cgagGGTACCttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
        tail: 100
      - arg: '--cut-on-right'
        ref: 'cacatatacctgccgttcactattatttagtgaaatgagatattatgatattttctgaattgtgattaaaaaggcaactttatgcccatgcaacagaaac'
        head: 100


    experiments: 

      455-lasso:
        filter-out: 'nuconfig/filter-out-r6k-psc101-937.fasta'
        pools:
          - '455_lasso_bc_pool_I'
          - '455_lasso_pool_I'
          - '455_lasso_pool_II'
          - '455_lasso_pool_III'
        plasmids:
          - 1064-937-lasso
        extract-barcode:
          1064-937-lasso:
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            knowncodes: 'nuconfig/recipient-barcodes.fasta' 
              # implement how to query known barcodes by sample
        pre-assembly:
          samlami: *samlami-internal
        assembly-method: msa
        post-assembly:
          cutoutdonor:
            itermae: 'nuconfig/itermae_payload_1064_lasso.yaml'
            target-fasta: 'nuconfig/lorenzo_orfs_fixed_230830.fasta'
            target-size: 'long'



barcode-barcode analysis
--------------------------------------------


Command run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    nextflow main.nf -c nextflow.config \
        -resume -ansi-log false -with-dag reports/dag.html \
        -profile common,blocky \
        --slurm_queue high \
        -params-file nuconfig/config-bc.yaml



BPS config yaml
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

yes



::

    pools:

      idt_op1_single_s1:
        plasmids:
        - 1071-937-op
        demux: |
          unclassified	idt_op1_s1

      idt_op1_single_s6:
        plasmids:
        - 1071-937-op
        demux: |
          unclassified	idt_op1_s6

      idt_op1_single_s7:
        plasmids:
        - 1071-937-op
        demux: |
          unclassified	idt_op1_s7

      idt_op1_single_s9:
        plasmids:
        - 1071-937-op
        demux: |
          unclassified	idt_op1_s9

      idt_op1_single_s8:
        plasmids:
        - 1071-937-op
        demux: |
          unclassified	idt_op1_s8

      idt_op1_single_s10:
        plasmids:
        - 1071-937-op
        demux: |
          unclassified	idt_op1_s10

      idt_op1_pool_I:
        plasmids:
        - 1071-937-op
        demux: |
          barcode02	idt_op1_s2
          barcode03	idt_op1_s3
          barcode04	idt_op1_s4
          barcode05	idt_op1_s5
          barcode07	idt_op1_s7
          barcode08	idt_op1_s8
          barcode09	idt_op1_s9
          barcode10	idt_op1_s10

      idt_op2_pool_I:
        plasmids:
        - 1071-937-op
        demux: |
          barcode10	idt_op2_s1
          barcode11	idt_op2_s2
          barcode12	idt_op2_s3
          barcode13	idt_op2_s4
          barcode14	idt_op2_s5
          barcode15	idt_op2_s6
          barcode16	idt_op2_s7
          barcode17	idt_op2_s8
          barcode18	idt_op2_s9
          barcode19	idt_op2_s10

      twist_op_s1-20:
        plasmids:
        - 1064-937-op
        demux: |
          barcode01	twist_op_s1
          barcode02	twist_op_s2
          barcode03	twist_op_s3
          barcode04	twist_op_s4
          barcode05	twist_op_s5
          barcode06	twist_op_s6
          barcode07	twist_op_s7
          barcode08	twist_op_s8
          barcode09	twist_op_s9
          barcode10	twist_op_s10
          barcode11	twist_op_s11
          barcode12	twist_op_s12
          barcode13	twist_op_s13
          barcode14	twist_op_s14
          barcode15	twist_op_s15
          barcode16	twist_op_s16
          barcode17	twist_op_s17
          barcode18	twist_op_s18
          barcode19	twist_op_s19
          barcode20	twist_op_s20

      twist_op_s1-35:
        plasmids:
        - 1064-937-op
        demux: |
          barcode01	twist_op_s1
          barcode02	twist_op_s2
          barcode03	twist_op_s3
          barcode04	twist_op_s4
          barcode05	twist_op_s5
          barcode06	twist_op_s6
          barcode07	twist_op_s7
          barcode08	twist_op_s8
          barcode09	twist_op_s9
          barcode10	twist_op_s10
          barcode11	twist_op_s11
          barcode12	twist_op_s12
          barcode13	twist_op_s13
          barcode14	twist_op_s14
          barcode15	twist_op_s15
          barcode16	twist_op_s16
          barcode17	twist_op_s17
          barcode18	twist_op_s18
          barcode19	twist_op_s19
          barcode20	twist_op_s20
          barcode21	twist_op_s21
          barcode22	twist_op_s22
          barcode23	twist_op_s23
          barcode24	twist_op_s24
          barcode25	twist_op_s25
          barcode26	twist_op_s26
          barcode27	twist_op_s27
          barcode28	twist_op_s28
          barcode29	twist_op_s29
          barcode30	twist_op_s30
          barcode31	twist_op_s31
          barcode32	twist_op_s32
          barcode33	twist_op_s33
          barcode34	twist_op_s34
          barcode35	twist_op_s35

      twist_op_s36-72:
        plasmids:
        - 1064-937-op
        demux: |
          barcode36	twist_op_s36
          barcode37	twist_op_s37
          barcode38	twist_op_s38
          barcode39	twist_op_s39
          barcode40	twist_op_s40
          barcode41	twist_op_s41
          barcode42	twist_op_s42
          barcode43	twist_op_s43
          barcode44	twist_op_s44
          barcode45	twist_op_s45
          barcode46	twist_op_s46
          barcode47	twist_op_s47
          barcode48	twist_op_s48
          barcode49	twist_op_s49
          barcode50	twist_op_s50
          barcode51	twist_op_s51
          barcode52	twist_op_s52
          barcode53	twist_op_s53
          barcode54	twist_op_s54
          barcode55	twist_op_s55
          barcode56	twist_op_s56
          barcode57	twist_op_s57
          barcode58	twist_op_s58
          barcode59	twist_op_s59
          barcode60	twist_op_s60
          barcode61	twist_op_s61
          barcode62	twist_op_s62
          barcode63	twist_op_s63
          barcode64	twist_op_s64
          barcode65	twist_op_s65
          barcode66	twist_op_s66
          barcode67	twist_op_s67
          barcode68	twist_op_s68
          barcode69	twist_op_s69
          barcode70	twist_op_s70
          barcode71	twist_op_s71
          barcode72	twist_op_s72

      twist_op_s1-30:
        plasmids:
          - 1064-937-op
        demux: |
          barcode25	twist_oligos_1
          barcode26	twist_oligos_2
          barcode27	twist_oligos_3
          barcode28	twist_oligos_4
          barcode29	twist_oligos_5
          barcode30	twist_oligos_6
          barcode31	twist_oligos_7
          barcode32	twist_oligos_8
          barcode33	twist_oligos_9
          barcode34	twist_oligos_10
          barcode35	twist_oligos_11
          barcode36	twist_oligos_12
          barcode37	twist_oligos_13
          barcode38	twist_oligos_14
          barcode39	twist_oligos_15
          barcode40	twist_oligos_16
          barcode73	twist_oligos_17
          barcode74	twist_oligos_18
          barcode75	twist_oligos_19
          barcode76	twist_oligos_20
          barcode77	twist_oligos_21
          barcode78	twist_oligos_22
          barcode79	twist_oligos_23
          barcode80	twist_oligos_24
          barcode81	twist_oligos_25
          barcode82	twist_oligos_26
          barcode83	twist_oligos_27
          barcode84	twist_oligos_28
          barcode85	twist_oligos_29
          barcode86	twist_oligos_30

    plasmids:

      1071-937-op:
        signature: 'tctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccCAGGCGGGCTCACCTCCGTGtggGCGGCCATggcgcgcc'

      1064-937-op:
        signature: 'gggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtggGCGGCCATCGAGGGCTAGAATTACCTACCGGCCTCCACCATGCCTGCGggcgcgcc'

    runs:

      20220204_1659_MN21180_AHY362_20fff81e:
        fastq: '/archive/l/seq_data/bps/20220204_1659_MN21180_AHY362_20fff81e/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-FLG001'
        pool: 'idt_op1_single_s1'
        medaka-model: 'r941_min_sup_g507'
        header-parse: 'guppy'

      20220210_2158_MN21180_AHV186_47bfe197:
        fastq: '/archive/l/seq_data/bps/20220210_2158_MN21180_AHV186_47bfe197/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-FLG001'
        pool: 'idt_op1_single_s6'
        medaka-model: 'r941_min_sup_g507'
        header-parse: 'guppy'

      20220216_2034_MN21180_ahy372_ed3cdbcb:
        fastq: '/archive/l/seq_data/bps/20220216_2034_MN21180_ahy372_ed3cdbcb/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-FLG001'
        pool: 'idt_op1_single_s7'
        medaka-model: 'r941_min_sup_g507'
        header-parse: 'guppy'

      20220221_1007_MN21180_aig649_5975f02f:
        fastq: '/archive/l/seq_data/bps/20220221_1007_MN21180_aig649_5975f02f/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-FLG001'
        pool: 'idt_op1_single_s8'
        medaka-model: 'r941_min_sup_g507'
        header-parse: 'guppy'

      20220217_1545_MN21180_aib485_85e73edd:
        fastq: '/archive/l/seq_data/bps/20220217_1545_MN21180_aib485_85e73edd/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-FLG001'
        pool: 'idt_op1_single_s9'
        medaka-model: 'r941_min_sup_g507'
        header-parse: 'guppy'

      20220222_1011_MN21180_aht720_3d126414:
        fastq: '/archive/l/seq_data/bps/20220222_1011_MN21180_aht720_3d126414/basecalling/output/*.fastq.gz'
        flowcell: 'FLO-FLG001'
        pool: 'idt_op1_single_s10'
        medaka-model: 'r941_min_sup_g507'
        header-parse: 'guppy'


      20220219_2042_MN21180_FAK29156_8060e67d:
        fastq: '/archive/l/seq_data/bps/20220219_2042_MN21180_FAK29156_8060e67d/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN106D'
        pool: 'idt_op1_pool_I'
        medaka-model: 'r941_min_sup_g507'
        header-parse: 'modified'

      20220220_1413_MN21180_FAK29156_1d6cb4e1:
        fastq: '/archive/l/seq_data/bps/20220220_1413_MN21180_FAK29156_1d6cb4e1/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN106D'
        pool: 'idt_op1_pool_I'
        medaka-model: 'r941_min_sup_g507'
        header-parse: 'modified'

      20220505_1256_MN21180_FAT28805_36dcfc86:
        fastq: '/archive/l/seq_data/bps/20220505_1256_MN21180_FAT28805_36dcfc86/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: 'idt_op2_pool_I'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20220506_1604_MN21180_FAT28805_e2260b2f:
        fastq: '/archive/l/seq_data/bps/20220506_1604_MN21180_FAT28805_e2260b2f/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: 'idt_op2_pool_I'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20220810_1419_MN24357_FAT11143_dc331389:
        fastq: '/archive/l/seq_data/bps/20220810_1419_MN24357_FAT11143_dc331389/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: 'twist_op_s1-20'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20220826_1442_MN24357_FAT60123_d777449a:
        fastq: '/archive/l/seq_data/bps/20220826_1442_MN24357_FAT60123_d777449a/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: 'twist_op_s1-35'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20220928_1634_MN21180_FAT60123_fc2fcf63:
        fastq: '/archive/l/seq_data/bps/20220928_1634_MN21180_FAT60123_fc2fcf63/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: 'twist_op_s36-72'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20220929_1648_MN21180_FAT60123_2c6fef82:
        fastq: '/archive/l/seq_data/bps/20220929_1648_MN21180_FAT60123_2c6fef82/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: 'twist_op_s36-72'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20220930_1046_MN21180_FAU29365_e667d692:
        fastq: '/archive/l/seq_data/bps/20220930_1046_MN21180_FAU29365_e667d692/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: 'twist_op_s36-72'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      20221003_1553_MN21180_FAU29365_de37c928: 
        fastq: '/archive/l/seq_data/bps/20221003_1553_MN21180_FAU29365_de37c928/basecalling/output/demuxed.fastq.gz'
        flowcell: 'FLO-MIN112'
        pool: 'twist_op_s36-72'
        medaka-model: 'r104_e81_sup_g610'
        header-parse: 'modified'

      plasmidsaurus_PAQ18854:
        fastq: '/archive/l/seq_data/bps/20230625_plasmidsaurus/Li_*/*fastq.gz'
        flowcell: 'FLO-PRO114M'
        pool: 'twist_op_s1-30'
        medaka-model: 'r1041_e82_400bps_sup_v4.2.0'
        header-parse: 'guppy'


    #####
    #####
    ##### TODO conceptually rewrite so there's a top level ontology category of
    ##### 'assembly', so all the pre-assembly and assembly steps can be pooled
    ##### across 'experiments', which maybe should just be 'analyses' then
    #####
    #####

    samlami-rec-bc: &samlami-rec-bc
      - arg: '--cut-on-left'
        ref: 'ctcaagcaaggttttcagtataatgttacatgcgtacacgcgtctgtacagaaaaaaaagaaaaatttgaaatataaataacgttcttaatactaacata'
        tail: 100
      - arg: '--cut-on-right'
        ref: 'ttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
        head: 100

    samlami-internal-1064: &samlami-internal-1064
      - arg: '--cut-on-left'
        ref: 'cgagGGTACCttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
        tail: 100
      - arg: '--cut-on-right'
        ref: 'cacatatacctgccgttcactattatttagtgaaatgagatattatgatattttctgaattgtgattaaaaaggcaactttatgcccatgcaacagaaac'
        head: 100

    samlami-internal-1071: &samlami-internal-1071
      - arg: '--cut-on-left'
        ref: 'cgagGGTACCttgccctctctcttcattcagggtcatgagaggcacgccattcaaggggagaagtgagatcggtaccGGGGCCACTAGGGACAGGATtgg'
        tail: 100
      - arg: '--cut-on-right'
        ref: 'atgtcacatctcgcagaactggttgccagtgcgaaggcggccattagccaggcgtcagatgttgccgcgttagataatgtgcgcgtcgaatatttgggta'
        head: 100

    experiments: 

      twist-oligos-bc:
        pools:
          - 'twist_op_s1-20'
          - 'twist_op_s1-35'
          - 'twist_op_s36-72'
          - 'twist_op_s1-30'
        filter-out: 'nuconfig/filter-out-r6k-psc101-937.fasta'
        plasmids:
        - 1064-937-op
        extract-barcode:
          1064-937-op:
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            knowncodes: 'nuconfig/recipient-barcodes.fasta' 
        pre-assembly:
          samlami: *samlami-internal-1064
        assembly-method: msa
        post-assembly:
          cutoutdonor:
            itermae: 'nuconfig/itermae_payload_1064.yaml'
            target-fasta: 'nuconfig/ref_1100_oligos_220215.fasta'
            target-size: 'short'

      idt-oligo-bc:
        pools: 
          - 'idt_op1_single_s1'
          - 'idt_op1_single_s6'
          - 'idt_op1_single_s7'
          - 'idt_op1_single_s8'
          - 'idt_op1_single_s9'
          - 'idt_op1_single_s10'
          - 'idt_op2_pool_I'
          - 'idt_op1_pool_I'
        filter-out: 'nuconfig/filter-out-r6k-psc101-937.fasta'
        plasmids:
          - 1071-937-op
        extract-barcode:
          1071-937-op:
            samlami: *samlami-rec-bc
            itermae: 'nuconfig/itermae_BC_937.yaml'
            knowncodes: 'nuconfig/recipient-barcodes.fasta' 
        pre-assembly:
          samlami: *samlami-internal-1071
        assembly-method: msa
        post-assembly:
          cutoutdonor:
            itermae: 'nuconfig/itermae_payload_1071.yaml'
            target-fasta: 'nuconfig/ref_1100_oligos_220215.fasta'
            target-size: 'short'









