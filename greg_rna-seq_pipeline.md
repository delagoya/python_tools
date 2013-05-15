# Greg's RNA-Seq pipeline

## (Work in progress)

1. run RUM/STAR
    * perhaps a wrapper script that takes in a dir to search for FASTQ sets
    * this can then write out a shell script for running all of the RUM/STAR jobs
1. filter out non-standard chromosomes
    * Split the files on a per-chromosome basis
1. split into uniques and non-uniques
    * add as an uption to above, same script
1. BAM location index
1. run mapping stats
    * likely a part of the above scripts, send to a CSV file.
1. resample uniques & non-uniques
    * resample pretty much done, need to run on multi-inputs of U & NU
    * should take mapping stats as input
    * should only need to define the target
1. concatenate
    * possibly this is already done.
1. make coverage plots
    * Current method is to re-run postprocess steps of RUM
        1. make junction calls
        1. sam2rum on sam files output from junctions script
        1. sort Unique
        1. sort NU
        1. coverage plots Unique
        1. coverage plots NU
1. merge coverage plots for each experimental group
1. get inferred internal exons
1. concatenate with annotated exons
1. quantify exons
1. get annotated introns
1. quantify introns
1. exons to spreadhseet
1. introns to spreadsheet
1. junctions to spreadsheet
1. make master spreadsheet
1. add annotation to spreadsheet
1. remove null features
1. coverage plots to bigwig
1. upload coverage and junctions plots to genome browser
    * either by sample or by condition
1. hierarchical clustering
1. differential expression

