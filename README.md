# Scripts for working with SAM/BAM 

These are a set of Python scripts for working with SAM/BAM files. 

## Dependencies

These scripts depend on python version >= 2.7.3 and the pysam and numpy packages. 

### PGFI Cluster setup

If you are on the PGFI cluster, you have the option of using the AnacondaCE python distribution as described below, or you can set up your initial working environment using `user modules` and `virtualenv`:


```shell
module load python-2.7.3
virtualenv $HOME/my_python  --system-site-packages
source $HOME/my_python/bin/activate
pip install pysam
```

To activate this environment for every session add these line to your `.bash_profile` or `.bashrc` file (whichever one you regularly use):

```shell
module load python-2.7.3
source $HOME/my_python/bin/activate
```


### AnacondaCE

By far the easiest way to get started is to use the AnacondaCE Python distribution from http://continuum.io/downloads.html

Once installed and activated, you should be able to use these scripts. 

Windows users may encounter some errors with pysam, in which case I wish you good luck and god speed, because there is nothing I can do about that.

## USAGE

### `sam_filter.py `

This script filters any read that does not align to a canonical chromosome. You can optionally write the mitochondrial chromosome alignments to another file. 

### `sam_sample.py`

This script downsamples a SAM/BAM file that is ordered by read ID to a close-to-exact count of reads. Required input is the number of total reads in the file and the desired target number.

### `pull_mapping_stats.py`

This script pulls the mapping statistics from a set of RUM alignments. Required is the file glob pattern to find the mapping_stats.txt file.



