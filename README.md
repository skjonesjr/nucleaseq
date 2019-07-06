# NucleaSeq

A package for the processing of NucleaSeq data, as described in the manuscript:

### Massively parallel kinetic profiling of natural and engineered CRISPR nucleases

**Stephen K. Jones Jr, John A. Hawkins, Nicole V. Johnson, Cheulhee Jung, Kuang Hu, James R. Rybarski, Janice S. Chen, Jennifer A. Doudna, William H. Press, and Ilya J. Finkelstein**

*bioRxiv* 

### Installation

The following instructions should work across platforms, except that installing virtualenv with apt-get is Ubuntu specific. For other platforms, install virtualenv appropriately if desired.

First, clone the repository to a local directory:

```
git clone https://github.com/hawkjo/nucleaseq.git
```

Optionally, you can install into a virtual environment (recommended):

```
sudo apt-get install -y virtualenv
cd nucleaseq
virtualenv envnucleaseq
. envnucleaseq/bin/activate
```

Now install dependencies:

```
pip install numpy==1.13.3 && pip install cython==0.24 && pip install -r requirements.txt
```

Note that some versions of pip require each line of the requirements.txt file be installed in a
separate call.

Next install the freebarcodes package in the same virtual environment, found at https://github.com/hawkjo/freebarcodes.

Finally, install nucleaseq with

```
python setup.py install
```

### Usage

```
Usage:
  nucleaseq setup_run <min_seq_len> <max_seq_len> <out_prefix> <sample_dirs>... [-v | -vv | -vvv]
  nucleaseq preprocess <uncut_or_cut> <targets_file> <target_name> <pamtarg_pos> <exploded_oligos_file> <max_primer_err> <max_bc_err> <read_names_by_seq_file> <read_names_by_sample_file> <out_prefix> <start> <inc> <large_inc> [--nprocs=<nprocs>] [-v | -vv | -vvv]
```

#### Coming soon:
Jupyter notebooks for convenient organization of output results and figures.
