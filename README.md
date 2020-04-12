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

### A note on memory

NucleaSeq data processing and analysis is memory intensive. All work for the original paper was
performed on a server with 256GB of RAM. Less is possible with reduced the parallelization, but
results in correspondingly slower runtimes. Memory usage also depends on library design, and
analysis of the original experiments required at least 64GB of RAM. For the first run, monitor
memory use carefully, and select the machine and parallelization selections accordingly.

### Usage

The NucleaSeq pipeline is a combination of command line invocations and Jupyter notebooks,
typically performed in the order below. The notebooks are located in the `notebooks` folder and
contain instructions for their use.

#### Library design
Notebook: `NucleaSeq Oligo Design.ipynb`

#### Raw data preprocessing
Command line:
```
  nucleaseq setup_run ...
  nucleaseq preprocess ...
```
Type `nucleaseq --help` for details.

#### Analysis
Notebooks, in the following order:
```
NucleaSeq Full Data Tabulation and Normalization Factors.ipynb
NucleaSeq Full Data Analysis.ipynb
NucleaSeq Cut Data Processing.ipynb
```

### CRISPR Model

The CRISPR Model described in the paper is available in this software package. Its use is explained
in the notebook `CRISPR Model Calculator Tutorial.ipynb`.
