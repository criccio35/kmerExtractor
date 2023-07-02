* Written by Camila Riccio, Mauricio Peñuela, Camilo Rocha and Jorge Finke
* Last update: 30/06/23 

kmerExtractor
=============

This is a Python3 implementation to computes the k-mer frequencies for windows within the chromosomes of an organism's sequence.

A genome sequence is composed of four basic types of nucleotides: adenine (A), cytosine (C), guanine (G), and thymine (T).
All possible nucleotide substrings of length  are called k-mers. The k-mer frequency can be used to unravel patterns within
genomic sequences. Usually these frequencies are calculated to characterize complete genomes. However, this approach can
provide a higher resolution of the genomic composition if frequencies are calculated for small windows within the sequence
of each chromosome. Additionally, this windowed resolution allows the frequency of k-mers to be used as input for machine
learning models that can help predict the behavior of other windowed properties within the chromosome.
The kmerExtractor module was designed to fulfill this purpose.

kmerExtractor takes as input a FASTA file with the nucleotide sequence of an organism separated by chromosomes.
It splits the chromosomes into windows of the desired size (e.g., 100 kb or 1 Mb).
For a given k, the kmerExtractor computes the frequency of the k-mers in the different windows within the chromosomes.
It returns a dataframe, and saves it as a csv file, indicating the chromosome, the window index, the start and end positions
of the window, and the frequencies of the corresponding k-mers.
A great advantage of this module is that it calculates and saves the frequencies of the k-mers as it reads the FASTA file,
optimizing space by avoiding storing the enrire sequence. So it is capable of handling very large genomes.

For visualization of the obtained results, the module can generate a lineplot for the frequency of a subset of k-mers along
a given chromosome, and a stacked barplot of the frequency of all k-mers discriminated by chromosomes.
Furthermore, it can generate the Chaos Game Representation (CGR) and the Frequency Chaos Game Representation (FCGR)
of the sequence, a chromosome, or a window within a chromosome.

Setup
------
Clone the repository::

  git clone git@github.com/criccio35/kmerExtractor


Requirements
------------
Install all of the modules listed in the Python requirements file into the project environment::

  pip install -r requirements.txt

How to use
----------

Make sure the kmerExtractor.py file is in the same folder
as the notebook or python script where you want to use the module.
Load the module as shown below::

  import kmerExtractor as kex

Find a complete example and the online documentation `here <https://criccio35.github.io/kmerExtractor/>`_
or in pdf `here <main/docs/kmerextractor.pdf>`_.


