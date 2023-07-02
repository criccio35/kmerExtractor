from itertools import product
import os
import re
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

def kmers_in_sequence(sequence,k):
    '''
    Calculate the frequency of k-mers in a given sequence.
    Only k-mers containing the four basic uppercase nucleotides (A, C, G, T) are accepted.
    Any substrings of length k that have any other character in their composition will be disregarded and thus not considered as k-mers.

    :param sequence: Nucleotide sequence.
    :type sequence: str
    :param k: Length of k-mers.
    :type k: int

    :return: Dictionary containing the k-mer frequencies. Keys are the k-mers and values are the corresponding frequencies.
    :rtype: dict
    '''
    kmers_list = [''.join(km) for km in product('ACGT', repeat = k)]
    kmers = {i:0 for i in kmers_list}
    key = ''
    for nt in sequence:
        if len(key)==k:
            if key in kmers:
                kmers[key] += 1
            key=key[1:]+nt
        else:
            key = key + nt
    # count last k-mer
    if (len(key) == k) and (key in kmers):
        kmers[key] += 1
    return kmers


def kmers_by_window(filepath, k, window_size=100_000, output_path=None):
    '''
    Calculate the frequency of k-mers by windows within each of the chromosome sequences.
    First, the entire sequence contained within the FASTA file is read.
    The sequences of each chromosome are then identified and separated.
    Subsequently, each chromosome is separated into windows of the specified size.
    Finally, the frequencies of the k-mers are calculated for each window.

    :param filepath: Location of the FASTA file containing the nucleotide sequence of an organism separated by chromosomes.
    :type filepath: str
    :param k: Length of k-mers.
    :type k: int
    :param window_size: Window size, in number of base pairs, to be used to count the frequency of k-mers within the chromosome.
                        Must be less than or equal to the length of the smallest chromosome.
                        Default to 100000.
    :type window_size: int
    :param output_path: Path to the directory where the output csv file will be saved.
                        If not specified, the output file will be saved in the same directory as the input file. Default to None
    :type output_path: str

    :return: Dataframe indicating the chromosome, the window index, the start and end positions of the window,
             and the frequencies of the corresponding k-mers.
    :rtype: pandas.DataFrame
    '''

    def hop(start, stop, step):
        for i in range(start, stop, step):
            yield i
        yield stop

    f_in = open(filepath, 'r')
    all_file = f_in.read()
    f_in.close()

    path, filename = os.path.split(filepath)
    file_noext = os.path.splitext(filename)[0]
    if output_path is not None:
        outputfile_path = output_path + '/CGRW_k{0}_{1}.csv'.format(k, file_noext)  # save in the specified directory
    else:
        outputfile_path = path + '/CGRW_k{0}_{1}.csv'.format(k, file_noext)  # save in the input directory

    f_out = open(outputfile_path, 'w')

    kmers_list = [''.join(km) for km in product('ACGT', repeat=k)]

    print('chromosome,window,start,end', end=',', file=f_out)
    print(*kmers_list, sep=',', file=f_out)

    seq_all = "".join(all_file.split("\n"))
    del all_file
    seq_bychr = re.split('>chr[1-9]*', seq_all)[1:]
    chr_names = re.findall('>chr[1-9]*', seq_all)
    chr_len = [len(seq) for seq in seq_bychr]
    del seq_all

    if window_size > min(chr_len):
        raise ValueError('Window size must be less than or equal to the length of the smallest chromosome.')

    for chr_ix in range(len(seq_bychr)):
        position = list(hop(0, len(seq_bychr[chr_ix]), window_size))
        window = 0
        for start, end in zip(position[:-1], position[1:]):
            chr_seq = seq_bychr[chr_ix][start:end]
            kmers_w = kmers_in_sequence(chr_seq, k)

            print(chr_names[chr_ix][1:] + ',' + str(window) + ',' + str(start + 1) + ',' + str(end), end=',',
                  file=f_out)
            print(*list(kmers_w.values()), sep=',', file=f_out)

            window += 1

    f_out.close()
    df = pd.read_csv(outputfile_path)

    return df


def kmers_by_window_opt(filepath, k, window_size=100_000, output_path=None):
    '''
    Calculate the frequency of k-mers by windows within each of the chromosome sequences.
    It calculates the k-mers, as it reads the FASTA file, optimizing space by avoiding storing the sequence.

    :param filepath: Location of the FASTA file containing the nucleotide sequence of an organism separated by chromosomes.
    :type filepath: str
    :param k: Length of k-mers.
    :type k: int
    :param window_size: Window size, in number of base pairs, to be used to count the frequency of k-mers within the chromosome.
                        Must be greater than the number of nucleotides per line in the input FASTA file. Default to 100000
    :type window_size: int
    :param output_path: Path to the directory where the output csv file will be saved.
                        If not specified, the output file will be saved in the same directory as the input file. Default to None
    :type output_path: str

    :return: Dataframe containing the k-mer frequencies for all windows of each chromosome.
    :rtype: pandas.DataFrame
    '''

    kmers_list = [''.join(km) for km in product('ACGT', repeat=k)]

    input_file = open(filepath, 'r')
    path, filename = os.path.split(filepath)
    file_noext = os.path.splitext(filename)[0]
    if output_path is not None:
        outputfile_path = output_path + '/CGRW_k{0}_{1}.csv'.format(k, file_noext)  # save in the specified directory
    else:
        outputfile_path = path + '/CGRW_k{0}_{1}.csv'.format(k, file_noext)  # save in the input directory
    f = open(outputfile_path, 'w')

    kmers = {i: 0 for i in kmers_list}
    #key = ''
    nts_in_line_prev = 0

    print('chromosome,window,start,end', end=',', file=f)
    print(*list(kmers.keys()), sep=',', file=f)

    for line in tqdm(input_file.readlines()):
        if line[0] == '>':  # new chromosome start
            window = 0
            count = 0
            start = 1
            end = 0
            chromosome = line[1:-1] # the first and last characters are '>' and '\n', respectively
            kmers = {i: 0 for i in kmers_list}
            key = ''

        else:  # if line[0] != '>':
            nts_in_line = len(line) - 1  # last character is '\n'
            count += nts_in_line

            if nts_in_line_prev == 0:
                nts_in_line_prev = nts_in_line

            if count <= window_size:
                stop = nts_in_line  # count k-mers until line end
            else:
                stop = nts_in_line - (count - window_size)  # count k-mers until window end
                count -= nts_in_line - stop

            for i in range(stop):
                nt = line[i].upper()
                end += 1
                if len(key) < k:
                    key = key + nt
                else:
                    if key in kmers:
                        kmers[key] += 1
                    key = key[1:] + nt

            if count >= window_size:  # window is over
                if key in kmers:
                    kmers[key] += 1  # add last k-mer in window
                # print window info and k-mer frequencies
                print(chromosome, window, start, end, *list(kmers.values()), sep=',', file=f)

                # new window info
                window += 1
                start = end + 1
                end = start - 1
                count = 0
                kmers = {i: 0 for i in kmers_list}
                key = ''

                # start new window counting in the rest of the line
                count += nts_in_line - stop
                for i in range(stop, len(line) - 1):
                    nt = line[i].upper()
                    end += 1
                    if len(key) < k:
                        key = key + nt
                    else:
                        if key in kmers:
                            kmers[key] += 1
                        key = key[1:] + nt

            if nts_in_line < nts_in_line_prev:  # chromosome is over but window is not
                if key in kmers:
                    kmers[key] += 1  # add last k-mer in window
                # print window info and k-mer frequencies
                print(chromosome, window, start, end, *list(kmers.values()), sep=',', file=f)

            nts_in_line_prev = nts_in_line

    if key in kmers:
        kmers[key] += 1  # last k-mer in window
    # print window info and k-mer frequencies
    print(chromosome, window, start, end, *list(kmers.values()), sep=',', file=f)

    f.close()
    df = pd.read_csv(outputfile_path)

    return df


def FCGR(kmers):
    '''
    Organize the k-mers and their associated values in the common structure of 
    Frequency Chaos Game Representation (FCGR), i.e., a matrix of dimensions N x N,
    where N = 4^(k/2); with the upper left corner equal to the k-mer of C, 
    the upper right corner the k-mer of G, the lower left corner the k-mer of A,
    and the lower right corner the k-mer of T.

    :param kmers: dictionary mapping the k-mers to a real value. All k-mers must be of the same length.
    :type kmers: dict

    :returns: FCGR matrix and the names of the k-mers in the matrix.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    '''
    k = len(list(kmers.keys())[0])
    N = 2**k
    counts = [[0 for _ in range(N)] for _ in range(N)]
    names = [[None for _ in range(N)] for _ in range(N)]
    maxx, maxy = N, N
    x, y = 0, 0
    for key,value in kmers.items():
        rkey = key[::-1]
        for char in rkey:
            if char == 'G': x += maxx >> 1     # maxx // 2
            elif char == 'A': y += maxy >> 1  # maxy // 2
            elif char == 'T' or char == 'U': 
                x += maxx >> 1
                y += maxy >> 1
            maxx >>= 1
            maxy >>= 1
        
        counts[y][x] = value  
        names[y][x] = key
        maxx, maxy = N, N
        x, y = 0, 0
 
    return np.array(counts), np.array(names)


def plot_kmers_across_windows(df, kmer_names, chromosome, figsize=(12,4), ax=None):
    '''
    Plot the frequencies of the given k-mers across the windows of a given chromosome.

    :param df: dataframe with the k-mer frequencies by window.
    :type df: pandas.DataFrame
    :param kmer_names: list of k-mer names to plot.
    :type kmer_names: list
    :param chromosome: name of the chromosome to plot.
    :type chromosome: str
    :param figsize: The size in inches of the figure to create. default to (12,4).
    :type figsize: tuple[int,int]
    :param ax: The Axes object containing the plot. If None, a new figure and axes is created. Default to None.
    :type ax: matplotlib.axes._axes.Axes

    :return: The Axes object containing the plot.
    :rtype: matplotlib.axes._axes.Axes
    '''
    window_size = df.loc[0,'end']
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)
    else:
        ax = ax

    for kmer in kmer_names:
        sns.lineplot(data=df[df.chromosome == chromosome],x='window',y=kmer, label=kmer, ax=ax)
    ax.set_title('Chromosome {}'.format(chromosome),fontsize=15)
    ax.set_xlabel('Chromosome position (x{} bp)'.format(window_size),fontsize=12)
    ax.set_ylabel('k-mer frequency',fontsize=12)
    ax.legend()

    return ax


def plot_kmers_freq_within_chromosomes(df, figsize=(15,6),ax=None):
    '''
    Plot the sum of the k-mer frequencies within each chromosome.

    :param df: dataframe with the k-mer frequencies by window.
    :type df: pandas.DataFrame
    :param figsize: The size in inches of the figure to create. default to (15,6).
    :type figsize: tuple[int,int]
    :param ax: Axes in which to draw the plot, otherwise use the currently-active Axes. Default to None.
    :type ax: matplotlib.axes.Axes

    :return: The Axes object containing the plot.
    :rtype: matplotlib.axes._axes.Axes
    '''
    df_agg = df.groupby('chromosome').agg({k:'sum' for k in df.columns[4:]}).T

    axi = df_agg.plot(kind='bar', stacked=True, figsize=figsize, ax=ax)

    axi.set_ylabel('k-mer frequency', fontsize=12)
    axi.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

    return axi


def plot_FCGR(df, chromosome=None, window=None, figsize=(7, 6), colormap='bwr', ax=None):
    '''
    Plot the Frequency Chaos Game Representation (FCGR) of the k-mers in the input dataframe.

    :param df: Dataframe with the k-mer frequencies by windows across chromosomes.
    :type df: pandas.DataFrame
    :param chromosome: Name of the chromosome to plot. If None, parameter window must be None,
                       plotting the sum of k-mers across all chromosomes. Default to None.
    :type chromosome: str
    :param window: Number of the window to plot. If None, plot the sum of all windows for each chromosome. Default to None.
    :type window: int
    :param figsize: The size in inches of the figure to create. Default to (7,6).
    :type figsize: tuple[int,int]
    :param colormap: Matplotlib colormap name. The mapping from data values to color space. Default to 'bwr'.
    :type colormap: str
    :param ax: Axes in which to draw the plot, otherwise use the currently-active Axes. Default to None.
    :type ax: matplotlib.axes.Axes

    :return: The Axes object containing the plot.
    :rtype: matplotlib.axes._axes.Axes
    '''
    kmer_names = df.columns[4:]
    if window is None:
        df = df.groupby('chromosome').agg({k: 'sum' for k in df.columns[4:]}).T
        if chromosome is None:
            vals = df.sum(1).values
        else:
            vals = df[chromosome].values
    else:
        if chromosome is None:
            raise ValueError('Parameter chromosome must be specified if window is not None.')
        df = df[(df.chromosome == chromosome) & (df.window == window)]
        vals = df[kmer_names].values[0]

    v_min = min(vals)
    v_max = max(vals)
    vals_norm = (vals - v_min) / (v_max - v_min + 1e-6)

    # dictionary with normalized values
    kmer2freq = {k: v for k, v in zip(kmer_names, vals_norm)}

    # FCGR matrices
    FCGR_vals, FCGR_names = FCGR(kmer2freq)

    # -------------------------------------------
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = ax

    sns.heatmap(FCGR_vals, annot=FCGR_names, fmt='',
                yticklabels=False,
                xticklabels=False,
                vmin=0, vmax=1,
                square=True,
                cmap=colormap, ax=ax)

    return ax


def CGR(sequence, d=1):
    '''
    Generate the X and Y coordinates of the Chaos Game Representation (CGR) for the input sequence.
    Purine nucleotides (A and G) are placed on the minor diagonal, and the pyrimidine nucleotides (C and T) on the main diagonal.
    The CGR is built in a square of dimensions 2d x 2d with center at the coordinates (0,0).
    Only the four basic uppercase nucleotides (A, C, G, T) are accepted. Any other character is ignored.

    :param sequence: Nucleotide sequence.
    :type sequence: str
    :param d: Half of the side of the square. Default to 1.
    :type d: int

    :return: The x and y coordinates of the CGR.
    :rtype: tuple[list[float],list[float]]
    '''
    ACGT = {'A': [0 - d, 0 - d], 'C': [0 - d, 0 + d], 'G': [0 + d, 0 + d], 'T': [0 + d, 0 - d]}

    X, Y = [0], [0]
    j = 0
    for i in tqdm(range(len(sequence))):
        nt = sequence[i]
        if nt in 'ACGT':
            x = (X[j] + ACGT[nt][0]) / 2
            y = (Y[j] + ACGT[nt][1]) / 2
            X.append(x)
            Y.append(y)
            j += 1

    return X, Y


def plot_CGR(sequence, d=1,markersize=1, ax=None,figsize=(6,6)):
    '''
    Plot the Chaos Game Representation (CGR) of the input sequence.
    Purine nucleotides (A and G) are placed on the minor diagonal, and the pyrimidine nucleotides (C and T) on the main diagonal.
    The CGR is built in a square of dimensions 2d x 2d with center at the coordinates (0,0).
    Only the four basic uppercase nucleotides (A, C, G, T) are accepted. Any other character is ignored.

    :param sequence: Nucleotide sequence.
    :type sequence: str
    :param d: Half of the side of the square in which the CGR is built. Default to 1.
    :type d: int

    :return: The Axes object containing the plot.
    :rtype: matplotlib.axes._axes.Axes
    '''

    X,Y = CGR(sequence, d)
    # -------------------------------------------
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = ax

    ax.plot(X, Y, 'ko', markersize=markersize)
    ax.set_xlim((-d, d))
    ax.set_ylim((-d, d))
    ax.set_xticks([])
    ax.set_yticks([])

    return ax
