import numpy as np
import pandas as pd
import os
import random
import re
import argparse
import openpyxl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from scipy.cluster.hierarchy import dendrogram, linkage
from openpyxl.styles import PatternFill
from collections import defaultdict
from multiprocessing import Pool


from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import plotly.graph_objects as go


def calculate_nucleotide_ratio(filepath, is_rna=False):
    """Calculate the nucleotide ratio of the given sequences and print the result.

    Parameters
    ----------
    filepath : string (/path/to/*.xlsx or .fasta)
        The file path of the sequence file. Supported formats are (.xlsx) and (.fasta).
    is_rna : bool (default: False)
        Set this to True if the input sequences are RNA sequences. Default is False for DNA sequences.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        Raised if an unsupported file format is provided.

    Examples
    --------
    >>> calculate_nucleotide_ratio('/path/to/sequences.fasta')
    A: 30.00%
    T: 20.00%
    G: 25.00%
    C: 25.00%

    >>> calculate_nucleotide_ratio('/path/to/sequences.xlsx', is_rna=True)
    A: 30.00%
    U: 20.00%
    G: 25.00%
    C: 25.00%
    """

    if filepath.endswith('.xlsx'):
        df = pd.read_excel(filepath) 
        sequences = df['Seq']
    elif filepath.endswith('.fasta'):
        sequences = []
        for record in SeqIO.parse(filepath, 'fasta'):  
            sequences.append(str(record.seq))
    else:
        raise ValueError('Unsupported file format. Only (.xlsx) and (.fasta) files are supported.')

    if is_rna:
        base_counts = {'A': 0, 'U': 0, 'G': 0, 'C': 0}  # Initialize base counts for RNA
    else:
        base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}  # Initialize base counts for DNA

    for sequence in sequences:
        for base in sequence:
            if base.upper() in base_counts:
                base_counts[base.upper()] += 1

    total_bases = sum(base_counts.values())
    base_ratios = {base: count / total_bases for base, count in base_counts.items()}

    for base, ratio in base_ratios.items():
        print(f'{base}: {ratio:.2%}')


def kmer_statistics(filepath, K, n_Kmers, output_filepath, weigh_sequences_by_count=False):
    """Perform k-mer statistics on the given sequences and write the result to an Excel file.
    If weigh_sequences_by_count is True, it will use the 'Count' column from the input .xlsx file or parse count from .fasta header to weigh the k-mer frequencies.

    Parameters
    ----------
    filepath : string (/path/to/*.xlsx or .fasta)
        The file path of the sequence file. Supported formats are (.xlsx) and (.fasta).
    K : int
        The length of the k-mers to consider.
    n_Kmers : int
        The number of top k-mers to include in the output.
    output_filepath : string (/path/to/*.xlsx)
        The file path to save the output Excel file.
    weigh_sequences_by_count : bool, optional
        If True, uses the 'Count' column from the input .xlsx file or parse count from .fasta header to weigh the k-mer frequencies.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        Raised if an unsupported file format is provided.

    Examples
    --------
    >>> kmer_statistics('/path/to/sequences.fasta', 3, 5, '/path/to/output.xlsx', weigh_sequences_by_count=True)
    """


    sequences = []
    counts = []
    original_counts = []  


    if filepath.endswith('.xlsx'):
        df = pd.read_excel(filepath)
        sequences = df['Seq'].tolist()
        original_counts = df['Count'].tolist()
        if weigh_sequences_by_count:
            counts = df['Count'].tolist()
        else:
            counts = [1] * len(sequences)
    elif filepath.endswith('.fasta'):
        with open(filepath, 'r') as fasta_file:
            for line in fasta_file:
                if line.startswith('>'):
                    parts = line.split('-')
                    count = int(parts[-1].strip()) if len(parts) > 1 else 1
                    counts.append(count if weigh_sequences_by_count else 1)
                    original_counts.append(count)
                else:
                    sequences.append(line.strip())
    else:
        raise ValueError('Unsupported file format. Only (.xlsx) and (.fasta) files are supported.')
    
    kmer_counts = defaultdict(int)

    for sequence, count in zip(sequences, counts):
        kmers_in_sequence = {sequence[i:i+K] for i in range(len(sequence) - K + 1)}
        for kmer in kmers_in_sequence:
            kmer_counts[kmer] += count

    sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
    sorted_kmers = sorted_kmers[:n_Kmers]

    kmer_sequences = {kmer: [] for kmer, _ in sorted_kmers}

    if filepath.endswith('.fasta'):
        for sequence, original_count in zip(sequences, original_counts):
            for kmer, _ in sorted_kmers:
                if kmer in sequence:
                    kmer_sequences[kmer].append((sequence, original_count))
    else:
        for sequence, original_count in zip(sequences, original_counts):
            for kmer, _ in sorted_kmers:
                if kmer in sequence:
                    kmer_sequences[kmer].append((sequence, original_count))


    wb = openpyxl.Workbook(write_only=False)  
    ws = wb.active
    ws.append(['Kmer', 'Frequency', 'Sequences'])

    for kmer, count in sorted_kmers:
        sequences_info = kmer_sequences[kmer]
        sequence_str = ", ".join([f"{seq}({cnt})" for seq, cnt in sequences_info])
        ws.append([kmer, count, sequence_str])

    wb.save(output_filepath)


    # Plotting the results
    html_colors = ['#CECCE5', '#D1C7AE', '#F18C25', '#ED6F6E', '#C6E3D8', '#4FBCA7', '#EDFB9A', '#E64825', '#FAECA8', '#B0D9A5']
    df = pd.read_excel(output_filepath)
    top10_df = df.head(10)

    df2 = pd.read_excel(filepath)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 10), dpi=100, gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.1}, sharex=True)

    bar_width = 0.37
    bars = ax1.bar(top10_df.index, top10_df['Frequency'], color=np.random.choice(html_colors, len(top10_df), replace=False), width=bar_width)

    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2, height, str(height),
                 ha='center', va='bottom', fontsize=10, fontweight='bold', color='black',
                 bbox=dict(facecolor='white', edgecolor='black', pad=3))

    ax1.set_xticks(range(len(top10_df)))
    ax1.set_xticklabels([])

    ax1.set_ylabel('Frequency', fontsize=14, fontweight='bold', labelpad=32)
    ax1.set_title('Kmer_statistics', fontsize=18, fontweight='bold')

    ax1.set_ylim(10, ax1.get_ylim()[1])

    ax1.tick_params(axis='y', labelsize=12)

    monospace_font = {'family': 'monospace', 'weight': 'bold', 'size': 'large'}
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=top10_df['Kmer'][i],
                                  markerfacecolor=bars[i].get_facecolor(), markersize=10) for i in range(len(top10_df))]
    ax1.legend(handles=legend_elements, bbox_to_anchor=(1.03, 1.016), loc='upper left', facecolor='white', edgecolor='black', shadow=True)

    line_values = []
    for kmer in top10_df['Kmer']:
        found_value = df2['Seq'].str.contains(kmer).idxmax() if df2['Seq'].str.contains(kmer).any() else np.nan
        line_values.append(found_value + 1 if pd.notna(found_value) else np.nan)

    for i, kmer in enumerate(top10_df['Kmer']):
        ax2.plot(i, line_values[i], marker='o', color=bars[i].get_facecolor(), markersize=8, linestyle='-')
    ax2.plot(range(len(top10_df)), line_values, linestyle='-', color='black')

    ax2.set_ylabel('Rank', fontsize=14, fontweight='bold', labelpad=35)
    ax2.set_ylim(50, 1000)
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))
    ax2.tick_params(axis='y', labelsize=12)

    ax2.set_xticks(range(len(top10_df)))
    ax2.set_xticklabels([])
    ax2.set_xlabel('Top10  K-mers', fontsize=14, fontweight='bold', labelpad=10)

    plt.tight_layout()
    plt.savefig(output_filepath.replace('.xlsx', '_plot.png'))  # Save plot as PNG
    plt.show()


def calculate_common_kmers(file_paths, output_path):
    """Calculate the common k-mers and their frequencies across multiple files and save the results to an Excel file.
    
    Parameters
    ----------
    file_paths : list of string
        The list of file paths of the k-mer statistics Excel files.
    output_path : string (/path/to/*.xlsx)
        The file path to save the output Excel file.
    
    Returns
    -------
    None
    
    Examples
    --------
    >>> calculate_common_kmers(['/path/to/file1.xlsx', '/path/to/file2.xlsx'], '/path/to/output.xlsx')
    """
    
    kmers = []
    freqs = []
    column_names = []
    for i, path in enumerate(file_paths):
        df = pd.read_excel(path)
        kmers.append(df['Kmer'].tolist())
        freqs.append(df['Frequency'].tolist())
        col_name = f'Counts in file{i+1}(Rank)'
        column_names.append(col_name)

    common_kmers = set(kmers[0]).intersection(*kmers[1:])

    # Calculate the frequencies and ranks for the common k-mers in each file
    counts = []
    for kmer in common_kmers:
        freq_rank_list = []
        for i in range(len(file_paths)):
            freq_index = kmers[i].index(kmer)
            freq_value = freqs[i][freq_index]
            rank = freq_index + 1
            freq_rank = f"{freq_value} ({rank})"
            freq_rank_list.append(freq_rank)
        total_freq = sum(freqs[j][kmers[j].index(kmer)] for j in range(len(file_paths)))
        counts.append([kmer] + freq_rank_list + [total_freq])

    column_names = ['Kmer'] + column_names + ['Total Count']
    df_counts = pd.DataFrame(counts, columns=column_names)
    df_counts = df_counts.sort_values(by=['Total Count'], ascending=False)
    df_counts = df_counts.reset_index(drop=True)
    df_counts = df_counts.drop(columns=['Total Count'])
    df_counts.to_excel(output_path, index=False)
    print(f'Successfully saved common_kmers to {output_path}.')


    # Plotting
    df = pd.read_excel(output_path, nrows=10)
    kmer_names = df.iloc[:, 0].tolist()  

    counts_list = []
    ranks_list = []
    for col in df.columns[1:]:
        counts = df[col].str.extract(r'(\d+)')[0].astype(int)
        ranks = df[col].str.extract(r'\((\d+)\)')[0].astype(int)
        counts_list.append(counts)
        ranks_list.append(ranks)

    html_colors = ['#ED6F6E', '#B0D9A5', '#CECCE5', '#F18C25', '#EDFB9A', '#FAECA8', '#E64825', '#D1C7AE', '#C6E3D8', '#4FBCA7']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9), dpi=100, sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    bar_width = 0.8 / len(counts_list)  
    index = range(len(kmer_names))

    for i, (counts, color) in enumerate(zip(counts_list, html_colors)):
        bars = ax1.bar([j - 0.4 + i * bar_width for j in index], counts, bar_width, label=f'Shared K-mer in file {i+1}', color=color)
        for rect in bars:
            height = rect.get_height()
            ax1.text(rect.get_x() + rect.get_width() / 2., height, '%d' % int(height),
                     ha='center', va='bottom', fontsize=10, fontweight='bold', color='black', bbox=dict(facecolor='white', edgecolor='black', pad=3))

    ax1.set_ylabel('Frequency', fontsize=14, fontweight='bold', labelpad=32)
    ax1.set_title('Calculate_common_Kmers', fontsize=18, fontweight='bold')
    ax1.tick_params(axis='y', labelsize=12)
    ax1.legend(fontsize=12, facecolor='white', edgecolor='black', shadow=True)
    ax1.set_ylim(0, ax1.get_ylim()[1])

    for i, (ranks, color) in enumerate(zip(ranks_list, html_colors)):
        ax2.plot(index, ranks, marker='o', linestyle='-', color=color, linewidth=2, label=f'File{i+1} Rank')

    ax2.set_xlabel('Shared  K-mers', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Rank', fontsize=14, fontweight='bold', labelpad=35)
    ax2.tick_params(axis='y', labelsize=12)
    ax2.yaxis.set_major_locator(MultipleLocator(5))
    ax2.set_ylim(0, 33)

    ax2.set_xticks([i for i in index])
    ax2.set_xticklabels(kmer_names, rotation=0, fontsize=20)
    ax1.set_xticklabels([])  
    ax2.legend().remove()

    plt.tight_layout()
    plt.show()


def calculate_common_seqs(file_paths, output_path):
    """Calculate the common sequences and their frequencies across multiple files, save the result to an Excel file, and visualize the results.

    Parameters
    ----------
    file_paths : list of string
        The list of file paths of the sequence statistics Excel files or fasta files.
    output_path : string
        The file path to save the output Excel file.

    Returns
    -------
    None

    Examples
    --------
    >>> calculate_common_seqs(['/path/to/file1.xlsx', '/path/to/file2.fasta'], '/path/to/output.xlsx')
    """
    seqs = []
    counts = []
    column_names = []
    sequences_dict_list = []
    row_indices = []

    for i, path in enumerate(file_paths):
        if path.endswith('.xlsx'):
            df = pd.read_excel(path)
            seqs.append(df['Seq'].tolist())
            counts.append(df['Count'].tolist())
            row_indices.append(df.index.tolist())
        elif path.endswith('.fasta'):
            sequences_dict = defaultdict(int)
            for idx, record in enumerate(SeqIO.parse(path, 'fasta')):
                freq = int(record.id.split('-')[1])
                sequences_dict[str(record.seq)] += freq
            seqs.append(list(sequences_dict.keys()))
            counts.append(list(sequences_dict.values()))
            sequences_dict_list.append(sequences_dict)
            row_indices.append(list(range(len(sequences_dict))))
        else:
            raise ValueError('Unsupported file format. Only (.xlsx) and (.fasta) files are supported.')

        col_name = f'Counts in file{i+1}'
        column_names.append(col_name)

    common_seqs = set(seqs[0]).intersection(*seqs[1:])

    seq_counts = []
    for seq in common_seqs:
        seq_count_list = []
        for i in range(len(file_paths)):
            if file_paths[i].endswith('.xlsx'):
                count_index = seqs[i].index(seq)
                count_value = counts[i][count_index]
                row_index = row_indices[i][count_index] + 1
                formatted_value = f"{count_value}({row_index})"
                seq_count_list.append(formatted_value)
            elif file_paths[i].endswith('.fasta'):
                sequences_dict = sequences_dict_list[i]
                count_value = sequences_dict[seq]
                row_index = list(sequences_dict.keys()).index(seq) + 1
                formatted_value = f"{count_value}({row_index})"
                seq_count_list.append(formatted_value)
        
        total_count = sum(counts[j][seqs[j].index(seq)] for j in range(len(file_paths)))
        seq_counts.append([seq] + seq_count_list + [total_count])

    column_names = ['Seq'] + column_names + ['Total Counts']
    df_seq_counts = pd.DataFrame(seq_counts, columns=column_names)
    df_seq_counts = df_seq_counts.sort_values(by=['Total Counts'], ascending=False)
    df_seq_counts = df_seq_counts.reset_index(drop=True)

    df_seq_counts.to_excel(output_path, index=False)
    print(f'Successfully saved common_seqs to {output_path}')


    df = pd.read_excel(output_path, nrows=10)

    for i in range(len(file_paths)):
        df[f'Counts in file{i+1}'] = df[f'Counts in file{i+1}'].astype(str)
        df[f'A_file{i+1}'] = df[f'Counts in file{i+1}'].str.extract(r'(\d+)\(').astype(int)
        df[f'B_file{i+1}'] = df[f'Counts in file{i+1}'].str.extract(r'\((\d+)\)').astype(int)

    x_labels = [f'{tuple(df[f"B_file{i+1}"][j] for i in range(len(file_paths)))}' for j in range(len(df))]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 9), dpi=100, sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    bar_width = 0.8 / len(file_paths)
    index = range(len(df))

    colors = ['#E76254', '#376795', '#4CAF50', '#FFC107', '#9C27B0', '#03A9F4', '#FF5722', '#795548']
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p']

    for i in range(len(file_paths)):
        ax1.bar([j - (0.4 - i * bar_width) for j in index], df[f'A_file{i+1}'], bar_width, label=f'Shared Seq in file{i+1}', color=colors[i % len(colors)])

    ax1.set_ylabel('Counts', fontsize=14, fontweight='bold')
    ax1.set_title('Calculate_common_seqs', fontsize=18, fontweight='bold')
    ax1.tick_params(axis='y', labelsize=12)
    ax1.legend(fontsize='large', facecolor='white', edgecolor='black', shadow=True)

    for i in range(len(file_paths)):
        for j, rect in enumerate(ax1.patches[i * len(df):(i + 1) * len(df)]):
            height = rect.get_height()
            ax1.text(rect.get_x() + rect.get_width() / 2., height, '%d' % int(height),
                     ha='center', va='bottom', fontsize=10, fontweight='bold', color='black', bbox=dict(facecolor='white', edgecolor='black', pad=3))

    for i in range(len(file_paths)):
        ax2.plot(index, df[f'B_file{i+1}'], marker=markers[i % len(markers)], linestyle='-', color=colors[i % len(colors)], linewidth=2, label=f'Ranks in file{i+1}')

    ax2.set_xlabel('Shared Sequences', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Ranks', fontsize=14, fontweight='bold', labelpad=32)
    ax2.tick_params(axis='y', labelsize=12)
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax2.set_ylim(0, 20)
    ax2.set_xticks([i for i in index])
    ax2.set_xticklabels(x_labels, rotation=0, fontsize=10)
    ax1.set_xticklabels([]) 
    ax2.legend().remove()

    plt.tight_layout()
    plt.show()


def kmer_concatenation(filepath, output_path, sheet_name):
    """Concatenate kmers based on certain rules and save the results to an Excel file.

    Parameters
    ----------
    filepath : string (/path/to/*.xlsx)
        The file path of the input Excel file.
    output_path : string (/path/to/*.xlsx)
        The file path to save the output Excel file.
    sheet_name : string
        The name of the sheet in the input Excel file.

    Returns
    -------
    None

    Examples
    --------
    >>> kmer_concatenation('/path/to/input.xlsx', '/path/to/output.xlsx', 'Sheet1')
    """

    df = pd.read_excel(filepath, sheet_name=sheet_name, nrows=50)
    kmer_list = df['Kmer'].tolist()
    concatenated_sequences = []
    
    for i in range(len(kmer_list)):
        for j in range(len(kmer_list)):
            if i != j:
                kmer1 = kmer_list[i]
                kmer2 = kmer_list[j]
                if kmer2.startswith(kmer1[1:]):
                    sequence = kmer1 + kmer2[len(kmer1) - 1:]
                    concatenated_sequences.append((sequence, kmer1, kmer2))
    
    print('The results of the first round of concatenating：')
    for sequence, kmer1, kmer2 in concatenated_sequences:
        print(sequence + '\t' + kmer1 + '+' + kmer2)
    
    result_df = pd.DataFrame({'Concatenating results': [sequence for sequence, _, _ in concatenated_sequences],
                              'Kmers for concatenating': [kmer1 + '+' + kmer2 for _, kmer1, kmer2 in concatenated_sequences]})
    
    # Perform subsequent rounds of concatenating until no more sequences can be concatenated
    round_num = 2
    while concatenated_sequences:
        print('-------------------------------')
        print('The results of the {} round of concatenating:'.format(round_num))
        new_kmer_list = []
        for sequence, kmer1, kmer2 in concatenated_sequences:
            new_kmer_list.append(sequence)
        
        # Reset the concatenated_sequences list and perform concatenation on the new_kmer_list
        concatenated_sequences = []
        for i in range(len(new_kmer_list)):
            for j in range(len(new_kmer_list)):
                if i != j:
                    kmer1 = new_kmer_list[i]
                    kmer2 = new_kmer_list[j]
                    if kmer2.startswith(kmer1[1:]):
                        sequence = kmer1 + kmer2[len(kmer1) - 1:]
                        concatenated_sequences.append((sequence, kmer1, kmer2))
        
        for sequence, kmer1, kmer2 in concatenated_sequences:
            print(sequence + '\t' + kmer1 + '+' + kmer2)
        
        round_df = pd.DataFrame({'Concatenating results': [sequence for sequence, _, _ in concatenated_sequences],
                                 'Kmers for concatenating': [kmer1 + '+' + kmer2 for _, kmer1, kmer2 in concatenated_sequences]})
        result_df = result_df._append(round_df, ignore_index=True)
        round_num += 1
        
        if len(concatenated_sequences) == 1:
            break
    
    with pd.ExcelWriter(output_path) as writer:
        result_df.to_excel(writer, index=False)


def levenshtein_distance(s, t):
    """Calculate the Levenshtein distance between two strings.

    Parameters
    ----------
    s : str
        The first string.
    t : str
        The second string.

    Returns
    -------
    int
        The Levenshtein distance between the two input strings.

    Examples
    --------
    >>> levenshtein_distance("ATGCTA", "AGTCTA")
    2
    """
    m, n = len(s), len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
        
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(dp[i][j - 1], dp[i - 1][j], dp[i - 1][j - 1])
                
    return dp[m][n]


def cluster_tree(filepath):
    """Generate a clustering dendrogram based on Levenshtein Edit Distances (LEDs) between sequences and/or affinity values.

    Parameters
    ----------
    filepath : str (/path/to/*.xlsx)
        The file path of the Excel file containing sequences and their corresponding IDs, along with optional 'Affinity' values.

    Returns
    -------
    None

    Examples
    --------
    >>> cluster_tree('/path/to/file.xlsx')
    """


    data = pd.read_excel(filepath)
    
    if 'Affinity' not in data.columns:  # Check if 'Affinity' column is present
        sequences = data['Seq'].tolist()
        labels = data['ID'].tolist()

        dist_matrix = np.zeros((len(sequences), len(sequences)))
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                dist = levenshtein_distance(sequences[i], sequences[j])
                dist_matrix[i][j] = dist
                dist_matrix[j][i] = dist

        linkage_matrix = linkage(dist_matrix, method='average')

        fig, ax = plt.subplots(figsize=(24, 60))
        dendrogram(linkage_matrix, labels=labels, orientation='right', leaf_font_size=24,color_threshold=50)
        ax.set_xlabel('LEDs',fontdict={'fontsize':22})
        ax.set_title('Clustering dendrogram based on Levenshtein Edit Distances(LEDs)', fontdict={'fontsize':22,'style':'italic'})
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=6)

        plt.show()
    
    else:
        sequences = data['Seq'].tolist()
        labels = data['ID'].tolist()
        affinities = [float(aff.strip('%')) for aff in data['Affinity']]  

        max_affinity = max(affinities)
        min_affinity = min(affinities)

        label_dict = {}
        affinity_list = []
        for i, label in enumerate(labels):
            label_dict[label] = affinities[i]
            affinity_list.append(affinities[i])

        max_affinity = max(affinities)
        min_affinity = min(affinities)
                
        label_dict = {}
        affinity_list = []  # Add a list to store the affinity values of all nodes
        for i, label in enumerate(labels):
            label_dict[label] = affinities[i]
            affinity_list.append(affinities[i])

        dist_matrix = np.zeros((len(sequences), len(sequences)))
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                dist = levenshtein_distance(sequences[i], sequences[j])
                dist_matrix[i][j] = dist
                dist_matrix[j][i] = dist

        linkage_matrix = linkage(dist_matrix, method='average')
        color_threshold = 16  
        fig, ax = plt.subplots(figsize=(13, 13))
        dendrogram(linkage_matrix, labels=labels, orientation='right', leaf_font_size=20, color_threshold=color_threshold)


        lbls = ax.get_ymajorticklabels()
        heat_map = np.zeros((len(labels), len(labels)))

        for i, lbl in enumerate(lbls):
            x = -4.75  
            y = lbl.get_position()[1] 
            label = lbl.get_text()
            if label in label_dict:
                aff = label_dict[label]
                bar_length = (aff - min_affinity) / (max_affinity - min_affinity) * 10  
                cmap = plt.cm.get_cmap('Greens')  
                norm = mcolors.Normalize(vmin=min_affinity, vmax=max_affinity * 1.1)  
                color = cmap(norm(aff)) 
                ax.barh(y, bar_length, left=0, height=5, color=color)  
                index = labels.index(label)
                heat_map[index, :] = (aff - min_affinity) / (max_affinity - min_affinity)  
            if lbl.get_text() in ['R15-10', 'R15-5', 'R15-13', 'R15-6', 'R15-56', 'R15-8', 'R15-2']:  # Replace with the IDs you want to change colors
                lbl.set_color('red')  

        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, orientation='vertical', shrink=1, pad=0.05)
        cbar.ax.tick_params(labelsize=20)
        cbar.set_label('Fluorescence Recovery(%)', fontsize=20)

        im = ax.imshow(heat_map, cmap='Greens', aspect='auto', extent=[-5.5, -4.5, 0, len(labels)], alpha=0.7, vmin=0, vmax=1)
        ax.set_xlabel('LED',fontdict={'fontsize':28})
        ax.set_title('Clustering dendrogram based on Levenshtein Edit Distance(LED)', fontdict={'fontsize':28,'style':'italic'})
        ax.tick_params(axis='x', labelsize=20)

        plt.xticks(np.arange(0, 45, 5))
        plt.xlim(0, 35)

        plt.show()


def kmer_evolution(files, evolution_count, output_file=None):
    """Perform evolution analysis on k-mers between multiple files and save the results to an Excel file.

    Parameters
    ----------
    files : list of string
        The list of file paths of the k-mer statistics Excel files.
    evolution_count : int, optional
        The maximum number of evolutions allowed for a k-mer to be considered as evolved. Default is 1.
    output_file : string (/path/to/*.xlsx), optional
        The file path to save the output Excel file. If not provided, the results will not be saved.

    Returns
    -------
    None

    Examples
    --------
    >>> kmer_evolution(['/path/to/file1.xlsx', '/path/to/file2.xlsx'], evolution_count=2, output_file='/path/to/output.xlsx')
    """

    all_results_df = pd.DataFrame()
    
    for i in range(len(files)-1):
        file1 = files[i]
        file2 = files[i+1]

        file1_name = os.path.splitext(os.path.basename(file1))[0]
        file2_name = os.path.splitext(os.path.basename(file2))[0]

        file1_data = pd.read_excel(file1)
        file1_kmers = file1_data["Kmer"].tolist()
        file1_freqs = file1_data["Frequency"].tolist()

        file2_data = pd.read_excel(file2)
        file2_kmers = file2_data["Kmer"].tolist()
        file2_freqs = file2_data["Frequency"].tolist()

        output_rows = []
        
        # Compare k-mers between the two files
        for kmer in file1_kmers:
            if len(kmer) > 1 and kmer not in file2_kmers:
                for j in range(len(file2_kmers)):
                    evolved_kmer = None
                    if (
                        len(file2_kmers[j]) > 1
                        and len(file2_kmers[j]) == len(kmer)
                        and file2_kmers[j] not in file1_kmers
                    ):
                        evolution_count_cur = 0
                        evolution_pos = []
                        for pos in range(len(kmer)):
                            if kmer[pos] != file2_kmers[j][pos]:
                                evolution_count_cur += 1
                                evolution_pos.append(pos)
                                if evolution_count_cur > evolution_count:
                                    break
                        if evolution_count_cur == evolution_count:
                            evolved_kmer = file2_kmers[j]
                            evolution_info = ", ".join([f"Position {pos+1}: {kmer[pos]}->{evolved_kmer[pos]}" for pos in evolution_pos])

                            output_rows.append([
                                "".join(kmer),
                                file1_freqs[file1_kmers.index(kmer)],
                                "".join(evolved_kmer),
                                file2_freqs[j],
                                evolution_info
                            ])

        results_df = pd.DataFrame(output_rows, columns=[
            f"Raw Kmer in {file1_name}",
            f"Raw Frequency in {file1_name}",
            f"Evolved Kmer in {file2_name}",
            f"Evolved Frequency in {file2_name}",
            f"Evolution Info {file1_name}~{file2_name}"
        ])
        results_df.drop_duplicates(subset=[f"Raw Kmer in {file1_name}"], inplace=True)
        results_df.reset_index(drop=True, inplace=True) 

        all_results_df = pd.concat([all_results_df, results_df], axis=1)


    if output_file:   
        writer = pd.ExcelWriter(output_file, engine="openpyxl") 
        all_results_df.to_excel(writer, index=False)
        
        worksheet = writer.sheets['Sheet1']
        
        for i in range(6, len(all_results_df.columns), 6):
            worksheet.insert_cols(i)
        writer._save()

    print("Done for kmer_evolution!")


def kmer_mergence(input_path, output_path):
    """Merge k-mers between adjacent groups in an input Excel file and save the results to an output Excel file.

    Parameters
    ----------
    input_path : string (/path/to/*.xlsx)
        The file path of the input Excel file containing k-mer statistics.
    output_path : string (/path/to/*.xlsx)
        The file path to save the output Excel file.

    Returns
    -------
    None

    Examples
    --------
    >>> kmer_mergence('/path/to/input.xlsx', '/path/to/output.xlsx')
    """
    df = pd.read_excel(input_path)

    group_columns = [df.columns[i:i+6] for i in range(0, df.shape[1], 6)]
    first_group = group_columns[0]
    last_group = group_columns[-1]

    first_column = first_group[0]
    last_column = last_group[2]

    merge_kmers = []
    first_freqs = []
    last_freqs = []

    for i in range(len(group_columns)-1):
        current_group = group_columns[i]
        next_group = group_columns[i+1]

        current_kmers = df[current_group[2]].tolist()
        next_kmers = df[next_group[0]].tolist()

        current_freqs = df[current_group[3]].tolist()
        next_freqs = df[next_group[1]].tolist()

        merged_kmers = list(set(current_kmers + next_kmers))
        merge_kmers.append(pd.Series(merged_kmers, name=f'Merged Kmer in {current_group[2]}~{next_group[0]}'))

        first_freqs.extend(current_freqs)
        last_freqs.extend(next_freqs)

    merge_series = pd.concat(merge_kmers).drop_duplicates()

    first_freq_series = pd.Series(first_freqs, name="Merged Frequency in First Group")
    last_freq_series = pd.Series(last_freqs, name="Merged Frequency in Last Group")

    merged_frequencies = []

    for index, row in merge_series.items():
        frequencies = []
        for i in range(len(group_columns)-1):
            current_group = group_columns[i]
            next_group = group_columns[i+1]

            current_freqs = df[df[current_group[2]] == row][current_group[3]].values
            next_freqs = df[df[next_group[0]] == row][next_group[1]].values

            if len(current_freqs) > 0:
                frequencies.extend(current_freqs)
            elif len(next_freqs) > 0:
                frequencies.extend(next_freqs)

        if frequencies:
            merged_frequencies.append(frequencies[0])  # Only keep the first frequency value found
        else:
            merged_frequencies.append(0)

    merge_series = merge_series.reset_index(drop=True).rename("Merged Kmer")
    merge_frequencies = pd.Series(merged_frequencies, name="Merged Frequency")

    merged_df = pd.concat([merge_series, merge_frequencies], axis=1)

    first_freq_col = df[first_group[1]]
    last_freq_col = df[last_group[3]]

    merged_df = merged_df.sort_values(by="Merged Frequency", ascending=False).reset_index(drop=True)
    new_df = pd.concat([df[first_column], first_freq_col, merged_df["Merged Kmer"], merged_df["Merged Frequency"], df[last_column], last_freq_col], axis=1)

    new_df.replace(0, np.nan, inplace=True)
    new_df.dropna(axis=0, how='all', inplace=True)
    new_df.reset_index(drop=True, inplace=True)
    new_df.to_excel(output_path, index=False)

    print("Done for kmer_mergence!")


def create_sankey_chart(file_path):
    """Create a Sankey chart based on data from an Excel file and display it.

    Parameters
    ----------
    file_path : string (/path/to/*.xlsx)
        The file path to the Excel file containing the data.

    Returns
    -------
    None

    Examples
    --------
    >>> create_sankey_chart('/path/to/data.xlsx')
    """
    random.seed(20000706)
    xls = pd.ExcelFile(file_path)
    sheet_names = xls.sheet_names

    dfs = []
    for sheet_name in sheet_names:
        df = xls.parse(sheet_name)
        dfs.append(df)
    for i in range(len(dfs)):
        exec(f"df{i+1} = dfs[i]")

    nodes = {}
    # Extract nodes from each DataFrame and store them in the nodes dictionary
    for i in range(1, len(dfs)+1):
        exec(f'nodes[{i}] = df{i}.iloc[:, 0].tolist()')
        if i == len(dfs):
            exec(f'nodes[{i+1}] = df{i}.iloc[:, 2].tolist()')
        else:
            exec(f'nodes[{i}] = df{i}.iloc[:, 0].tolist()')

    keys = sorted(nodes.keys())
    # Create individual lists of nodes using dynamic variable names
    for i, key in enumerate(keys):
        exec(f"nodes{i+1} = nodes[{key}]")

    # Add extra nodes to the existing nodes using dynamic variable names
    for i in range(len(dfs), 2, -1):
        exec(f'nodes{i}_extra = list(set(nodes{i+1}).difference(set(nodes{i-1})))')
        exec(f'nodes{i} += nodes{i}_extra')

    node_colors = {}
    all_nodes = []
    for i in range(1, len(dfs) + 2):
        all_nodes.extend(eval(f"nodes{i}"))


    html_colors = ['#B4DAFE', '#8EB1D8', '#6690B8', '#436490', '#1F4A62', '#E76254', '#EE8A47',
                   '#F7AA58', '#FFD06F', '#FFE6B7', '#AADCE0', '#72BCD5', '#528FAD', '#376795',
                   '#1E466E', '#FFD47F', '#F7C1CF', '#7B92C7', '#ADD9EE', '#F58383', '#FBB46F',
                   '#FAEE85', '#A8CAE8', '#FCBB44', '#F1766D', '#7A70B5', '#839DD1', '#FAECA8',
                   '#F18C25', '#FDD5C0', '#DB614F', '#F6C0CC', '#E64825', '#CAC0E1', '#715EA9',
                   '#ABDAEC', '#6A9ACE', '#97D1A0', '#1E803D', '#FCDFBE', '#F3DAC0', '#F7C4C1',
                   '#E3C6E0', '#CECCE5', '#C3E2EC', '#BCD1BC', '#DBEDC5', '#FDD379', '#F7C2CD',
                   '#A6DAEF', '#B0D9A5', '#E68D3D', '#E26472', '#6270B7', '#077535', '#F1AEA7',
                   '#ED6F6E', '#9D9ECD', '#5867AF', '#FCBB44', '#F1766D', '#7A70B5', '#839DD1',
                   '#FAC074', '#F28147', '#9FD4AE', '#68BD48', '#0D8B43', '#5560AC', '#F49E39',
                   '#E7483D', '#918AC2', '#8FC751', '#DF8D44', '#A6519E', '#6565AE', '#357337',
                   '#ECA8A9', '#74AED4', '#D3E2B7', '#CFAFD4', '#F7C97E', '#67ADB7', '#F5E1D8',
                   '#E4A6BD', '#F3D8E1', '#AFACB7', '#F3A17C', '#FAC45A', '#78A040', '#36600E',
                   '#E6E0B0', '#C2C1E0', '#D0E4EF', '#F7CCC6', '#D4B6D8', '#AA3E53', '#AEC7E8',
                   '#D89C7C', '#CF9198', '#E7CFD4', '#FF9896', '#1A3B30', '#509579', '#99BDBD',
                   '#C6E3D8', '#5A93BD', '#957064', '#C1AEA1', '#C4C7B4', '#EDC6BC', '#FFDA66',
                   '#48597E', '#8D9FBB', '#C1C2D2', '#C6D7E3', '#E0E7F5', '#544559', '#B5AABD',
                   '#AABDC4', '#D1C7AE', '#E5D9DA', '#6E746A', '#8D958F', '#C8CCC1', '#E7DBC1',
                   '#E3E5D7', '#EB6468', '#F2A7A2', '#F4CFD6', '#FCE694', '#C3E4F5', '#6A88C2',
                   '#EDFB9A', '#E6D7C4', '#FF999A', '#C6C5E5', '#CEE6C4', '#E6C4E5', '#F0DB63',
                   '#4FBCA7', '#BB4130', '#BB8532', '#B97865', '#333333', '#D62728']
    random.shuffle(html_colors)

    color_index = 0  

    for node in all_nodes:
        if color_index >= len(html_colors):
            color_index = 0
            random.shuffle(html_colors)
        node_colors[node] = html_colors[color_index]
        color_index += 1

    links = []

    # Create links between nodes based on the data in the DataFrames
    if len(nodes) == 4:
        for i, node in enumerate(dfs[0].iloc[:, 0]):
            link = dict(source=nodes[1].index(node),
                        target=len(nodes[1]) + nodes[2].index(dfs[0].iloc[i, 2]),
                        value=dfs[0].iloc[i, 1],
                        color=node_colors[node])
            links.append(link)

        for i, node in enumerate(dfs[1].iloc[:, 0]):
            link = dict(source=len(nodes[1]) + nodes[2].index(node),
                        target=len(nodes[1]) + len(nodes[2]) + nodes[3].index(dfs[1].iloc[i, 2]),
                        value=dfs[1].iloc[i, 1],
                        color=node_colors[node])
            links.append(link)

    else:
    #For a general case where there are more than three DataFrames
        # for i, node in enumerate(dfs[0].iloc[:, 0]):
        #     link = dict(source=nodes[1].index(node),
        #                 target=len(nodes[1]) + nodes[2].index(dfs[0].iloc[i, 2]),
        #                 value=dfs[0].iloc[i, 1],
        #                 color=node_colors[node])
        #     links.append(link)

        # for i, node in enumerate(dfs[1].iloc[:, 0]):
        #     link = dict(source=len(nodes[1]) + nodes[2].index(node),
        #                 target=len(nodes[1]) + len(nodes[2]) + nodes[3].index(dfs[1].iloc[i, 2]),
        #                 value=dfs[1].iloc[i, 1],
        #                 color=node_colors[node])
        #     links.append(link)

        # for i, node in enumerate(dfs[2].iloc[:, 0]):
        #     link = dict(source=len(nodes[1]) + len(nodes[2]) + nodes[3].index(node),
        #                 target=len(nodes[1]) + len(nodes[2]) + len(nodes[3]) + nodes[4].index(dfs[2].iloc[i, 2]),
        #                 value=dfs[2].iloc[i, 3],
        #                 color=node_colors[node])
        #     links.append(link)

        # for i, node in enumerate(dfs[3].iloc[:, 0]):
        #     link = dict(source=len(nodes[1]) + len(nodes[2]) + len(nodes[3]) + nodes[4].index(node),
        #                 target=len(nodes[1]) + len(nodes[2]) + len(nodes[3]) + len(nodes[4]) + nodes[5].index(dfs[3].iloc[i, 2]),
        #                 value=dfs[3].iloc[i, 3],
        #                 color=node_colors[node])
        #     links.append(link)



        for i, df in enumerate(dfs):
            for j, node in enumerate(df.iloc[:, 0]):
                source_index = sum(len(nodes[k]) for k in range(1, i+1)) + nodes[i+1].index(node)
                target_index = sum(len(nodes[k]) for k in range(1, i+2)) + nodes[i+2].index(df.iloc[j, 2])
                link = dict(source=source_index,
                            target=target_index,
                            value=df.iloc[j, 1],
                            color=node_colors[node])
                links.append(link)



    node_labels = [node for node in all_nodes]  
    node_colors = [node_colors[node] for node in all_nodes]  
    fig = go.Figure(data=[go.Sankey(
        node=dict(label=node_labels, color=node_colors),
        link=dict(source=[link['source'] for link in links],
                  target=[link['target'] for link in links],
                  value=[link['value'] for link in links],
                  color=[link['color'] for link in links]))
    ])
    fig.update_layout(title_text="Sankey Diagrams(K=7, evolution=2)", width=1600, height=1200)

    fig.show()
    
    
##################################################################################################################
##################################################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = 'APV-Sankey.py',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = ("""

                       
    ___    ____ _    __     _____             __             
   /   |  / __ \ |  / /    / ___/____ _____  / /_____  __  __
  / /| | / /_/ / | / /_____\__ \/ __ `/ __ \/ //_/ _ \/ / / /
 / ___ |/ ____/| |/ /_____/__/ / /_/ / / / / ,< /  __/ /_/ / 
/_/  |_/_/     |___/     /____/\__,_/_/ /_/_/|_|\___/\__, /  
                                                    /____/   

                                          
The program (NAME:APV-Sankey) is a bioinformatics aptamers analysis toolbox for analyzing sequence data. \n
It includes functions for calculating base ratios, kmer statistics, identifying common kmers and sequences, \n
concatenating k-mers, clustering, performing kmer evolution analysis, creating sankey diagrams, and more. \n
It also provides the functionality to create a Sankey chart based on the input data. \n
The program accepts command-line arguments to specify the input files, parameters, and output file paths. \n
Users can input sequence files in various formats and customize analysis parameters to gain insights into their data.
                      """))

    parser.add_argument('--input', help='Path to the sequence file (.xlsx or .fasta)')
    parser.add_argument('--rna', action='store_true', help='Flag indicating whether the sequences are RNAs')
    parser.add_argument('--k', type=int, default=8, help='Kmer length for Kmer statistics (default: 8)')
    parser.add_argument('--n_Kmers', type=int, default=200, help='Number of top Kmers to retrieve (default: 200)')
    parser.add_argument('--output', type=str, help='Output filepath for Kmer statistics')
    parser.add_argument('--common_kmers_input', nargs='+', help='Common_kmers_input file paths')
    parser.add_argument('--common_kmers_output', help='Output filepath for common kmers')
    parser.add_argument('--common_seqs_input', nargs='+', help='Common_seqs_input file paths')
    parser.add_argument('--common_seqs_output', help='Output filepath for common seqs')
    parser.add_argument('--kmer_concatenation_input', help='Kmer_concatenation_input file paths')
    parser.add_argument('--kmer_concatenation_output', help='Output filepath for kmer_concatenation')
    parser.add_argument('--cluster_filepath', help='Path for clustering')
    parser.add_argument("--kmer_evolution_input", nargs="+", help="Kmer_evolution_input file paths")
    parser.add_argument("--evolution-count", type=int, default=2, help="Number of allowed evolutions in a kmer(default: 1)")
    parser.add_argument("--kmer_evolution_output", help="Output filepath for kmer_evolution")
    parser.add_argument('--kmer_mergence_input', type=str, help='Kmer_mergence_input file paths')
    parser.add_argument('--kmer_mergence_output', type=str, help='Output filepath for kmer_mergence')
    parser.add_argument('--sankey_input', type=str, help='Path for sankey diagrams')
    args = parser.parse_args()


    calculate_nucleotide_ratio(args.input, is_rna=args.rna)
    print("Calculate_nucleotide_ratio completed.")
    kmer_statistics(args.input, args.k, args.n_Kmers, args.output)
    print("Kmer_statistics completed.")

    calculate_common_kmers(args.common_kmers_input, args.common_kmers_output)
    calculate_common_seqs(args.common_seqs_input, args.common_seqs_output)

    kmer_concatenation(args.kmer_concatenation_input, args.kmer_concatenation_output, sheet_name='Sheet')
    cluster_tree(args.cluster_filepath)
    
    kmer_evolution(args.kmer_evolution_input, args.evolution_count, args.kmer_evolution_output)
    kmer_mergence(args.kmer_mergence_input, args.kmer_mergence_output)
    create_sankey_chart(args.sankey_input)