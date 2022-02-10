import pandas as pd
from pandas import DataFrame
from tqdm.notebook import tqdm

from helpers.UTR import get_UTR_seq_nondash


class GeneUTRSequenceNotFoundError(Exception):
    pass


def _get_longest(sequences):
    sequences_nondash = [get_UTR_seq_nondash(seq) for seq in sequences]

    longest_seq = sequences_nondash[0]
    for seq in sequences_nondash:
        if len(seq) > len(longest_seq):
            longest_seq = seq

    return longest_seq


def _get_UTR_seq_and_gene_symbol_from_ENSMUSG(
        UTR_seq_data: DataFrame,
        ENSMUSG_ID: str,
        return_nondash=True,
        ignore_conversation=True,
        ignore_gene_version=True,
) -> [str, str]:
    """
    Returns the UTR sequence and gene symbol.

    Parameters
    ----------
        UTR_seq_data : <DataFrame>
            Given UTR sequence data.

        ENSMUSG_ID : <str>
            Given Gene id in ENSMUSG format.

        return_nondash : <bool>
            Return the nondash form of sequence if this parameter is True.

        ignore_conversation : <bool>
            Ignores the capitalization of the sequence letters (i.e. upper/lower case), by setting all characters
            into upper case.

        ignore_gene_version : <bool>
            The ENSMUSG have a number indicating the version of the gene. If this parameter set to True, we ignore this
            version.
            E.g.
                ENSMUSG00000035401.8 → ENSMUSG00000035401

    Returns
    -------
        UTR_seq : <str>
            The sequence associated with given gene.
    """
    if ignore_gene_version:
        query = UTR_seq_data[UTR_seq_data["Gene ID"].apply(lambda x: x.split(".")[0]) == ENSMUSG_ID]
    else:
        query = UTR_seq_data[UTR_seq_data["Gene ID"] == ENSMUSG_ID]

    try:
        [UTR_seq] = query["UTR sequence"]
    except ValueError:
        if len(query["UTR sequence"]) == 0:
            # print(f"Gene {ENSMUSG_ID} is not found.")
            raise GeneUTRSequenceNotFoundError

        UTR_seq = _get_longest(query["UTR sequence"])

    try:
        [gene_symbol] = query["Gene Symbol"].unique()
    except:
        print("number of unique values is not 1.")
        print(query["Gene Symbol"])
        raise

    if ignore_conversation:
        UTR_seq = UTR_seq.upper()

    if return_nondash:
        UTR_seq = get_UTR_seq_nondash(UTR_seq)

    return UTR_seq, gene_symbol


def find_motifs(motifs, sequence):
    motif_to_counts = {}
    for motif in motifs:
        motif_to_counts[motif] = sequence.count(motif)

    return motif_to_counts


def get_table_entry(
        input_gene,
        UTR_seq_data,
        motifs,
):

    sequence, gene_symbol = _get_UTR_seq_and_gene_symbol_from_ENSMUSG(
        UTR_seq_data=UTR_seq_data,
        ENSMUSG_ID=input_gene,
        return_nondash=True
    )

    entry = [
        (
            gene_symbol,
            input_gene,
            sequence,
            len(sequence)
        )
    ]

    entry_data = pd.DataFrame(
        entry,
        columns=[
            "GENE_SYMBOL", "GENE_ENSMUSG", "SEQUENCE", "SEQUENCE_LENGTH"
        ]
    )

    motif_counts = find_motifs(motifs, sequence)

    for motif, count in motif_counts.items():
        entry_data[motif] = count
        entry_data[f"{motif}_normalized"] = count / len(sequence)

    return entry_data


def construct_table(
        input_genes,
        UTR_seq_data,
        motifs
):
    entries = []  # dataframes
    genes_UTR_seq_not_found = []
    for input_gene in tqdm(input_genes):
        try:
            entry = get_table_entry(
                input_gene=input_gene,
                UTR_seq_data=UTR_seq_data,
                motifs=motifs
            )
            entries.append(entry)

        except GeneUTRSequenceNotFoundError:
            genes_UTR_seq_not_found.append(input_gene)
            continue

    constructed_table = pd.concat(entries, ignore_index=True)

    print(f"Number of genes whose UTR sequence not found: {len(genes_UTR_seq_not_found)}")

    return constructed_table, genes_UTR_seq_not_found


def read_genes(file_name, unique=False):
    with open(file_name) as file:
        genes = [line.strip() for line in file.readlines()]

    if unique:
        return sorted(set(genes))

    else:
        return genes
