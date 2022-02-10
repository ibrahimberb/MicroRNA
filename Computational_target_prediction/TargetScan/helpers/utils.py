import pandas as pd
from helpers.etc import get_site_positions
from helpers.UTR import get_UTR_seq_from_ENSMUSG


def find_motifs(motifs, sequence):
    motif_to_counts = {}
    for motif in motifs:
        motif_to_counts[motif] = sequence.count(motif)

    return motif_to_counts


def get_table_entry(input_gene, pred_target_info_data, UTR_seq_data, gene_mapping, representative_transcripts, motifs):
    query = pred_target_info_data[
        pred_target_info_data["Gene ID"].apply(lambda x: x.split(".")[0]) == input_gene
        ]

    site_positions = get_site_positions(query)

    entry = [
        (
            gene_mapping[input_gene],
            input_gene,
            len(query),
            site_positions,
        )
    ]

    entry_data = pd.DataFrame(
        entry, columns=["GENE_NAME", "GENE_ENSMUSG", "#_SITE", "SITE_POSITIONS"]
    )

    sequence = get_UTR_seq_from_ENSMUSG(
        UTR_seq_data=UTR_seq_data,
        ENSMUSG_ID=input_gene,
        representative_transcripts=representative_transcripts,
        return_nondash=True
    )

    motif_counts = find_motifs(motifs, sequence)

    for motif, count in motif_counts.items():
        entry_data[motif] = count

    return entry_data


def construct_table(input_genes, pred_target_info_data, UTR_seq_data, gene_mapping, representative_transcripts, motifs):
    entries = []  # dataframes
    for input_gene in input_genes:
        entry = get_table_entry(
            input_gene=input_gene,
            pred_target_info_data=pred_target_info_data,
            UTR_seq_data=UTR_seq_data,
            gene_mapping=gene_mapping,
            representative_transcripts=representative_transcripts,
            motifs=motifs
        )
        entries.append(entry)

    constructed_table = pd.concat(entries, ignore_index=True)

    return constructed_table


def read_genes(file_name):
    with open(file_name) as file:
        genes = [line.strip() for line in file.readlines()]

    return genes
