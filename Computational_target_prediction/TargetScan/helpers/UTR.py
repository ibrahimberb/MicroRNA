from typing import List

from pandas import DataFrame


def get_UTR_seq_from_ENSMUST(UTR_seq_data: DataFrame, ENSMUST_ID: str) -> str:
    """
    Returns the UTR Sequence.
    """
    query = UTR_seq_data[UTR_seq_data["Refseq ID"] == ENSMUST_ID]
    [UTR_seq] = query["UTR sequence"]
    return UTR_seq


def get_UTR_seq_from_ENSMUSG(
        UTR_seq_data: DataFrame,
        ENSMUSG_ID: str,
        representative_transcripts: List[str],
        return_nondash=True,
        ignore_conversation=True,
        ignore_gene_version=True,
) -> str:
    """
    Returns the UTR sequence.

    Parameters
    ----------
        UTR_seq_data : <DataFrame>
            Given UTR sequence data.

        ENSMUSG_ID : <str>
            Given Gene id in ENSMUSG format.

        representative_transcripts : <List[str]>
            Ensures that there will be one transcript per gene, selected for being the most prevalent.
            This will be used in the case of multiple transcripts occur for a given gene to select desired transcript.

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
        multiple_transcripts = [transcript.split(".")[0] for transcript in query["Refseq ID"].tolist()]
        intersection_transcripts = set(representative_transcripts).intersection(multiple_transcripts)
        [selected_transcript] = intersection_transcripts
        query = UTR_seq_data[UTR_seq_data["Refseq ID"].str.contains(selected_transcript)]
        [UTR_seq] = query["UTR sequence"]

    if ignore_conversation:
        UTR_seq = UTR_seq.upper()

    if return_nondash:
        return get_UTR_seq_nondash(UTR_seq)
    else:
        return UTR_seq


def get_UTR_seq_nondash(UTR_seq: str):
    if not isinstance(UTR_seq, str):
        raise ValueError("Expected type <str>")
    UTR_seq_nondash = ''.join(
        [char for char in UTR_seq if char != "-"]
    )
    return UTR_seq_nondash
