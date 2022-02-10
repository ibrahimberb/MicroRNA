from datetime import datetime
from typing import List
import logging
import os.path as op


def extract_genes(genes: List[str], file_name, overwrite=False, file_extension="txt"):
    file_date = datetime.today().strftime('%Y-%m-%d')
    file_name = f'{file_name}_{file_date}.{file_extension}'

    # Ensure the file is not exists before creating to prevent overwriting.
    if op.isfile(file_name) and not overwrite:
        print(f"File {file_name} is already exist.\n"
              "To overwrite existing file, use `overwrite=True`.")
        raise FileExistsError

    else:
        with open(file_name, "w") as file:
            for gene in genes:
                file.write(f"{gene}\n")

        print(f"Genes are exported successfully into file {file_name}")
