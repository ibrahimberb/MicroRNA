import pandas as pd

DATA_PATH = "Predicted_Target_Locations.default_predictions.mm10.bed"

df = pd.read_csv(DATA_PATH, sep='\t', comment='t', header=None)
header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
df.columns = header[:len(df.columns)]

print(df)
