require(biomaRt)

# genes_short <- c('ENSMUSG00000031201', 'ENSMUSG00000017146',
#                  'ENSMUSG00000041147', 'ENSMUSG00000034218', 'ENSMUSG00000059552')

genes <- (readLines("target_genes_mouse_miRNA_223_2022-01-10.txt"))

ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

annot <- getBM(
  attributes = c(
    'hgnc_symbol',
    'external_gene_name',
    'ensembl_gene_id',
    'entrezgene_id',
    'gene_biotype'),
  filters = 'external_gene_name',
  values = genes,
  mart = ensembl)

annot <- merge(
  x = as.data.frame(genes),
  y =  annot,
  by.y = 'external_gene_name',
  all.x = T,
  by.x = 'genes')

annot

write.csv(x=annot, file="ENSMUS_converted.csv")
