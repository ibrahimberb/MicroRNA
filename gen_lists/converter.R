require(biomaRt)

myfiles <- lapply(
  list.files(pattern = "*.txt"), read.delim
)

genes <- (readLines("target_genes_mouse_miRNA_223_2022-01-10.txt"))


ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

annot <- getBM(
  attributes = c(
    'mgi_symbol',
    'external_gene_name',
    'ensembl_gene_id',
    'entrezgene_id',
    'gene_biotype'),
  filters = 'ensembl_gene_id',
  values = genes,
  mart = ensembl)

annot <- merge(
  x = as.data.frame(genes),
  y = annot,
  by.y = 'ensembl_gene_id',
  all.x = T,
  by.x = 'genes')

annot