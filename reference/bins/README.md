# QDNAseq bins

QDNAseq bin annotations for T2T-CHM13 v2.0 at several bin sizes (`T2T.<size>.SR50.rda`, 50 bp reads).

The pipeline only uses **`T2T.500kbp.SR50.rda`**, set by `references.qdnaseq_bins` in config/samples.yaml. The other sizes are kept for future development, for example CNV plots at a different resolution.

Changing `qdnaseq_bins` alone is not enough to switch bin size:
- Output file names contain `_500kbp`.
- The CNV gene annotation file (`references.cnv_gene_bed`) holds bin numbers for 500 kb bins and would need to be regenerated for the new size.
