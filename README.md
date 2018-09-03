
# Enrichment test for gene set over-representation for de-novo mutations associated with neurodevelopmental disorders based on denovo-db database

To run enrichment, install addR - contains all relevant functions, including dependencies for dnmr()
Data visualisation via Heat() - a wrapper for heatmap.2 function (included in addR), values='pval' converts P-values of enrichment to -log10(P-value) & adds R style significance values e.g. *** <0.0001

```
devtools::install_github("ks4471/addR")
library(adds)
```
Input is expected as a list of named gene sets, enrichment performed for each gene set. Example gene list available to download from dropbox https://www.dropbox.com/s/y2szo7ywloleh92/example.gene.list.Rdata?dl=0
```
input_gene_sets=list(gene_set1=c("ABCA4","ABCG4","ACOT7","ACSM4"),gene_set_2=c("AASDHPPT","ABCE1","ABHD13","ABRAXAS2"))
dnm_enrich=dnmr(input_gene_sets)
Heat(as.matrix(dnm_enrich$pval),values='pval')
```




# data from denovo-db:
Turner, T. N. et al. denovo-db: a compendium of human de novo variants. Nucleic Acids Res. 45, D804–D811 (2017).

# used in:
1.  Delahaye-Duriez, A. et al. Rare and common epilepsies converge on a shared gene regulatory network providing opportunities for novel antiepileptic drug discovery. Genome Biol. 17, 245 (2016).
2.  Johnson, M. R. et al. Systems genetics identifies a convergent gene network for cognition and neurodevelopmental disease. Nat. Neurosci. 19, 223–232 (2015).

