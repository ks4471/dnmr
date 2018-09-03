
# Enrichment test for gene set over-representation for de-novo mutations associated with neurodevelopmental disorders based on denovo-db database
`
input_gene_sets=list(gene_set1=c("ABCA4","ABCG4","ACOT7","ACSM4"),gene_set_2=c("AASDHPPT","ABCE1","ABHD13","ABRAXAS2"))

dnm_enrich=dnmr(input_gene_sets)
`
wrapper for heatmap.2 function, values='pval' converts P-values of enrichment to -log10(Pvalue) & adds R style significance values e.g. *** <0.0001
`
Heat(as.matrix(dnm_enrich$pval),values='pval')
`
# data from denovo-db:
Turner, T. N. et al. denovo-db: a compendium of human de novo variants. Nucleic Acids Res. 45, D804–D811 (2017).

# test used in:
1.  Delahaye-Duriez, A. et al. Rare and common epilepsies converge on a shared gene regulatory network providing opportunities for novel antiepileptic drug discovery. Genome Biol. 17, 245 (2016).
2.  Johnson, M. R. et al. Systems genetics identifies a convergent gene network for cognition and neurodevelopmental disease. Nat. Neurosci. 19, 223–232 (2015).



