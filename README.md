# Summary
The aim here is to create a series of plots a la [King et al 2013](https://www.nature.com/articles/nature12504) 
and [Mabb et al 2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156439) 
showing how gene expression changes as a function of gene length. 

Here we will use two datasets:

1. First we will use data from [Mabb et al 2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156439) 
to demonstrate how this new code reproduces previously published Figure 2B (the red line comparing WT Topotecan vs WT Vehicle). 
    * These results can be viewed [here](https://github.com/jeremymsimon/Zylka_geneLength_binnedPlots/blob/main/GSE79951_lineplots.md)

2. Then we will use data from [Jimenez et al 2020](https://molecularautism.biomedcentral.com/articles/10.1186/s13229-020-00369-8) 
to explore whether gene expression changes coincide with the age of the animal and whether those effects are exacerbated in Chd8 heterozygotes.
    * These results can be viewed [here](https://github.com/jeremymsimon/Zylka_geneLength_binnedPlots/blob/main/GSE142208_lineplots.md)


For each dataset, we retrieve the data directly from GEO as published, gene lengths as annotated by [GENCODE](https://www.gencodegenes.org/mouse),
and bin genes into groups of 200 based on their length. We aggregate biological replicate expression values by taking their mean per gene and per condition.
Then within each bin, we compute and plot the mean log2-fold-change in expression relative to the mean gene length, plotted on a log10 scale. 

Note that in the future, in any context where we've already done a differential expression test with DESeq2, we should plug in the shrunken log2-fold-changes from that test rather than computing them by hand as we did here.
