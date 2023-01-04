Exploration of gene expression changes in the cortex following topotecan
treatment
================

# Summary

The aim here is to create a series of plots a la [King et al
2013](https://www.nature.com/articles/nature12504) and [Mabb et al
2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156439)
showing how gene expression changes as a function of gene length.

Here we will use data from [Mabb et al
2016](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156439)
to demonstrate how this new code reproduces previously published Figure
2B (the red line comparing WT Topotecan vs WT Vehicle).

We retrieve the data directly from GEO as published, gene lengths as
annotated by [GENCODE
vM25](https://www.gencodegenes.org/mouse/release_M25.html), and bin
genes into groups of 200 based on their length. We aggregate biological
replicate expression values by taking their mean per gene and per
condition. Then within each bin, we compute and plot the mean
log2-fold-change (topotecan vs vehicle) in expression relative to the
mean gene length, plotted on a log10 scale.

# Load libraries

``` r
library(DESeq2)
library(tidyverse)
library(plyranges)
library(ggh4x)
library(patchwork)
```

# Import in matrix from earlier paper directly from GEO

Here since the data supplied are raw counts, we will briefly read them
into DESeq2 to VST-normalize them Note, this is not exactly the same
normalization as employed in the publication

``` r
GSE79951 <- read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79951/suppl/GSE79951_WT_Top1cKO_Veh_Topot_rawCounts_matrix.txt.gz") %>%
    dplyr::select(-contains("cKO")) %>%
    dplyr::rename("WT_Topotecan_Rep1" = `WT-Topot-Rep1`,
                  "WT_Topotecan_Rep2" = `WT-Topot-Rep2`,
                  "WT_Topotecan_Rep3" = `WT-Topot-Rep3`,
                  "WT_Vehicle_Rep1" = `WT-Veh-Rep1`,
                  "WT_Vehicle_Rep2" = `WT-Veh-Rep2`,
                  "WT_Vehicle_Rep3" = `WT-Veh-Rep3`
        )
```

    ## Rows: 23893 Columns: 13
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (1): Gene
    ## dbl (12): Top1-cKO-Topot-Rep1, Top1-cKO-Topot-Rep2, Top1-cKO-Topot-Rep3, Top1-cKO-Veh-Rep1, Top1-cKO-Veh-Rep2, Top1-cKO-Veh-Rep3, W...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
rep <- factor(str_sub(colnames(GSE79951)[-1],-4,-1))
treatment <- factor(str_replace_all(colnames(GSE79951)[-1],"(.+)_(.+)_Rep.","\\2"))

coldata <- data.frame(row.names=colnames(GSE79951)[-1], rep, treatment)

dds <- DESeqDataSetFromMatrix(countData = as.data.frame(column_to_rownames(GSE79951,var="Gene")), colData = coldata, design = ~treatment)
```

    ## converting counts to integer mode

``` r
vsd <- vst(dds)
```

# Now aggregate biological replicates

``` r
GSE79951.summarized <- rownames_to_column(as.data.frame(assay(vsd)),var="Gene") %>%
    as_tibble() %>% 
    pivot_longer(cols=!Gene,names_to="Sample",values_to="Expression") %>%
    separate(Sample,c("Genotype","Treatment","Rep"),sep="_",remove=F) %>%
    dplyr::select(-Genotype) %>%
    group_by(Treatment,Gene) %>%
    summarize(Mean = mean(Expression))  
```

    ## `summarise()` has grouped output by 'Treatment'. You can override using the `.groups` argument.

# Download gene annotation and compute gene lengths, taking max length per gene symbol

Original paper used mm9 annotations, so we’ll use GENCODE M1

``` r
options(timeout = max(300, getOption("timeout")))
gtf <- read_gff("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz")
lengths <- gtf %>%
    group_by(gene_name) %>% 
    summarize(Max = max(width)) %>% 
    as_tibble() %>%
    dplyr::rename("Gene" = gene_name, "Length" = Max)
```

# Restrict lengths and expression dataset to common set of gene symbols

``` r
intersect <- intersect(lengths$Gene,GSE79951.summarized$Gene)
GSE79951.summarized.inter <- GSE79951.summarized %>%
    dplyr::filter(Gene %in% intersect)

lengths.inter <- lengths %>%
    dplyr::filter(Gene %in% intersect)
```

# Bin genes into groups of 200 by length

Uses solution from [this post](https://stackoverflow.com/a/27626007)

``` r
chunk <- function(x, n) (mapply(function(a, b) (x[a:b]), seq.int(from=1, to=length(x), by=n), pmin(seq.int(from=1, to=length(x), by=n)+(n-1), length(x)), SIMPLIFY=FALSE))
z <- chunk(lengths.inter$Gene[order(lengths.inter$Length,decreasing=T)],200)
names(z) <- paste0("Chunk",seq(1:length(z)))

genes.chunked <- stack(z) %>% 
    as_tibble() %>%
    dplyr::rename("Gene" = values, "Chunk" = ind)
```

# Plot samples comparing Topotecan vs vehicle

``` r
inner_join(GSE79951.summarized.inter,lengths.inter,by="Gene") %>%
    inner_join(genes.chunked,by="Gene") %>%
    mutate(LogLength = log10(Length+1)) %>%
    pivot_wider(names_from=Treatment,values_from=Mean) %>%

    ## Compute ratios between Topotecan and vehicle
    mutate(Ratio = Topotecan - Vehicle) %>%

    # Compute average ratio for each gene chunk 
    dplyr::select(Gene,Length,Chunk,LogLength,Ratio) %>%
    group_by(Chunk) %>%
    summarize(MeanRatio = mean(Ratio,na.rm=T),MeanLength = mean(LogLength,na.rm=T),n = dplyr::n()) %>%
    dplyr::filter(MeanLength > 3) %>%
    ggplot(aes(x=MeanLength,y=MeanRatio)) +
    geom_point() +
    geom_smooth() +
    coord_cartesian(ylim=c(-2, 2)) +
    xlab("log10 average gene length in bin") +
    ylab("log2 ratio, topotecan/vehicle")
```

    ## `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

![](GSE79951_lineplots_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
