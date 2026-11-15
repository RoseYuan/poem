# Details about each evaluation metrics

## Partition-based metrics

| Min_level | Metric | Calculation |
|:---|:---|:---|
| dataset | Rand Index (RI) | $`\frac{a+d}{n(n-1)/2}`$; the ratio of the sum of true positive and true negative pairs to the total number of object pairs. |
| class/cluster | Wallace Homogeneity (WH) | $`\frac{a}{a+c}`$; the ratio of the true positive pairs to the total number of object pairs that are in the same cluster in $`P`$. |
| class/cluster | Wallace Completeness (WC) | $`\frac{a}{a+b}`$; the ratio of the true positive pairs to the total number of object pairs that are in the same classes in $`G`$. |
| dataset | Adjusted Rand Index (ARI) | $`\frac{\text{RI}-\mathrm{E}(\text{RI})}{1-\mathrm{E}(\text{RI})} = \frac{2(ad-bc)}{(a+b)(b+d)+(a+c)(c+d)}`$; adjusting RI by accounting for the expected similarity of all pairings due to chance using the Permutation Model for clusterings. ARI is the harmonic mean of AWH and AWC. |
| dataset | Normalized Class Size Rand Index (NCR) | A normalized version of RI, where each concordance quantities are divided by the maximum possible concordance values for their respective class. |
| dataset | Mutual Information (MI) | $`H(G) - H(G|P)`$; the difference between the shannon entropy of $`G`$ and the conditional entropy of $`G`$ given $`P`$. |
| class/cluster | Adjusted Wallace Homogeneity (AWH), Adjusted Wallace Completeness (AWC), and Adjusted Mutual Information (AMI) | Chance adjusted version of WH, WC and MI, respectively. For a metric M, the chance adjusted version of it is $`\frac{\text{M}-\mathrm{E}(\text{M})}{1-\mathrm{E}(\text{M})}`$. |
| dataset | (Entropy-based) Homogeneity (EH) | $`1-\frac{H(G|P)}{H(G)}`$ if $`H(G,P)\neq0`$, $`1`$ otherwise; the ratio of MI to the individual entropy of $`G`$. |
| dataset | (Entropy-based) Completeness (EC) | $`1-\frac{H(P|G)}{H(P)}`$ if $`H(P,G)\neq0`$, $`1`$ otherwise; the ratio of MI to the individual entropy of $`P`$. |
| class/cluster | V Measure (VM) | $`\frac{(1+\beta)\times\text{EH}\times\text{EC}}{\beta\times\text{EH}+\text{EC}}`$; the harmonic mean between EH and EC. It is identical to normalized mutual information (NMI) when arithmetic mean is used for averaging in NMI calculation. |
| class/cluster | (weighted average) F Measure (wFM) | Here we calculate weighted F1-score, where the weights are based on the sizes of classes. |

Partition-based metrics. The notation used is common throughout the
table: consider comparing the predicted partition $`P`$ to the
ground-truth partition $`G`$; $`a`$ is the number of pairs that are in
the same group both in $`P`$ and $`G`$; $`b`$ is the number of pairs
that are in the same class in $`G`$ but in different clusters in $`P`$;
$`c`$ is the number of pairs that are in different classes in $`G`$ but
in the same cluster in $`P`$; $`d`$ is the number of pairs that are in
different groups both in $`P`$ and $`G`$; $`n`$ is the total number of
objects; $`E`$ is the expectation operator; $`H(⋅)`$ is the Shannon
entropy; \$eta\$ is the ratio of weight attributed to homogeneity vs
completeness; the expactation value of RI, WH, and WC is calculated when
assuming a generalized hypergeometric model.

## Embedding-based metrics

| Min_level | Metric | Calculation |
|:---|:---|:---|
| dataset | Silhouette score | $`\frac{n-m}{\text{max}(m, n)}`$, where $`n`$ is the mean distance between a sample and the nearest class that the sample is not a part of, and $`m`$ is the mean intra-class distance. |
| dataset | Composed Density between and within Clusters (CDbw) | The CDbw index consists of three main components: cohesion, compactness, and separation between clusters. It uses multiple representative points selected from each cluster to calculate intra-cluster density and between-cluster distances, reflecting the geometry of the clusters and capturing changes in intra-cluster density. |
| dataset | Density Based Clustering Validation index (DBCV) | A density-based index that computes the least dense region inside a cluster and the most dense region between the clusters, to measure the within and between cluster density connectedness of clusters. |

Embedding-based metrics.

## Graph-based metrics

| Min_level | Metric | Calculation |
|:---|:---|:---|
| dataset | Modularity | For a given graph partition, it quantifies the number of edges within communities relative to what would be expected by random chance. $`Q = \frac{1}{2m} \sum_{ij} \left( A_{ij} - \gamma \frac{k_i k_j}{2m} \right) \delta(c_i, c_j)`$, where $`m`$ is the number of edges, $`A`$ is the adjacency matrix of the graph, $`k_i`$ is the (weighted) degree of $`i`$, $`\gamma`$ is the resolution parameter, and $`\delta(c_i, c_j)`$ is $`1`$ if $`i`$ and $`j`$ are in the same community else $`0`$. |
| element | Local Inverse Simpson’s Index (LISI) | For a given node in a weighted kNN graph, the expected number of nodes needed to be sampled before two nodes are drawn from the same classes within its neighborhood. |
| element | Neighborhood Purity (NP) | For each node in a graph, the proportion of its neighborhood that is of the same class as it. |
| element | Proportion of Weakly Connected (PWC) | For a given community in a graph, the proportion of nodes that have more connections to the outside of the community than the inside of the community. |
| element | Cohesion | The minimum number of nodes that must be removed to split a graph. |
| class/cluster | Adhesion | The minimum number of edges that must be removed to split a graph. |
| class/cluster | Adjusted Mean Shortest Path (AMSP) | A measure of the disconnectness and spread of the subgraph connecting elements of a given class. If the graph subclass is disconnected, the mean shortest path of each connected subgraph $`m`$ are summed. $`\frac{\sum_{i} (1+m_i)}{\sqrt{N}}`$, where $`m`$ is the mean shortest path and $`N`$ is the number of nodes of the given class. Note that the normalization for size is only approximative, and only applicable for kNN graphs. |
| class/cluster | Neighborhood Class Enrichment (NCE) | The log2 fold-enrichment (i.e. over-representation) of the node’s class among its nearest neighbors, over the expected given its relative abundance. |

Graph-based metrics.

## Metrics for spatial clusterings

| Min_level | Metric | Calculation |
|:---|:---|:---|
| class/cluster | Percentage of Abnormal Spots (PAS) | PAS measures the percentage of abnormal spots, which is defined as spots with a spatial domain label differing from more than half of its nearest neighbors. |
| class/cluster | Spatial Chaos Score (CHAOS) | CHAOS is the mean length of the graph edges in the 1-nearest neighbor (1NN) graph for each domain averaged across domains. |
| element | Entropy-based Local indicator of Spatial Association (ELSA) | For a site $`i`$, $`E_i = E_{ai} \times E_{ci}`$, where $`E_{ai}`$ summarizes the dissimilarity between site $`i`$ and the neighbouring sites, and $`E_{ci}`$ quantifies the diversity of the categories within the neighbourhood of site $`i`$. |
| dataset | Spatial RI, ARI, WH, WC, AWH, and AWC | Spatial versions of the pair-sorting indices, based on fuzzy versions of the metrics. Specifically, we use the Normalized Degree of Concordance (NDC, see Hullermeier et al., 2012) and the Adjusted Concordance Index (ACI, see D’Errico et al., 2021) as fuzzy versions of RI and ARI respectively, and developed fuzzy versions of the other metrics using the same logic. In the spatial context, we first make a fuzzy version of the true labels based on the spatial neighborhood, and then track the maximum pair concordance between the predicted labels and either the hard or fuzzy ground truth. |
| element | Spot-wise Pair Concordance (SPC) | The proportion, for each spot, of the pairs it forms with all other spots that are concordant (i.e. in the same partition or not in both) across the clustering and ground truth. This value will be the same for all spots that share the same combination of cluster and class, and is especially useful for visualization. A variant of this can be computed that ignores negative pairs (i.e. that are discordant in both the clustering and ground truth). When negative pairs are included, the average of SPC equals to the Rand Index. |
| element | Spatial SPC | Like the non-spatial Spot-wise Pair Concordance, with the difference that the clustering is evaluated against both a ‘hard’ and ‘fuzzy’ version of the ground truth, as for the computation of the Spatial versions of the pair-sorting indices. |
| dataset | Spatial Set Matching Accuracy | An accuracy that downweights misclassifications based on the spatial neighborhood. Instead of counting as zero in the accuracy computation, the misclassified node counts as the proportion of its spatial neighborhood that is of node’s predicted class. |

Metrics for spatial clusterings.

## Session info

``` r
sessionInfo()
```

    ## R version 4.4.2 (2024-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Zurich
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] BiocStyle_2.32.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.6.5         svglite_2.1.3       cli_3.6.3          
    ##  [4] knitr_1.48          rlang_1.1.4         xfun_0.46          
    ##  [7] stringi_1.8.4       textshaping_0.3.6   jsonlite_1.8.8     
    ## [10] glue_1.8.0          colorspace_2.1-1    htmltools_0.5.8.1  
    ## [13] ragg_1.3.2          sass_0.4.9          scales_1.3.0       
    ## [16] rmarkdown_2.27      munsell_0.5.1       evaluate_0.24.0    
    ## [19] jquerylib_0.1.4     kableExtra_1.4.0    fastmap_1.2.0      
    ## [22] yaml_2.3.10         lifecycle_1.0.4     bookdown_0.40      
    ## [25] stringr_1.5.1       BiocManager_1.30.23 compiler_4.4.2     
    ## [28] fs_1.6.4            htmlwidgets_1.6.4   rstudioapi_0.16.0  
    ## [31] systemfonts_1.1.0   digest_0.6.36       viridisLite_0.4.2  
    ## [34] R6_2.5.1            magrittr_2.0.3      bslib_0.8.0        
    ## [37] tools_4.4.2         xml2_1.3.6          pkgdown_2.1.1      
    ## [40] cachem_1.1.0        desc_1.4.3
