## code to prepare `metric_info` dataset goes here
metric_info <- data.frame(Name=c("SW",
                                 
                                 "meanSW", "minSW", "pnSW", "dbcv",
                                 
                                 "meanSW", "meanClassSW", "pnSW", "minClassSW", 
                                 "cdbw", "cohesion (cdbw)", "compactness", 
                                 "sep", "dbcv",
                                 
                                 "SI","ISI","NP","NCE",
                                 
                                 "SI","ISI","NP","NCE","AMSP","PWC","adhesion",
                                 "cohesion (graph)",
                                 
                                 "SI","ISI","NP","AMSP","PWC","NCE", 
                                 "adhesion","cohesion (graph)",
                                 
                                 "SPC", "ASPC",
                                 
                                 "WC","WH","AWC","AWH","FM",
                                 
                                 "RI","WC","WH","ARI","AWC","AWH","NCR","FM",
                                 "MI","AMI","VI","EH","EC","VM","VDM","MHM",
                                 "MMM","Mirkin","Accuracy",
                                 
                                 "fuzzySPC",
                                 
                                 "fuzzyWH", "fuzzyAWH", "fuzzyWC", "fuzzyAWC",
                                 
                                 "fuzzyRI", "fuzzyARI", "fuzzyWH", "fuzzyAWH", 
                                 "fuzzyWC", "fuzzyAWC",
                                 
                                 "nsSPC", "NPC",
                                 
                                 "nsWH","nsAWH", 
                                 "nsWC","nsAWC",
                                 
                                 "nsRI","nsARI","nsWH",
                                 "nsAWH", "nsWC","nsAWC",
                                 "nsAccuracy", "SpatialARI", "SpatialRI", 
                                 "SpatialSPC", 
                                 
                                 "PAS", "ELSA",
                                 
                                 "CHAOS", "PAS", "ELSA",
                                 
                                 "CHAOS", "PAS", "ELSA", "MPC", "PC", "PE"
                                 ), 
                          Levels=c("element",
                                   rep("class",4),
                                   rep("dataset",9),
                    
                                   rep("element",4),
                                   rep("class",8),
                                   rep("dataset",8),
                                   
                                   rep("element",2),
                                   rep("class",5),
                                   rep("dataset",19),
                                   
                                   rep("element",1),
                                   rep("class",4),
                                   rep("dataset",6),
                                   
                                   rep("element",2),
                                   rep("class",4),
                                   rep("dataset",9),
                                   
                                   rep("element",3),
                                   rep("class",3),
                                   rep("dataset",6)
                                        ), 
                          MainFunction=c(
                            rep("getEmbeddingMetrics",14),
                            rep("getGraphMetrics",20),
                            rep("getPartitionMetrics",26),
                            rep("getFuzzyPartitionMetrics",11),
                            rep("getSpatialExternalMetrics",16),
                            rep("getSpatialInternalMetrics",11)
                            ))

metric_info <- metric_info %>% group_by(Name) %>% summarize(
  Levels = paste(unique(Levels), collapse = ", "),
  MainFunction = paste(unique(MainFunction), collapse = ", ")
)




tmp <- data.frame(metric=c(
           "Adjusted Mutual Information",
           "Adjusted Mean Shortest Path",
           "Adjusted Rand Index",
           "Adjusted Spot-wise Pair Concordance",
           "Adjusted Wallace Completeness",
           "Adjusted Wallace Homogeneity",
           "Set-matching Accuracy",
           "Spatial Chaos Score",
           "(Entropy-based) Completeness",
           "(Entropy-based) Homogeneity",
           "Entropy-based Local Indicator of Spatial Association (ELSA score)",
           "F-measure/weighted average F1 score",
           "Inverse Simpson’s Index",
           "Meila-Heckerman Measure",
           "Mutual Information",
           "Maximum-Match Measure",
           "Modified Partition Coefficient",
           "Mirkin Metric",
           "Neighborhood Class Enrichment",
           "Normalized Class Size Rand Index",
           "Neighborhood Purity",
           "Neighboring Pair Concordance for Spatial Clustering",
           "Proportion of Abnormal Spots (PAS score)",
           "Partition Coefficient",
           "Partition Entropy",
           "Proportion of Weakly Connected",
           "Rand Index",
           "Simpson’s Index",
           "Spot-wise Pair Concordance",
           "Silhouette Width",
           "Adjusted Rand Index weighted by distance for Spatial Clustering",
           "Rand Index weighted by distance for Spatial Clustering",
           "Spot-wise Pair Concordance weighted by distance for Spatial Clustering",
           "Van Dongen Measure",
           "Variation of Information",
           "V-measure",
           "Wallace Completeness",
           "Wallace Homogeneity",
           "Adhesion of a graph",
           "Composed Density between and within Clusters (CDbw)",
           "CDbw Cohesion",
           "Cohesion of a graph",
           "CDbw Compactness",
           "Density Based Clustering Validation Index (DBCV)",
           "Fuzzy version of Adjusted Rand Index", 
           "Fuzzy version of Adjusted Wallace Completeness",
           "Fuzzy version of Adjusted Wallace Homogeneity", 
           "Fuzzy version of Rand Index", 
           "Fuzzy version of Spot-wise Pair Concordance",
           "Fuzzy version of Wallace Completeness", 
           "Fuzzy version of Wallace Homogeneity", 
           "Mean of Class-level Silhouette Width",
           "Mean Silhouette Width",
           "Minimal of Class-level Silhouette Width",
           "Minimal Silhouette Width", 
           "Neighborhood smoothed Adjusted Rand Index for Spatial Clustering",
           "Neighborhood smoothed Adjusted Wallace Completeness for Spatial Clustering",
           "Neighborhood smoothed Adjusted Wallace Homogeneity for Spatial Clustering", 
           "Neighborhood smoothed Set-matching Accuracy for Spatial Clustering",
           "Neighborhood smoothed Rand Index for Spatial Clustering",
           "Neighborhood smoothed Spot-wise Pair Concordance for Spatial Clustering",
           "Neighborhood smoothed Wallace Completeness for Spatial Clustering",
           "Neighborhood smoothed Wallace Homogeneity for Spatial Clustering",
           "Proportion of Negative Silhouette Width", 
           "CDbw seperation"
                        
))

metric_info$Description <- tmp$metric                             
metric_info <- data.frame(metric_info)                            
usethis::use_data(metric_info, overwrite = TRUE)
