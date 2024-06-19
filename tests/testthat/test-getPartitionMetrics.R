#################
# generate 3 test datasets in "Warrens, Matthijs J., et, al. Journal of Classification 39.3 (2022): 487-509".
#################

library(dplyr)
total_instances <- 56
# dataset 1 (in Table 6)
reference_partition <- c(rep("U1", 20), rep("U2", 20), rep("U3", 8), rep("U4", 8))
trial_partition <- c(rep("Z1", 20), rep("Z2", 20), rep("Z3", 4), rep("Z4", 4), rep("Z3", 4), rep("Z4", 4))
dataset1 <- data.frame(
  Instance = 1:total_instances,
  Reference_Partition = reference_partition,
  Trial_Partition = trial_partition
)

# dataset 2 (in Table 8)
reference_partition <- c(
  rep("U1", 20),  # U1 has 20 instances
  rep("U2", 20),  # U2 has 20 instances
  rep("U3", 8),   # U3 has 8 instances
  rep("U4", 8)    # U4 has 8 instances
)
trial_partition <- c(
  rep("Z1", 10), rep("Z2", 10),  # U1 split between Z1 and Z2
  rep("Z1", 10), rep("Z2", 10),  # U2 split between Z1 and Z2
  rep("Z3", 8),                  # U3 goes to Z3
  rep("Z4", 8)                   # U4 goes to Z4
)
dataset2 <- data.frame(
  Instance = 1:total_instances,
  Reference_Partition = reference_partition,
  Trial_Partition = trial_partition
)

# dataset 3 (in Table 1)
# Define the total number of instances
total_instances <- 336

# Create vectors for the reference and trial partitions
reference_partition <- c(
  rep("cp", 5),  rep("im", 8),  rep("imL", 0), rep("imS", 1), rep("imU", 0), rep("om", 2), rep("omL", 0), rep("pp", 46),  # Z1
  rep("cp", 0),  rep("im", 0),  rep("imL", 1), rep("imS", 0), rep("imU", 0), rep("om", 18), rep("omL", 5), rep("pp", 1),  # Z2
  rep("cp", 137), rep("im", 1), rep("imL", 0), rep("imS", 0), rep("imU", 1), rep("om", 0),  rep("omL", 0), rep("pp", 4),  # Z3
  rep("cp", 1),  rep("im", 68), rep("imL", 1), rep("imS", 1), rep("imU", 34), rep("om", 0),  rep("omL", 0), rep("pp", 1)   # Z4
)

trial_partition <- c(
  rep("Z1", 62),
  rep("Z2", 25),
  rep("Z3", 143),
  rep("Z4", 106)
)
dataset3 <- data.frame(
  Instance = 1:total_instances,
  Reference_Partition = as.factor(reference_partition),
  Trial_Partition = as.factor(trial_partition)
)

#################
# Calculate all the class level partition metrics
#################
# Check the matching table
# print(table(dataset3$Reference_Partition, dataset3$Trial_Partition))

# Check the results
test_that("dataset1: RI,ARI,WC,AWC,WH,AWH match the expected values", {
  res <- getPartitionMetrics(dataset1$Reference_Partition, dataset1$Trial_Partition)
  expect_equal(round(as.numeric(res), 2), c(0.96,0.93,0.93,0.90,0.90,0.90))
})
test_that("dataset2: RI,ARI,WC,AWC,WH,AWH match the expected values", {
  res <- getPartitionMetrics(dataset2$Reference_Partition, dataset2$Trial_Partition)
  expect_equal(round(as.numeric(res), 2), c(0.74,0.54,0.54,0.36,0.36,0.36))
})
test_that("dataset3: RI,ARI,WC,AWC,WH,AWH match the expected values", {
  res <- getPartitionMetrics(dataset3$Reference_Partition, dataset3$Trial_Partition)
  expect_equal(round(as.numeric(res), 2), c(0.89,0.88,0.75,0.73,0.83,0.65))
})