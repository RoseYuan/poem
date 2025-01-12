## code to prepare `noisy_moon` dataset goes here
library(dplyr)
moon_maker <- function(n = 50, noise = 0.1,
                       x_center = 0, y_center = 0, radius = 1, seed=42) {
  set.seed(seed)
  moon <- tibble(
    i = seq_len(n),
    x = x_center + radius * cos(pi * i/n) + rnorm(n, 0, sd = noise * radius),
    y = y_center + radius * sin(pi * i/n) + rnorm(n, 0, sd = noise * radius),
  )
  return(data.frame(moon))
}
# Combine the data into a data frame for easier manipulation
df1 <- moon_maker(x_center = 1, y_center = -1)
df1$label = 1
df2 <- moon_maker(x_center = 2, y_center = 0.5)
df2$label = 2
df2$y <- - df2$y
data <- rbind(df1, df2)
data$label <- as.factor(data$label)
noisy_moon <- data[, c("x", "y", "label")]


# Perform K-means clustering
kmeans_result <- kmeans(data, centers = 2)
kmeans_labels <- kmeans_result$cluster
noisy_moon$kmeans_label <- factor(kmeans_labels)

# Perform HDBSCAN clustering
library(dbscan)
hdbscan_result <- hdbscan(as.matrix(data[,c("x", "y")]), minPts = 5)
hdbscan_labels <- hdbscan_result$cluster
noisy_moon$hdbscan_label <- factor(hdbscan_labels)

# Plot the dataset
ggplot(noisy_moon, aes(x = x, y = y, color=label)) +
  geom_point() +
  labs(title = "Noisy Moons Dataset")

# Plot the K-means clustering result
ggplot(noisy_moon, aes(x = x, y = y, color = kmeans_label)) +
  geom_point() +
  labs(title = "K-means Clustering", color = "Cluster")

# Plot the HDBSCAN clustering result
ggplot(noisy_moon, aes(x = x, y = y, color = hdbscan_label)) +
  geom_point() +
  labs(title = "HDBSCAN Clustering", color = "Cluster")

usethis::use_data(noisy_moon, overwrite = TRUE)
