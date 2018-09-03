library(tidyverse)

n_pop <- 3
pop_sizes <- c(1000, 1000, 1000, 1000) # sizes of the different populations
pop_means <- c(1,2,3,4) # means of the different populations
pop_sds <- c(0.5, 0.5, 0.5, 0.5) # sds of the different populations
pop_params <- list(pop_means, pop_sds, pop_sizes)

n_markers <- 20 # number of markers
n_cells <- sum(pop_sizes) # total number of cells

create_data <- function(mean, sd, size){
  matrix(rnorm(size * n_markers, mean, sd), ncol = n_markers)
}

sim_data <- do.call(rbind, pmap(pop_params, ~create_data(..1, ..2, ..3)))

sim_data <- asinh((sim_data^10)/5)

colnames(sim_data) <- paste0("Marker", seq_len(n_markers))

sample_ids <- rep(LETTERS[1:2], c(2000, 2000))
bio_ids <- rep(c("T", "B", "T", "B"), each = 1000)

pca <- prcomp(sim_data)$x[,1:2]
pca <- data.frame(pca, sample = factor(sample_ids), bio = bio_ids)
pca

ggplot(pca, aes(x = PC1, y = PC2, col = sample, shape = bio)) + 
  geom_point() + 
  theme_bw()

M <- ruv::replicate.matrix(bio_ids)

norm_data <- RUVIII(Y = sim_data, M = M, ctl = 1:n_markers, k = 1)$new_Y

pca <- prcomp(norm_data)$x[,1:2]
pca <- data.frame(pca, sample = factor(sample_ids), bio = bio_ids)
pca

ggplot(pca, aes(x = PC1, y = PC2, col = sample, shape = bio)) + 
  geom_point() + 
  theme_bw()

