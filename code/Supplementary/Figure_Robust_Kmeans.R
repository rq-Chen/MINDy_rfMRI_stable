## Setup

library(readr)
library(dplyr)
library(here)
library(wrsk)

# IO
mdl_name <- "HCP_Rest_FIX_Simple_Mdl200_sess"
mtf_csv <- paste0("motifs_", mdl_name, ".csv")
kmeans_csv <- paste0("kmeans_", mdl_name, ".csv")

# Load data
dat <- read_csv(here("data", mtf_csv))
kmeans_idx <- read_csv(here("data", kmeans_csv))

# Number of clusters and sparsity parameters to try
k_list <- seq(2, 8)
s_list <- seq(1.1, sqrt(ncol(dat)), 1)

# Parallelization
n_perm <- 10
# Use one core on windows
cores <- ifelse(
  .Platform$OS.type == "windows", 1, min(c(detectCores() - 2, n_perm))
)


## Overwrite wrsk() with warnings suppressed

# wrskGap() contains a tryCatch() call to wrsk() which
# returns "warning" if it encounters a warning, which does
# happen. This is a workaround to suppress the warnings.
source(here("code", "Supplementary", "my_wrskGap.R"))


## Robust sparse kmeans for all K

opt_gap <- rep(-Inf, length(k_list))
opt_s <- rep(-Inf, length(k_list))
solutions <- list()
for (i in seq_along(k_list)) {
  k <- k_list[i]
  res <- my_wrskGap(dat, k, s_list, npermute = n_perm, cores)
  tmp <- which.max(res$gap)
  tmp <- which(res$gap[1:tmp] > (res$gap[tmp] - res$se[tmp]))[1]
  opt_s[i] <- res$s[tmp]
  opt_gap[i] <- res$gap[tmp]
  solutions[[i]] <- res$resFinal[[tmp]]

  # Save after each K to avoid losing progress
  res <- list(
    opt_gap = opt_gap,
    opt_s = opt_s,
    solutions = solutions
  )
  saveRDS(res, here("data", paste0("robust_kmeans_", mdl_name, ".rds")))
}


## For k = 4

res <- my_wrskGap(dat, 4, s_list, npermute = n_perm, cores)
saveRDS(res, here("data", paste0("robust_kmeans_k4_", mdl_name, ".rds")))


# Skip the rest
q()


## Plots

library(ggplot2)
library(purrr)

# Read results
k4 <- readRDS(here("data", paste0("robust_kmeans_k4_", mdl_name, ".rds")))
tmp <- which.max(k4$gap)
tmp <- which(k4$gap[1:tmp] > (k4$gap[tmp] - k4$se[tmp]))[1]
opt_res <- k4$resFinal[[tmp]]

# Combine the obsweights in all solutions
obsweights <- k4$resFinal %>%
  imap(function(x, idx) data.frame(w = x$obsweights, s = k4$s[[idx]])) %>%
  list_rbind() %>%
  tibble() %>%
  mutate(s = factor(s)) %>%
  group_by(s)

# Generate a density plot for obsweights
options(repr.plot.width = 6, repr.plot.height = 6)
figdir <- here("figures", mdl_name)
if (!dir.exists(figdir)) dir.create(figdir)
obsweights %>%
  ggplot(aes(x = w, group = s, color = s)) +
  geom_density() +
  scale_color_discrete_sequential(name = "Sparsity") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", ) +
  labs(x = "Observed weights", y = "Density") +
  theme_minimal() +
  theme(text = element_text(size = 14))

ggsave(here(figdir, paste0("robust_kmeans_k4_obsweights_",
                           mdl_name, ".png")), dpi = 300)

# Save cluster centers
new_dat <- dat %>%
  mutate(obsweight = opt_res$obsweights) %>%
  mutate(cluster = factor(opt_res$outclusters)) %>%
  group_by(cluster)
new_dat %>%
  summarise(across(!obsweight, ~ weighted.mean(.x, obsweight))) %>%
  write_csv(here("data", paste0("robust_kmeans_k4_centers_", mdl_name, ".csv")))
