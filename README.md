# MotifFinder: Find enriched motifs in a set of DNA sequences

By Simon Myers & Nicolas Altemose

This is an R package for finding enriched motifs in a set of DNA sequences using an iterative Gibbs sampler described in [Altemose et al. eLife 2017](https://elifesciences.org/articles/28383).

This code is ported from the original at [https://github.com/altemose/PRDM9-map](https://github.com/altemose/PRDM9-map)

If you use this program, please cite [Altemose et al. eLife 2017](https://elifesciences.org/articles/28383).

This is free software shared in the hope it may be of use; no warranty is given or implied.

## Installation & Usage

```R
# install.packages("devtools")
devtools::install_github("myersgroup/MotifFinder")

# simulate set of sequences enriched for a motif
set.seed(42)
simulated_sequences <- simulate_sequences(motif="ATGCATGA")

# run MotifFinder
motif_found <- findamotif(simulated_sequences, len=7)

# visualise the motif found
seqLogo::seqLogo(get_PWM(motif_found))
```

![](vignettes/vignette_files/figure-markdown_github/unnamed-chunk-4-1.png)
