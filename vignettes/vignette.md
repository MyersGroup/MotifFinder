MotifFinder: Vignette (Tutorial)
================
Daniel Wells
2017-12-09

Simulate Data
-------------

We create 300 DNA sequences of length 200. Then we add the motif "ATGCATGA" to ~10% of the sequences at the center.

``` r
# create random sequences
sequences <- character(length = 300)
for (i in seq_along(sequences)){
  sequences[i] <- rep(paste(sample(c("A","T","G","C"), 200, T), collapse = ''))
}

# add enriched motif at position 101 to 108 of ~10% of sequences
for (i in seq(from=1, to=length(sequences), by=10)){
  substr(sequences[i],101,108) <- "ATGCATGA"
}

str(sequences)
```

    ##  chr [1:300] "GCGAAGCGGGCGAACGCCCGGGGTAGTGCTAGGAATCATCACTTTACCAGTGTACGATGCGTGAATTTTTGTCTTCTTTACGAAATGTGTAATGCGAATCATGCATGATAG"| __truncated__ ...

Run MotifFinder
===============

We run MotifFinder with a length slightly shorter than the known motif and uniform scores.

``` r
motif_found <- findamotif(sequences, len=7, scores=rep(1,300))
```

Plot the Motif(s) Found
-----------------------

We can see that we have recovered the motif.

``` r
scorematset=exp(motif_found$scoremat)
compmat=scorematset[,c(4:1)]

scorematset=scorematset/rowSums(scorematset)
compmat=compmat/rowSums(compmat)

pwm = t(scorematset)
compmat=compmat[nrow(compmat):1,]

seqLogo::seqLogo(t(compmat))
```

![](vignette_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
seqLogo::seqLogo(pwm)
```

![](vignette_files/figure-markdown_github/unnamed-chunk-4-2.png)
