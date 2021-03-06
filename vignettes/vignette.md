MotifFinder: Vignette (Tutorial)
================
Daniel Wells
2019-03-01

## Simulate Data

We create 300 DNA sequences of length 200 with the motif “ATgTT\_GtCC”
around the center of 50% of the sequences.

``` r
library(MotifFinder)
set.seed(42)
simulated_sequences <- simulate_sequences(motif="ATgTT_GtCC")
```

    ## [1] "motif start position is 96"

``` r
str(simulated_sequences)
```

    ## List of 6
    ##  $ seqs          : Named chr [1:300] "CCTCTTTTGTATGGCCAAGGGGCGTAAAACTTAATCTCCGCTTTGCAGACTAGTGCGATTCGTGCGGCCGCTGGACCCAAGCATCCAGGCTCCCGGATTCCCTTACAAAAA"| __truncated__ "TCGGACTCCGTTGTGCAAACGTTTGCTGATAAGCCCGACCAAGCACTTCAGGAGTTAGTAAGACACATTGGCGAGCATAAGGCAGTTAGTAATAGGGACGAACATAAACGT"| __truncated__ "AGGTCACTATAGAGAATCTAATCCGATTGGGGTACATAATGCTCTTACGATTGGTTTGCCTACAAGTATCTACCTTGCACAGCACTCAAACTCTTACTAAAAGCACGACTA"| __truncated__ "TGTAACTCAGTACATCTACCGGGGTCTATACAAGCCTTCGGGTCATTATTAGGTAGGCCGGCCAATGGGACAAAGAATCTCTAACGCAGGCAATCGCTCGAGGCCTGGTCC"| __truncated__ ...
    ##   ..- attr(*, "names")= chr [1:300] "1" "2" "3" "4" ...
    ##  $ whichreg      : int [1:150] 226 94 252 273 147 117 210 220 21 251 ...
    ##  $ whichpos_s1   : int [1:150] 92 97 90 105 98 95 105 86 106 92 ...
    ##  $ whichpos      : num [1:150] 92 95 90 87 98 97 105 86 86 100 ...
    ##  $ whichrevstrand: int [1:150] 21 135 77 17 291 253 280 20 109 284 ...
    ##  $ truemotif     : chr "ATgTT_GtCC"

# Run MotifFinder

We run MotifFinder with a length slightly shorter than the known motif
length.

``` r
motif_found <- findamotif(simulated_sequences$seqs, len=7, stranded_prior = T, seed = 42)
```

## Plot the Motif(s) Found

We can see that we have recovered the motif.

``` r
library(ggseqlogo)
ggseqlogo(get_PWM(motif_found))
```

![](vignette_files/figure-gfm/plotlogo-1.png)<!-- -->

``` r
ggseqlogo(get_PWM(motif_found, complement=TRUE))
```

![](vignette_files/figure-gfm/plotlogo-2.png)<!-- -->

## Where is the motif

We can also check the location of the motifs found as well as the
location it thinks the motifs are.

``` r
plot_motif_location(motif_found)
```

![](vignette_files/figure-gfm/plotlocation-1.png)<!-- -->

``` r
plot(motif_found$prior)
```

![](vignette_files/figure-gfm/plotlocation-2.png)<!-- -->

## Check Convergence

It’s a good idea to plot the inferred probability the motif being in
each sequence (alpha parameter) for each iteration to check convergence.

``` r
plot(motif_found$alphas, ylim=c(0,1))
```

![](vignette_files/figure-gfm/plotalpha-1.png)<!-- -->

Other helper functions include downloading PWMs from the Jaspar or
Hocomoco databases.

``` r
Alx1 <- download_PWM("ALX1_MOUSE.H11MO.0.B")
str(Alx1)
```

    ## List of 2
    ##  $ pwm : num [1:12, 1:4] 0.3795 0.1211 0.6119 0.6785 0.0768 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:4] "A" "C" "G" "T"
    ##  $ name: chr "ALX1_MOUSE.H11MO.0.B"

## Extracting DNA from the genome

To use real DNA instead of simulating sequences you will need to
download the genome and extract regions that are specified using a BED
file.

``` bash
wget ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
samtools faidx motifs/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

bedtools getfasta -s -fi Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa -bed example_locations.bed -name > example_locations.fasta
```

We can then load these sequences using the helper function
load\_sequences

``` r
genomic_sequences <- load_sequences("example_locations.fasta")
```
