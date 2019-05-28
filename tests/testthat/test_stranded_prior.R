
context("Stranded Prior Works")

set.seed(42)
simulated_sequences <- simulate_sequences("ATgTT_GtCC", motif_position = 50, jitter = 0)

# run MotifFinder
motif_found <- findamotif(simulated_sequences$seqs, len=7, seed=258442, stranded_prior = T, motif_seed = "modal")

test_that("Prior is highest in correct bins", {
  expect_gte(min(rank(motif_found$prior)[8]), 10)
  expect_lt(motif_found$prior[3], 0.1)
})

verbosity0 <- capture_output(findamotif(simulated_sequences$seqs, len=7, seed=258442, nits = 3, verbosity = 0))
verbosity1 <- capture_output(findamotif(simulated_sequences$seqs, len=7, seed=258442, nits = 3, verbosity = 1))
verbosity2 <- capture_output(findamotif(simulated_sequences$seqs, len=7, seed=258442, nits = 3, verbosity = 2))
verbosity3 <- capture_output(findamotif(simulated_sequences$seqs, len=7, seed=258442, nits = 3, verbosity = 3))

test_that("Verbosity Flag works", {
  expect_equal(nchar(verbosity0), 0)
  expect_gt(nchar(verbosity1), 0)
  expect_gt(nchar(verbosity2), nchar(verbosity1))
  expect_gt(nchar(verbosity3), nchar(verbosity2))
})

examples <- extract_matches(motif_found, orderbyregprob = T)

test_that("Extraction nhar=motif length", {
  expect_equal(unique(nchar(extract_matches(motif_found))), 10)
})


test_that("Warning when no motif exists", {
  set.seed(42)
  simulated_sequences <- simulate_sequences("_", sequence_length = 100)

  expect_warning(findamotif(simulated_sequences$seqs, len=7, seed=258443, stranded_prior = T, motif_seed = "random"), "No motifs remaining")
})

test_that("Warning when motif too long", {
simulated_sequences <- simulate_sequences("AGCAGCTAGCTAGCTAAGCATCAGCGAGCAGCCACAGCACAGCATCAGCTAGTCGATATA", sequence_length = 100)

expect_warning(findamotif(simulated_sequences$seqs, len=7, seed=258442, motif_seed = "modal"), "No motifs remaining")
})

