
context("Sequence Simulation Checks")

test_that("Sequences are correct length and number", {
  expect_length(simulate_sequences("AAA", number_sequences = 10)$seqs, 10)
  expect_setequal(nchar(simulate_sequences("AAA", sequence_length = 50)$seqs), 50)
})



test_that("Sequences are only ACTG", {
  expect_true(grepl("^[ACTG]+$",paste0(simulate_sequences(motif="ACGT")$seqs, collapse = "")))
})

test_that("Enrichment is proportion containing motif when highprob=1 & All caps letters", {
  expect_equal(sum(grepl("GGGGGGGGGG",
                         simulate_sequences(motif="GGGGGGGGGG",
                                            highprob=1,
                                            enrichment=0.5,
                                            number_sequences=100, randomstrand = F)$seqs, fixed = T)),
               expected = 50, tolerance=1, scale=1)

  expect_equal(sum(grepl("GGGGGGGGGG",
                         simulate_sequences(motif="GGGGGGGGGG",
                                            highprob=1,
                                            enrichment=1,
                                            number_sequences=100, randomstrand = F)$seqs, fixed = T)),
               expected = 100, tolerance=1, scale=1)
})

test_that("At least one motif is inserted", {
  expect_error(simulate_sequences(motif="GGG",enrichment = 0.014, number_sequences = 100)$seqs)
})
