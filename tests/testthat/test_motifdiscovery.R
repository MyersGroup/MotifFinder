
context("Motif Discovery Works")

set.seed(42)
simulated_sequences <- simulate_sequences("ATgTT_GtCC", number_sequences = 300, sequence_length = 200,
                                          motif_position = NULL, enrichment = 0.5, jitter = 10, highprob = 0.85,
                                          lowprob = 0.6)

test_that("Simulated sequence lengths are correct", {
  expect_true(all(nchar(simulated_sequences$seqs)==200))
  expect_equal(length(simulated_sequences$seqs),300)
})


test_that("Simulated sequence is DNA only", {
  expect_equal(sort(unique(strsplit(paste0(simulated_sequences$seqs, collapse = ""),"")[[1]])),
               c("A","C","G","T"))
})

test_that("Simulated sequence indexes are within range", {
  expect_true(all(simulated_sequences$whichreg>0))
  expect_true(all(simulated_sequences$whichreg<=300))

  expect_true(all(simulated_sequences$whichpos>0))
  expect_true(all(simulated_sequences$whichpos<=300))

  expect_true(all(simulated_sequences$whicgpos_s1>0))
  expect_true(all(simulated_sequences$whicgpos_s1<=300))

  expect_true(all(simulated_sequences$whicgrevstrand>0))
  expect_true(all(simulated_sequences$whicgrevstrand<=300))

})

# run MotifFinder
motif_found <- findamotif(simulated_sequences$seqs, len=7, seed=258442)

test_that("Ranges looks ok", {
  expect_lte(min(motif_found$prior), 1)
  expect_gte(min(motif_found$prior), 0)
  expect_true(all(rowSums(exp(motif_found$scoremat)) == 1))
  expect_lte(min(motif_found$regprobs), 1)
  expect_gte(min(motif_found$regprobs), 0)
  expect_true(all(motif_found$beststrand %in% c(0,1)))
})


test_that("Precision and recall are ok", {
  conf_mat <- cfm(motif_found, simulated_sequences, complement = T)$overall
  expect_gt((conf_mat/colSums(conf_mat)[1])[1,1], 0.5)
  expect_lt((conf_mat/rowSums(conf_mat)[1])[1,2], 0.5)
})


test_that("PWM extraction works", {
  expect_equal(dim(get_PWM(motif_found)), c(4,10))
  expect_equal(dim(get_PWM(motif_found, complement = T)), c(4,10))
  expect_equal(rownames(get_PWM(motif_found)), c("A", "C", "G", "T"))
  expect_gt(min(get_PWM(motif_found)), 0)
  expect_lt(max(get_PWM(motif_found)), 1)
})

test_that("pcm2pwm works", {
  pcm <- t(rbind(matrix(c(rep(10,8),rep(100,8)),nrow=2),
               matrix(c(rep(20,7),0,rep(200,7),0),nrow=2)))
  pwm <- pcm2pwm(pcm)
  expect_gt(min(pwm), 0)
  expect_lt(max(pwm), 1)
  expect_setequal(rowSums(pwm), 1)
  expect_equal(nrow(pwm), 8)
  expect_equal(ncol(pwm), 4)
})

test_that("alpha in expected range (given seed)", {
  expect_lt(motif_found$alpha, 0.7)
  expect_gt(motif_found$alpha, 0.6)
})

test_that("Prior is highest in center", {
  expect_gte(min(rank(motif_found$prior)[5:6]), 9)
  expect_gt(sum(motif_found$prior[c(5,6)]),0.6)
})

test_that("Motif is found", {
  expect_equal(pwm2text(motif_found$scoremat, threshold = 0.7), "GGaC_AacAT")
})

test_that("Can find Motif location when updatemot=0", {
  motif_locs <- getmotifs(motif_found$scoremat,
            dimvec=nrow(motif_found$scoremat),
            seqs=simulated_sequences$seqs,
            updatemot = 0,
            maxwidth=2200,
            ourprior=rep(0.1,10),
            seed=258442)

  expect_equal(pwm2text(motif_locs$scoremat, threshold = 0.7), "GGaC_AacAT")

  expect_equal(mean(motif_found$whichpos), expected = 100, tolerance=5, scale=1)
})

# Jaspar server too slow error 503
# test_that("Download PWM from jaspar works", {
#   pwm <- download_PWM("MA0506.1")
#   expect_equal(ncol(pwm$pwm), 4)
#   expect_equal(colnames(pwm$pwm), c("A", "C", "G", "T"))
#   expect_gt(nchar(pwm$name), 1)
#   expect_gt(min(pwm$pwm), 0)
#   expect_lt(max(pwm$pwm), 1)
# })

test_that("Download PWM from hocomoco works", {
  pwm <- download_PWM("ALX1_MOUSE.H11MO.0.B")
  expect_equal(ncol(pwm$pwm), 4)
  expect_equal(colnames(pwm$pwm), c("A", "C", "G", "T"))
  expect_gt(nchar(pwm$name), 1)
  expect_gt(min(pwm$pwm), 0)
  expect_lt(max(pwm$pwm), 1)
})

test_that("Can plot motif location", {
  expect_true(is.ggplot(plot_motif_location(motif_found)))
  expect_true(is.ggplot(plot_motif_location(motif_found, linepos = 50, top_n = 10)))
  expect_equal(class(plot_motif_location(motif_found)$layers[[1]]$geom)[1], "GeomSegment")
})
