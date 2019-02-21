
context("Motif Discovery Works")

set.seed(42)
simulated_sequences <- simulate_sequences("ATgTT_GtCC", number_sequences = 300, sequence_length = 200,
                                          motif_position = NULL, enrichment = 0.5, jitter = 10, highprob = 0.85,
                                          lowprob = 0.6)

test_that("Simulated sequences are exactly the same, given seed & parameters above", {
  expect_known_hash(simulated_sequences,"a19d5c9005")
})

# run MotifFinder
motif_found <- findamotif(simulated_sequences, len=7, seed=258442)

test_that("results are exactly the same, given seed's above", {
  expect_known_hash(motif_found$prior,"fb622462d8")
  expect_known_hash(motif_found$alphas,"e88bde2511")
  expect_known_hash(motif_found$scoremat,"1d39f6fb9b")
  expect_known_hash(digest::digest(signif(motif_found$regprobs)),"4074286e8e")
  expect_known_hash(motif_found$whichpos,"df95dad424")
  expect_known_hash(motif_found$beststrand,"ae9dcdf14b")
})

test_that("PWM extraction works", {
  expect_known_hash(get_PWM(motif_found),"fe855048f2")
  expect_equal(dim(get_PWM(motif_found)), c(4,10))
  expect_equal(rownames(get_PWM(motif_found)), c("A", "C", "G", "T"))
  expect_gt(min(get_PWM(motif_found)), 0)
  expect_lt(max(get_PWM(motif_found)), 1)
})

test_that("PWM extraction works", {
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
  expect_lt(motif_found$alpha, 0.8)
  expect_gt(motif_found$alpha, 0.7)
})

test_that("Prior is highest in center", {
  expect_equal(rank(motif_found$prior)[5:6], c(9,10))
})

test_that("Motif is found", {
  expect_equal(pwm2text(motif_found$scoremat, threshold = 0.7), "ATgTT_GtCC")
})

test_that("Can find Motif location when updatemot=0", {
  motif_locs <- getmotifs(motif_found$scoremat,
            dimvec=nrow(motif_found$scoremat),
            seqs=simulated_sequences,
            updatemot = 0,
            maxwidth=2200,
            ourprior=rep(0.1,10),
            seed=258442)

  expect_known_hash(motif_locs,"c74a705ec0")
  expect_known_hash(motif_locs$whichpos,"5115dc46ac")

  expect_equal(mean(motif_found$whichpos), expected = 100, tolerance=5, scale=1)
  expect_equal(mean(motif_found$bestmatch), expected = 100, tolerance=5, scale=1)
})

test_that("Download PWM from jaspar works", {
  pwm <- download_PWM("MA0506.1")
  expect_equal(ncol(pwm$pwm), 4)
  expect_equal(colnames(pwm$pwm), c("A", "C", "G", "T"))
  expect_gt(nchar(pwm$name), 1)
  expect_gt(min(pwm$pwm), 0)
  expect_lt(max(pwm$pwm), 1)
})

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
  expect_equal(class(plot_motif_location(motif_found)$layers[[1]]$geom)[1], "GeomSegment")
})
