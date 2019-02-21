#' Simulate DNA sequences with an enriched motif
#'
#' @param motif character; string containing the motif sequence, can be upper or lower case ACGT or _
#' @param number_sequences numeric; number of sequences to generate
#' @param sequence_length numeric; length in bases of the sequences to be generated
#' @param motif_position numeric; start position of the motif in the sequences,
#' default is half the sequence length - half the motif width
#' @param enrichment numeric; what fraction of sequences should include the motif
#' @param jitter integer; how many base pairs to randomly +/- shift the position of the motif
#' @param highprob numeric; probability (between 0 and 1) of choosing a base when the motif has a capital/uppercase letter
#' @param lowprob numeric; probability (between 0 and 1) of choosing a base when the motif has a lowercase letter
#'
#' @return a character vector of random DNA sequences with an enriched motif
#'
#' @export

simulate_sequences <- function(motif, number_sequences=300, sequence_length=200,
                                       motif_position=NULL, enrichment=0.5, jitter=10, highprob=0.85, lowprob=0.6){

  # check sensible input
  stopifnot(nchar(motif) < sequence_length)
  stopifnot(all(strsplit(motif, split = "")[[1]] %in% c("A","T","G","C","a","t","c","g","_")))

  stopifnot(motif_position + nchar(motif) < sequence_length)
  stopifnot(motif_position >= 0)

  stopifnot(round(number_sequences * enrichment) > 1)
  stopifnot(enrichment <= 1)


  # default position to center of sequences
  if(is.null(motif_position)){
    motif_position <- round((sequence_length+1)/2 - nchar(motif)/2)
    print(paste("motif start position is", motif_position))
  }

  DNA <- c("A","T","G","C")

  # create random sequences
  example_sequences <- character(length = number_sequences)
  for (i in seq_along(example_sequences)){
    example_sequences[i] <- rep(paste(sample(DNA, sequence_length, T), collapse = ''))
  }

  # function to sample nucleotide with biased probability
  DNA_withprobability <- function(nucleotide, probability){
    nucleotide <- toupper(nucleotide)
    prob_vector <- numeric(length = 4)
    prob_vector[which(DNA == nucleotide)] <- probability
    prob_vector[which(DNA != nucleotide)] <- (1-probability)/3
    sample(DNA, 1, prob = prob_vector)
  }

  # function to generate DNA sequence according to motif code given
  generate_motif_string <- function(motif){
    motif_string <- character()
    for(i in 1:nchar(motif)){

      letter_code <- substr(motif,i,i)

      if(letter_code=="_"){ # equal probability
        motif_string <- c(motif_string, sample(DNA,1))

      }else if(toupper(letter_code)==letter_code){
        motif_string <- c(motif_string, DNA_withprobability(letter_code, highprob))

      }else{
        motif_string <- c(motif_string, DNA_withprobability(letter_code, lowprob))

      }
    }
    return(paste(motif_string, collapse = ""))
  }

  # add enriched motif
  for (i in sample(1:number_sequences, round(number_sequences * enrichment))){

    # change motif position slightly
    # careful not to sample negative and leave enough space for motif on right
    # if jitter=0, then n:n becomes a single number and sample samples from 1:n, so skip
    if(jitter>0){
      motif_position_random <- sample(max(0,(motif_position-jitter)):min((motif_position+jitter), sequence_length-nchar(motif)),1)
    }else{
      motif_position_random <- motif_position
    }

    substr(example_sequences[i], motif_position_random, motif_position_random+nchar(motif)-1) <- generate_motif_string(motif)
  }

  }

  names(example_sequences) <- seq_along(example_sequences)

  return(example_sequences)

}
