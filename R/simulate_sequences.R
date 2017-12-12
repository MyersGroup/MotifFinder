#' Simulate DNA sequences with an enriched motif
#'
#' @param motif character; string containing the motif sequence
#' @param number_sequences numeric; number of sequences to generate
#' @param sequence_length numeric; length in bases of the sequences to be generated
#' @param motif_position numeric; start position of the motif in the sequences
#' @param enrichment numeric; what fraction of sequences should include the motif
#'
#' @return a character vector of random DNA sequences with an enriched motif
#'
#' @export

simulate_sequences <- function(motif, number_sequences=300, sequence_length=200,
                                       motif_position=NULL, enrichment=0.2){

  # check sensible input
  stopifnot(nchar(motif) < sequence_length)
  stopifnot(motif_position + nchar(motif) < sequence_length)
  stopifnot(round(number_sequences * enrichment) > 1)

  # default position to center of sequences
  if(is.null(motif_position)){
    motif_position <- round(sequence_length / 2)
  }

  # create random sequences
  example_sequences <- character(length = number_sequences)
  for (i in seq_along(example_sequences)){
    example_sequences[i] <- rep(paste(sample(c("A","T","G","C"), sequence_length, T), collapse = ''))
  }

  # add enriched motif
  for (i in sample(1:number_sequences, round(number_sequences * enrichment))){
    substr(example_sequences[i], motif_position+1, motif_position+nchar(motif)) <- motif
  }

  return(example_sequences)

}
