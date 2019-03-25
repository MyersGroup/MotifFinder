#' Calculate Position-Weight Matrix
#'
#' @param motif the output of findamotif()
#' @param complement logical; if TRUE the reverse complement PWM is returned
#'
#' @return matrix; a 4xN position-weight matrix where N is the length of the motif
#'
#' @export

get_PWM <- function(motif, complement=FALSE){

  pwm = exp(motif$scoremat)

  if(complement==TRUE){
    pwm = pwm[,c(4:1)]
    }

  pwm = pwm / rowSums(pwm)

  if(complement==TRUE){
    pwm = pwm[nrow(pwm):1,]
    }

  pwm = t(pwm)

  rownames(pwm) <- c("A","C","G","T")

  return(pwm)

}

#' Export Position-Weight Matrix to file
#'
#' @param pwm matrix; a PWM e.g. the output of get_PWM() OR the output of findamotif()
#' @param name character; a name for the motif to include in the file header
#' @param file character; path & filename to save file to e.g. "data/motif1.txt"
#' @param format character; the file format to use, either "uniprobe" or "meme"
#' @param complement logical; if input is a list from findamotif() should
#' the motif be the reverse complment
#'
#' @return NULL
#'
#' @export

export_PWM <- function(pwm, name, file, format="meme", complement=FALSE){

  # check if pwm is raw pwm or output of findmotif()
  if(typeof(pwm)=="list"){
    pwm <- get_PWM(pwm, complement=complement)
  }

  if(format=="uniprobe"){
  rownames(pwm) <- c("A:","C:","G:","T:")
  write.table(name, file,row.names=F, col.names = F, quote = F) # write header
  write.table(pwm, file, row.names=T, col.names = F, quote = F, append = T)
  }else{

    meme_header <- paste0("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n\nMOTIF ",
                          name,"\n\nletter-probability matrix: alength= 4 w= ",ncol(pwm)," nsites= 20 E= 0")

    write.table(meme_header, file,row.names=F, col.names = F, quote = F) # write header
    write.table(format(t(pwm), digits=5), file, row.names=F, col.names = F, quote = F, append = T)
  }
  message(paste("Motif exported to",file))

}

#' Download PWM from Jaspar or Hocomoco
#'
#' @param id string; The ID of the motif to download
#' @param pseudocount numeric; value of pseudocount to add to every entry of the PWM (to avoid 0 or 1 counts)
#'
#' @return A list of two items, a PWM matrix and a name
#'
#' @examples
#' download_PWM("ALX1_MOUSE.H11MO.0.B")
#'
#' @import jsonlite
#' @export
#'

download_PWM <- function(id, pseudocount=NULL){

  if(!grepl(".H11MO.",id)){ #database=="jaspar"

    motif <- jsonlite::fromJSON(paste0("http://jaspar.genereg.net/api/v1/matrix/",id,".json"))
    pwm <- matrix(unlist(motif$pfm), ncol = 4)
    colnames(pwm) <- names(motif$pfm)
    pwm <- pwm[,c("A", "C", "G", "T")]
    pwm <- pcm2pwm(pwm, pseudocount=pseudocount)
    return(list(pwm=pwm, name=motif$name))

  }else if(grepl(".H11MO.",id)){ # database=="hocomoco"

    motif <- jsonlite::fromJSON(paste0("http://hocomoco11.autosome.ru/motif/",id,"/pcm.json"))
    colnames(motif) <- c("A", "C", "G", "T")
    pwm <- pcm2pwm(motif, pseudocount=pseudocount)
    return(list(pwm=pwm, name=id))
  }

}

#' Convert PCM to PWM
#'
#' @param pcm matrix; n by 4 matrix of PCM / PFM
#' @param pseudocount numeric; value of pseudocount to add to every entry of the PWM (to avoid 0 or 1 counts),
#' or vector of pseudocounts to add to each nt position
#'
#' @return A position weight matrix
#'
#' @examples
#' pcm <- t(rbind(matrix(c(rep(10,8),rep(100,8)),nrow=2),
#'                matrix(c(rep(20,7),0,rep(200,7),0),nrow=2)))
#' pcm2pwm(pcm)
#'
#' @export
#'

pcm2pwm <- function(pcm, pseudocount=NULL){
  stopifnot(ncol(pcm)==4)
  if(is.null(pseudocount)){
    pseudocount <- sapply(log(rowSums(pcm)), function(x) max(2,x))
  }
  pcm <- pcm + pseudocount
  pcm <- pcm / rowSums(pcm)
  return(pcm)
}

#' Check if a motif is present within another motif
#'
#' @param submotif matrix; motif to search for, format: position weight matrix, n by 4 (A, C, G, T), (output of pcm2pwm)
#' (default finds CpG dinucleotides)
#' @param motif matrix; position weight matrix to check in
#' @param ispcm logical; if TRUE will convert pcm to pwm (default: FALSE)
#' @param probabilistic logical; if FALSE the sum of the convolution is taken,
#' if TRUE (default), the product of the sum of each row is returned.
#' TRUE will only work if the submotif is itself a pwm rather than a generic kernal.
#'
#' @return score of submotif being present in motif
#'
#'
submotif <- function(motif, submotif=rbind(c(0,1,0,0),c(0,0,1,0)), ispcm=FALSE, probabilistic=TRUE){
  if(ispcm){
    motif <- pcm2pwm(t(motif))
  }

  if(nrow(motif) < nrow(submotif)){
    stop("The submotif must be shorter than the motif")
  }

  shift <- nrow(submotif)-1

  if(probabilistic){
    max(sapply(1:(nrow(motif)-shift), function(x) prod(rowSums(motif[x:(x+shift),] * submotif))))
  }else{
    max(sapply(1:(nrow(motif)-shift), function(x) sum(motif[x:(x+shift),] * submotif))) / nrow(submotif)
  }

}


#' Convert PWM matrix to text (character string)
#'
#' @param pwm matrix; n by 4 matrix of log scale numeric values
#' @param threshold numeric; value between 0.5 and 1 to determine when to display capital vs lowercase letter.
#'
#' @return A charachter string
#'
#' @examples
#' # pwm2text(position_weight_matrix)
#'
#' @import data.table jsonlite
#' @export
#'

pwm2text <- function(pwm, threshold=0.7){

  if(ncol(pwm)!=4){stop("The PWM shuold be a n by 4 matrix")}
  if(threshold<0.5){stop("Threshold should be >= 0.5 in order to chose a single letter at each position.")}

  # get column index of best nucleotide (with >0.5)
  tmp <- as.numeric(apply(pwm, 1, function(x) which(x == max(x) & x>log(0.5))[1] ))

  # which nt have a value lower than the threshold
  which_small <- pwm[cbind(1:nrow(pwm),tmp)] < log(threshold)
  which_small[is.na(which_small)] <- FALSE

  # add 4 to the index so they get small letters
  tmp[which_small] <-  tmp[which_small] + 4

  dna <- c("A","C","G","T","a","c","g","t")

  tmp <- dna[tmp]
  tmp[is.na(tmp)] <- "_"
  tmp <- paste(tmp, collapse="")
  return(tmp)

}

#' Plot location of motifs in input sequences
#'
#' @param found_motif list; The output of getmotifs()
#' @param linepos integer; position of line in plot to mark e.g. the TSS
#' @param top_n integer; how many input sequences to plot
#' @param linesize numeric; thickness of lines representing motif footprints
#'
#' @return A ggplot2 object
#'
#' @import data.table ggplot2 RColorBrewer
#' @export
#'

plot_motif_location <- function(found_motif, linepos=NULL, top_n=NULL, linesize=1){

  tmp <- found_motif$dt[!is.na(whichpos)]

  for (i in seq_along(found_motif$scorematdim)){
    tmp[whichmotif==i, motifend := whichpos + found_motif$scorematdim[i] - 1]
  }

  tmp[,whichstrand := as.factor(whichstrand)]
  tmp[,whichmotif := factor(whichmotif)]
  if(!is.null(names(found_motif$scorematdim))){
    levels(tmp$whichmotif) <- names(found_motif$scorematdim)
  }
  tmp[,maxregprob := max(regprob),by=sequence]

  if(!is.null(top_n)){
    tmp <- tmp[order(-maxregprob)][1:top_n]
  }

  tmp$seqID <- as.numeric(factor(tmp$sequence, levels=unique(tmp[order(maxregprob)]$sequence)))

  if(length(found_motif$scorematdim)>1){
    plot_glob <- ggplot(tmp, aes(whichpos, seqID, colour=whichmotif))
  }else{
    plot_glob <- ggplot(tmp, aes(whichpos, seqID, colour=whichstrand))
  }

  if(!is.null(linepos)){
    plot_glob <- plot_glob + geom_vline(xintercept = linepos)
  }

  plot_glob <- plot_glob +
    geom_segment(aes(x = whichpos, y = seqID, xend = motifend, yend = seqID), size=linesize) +
    scale_color_brewer(palette = "Set1") +
    labs(x="Position in Sequence (bp)", y="Sequence / Gene") +
    theme(legend.position = "bottom")

  return(plot_glob)
}

#' Mask occurrences of motif
#'
#' @param found_motif list; The output of getmotifs()
#' @param motif integer; Which motif to mask if there is more than one in the found_motif list.
#'
#' @return Same as output of getmotifs() but with motifs removed.
#'
#' @export
#'

mask_motif <- function(found_motif,motif=1){
  str_sub(found_motif$seqs[found_motif$whichregs],
          start = found_motif$whichpos,
          end = found_motif$whichpos + found_motif$scorematdim[motif] - 1) <- paste(rep("N", found_motif$scorematdim[motif]), collapse = "")
  return(found_motif)
}

#' Export character vector to file in FASTA format
#'
#' @param sequences named character vector; Sequences to be exported
#' @param file string; Filename and path for sequences to be saved in
#'
#' @return Nothing, file written to disk.
#'
#' @examples
#'
#' dna <- c(seqA="ATGCTAG",seqB="ATCGATGTT",seqC="TCGATCGAT")
#' export_FASTA(dna, "dna.fasta")
#' @export
#'

export_FASTA <- function(sequences, file){
  names(sequences) <- paste0(">",names(sequences))
  write.table(c(rbind(names(sequences), sequences)), file, row.names = FALSE, col.names = FALSE, quote = FALSE)
}


#' Load DNA sequences from a FASTA file
#'
#' @param sm.rm logical; should soft masked bases (lower case) be converted to N
#' or capitalised, default=TRUE
#' @param file string; Filename and path of FASTA file of sequences
#' @param flank integer; how many bases either side of the center should be loaded
#' @param proportion_N numeric; fraction of N bases allowed before the
#' sequence is removed
#'
#' @return charachter vector of DNA sequences
#'
#' @import seqinr stringr
#'
#' @export
#'
load_sequences <- function(file, sm.rm=T, flank=150, proportion_N=0.02){
  seqs <- seqinr::read.fasta(file, forceDNAtolower = FALSE, as.string = TRUE)
  if(sm.rm){
    seqs <-  gsub("a|t|c|g|n","N",as.character(seqs))
  }else{
    seqs <-  toupper(as.character(seqs))
  }

  names(seqs) <- 1:length(seqs)

  seq_length <- nchar(seqs[1])

  seqs <- substr(seqs, as.integer(seq_length/2) - flank, as.integer(seq_length/2) + flank)

  number_masked <- stringr::str_count(seqs, "N")

  return(seqs[number_masked <= as.integer(seq_length * proportion_N)])
}


#' Extract Motif Matches
#'
#' @param found_motif list; The output of getmotifs()
#' @param motif integer; Which motif to mask if there is more than one in the found_motif list.
#' @param orderbyregprob logical; should the extracted motifs be ordered from regprob high to low
#'
#' @return charachter vector of DNA sequences matching the motif
#'
#' @import stringr stringi
#'
#' @export
#'
extract_matches <- function(found_motif, motif=1, orderbyregprob=F){
  matches <- stringr::str_sub(found_motif$seqs[found_motif$whichregs],
                   start = found_motif$whichpos,
                   end = found_motif$whichpos + found_motif$scorematdim[found_motif$whichmot] -1)[found_motif$whichmot==motif]
  matches[found_motif$whichstrand==0] <- stringi::stri_reverse(chartr("acgtACGT", "tgcaTGCA", matches[found_motif$whichstrand==0]))

  names(matches) <- found_motif$whichstrand

  if(orderbyregprob){
    matches[order(found_motif$regprob[found_motif$whichregs],decreasing = T)]
  }

  return(matches)
}


#' Plot Sequence logo denovo motif vs Tomtom match
#'
#' @param query_motif list object output from findamotif()
#' @param tomtom_match data.frame; result from reading in the text output of Tomtom
#' columns required: "Target.ID", "Orientation", "Optimal.offset"
#' @param yaxis logical; If FALSE yaxis values and labels are hidden
#' @param titles character vector; titles to use for each logo if not the default
#'
#' @return A plot, ggplot2 object
#'
#' @import ggseqlogo ggplot2
#'
#' @export
#'

plot_tomtom_match <- function(query_motif=NULL, tomtom_match=NULL, titles=NULL, yaxis=TRUE){

  # Take first entry of data.frame
  i=1

  # Download db motif
  db_motif <- download_PWM(tomtom_match[i]$Target_ID)
  db_pwm <- t(db_motif$pwm)

  if(tomtom_match[i]$Orientation=="-"){
    reverse_c <- TRUE
  }else{
    reverse_c <- FALSE
  }

  # get denovo motif pwm
  dn_pwm <- get_PWM(query_motif, reverse_c)

  # initialise padding to 0
  db_padd <- dn_padd <- db_padd_e <- dn_padd_e  <- 0

  offset <- tomtom_match[i]$Optimal_offset

  # calculate padding
  if(!reverse_c){
    if(offset < 0){
      db_padd <- abs(offset)
    }else{
      dn_padd <- abs(offset)
    }
  }else{ # if reverse complement NB offset is if the jaspar is reverse complmented

    reverse_offset <- max(ncol(db_pwm),ncol(dn_pwm)) -
      min(ncol(db_pwm),ncol(dn_pwm)) - abs(offset)

    if(abs(offset) + ncol(db_pwm) < ncol(dn_pwm)){
      db_padd <- abs(reverse_offset)
    }else{
      dn_padd <- abs(reverse_offset)
    }
  }

  z <- function(n){matrix(0, nrow=4, ncol=n)}

  db_pwm <- cbind(z(db_padd),db_pwm)
  dn_pwm <- cbind(z(dn_padd),dn_pwm)

  # Padd motifs to equal length
  length_diff <- ncol(db_pwm)-ncol(dn_pwm)
  if(length_diff<0){
    db_pwm <- cbind(db_pwm,z(abs(length_diff)))
  }else{
    dn_pwm <- cbind(dn_pwm,z(abs(length_diff)))
  }

  pwm_list <- list(db_pwm,dn_pwm)

  if(!is.null(titles)){
    db_title <- titles[i]
  }else if(db_motif$name==tomtom_match[i]$Target_ID){
    db_title <- strsplit(db_motif$name,"_")[[1]][[1]]
  }else{
    db_title <- paste(db_motif$name, tomtom_match[i]$Target_ID)
  }

  names(pwm_list) <- c("Tomtom Match","Denovo")
  #paste("MotifFinder Denovo",tomtom_match[i]$Query_ID,tomtom_match[i]$side)

  #print(ggseqlogo(pwm_list,ncol=1,scales = "free_y"))

  p <- ggseqlogo(pwm_list, ncol=1, scales = "free_y") +
    facet_wrap(~seq_group,
               scales = "free_y",
               ncol = 1,
               strip.position = "right")+
    ggtitle(db_title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.2, size=12),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          strip.background = element_blank())

  #size=18, strip.text = element_text(size = 8),

  if(!yaxis){
    p <- p + theme(axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   strip.text.y = element_blank())
  }

  return(p)

}


#' Draw a sparkline bar graph with unicode block characters
#'
#' Rendered using
#' [block elements](https://en.wikipedia.org/wiki/Block_Elements).
#' In most common fixed width fonts these are rendered wider than regular
#' characters which means they are not suitable if you need precise alignment.
#' Based on the function in the pillar package.
#'
#' @param x A numeric vector between 0 and 1
#' @param safe Nominally there are 8 block elements from 1/8 height to full
#'   height (8/8). However, the half-height and full-height blocks appear
#'   to be rendered inconsistently (possibly due to font substitution).
#' @examples
#' \dontrun{
#' x <- seq(0, 1, length = 6)
#' spark_bar(x)
#' spark_bar(sample(x))
#'
#' # This might work if you're lucky
#' spark_bar(seq(0, 1, length = 8), safe = FALSE)
#'
#' spark_bar(c(0, NA, 0.5, NA, 1))
#' }
#' @export

spark_bar <- function(x, safe = TRUE) {
  stopifnot(is.numeric(x))

  bars <- vapply(0x2581:0x2588, intToUtf8, character(1))
  if (safe) {
    bars <- bars[-c(4, 8)]
  }

  factor <- cut(
    x,
    breaks = seq(0, 1, length.out = length(bars) + 1),
    labels = bars,
    include.lowest = TRUE
  )
  chars <- as.character(factor)
  chars[is.na(chars)] <- bars[length(bars)]

  structure(paste0(chars, collapse = ""), class = "spark")
}


#' Confusion Matrix
#'
#' @param found_motif list; The output of getmotifs() or findamotif()
#' @param simulated_sequences list; Output of simulate_sequences
#' @param complement logical; Did findamotif return the complemented version of the simulated motif
#' Change this if you get 0 true positives
#'
#' @return List of confusion matricies, one for overall and one each for each strand
#'
#' @details Note that for a positive match, all of the sequence, strand and position of the motif.
#' A mismatch in any of these will exclude a predicted +ve true +ve. However the total number of true -ve does not
#' count each position as one negative but each sequence which does not contain an added motif.
#'
#' Also note that there could be matches to your simulated motif by chance in random sequence, and as the motif
#' generation is probabilistic, some sequences with the motif added will not be a good match to the consensus motif.
#'
#' @export
#'

cfm <- function(motif_found, simulated_sequences, complement=F){

  simulated_sequences$strand <- (!1:length(simulated_sequences$seqs) %in% simulated_sequences$whichrevstrand)*1

  if(complement){
    motif_found$whichstrand[motif_found$whichstrand==0] <- -1
    motif_found$whichstrand[motif_found$whichstrand==1] <- 0
    motif_found$whichstrand[motif_found$whichstrand==-1] <- 1
  }

  reg_pos_s <- paste0(simulated_sequences$whichreg,"_",
                      simulated_sequences$whichpos,"_",
                      simulated_sequences$strand[simulated_sequences$whichreg])

  pre_pos_s <- paste0(motif_found$whichregs,"_",
                      motif_found$whichpos,"_",
                      motif_found$whichstrand)

  trueP <- sum(pre_pos_s %in% reg_pos_s)

  falseP <- sum(!pre_pos_s %in% reg_pos_s)

  totalP <- length(reg_pos_s)

  falseN = totalP - trueP

  trueN <- length(simulated_sequences$seqs)*(nchar(simulated_sequences$seqs[[1]])-nchar(simulated_sequences$truemotif))*2 - totalP - falseP


  #####


  trueP_s1 <- sum(pre_pos_s
                  %in%
                    reg_pos_s[simulated_sequences$strand[simulated_sequences$whichreg]==1])

  totalP_s1 <- sum(simulated_sequences$strand[simulated_sequences$whichreg]==1) # == sum(grepl("_1$",reg_pos_s))

  falseN_s1 = totalP_s1 - trueP_s1 # we do make predictions of motif for some of these sequences, but not at the right position

  falseP_s1 <- sum(!motif_found$whichregs %in% simulated_sequences$whichrevstrand) - trueP_s1  # i.e. total predicted + that are also strand==1 (!=0), - correct

  trueN_s1 <- sum(simulated_sequences$strand==1)*(nchar(simulated_sequences$seqs[[1]])-nchar(simulated_sequences$truemotif))*2 - totalP_s1 - falseP_s1 # should == total pred -ve, - FalseN_s1

  #####

  pred_neg <- c(1:1000)[!1:1000 %in% motif_found$whichregs]
  pred_neg_s1 <- pred_neg[!pred_neg %in% simulated_sequences$whichrevstrand]
  pred_pos_s1 <- motif_found$whichregs[!motif_found$whichregs %in% simulated_sequences$whichrevstrand]

  conf_mat_dimnames <- list(c("Pred +","Pred -"),c("True +","True -"))

  s1=matrix(c(trueP_s1, falseN_s1, falseP_s1, trueN_s1), ncol=2, dimnames = conf_mat_dimnames)
  overall=matrix(c(trueP, falseN, falseP, trueN), ncol=2, dimnames = conf_mat_dimnames)
  s0 = overall-s1

  return(list("overall"=overall,
              "s1"=s1,
              "s0"=s0))
}
