#' Calculate Position-Weight Matrix
#'
#' @param motif the output of findamotif()
#' @param complment logical; if TRUE the reverse complement PWM is returned
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
#' download_PWM("MA0506.1")
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

#' Convert PWM matrix to text (character string)
#'
#' @param pwm matrix; n by 4 matrix of log scale numeric values
#' @param treshold numeric; value between 0.5 and 1 to determine when to display capital vs lowercase letter.
#'
#' @return A charachter string
#'
#' @examples
#' pwm2text(position_weight_matrix)
#'
#' @import data.table jsonlite
#' @export
#'

pwm2text <- function(pwm, threshold=0.7){

  if(ncol(pwm)!=4){stop("The PWM shuold be a n by 4 matrix")}
  if(threshold<0.5){stop("Threshold should be >= 0.5 in order to chose a single letter at each position.")}

  # get column index of best nucleotide (with >0.5)
  tmp <- as.numeric(apply(pwm, 1,function(x) which(x == max(x) & x>log(0.5))))

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
#'
#' @return A ggplot2 object
#'
#' @import data.table ggplot2 RColorBrewer
#' @export
#'

plot_motif_location <- function(found_motif, linepos=NULL, top_n=NULL){

  tmp <- found_motif$dt[!is.na(whichpos)]

  for (i in seq_along(found_motif$scorematdim)){
    tmp[whichmotif==i, motifend := whichpos + found_motif$scorematdim[i]]
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
    transparency_alpha <- 0.5
  }else{
    plot_glob <- ggplot(tmp, aes(whichpos, seqID, colour=whichstrand))
    transparency_alpha <- 1
  }

  if(!is.null(linepos)){
    plot_glob <- plot_glob + geom_vline(xintercept = linepos)
  }

  plot_glob <- plot_glob +
    geom_segment(aes(x = whichpos, y = seqID, xend = motifend, yend = seqID), size=1, alpha=transparency_alpha) +
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
