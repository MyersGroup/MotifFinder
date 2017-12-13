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