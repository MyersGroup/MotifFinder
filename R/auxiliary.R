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
#' @param database string; Either "jaspar" or "hocomoco".
#' @param pseudocount numeric; value of pseudocount to add to every entry of the PWM (to avoid 0 or 1 counts) (Applied for jaspar only)
#'
#' @return A PWM matrix
#'
#' @examples
#' download_PWM("ALX1_MOUSE.H11MO.0.B", database="hocomoco")
#'
#' @import data.table jsonlite
#' @export
#'

download_PWM <- function(id, database="jaspar", pseudocount=5){

  if(database=="jaspar"){

    pwm <- fread(paste0("http://jaspar.genereg.net/api/v1/matrix/",id,".pfm"))
    pwm <- as.matrix(pwm + pseudocount) / colSums(pwm + pseudocount)
    rownames(pwm) <- c("A", "C", "G", "T")
    pwm <- t(pwm)

  }else if(database=="hocomoco"){

    pwm <- jsonlite::fromJSON(paste0("http://hocomoco11.autosome.ru/motif/",id,"/pwm.json"))
    pwm <- as.matrix(pwm)
    colnames(pwm) <- c("A", "C", "G", "T")
    pwm <- exp(pwm)

  }else{
    warning("database not recognised should be either 'jaspar' or 'hocomoco'")
  }

  return(pwm)

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
