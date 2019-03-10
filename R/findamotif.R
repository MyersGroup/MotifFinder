#' Find the single most enriched motif in a set of DNA sequences
#'
#' @param seqs a vector of strings giving the DNA sequences in which to find a motif
#' @param len length of motif to find (min=4)
#' @param scores a set of regional scores giving weights; e.g. ChIP-Seq enrichment values
#' @param nits number of iterations used for motif refinement
#' @param plen a parameter setting the geometric prior on how long each motif found should be. plen=0.05 corresponds to a mean length of 20bp and is the default. Setting plen large penalises longer motifs more
#' @param n_for_refine the top n_for_refine scoring regions only are used for motif refinement
#' @param prior a vector of length 10 probabilities giving the initial probability of a motif being found across different parts of the sequence from start:end. If left unspecified the initial prior is set at uniform and the algorithm tries to learn where motifs are, e.g. if they are centrally enriched.
#' @param updateprior a flag - should the algorithm update (learn) the prior on where the motifs occur within the DNA sequences(default is 1)?
#' @param seed integer; seed for random number generation, set this for exactly reproducible results.
#' @param verbosity integer; How verbose should this function be? 0=silent, 3=everything.
#' @param motif_blacklist charachter vector; motifs not to use as seed motif
#' @param motif_rank integer; which rank of seed motif to use (1st seed motif, 2nd etc.)
#' @param range integer; range around center to check for central enrichment
#' @param motif_seed string; "central", "modal", "random", or a string e.g. "ACGTGAC"
#'
#' @details
#' This function identifies a single PWM from an iterative Gibbs sampler described in Altemose et al. eLife 2017. Function 2 can refine multiple motifs further, jointly.
#'
#' The user must input a set of DNA sequences, a score for each sequence (e.g. an enrichment value or any other score), and a length for an initial motif (e.g. 8 bp) used to seed the algorithm.
#'
#' There are additional optional parameters.
#'
#' The program outputs a list of results, including information on the inferred PWM (i.e. motif found), as well as a probabilistic output of which regions contain this motif, and posterior distributions of the other parameters
#'
#'
#' @return List item with the following items:\cr
#' Details of input data given:
#' * seqs: the vector of input sequences used for finding motifs within
#' * trimmedseqs: the vector of input sequences used for finding motifs within, after trimming to shorten long input sequences
#'
#' Details of overall fitted model:
#' * scoremat: a matrix giving the pwm (log-scale) for the identified motif after iteration
#' * scorematdim: the length of the identified motif, and scoremat is of dimension scorematdimx4
#' * prior: a vector of length 10 probabilities giving the inferred probability of a motif being found across different parts of the sequence from start to end.
#' * alpha: a vector of probabilities giving the inferred probability of the motif being found within each input region
#' * bindmat: a version of scoremat accounting for the background sequence composition
#' * background is the inferred background model
#'
#' Details of output for given data:
#' * regprobs, regprob are in this case identical vectors giving the probability of the motif occurring in each given input sequence
#' * bestpos is a vector giving the best match to the motif in each given input sequence
#' * whichregs is a vector showing which input sequences had motifs identified in the final round of sampling of the Gibbs sampler
#' * whichpos: for motifs identified in regions described in whichreg, the start positions of motifs identified in the final round of sampling of the Gibbs sampler
#' * whichmot: not needed in this case
#' * whichstrand: for motifs identified in regions described in whichreg, the strand associated with motifs identified in the final round of sampling of the Gibbs sampler, relative to the input sequence
#'
#' @export
#' @import gtools

findamotif=function(seqs,len,scores=NULL,nits=50,scoring_its=5,n_for_refine=1000,prior=NULL,updateprior=1,plen=0.9,seed=NULL,verbosity=1, motif_rank=1,motif_blacklist=NULL,range=50,stranded_prior=F, motif_seed="central",conv_t=0, conv_n=200){

  if (is.null(seed)){
    seed <- sample.int(2^20, 1)
  }

  if (!is.na(seed)){
    set.seed(seed)
  }

  if(is.null(scores)){
    scores <- rep(1,length(seqs))
  }


  if(n_for_refine>length(seqs)) n_for_refine=length(seqs)
  regs=seqs

  if(motif_seed=="central"){
  if(verbosity>=3) print("Concatenating sequences....")
  seqs=paste(seqs,collapse="")
  if(verbosity>=3) print("....done")

  if(verbosity>=3) print("Replacing DNA letters with numbers....")
  seqv=as.vector(unlist(strsplit(seqs,"")))
  seqv = as.integer(chartr("ACGTacgtN","012345678",seqv))
  if(verbosity>=3) print("....done")

  if(verbosity>=3) print("Indexing....")
  nonrep=rep(0,(length(seqv)-len+1))
  for(i in 1:len){
    if(verbosity>=3) print(i)
    z=seqv[i:(i+length(nonrep)-1)]
    nonrep=nonrep+4^(len-i)*z
    nonrep[z>3]=-Inf
  }
  nonrep2=nonrep[!is.infinite(nonrep)]+1
  if(verbosity>=3) print("....done")

  if(verbosity>=3) print("Counting motifs....")
  res=1:(4^len)*0
  for(i in 1:length(nonrep2)){
    if(!i%%100000) if(verbosity>=3) print(i)
    res[nonrep2[i]]=res[nonrep2[i]]+1
  }
  bases=c("A","C","G","T")
  ournames=""
  for(i in 1:len) ournames=c(paste(bases[1],ournames,sep=""),paste(bases[2],ournames,sep=""),paste(bases[3],ournames,sep=""),paste(bases[4],ournames,sep=""))
  names(res)=ournames
  ournamesc=""
  for(i in 1:len) ournamesc=c(paste(ournamesc,bases[4],sep=""),paste(ournamesc,bases[3],sep=""),paste(ournamesc,bases[2],sep=""),paste(ournamesc,bases[1],sep=""))
  restot=res+res[ournamesc]
  if(verbosity>=3) print("....done")

  if(verbosity>=3) print("Finding potential starts....")
  lookups=order(-restot)[1:200]
  restot=restot[lookups]
  seqs=ournames[lookups]
  seqsc=ournamesc[lookups]
  seqs=seqs[restot>10]
  seqsc=seqsc[restot>10]
  if(!length(seqs)){
    if(verbosity>=3) print("No motif to start from is in more than 10 sequences")
    return(0)
  }
  pos=1
  while(pos<length(seqs)){
    check1=seqs[pos]
    check2=seqsc[pos]
    if(check1!=check2) {
      seqs=seqs[seqs!=check2]
      seqsc=seqsc[seqsc!=check1]
    }
    pos=pos+1
  }
  if(verbosity>=3) print("....done")

  if(verbosity>=3) print("Checking for central enrichment....")
  mids=nchar(regs)/2
  reg2=regs[nchar(regs)>=mids+range & mids-range>=1]
  midregs=substring(reg2,mids-range,mids+range)

  enrich=1:length(seqs)
  excess=enrich
  for(i in 1:length(seqs)){
    if(verbosity>=3) print(i)
    mot=seqs[i]
    motc=seqsc[i]
    ourset=unique(c(grep(mot, midregs, fixed = T),grep(motc, midregs, fixed = T)))
    ourset2=unique(c(grep(mot, reg2, fixed = T),grep(motc, reg2, fixed = T)))
    t2=sum(nchar(reg2)-len+1)
    t1=sum(nchar(midregs)-len+1)
    enrich[i]=length(ourset)/t1/(length(ourset2)-length(ourset))*(t2-t1)
    excess[i]=length(ourset)-(length(ourset2)-length(ourset))*t1/(t2-t1)
  }

  #####now take most enriched motif in top 50, and look at single-base alterations and extensions for their impact
  #####use pwm for motif refinement
  #####automatically obtain sequences without the motif


  excess=excess[enrich>1]
  seqs=seqs[enrich>1]
  seqsc=seqsc[enrich>1]
  enrich=enrich[enrich>1]
  if(!length(enrich)) {
    if(verbosity>=3) print("No motif to start from is centrally enriched")
    return(0)
  }

  #return(seqs)

  if(verbosity>=1) print("Top 5 start motifs:")
  if(verbosity>=1) print(seqs[order(-excess)[1:5]])

  mot = seqs[order(-excess)][!seqs[order(-excess)] %in% motif_blacklist][motif_rank]

  }else if(motif_seed=="random"){

    mot = paste0(sample(c("A","C","T","G"), T, size = len), collapse = "")

  }else if(motif_seed=="modal"){  # use modal seeding algo instead

    seeds <- do.call(paste0,expand.grid(rep(list(c("A","C","T","G")),len), stringsAsFactors = F))

    seeds_with_revcomp <- unique(sapply(seeds, function(x) paste0(sort(c(stri_reverse(chartr("ACGT", "TGCA", x)),x)), collapse = "|")))

    seed_count <- sapply(seeds_with_revcomp, function(x) sum(grepl(x,regs)))

    #seed_count <- str_count(seqs_concat,seeds_with_revcomp)
    #seed_count <- str_count(seqs_concat,seeds)
    #names(seed_count) <- seeds#_with_revcomp

    if(verbosity>=1) print("Top 5 start motifs:")
    print(head(sort(seed_count,T),5))
    mot = unlist(strsplit(names(sort(seed_count,T)),"|", fixed = T))[c(T,F)]
    mot = mot[!mot %in% motif_blacklist][motif_rank]
  }else{
    mot = motif_seed
  }

  if(verbosity>=1) print("Chose start motif:")
  if(verbosity>=1) print(mot)

  chosen_seed = mot

  if(verbosity>=3) print("Initialising....")
  mot=as.vector(unlist(strsplit(mot,"")))
  mot[mot=="A"]=1
  mot[mot=="C"]=2
  mot[mot=="G"]=3
  mot[mot=="T"]=4
  mot=as.double(as.vector(unlist(mot)))

  pwmstart=matrix(0.1,nrow=10,ncol=4)
  for(i in 1:len) pwmstart[i,mot[i]]=1
  pwmstart=pwmstart/rowSums(pwmstart)
  logpwm=log(pwmstart)

  seqtemp=regs[order(-scores)][1:n_for_refine]
  if(verbosity>=3) print("....done")

  if(verbosity>=3) print("Attempting to refine....")
  z=getmotifs(logpwm,length(logpwm[,1]),seqtemp,maxwidth=max(nchar(seqtemp)),alpha=0.5,incprob=0.99999,maxits=nits,plen=plen,updatemot=1,updatealpha=1,ourprior=prior,bg=-1,updateprior=updateprior,seed=NA,verbosity=verbosity,stranded_prior=stranded_prior, conv_n=conv_n, conv_t=conv_t)
  if(verbosity>=3) print("....done")

  if(is.null(z)){
    return(NULL)
  }

  if(verbosity>=1) print("Scoring regions (fixed motif) ....")
  z2=getmotifs(z$scoremat,z$scorematdim,regs,maxwidth=max(nchar(regs)),alpha=z$alpha,incprob=0.99999,maxits=scoring_its,plen=0.2,updatemot=0,updatealpha=1,ourprior=z$prior,updateprior=1,bg=-1,seed=NA,verbosity=verbosity,stranded_prior=stranded_prior)
  if(verbosity>=3) print("....done")

  z2$alphas <- z$alphas

  z2$seed <- as.integer(seed)

  z2$motif_seed <- chosen_seed

  return(z2)

}
