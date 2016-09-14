ripple_all<-function(input.seq,ws=4,LOD=3,tol=10E-2, phasing.cores = 4,
                     ripple.cores = 4, start = 1) {
  ## checking for correct objects
  if(!any(class(input.seq)=="sequence")) {
    stop(deparse(substitute(input.seq)),
         " is not an object of class 'sequence'")
  }
  if(ws < 2) stop("ws must be greater than or equal to 2")
  if(ws > 5) warning("WARNING: this operation may take a VERY long time\n\n")
  len <- length(input.seq$seq.num)
  ## computations unnecessary in this case
  if (len <= ws) stop("Length of sequence ",
                      deparse(substitute(input.seq)),
                      " is smaller than ws. You can use the ",
                      "compare function instead")

  ## allocate variables
  rf.init <- rep(NA,len-1)
  phase <- rep(NA,len-1)
  tot <- prod(1:ws)
  best.ord.phase <- matrix(NA,tot,len-1)
  best.ord.like <- best.ord.LOD <- rep(-Inf,tot)
  all.data <- list()

  ## gather two-point information
  list.init <- phases(input.seq)

  #### first position
  p <- start
  message("...", input.seq$seq.num[p-1], "|",
          paste(input.seq$seq.num[p:(p+ws-1)], collapse = "-"),"|",
          input.seq$seq.num[p+ws],"...")
  all.ord <- t(apply(perm.tot(input.seq$seq.num[p:(p+ws-1)]),1,function(x){
    return(c(head(input.seq$seq.num,p-1),
             x,tail(input.seq$seq.num,-p-ws+1)))
    }))

  poss <- mclapply(1:nrow(all.ord), mc.allow.recursive = TRUE,
                   mc.cores = ripple.cores, function(i){
                     message("Trying order ",i," of ",nrow(all.ord),
                             " for star position ",p)
                     mp <- list(seq.like = -Inf)
                     tryCatch({
                       mp <- map(make.seq(get(input.seq$twopt), all.ord[i,],
                                          twopt = input.seq$twopt),
                                 phase.cores = phasing.cores)
                     }, error = function(e){},
                     finally = {
                       return(mp)
                     })
                   })
  best <- which.max(sapply(poss,"[[","seq.like"))

  return(poss[[best]])
}

ripple_rand<-function(input.seq,ws=4,LOD=3,tol=10E-2, phasing.cores = 4,
                     ripple.cores = 4, start = 1, n = NULL, pref = "neutral") {
  ## checking for correct objects
  if(!any(class(input.seq)=="sequence")) {
    stop(deparse(substitute(input.seq)),
         " is not an object of class 'sequence'")
  }
  if(ws < 2) stop("ws must be greater than or equal to 2")
  if(ws > 5) warning("WARNING: this operation may take a VERY long time\n\n")
  if(is.null(input.seq$seq.like)) stop("You need to run map() at least once ",
                                       "before attempting to ripple.")
  if(is.null(n))
  {
    n <- prod(ws:1)/2
  }
  len <- length(input.seq$seq.num)
  ## computations unnecessary in this case
  if (len <= ws) stop("Length of sequence ",
                      deparse(substitute(input.seq)),
                      " is smaller than ws. You can use the ",
                      "compare function instead")

  ## allocate variables
  rf.init <- rep(NA,len-1)
  phase <- rep(NA,len-1)
  tot <- prod(1:ws)
  best.ord.phase <- matrix(NA,tot,len-1)
  best.ord.like <- best.ord.LOD <- rep(-Inf,tot)
  all.data <- list()

  ## gather two-point information
  list.init <- phases(input.seq)

  #### first position
  p <- start
  message("...", input.seq$seq.num[p-1], "|",
          paste(input.seq$seq.num[p:(p+ws-1)], collapse = "-"),"|",
          input.seq$seq.num[p+ws],"...")

  ref <- input.seq$seq.num[p:(p+ws-1)]

  probs <- drop(cor(ref,t(perm.tot(input.seq$seq.num[p:(p+ws-1)]))))

  probs <- (probs - min(probs)) / (max(probs) - min(probs))
  probs[probs == 1] <- 0

  all.ord <- t(apply(perm.tot(input.seq$seq.num[p:(p+ws-1)]),1,function(x){
    return(c(head(input.seq$seq.num,p-1),
             x,tail(input.seq$seq.num,-p-ws+1)))
  }))
  if(prefs == "similar")
  {
    all.ord <- all.ord[sample(1:nrow(all.ord),n,FALSE,probs),]
  }
  if(prefs == "dissimilar")
  {
    all.ord <- all.ord[sample(1:nrow(all.ord),n,FALSE,-probs),]
  }
  if(prefs == "neutral")
  {
    all.ord <- all.ord[sample(1:nrow(all.ord),n,FALSE),]
  }
  poss <- mclapply(1:nrow(all.ord), mc.allow.recursive = TRUE,
                   mc.cores = ripple.cores, function(i){
                     message("Trying order ",i," of ",nrow(all.ord),
                             " for star position ",p)
                     mp <- list(seq.like = -Inf)
                     tryCatch({
                       mp <- map(make.seq(get(input.seq$twopt), all.ord[i,],
                                          twopt = input.seq$twopt),
                                 phase.cores = phasing.cores)
                     }, error = function(e){},
                     finally = {
                       return(mp)
                     })
                   })
  best <- which.max(sapply(poss,"[[","seq.like"))
  if(poss[[best]]$seq.like > input.seq$seq.like)
  {
    return(poss[[best]])
  } else {
    return(input.seq)
  }
}

ripple_one <- function(input.seq,ws=4,LOD=3,tol=10E-2, phasing.cores = 4,
                       ripple.cores = 4, start = 1)
{
  if(!any(class(input.seq)=="sequence")) {
    stop(deparse(substitute(input.seq)),
         " is not an object of class 'sequence'")
  }
  if(ws < 2) stop("ws must be greater than or equal to 2")
  if(ws > 5) warning("WARNING: this operation may take a VERY long time\n\n")
  if(is.null(input.seq$seq.like)) stop("You need to run map() at least once ",
                                       "before attempting to ripple.")

  len <- length(input.seq$seq.num)
  ## computations unnecessary in this case
  if (len <= ws) stop("Length of sequence ",
                      deparse(substitute(input.seq)),
                      " is smaller than ws. You can use the ",
                      "compare function instead")

  ## allocate variables
  rf.init <- rep(NA,len-1)
  phase <- rep(NA,len-1)
  tot <- prod(1:ws)
  best.ord.phase <- matrix(NA,tot,len-1)
  best.ord.like <- best.ord.LOD <- rep(-Inf,tot)
  all.data <- list()

  ## gather two-point information
  list.init <- phases(input.seq)

  #### first position
  p <- start
  message("...", input.seq$seq.num[p-1], "|",
          paste(input.seq$seq.num[p:(p+ws-1)], collapse = "-"),"|",
          input.seq$seq.num[p+ws],"...")

  all.ord <- matrix(NA,(sum((ws-1):1) + 1) * 2,ws)
  all.ord[1,] <- input.seq$seq.num[p:(p+ws-1)]

  r <- 2
  for(i in 1:(ws - 1))
  {
    for(j in ws:(i+1))
    {
      all.ord[r,] <- all.ord[1,]
      tmp <- all.ord[r,i]
      all.ord[r,i] <- all.ord[r,j]
      all.ord[r,j] <- tmp
      r <- r+1
    }
  }
  for(r in (sum((ws-1):1) + 2):((sum((ws-1):1) + 1) * 2))
  {
    all.ord[r,] <- rev(all.ord[r - sum((ws-1):1) + 1,])
  }
  all.ord <- all.ord[!duplicated(apply(all.ord,1,paste,collapse="-")),]

  all.ord <- t(apply(all.ord,1,function(x){
    return(c(head(input.seq$seq.num,p-1),
             x,tail(input.seq$seq.num,-p-ws+1)))
  }))

  poss <- mclapply(1:nrow(all.ord), mc.allow.recursive = TRUE,
                   mc.cores = ripple.cores, function(i){
                     message("Trying order ",i," of ",nrow(all.ord),
                             " for star position ",p)
                     mp <- list(seq.like = -Inf)
                     tryCatch({
                       mp <- map(make.seq(get(input.seq$twopt), all.ord[i,],
                                          twopt = input.seq$twopt),
                                 phase.cores = phasing.cores)
                     }, error = function(e){},
                     finally = {
                       return(mp)
                     })
                   })
  best <- which.max(sapply(poss,"[[","seq.like"))
  if(poss[[best]]$seq.like > input.seq$seq.like)
  {
    return(poss[[best]])
  } else {
    return(input.seq)
  }
}

ripple_ord <- function(input.seq,ws=4,LOD=3,tol=10E-2, phasing.cores = 4,
                       ripple.cores = 4, method = "one", n = NULL,
                       pref = "neutral", start = 1)
{
  LG <- input.seq
  if(start + ws > length(input.seq$seq.num)) return(LG)
  if(method == "all"){
    for(i in start:(length(input.seq$seq.num) - ws))
    {
      LG <- ripple_all(LG,ws,LOD,tol,phasing.cores,ripple.cores,i)
    }
  }
  if(method == "random")
  {
    for(i in start:(length(input.seq$seq.num) - ws))
    {
      LG <- ripple_rand(LG,ws,LOD,tol,phasing.cores,ripple.cores,i,n,pref)
    }
  }
  if(method == "one")
  {
    for(i in start:(length(input.seq$seq.num) - ws))
    {
      LG <- ripple_one(LG,ws,LOD,tol,phasing.cores,ripple.cores,i)
    }
  }
  return(LG)
}

