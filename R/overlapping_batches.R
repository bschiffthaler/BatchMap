generate_overlapping_batches <- function(input.seq, size = 50, overlap = 15,
                                         silent = FALSE)
{
  start <- 1
  end <- size
  current <- 1
  res <- list()
  while(end <= length(input.seq$seq.num))
  {
    res[[current]] <- input.seq$seq.num[start:end]
    current <- current + 1
    start <- end - overlap
    if(end == length(input.seq$seq.num)) break
    end <- end + size - overlap - 1
    if(end > length(input.seq$seq.num)) end <- length(input.seq$seq.num)
  }
  sizes <- unlist(lapply(res, length))
  if(length(sizes) < 2 & ! silent)
  {
    warning("You should at least have two overlapping batches.",
            " Reconsider the size parameter.")
  }
  if(any(sizes/size > 1.25) & ! silent)
  {
    warning("One group is 25% bigger than the group size. ",
            "Consider adjusting parameters.")
  }
  return(res)
}

pick_batch_sizes <- function(input.seq, size = 50, overlap = 15, around = 5)
{
  test.sizes <- c(size, (size - 1):(size - around), (size + 1):(size + around))
  all.batches <- lapply(test.sizes, function(s){
    generate_overlapping_batches(input.seq, s, overlap, silent = TRUE)
  })
  x <- unlist(lapply(all.batches, function(f){
    ran <- range(unlist(lapply(f,length)))
    ran[2] - ran[1]
  }))
  x <- which(x == min(x))
  test.sizes[x[length(x)]] #prefer larger maps
}

map_overlapping_batches <- function(input.seq, size = 50, overlap = 10,
                        fun.order = NULL, phase.cores = 4,
                        ripple.cores = 1, verbosity = NULL, max.dist = Inf,
                        ws = 4, increase.every = 4, max.tries = 10,
                        min.tries = 0, ...)
{
  batches <- generate_overlapping_batches(input.seq, size, overlap)
  if("batch" %in% verbosity)
  {
    message("Have ", length(batches), " batches.")
    message("The number of markers in the final batch is: ",
            length(batches[[length(batches)]]))
    message("Prcoessing batch 1...")
  }
  LGs <- list()
  LG <- map(make.seq(get(input.seq$twopt), batches[[1]],
                     twopt = input.seq$twopt), phase.cores = phase.cores,
            verbosity = verbosity)
  round <- 1
  increment <- 0
  while(any(kosambi(LG$seq.rf) > max.dist) | round <= min.tries)
  {
    if(round > max.tries)
    {
      warning("Algorithm could not solve gaps in batch 1.")
      break
    }
    print(LG)
    if(round %% increase.every == 0) increment <- increment + 1
    LG <-  fun.order(LG, ripple.cores = ripple.cores, start = 1,
                     verbosity = verbosity, batches = batches,
                     ws = ws + increment, ...)
    round <- round + 1
  }
  LGs[[1]] <- LG
  for(i in 2:length(batches))
  {
    if("batch" %in% verbosity)
    {
      print(LGs[[i - 1]])
      message("Processing batch ",i,"...")
    }
    seeds <- tail(LGs[[i - 1]]$seq.phases, overlap)
    batches[[i]][1:(overlap+1)] <- tail(LGs[[i - 1]]$seq.num, overlap + 1)
    LG <- seeded.map(make.seq(get(input.seq$twopt),
                              batches[[i]],
                              twopt = input.seq$twopt),
                     verbosity = verbosity,
                     seeds = seeds)
    if(! is.null(fun.order ))
    {
      round <- 1
      increment <- 0
      while(any(kosambi(LG$seq.rf) > max.dist) | round <= min.tries )
      {
        if(round > max.tries)
        {
          warning("Algorithm could not solve gaps in batch ", i)
          break
        }
        if(round %% increase.every == 0) increment <- increment + 1
        LG <-  fun.order(LG, ripple.cores = ripple.cores, start=overlap+2,
                         verbosity = verbosity, batches = batches,
                         ws = ws + increment, ...)
        round <- round + 1
      }
    }
    LGs[[i]] <- LG
  }
  final.seq <- LGs[[1]]$seq.num
  final.phase <- LGs[[1]]$seq.phases
  for(i in 2:length(batches))
  {
    start <- length(final.seq) - overlap
    final.seq[start:length(final.seq)] <- head(LGs[[i]]$seq.num, overlap + 1)
    final.seq <- c(final.seq,
                   LGs[[i]]$seq.num[(overlap + 2):length(LGs[[i]]$seq.num)])
    start <- length(final.phase) - overlap + 1
    final.phase[start:length(final.phase)] <- head(LGs[[i]]$seq.phases, overlap)
    final.phase <- c(final.phase,
                   LGs[[i]]$seq.phases[(overlap + 1):length(LGs[[i]]$seq.phases)])
  }
  if("batch" %in% verbosity)
  {
    message("Final call to map...")
  }

  mp <- map(make.seq(get(input.seq$twopt), final.seq, final.phase, input.seq$twopt),
            verbosity = verbosity)
  return(list(Map = mp, LGs = LGs))
}
