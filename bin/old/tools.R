options("stringsAsFactors" = F)
options("scipen" = 20)
require(stringi)
require(data.table)
require(Biostrings)
require(matrixStats)
# require(GenomicRanges)
require(parallel)
print("note: strings will not be converted as factors")

bigObjects = function(cutoff = 2^20) {
  obj = sort(sapply(ls(globalenv()), function(x){object.size(get(x))}), decreasing = T)
  obj[obj > cutoff]/2^20
}

roundTo = function(x, multiple = 1) {
    round(x/multiple)*multiple
}


ceilingTo = function(x, multiple = 1) {
  ceiling(x/multiple)*multiple
}


floorTo = function(x, multiple = 1) {
  floor(x/multiple)*multiple
}


rescale <- function(x, range) {
  ampl <- max(range) - min(range)
  amplx <- diff(range(x, na.rm= T))
  xb <- x * ampl / amplx
  xb - min(xb, na.rm = T) + min(range)
}

mids <- function(v) {
  i <- 2:length(v)
  (v[i - 1L] + v[i]) / 2L
}

ranks <- function(...) {
  ord <- order(...)
  ord[ord] <- 1:length(ord)
  ord
}

prodNA <- function(x) {
  p <- prod(x, na.rm = T)
  if (all(is.na(x))) p <- as.numeric(NA)
  p
}

# T for numbers that do not respect the ascendign order in a vector
respectsOrder = function(x, decreasing = F) {
  x = x*2 + rep(0:1, length.out = length(x))
  m = match(x, sort(x, decreasing = decreasing))
  c(T, diff(m) < 0)
}


# makes a histogram, but avoids biases due to edge effects
hist2 <- function(v, ..., returnObject = F) {
  h <- hist(v, plot = F)
  newbreaks <- rescale(h$breaks, range(v))
  h <- (hist(v, breaks = newbreaks, ...))
  if (returnObject) invisible(h) else invisible(NULL)
}

pairMatrix <- function(rows, cols, values, default = NA, symetric = T) {
  nams <- unique(c(rows, cols))
  mat <- matrix(default, ncol = length(nams), nrow = length(nams))
  colnames(mat) <- nams
  rownames(mat) <- nams
  mat[cbind(rows, cols)] <- values
  if (symetric) mat[cbind(cols, rows)] <- values
  mat
}


# takes a data.table with a least to columns and makes it a named vector. names
# and values represent column numbers
namedVector <- function(dt, names = 1, values = 2) {
  setnames(dt, c(names, values), c("names", "values"))
  values <- dt$values
  names(values) <- dt$names
  values
}

toOpposite <- function(f) {
  x <- 1L * f
  x[x == 0] <- -1L
  x
}


# rearranges a list whose names are integers so that the names of elements of
# the list are placed at equivalent indices of a new list.
reList <- function(spl, len = NA) {
  indices <- as.integer(names(spl))
  if (any(is.na(indices))) {
    return(spl)
    stop("some names are not convertible to integers")
  }
  m <- max(indices)
  if (!is.na(len) & len >= m) m <- as.integer(len)
  temp <- vector(mode = "list", m)
  temp[indices] <- spl
  temp
}


# cumulated sum that resets itself each time a new level of "fact" starts. No
# level of "fact" must be intertwined with others
cumsums <- function(v, fact) {
  cs <- cumsum(v)
  last <- which(diff(fact) != 0L)
  last <- last[-length(last)]
  diffs <- diff(c(0L, cs[last]))
  v[last + 1L] <- v[last + 1L] - diffs
  cumsum(v)
}


# returns a cumulated sum but not started at the first element of a vector
cumsumStart <- function(vect, start = 0) {
  cs <- cumsum(vect)
  cs <- c(start, cum[-length(cum)] + start)
  setNames(cs, names(vect))
}


fromClosest <- function(n, start, end) {
  f <- end < start
  temp <- start
  start[f] <- end[f]
  end[f] <- temp[f]
  pmin(n - start, end - n)
}


# paste columns of a matrix or dataframe together
pasteCols <- function(mat, sep = "") {
  mat <- as.data.frame(mat)
  do.call(stri_c, c(mat, sep = sep))
}


lastChars <- function(string, n) {
  stri_sub(string, nchar(string) - n + 1L, length = n)
}


substrFrom <- function(string, fromStart, fromEnd = 0L) {
  stri_sub(string, fromStart, stri_length(string) - fromEnd)
}


# turns a string vector into integers
toInteger <- function(string, sorted = F, unique = T) {
  if (unique) u <- unique(string) else u <- string
  if (sorted) u <- sort(u)
  if(is.numeric(string)) {
    return(match(string, u))
  }
  return(chmatch(string, u))
}


tableMatch <- function(dt1, dt2, nomatch = as.integer(NA)) {
  if (!is.data.table(dt1)) dt1 <- as.data.table(dt1)
  if (!is.data.table(dt2)) dt2 <- as.data.table(dt2)
  nam <- stri_c("V", 1:ncol(dt1))
  setnames(dt1, nam)
  setnames(dt2, nam)
  dt2[, c("N1", "N2") := .(0L, 1:.N)]
  dt1[, c("N1", "N2") := .(1:.N, as.integer(nomatch))]
  dt <- rbind(dt2, dt1)
  dt <- dt[, .(N1, N2 = N2[1]), by = nam]
  m <- integer(nrow(dt1))
  f <- dt$N1 > 0L
  m[dt[f, N1]] <- dt[f, N2]
  m
}


equalClasses <- function(x, n) {
  f <- .bincode(x, quantile(x, seq(0, 1, length.out = n + 1), na.rm = T), include.lowest = T)
  f
}


splitEqual <- function(x, f = NA, n, ...) {
  l <- length(x)
  if (is.data.frame(x)) {
    l <- nrow(x)
  }
  if (all(is.na(f))) {
    f <- cut(1:l, n)
  } else {
    f <- equalClasses(f, n)
  }
  split(x, f, ...)
}

toTable <- function(var1, var2, values, missing = NA, makeSymetric = T) {
  nams <- unique(c(var1, var2))
  mat <- matrix(missing, ncol = length(nams), nrow = length(nams))
  rownames(mat) <- nams
  colnames(mat) <- nams
  mat[cbind(var1, var2)] <- values
  if (makeSymetric) mat[cbind(var2, var1)] <- values
  mat
}

freadLines <- function(file, ...) {
  fread(file, header = F, sep = "\n", ...)$V1
}

xyDistance <- function(x0, y0, x1, y1) {
  dX <- x1 - x0
  dY <- y1 - y0
  l <- sqrt(dX^2 + dY^2)
}


# converts contig positions into scaffold positions. The source dataframe
# specifies the contig and the position (columns 1 and 2) and the description
# describes the scaffold, the contig and its start position (in this order). If
# some contigs names are missing in scaffolds that have just one contig, contigs
# must be described by a prefix and a number, with no separator between
contigToScaffoldPosition <- function(source, description, sep = "_", prefix = "contig") {
  scaffCont <- testSplit(source[, 1], sep, 1:2)
  miss <- scaffCont[, 2] == ""
  if (any(miss)) {
    contig1 <- stri_c(prefix, "1")
    source[miss, 1] <- stri_c(source[miss, 1], contig1, sep = sep)
  }

  # so that contigs in the description have unique names
  pasted <- stri_c(description[, 1], description[, 2], sep = sep)
  indices <- chmatch(source[, 1], pasted)
  pos <- source[, 2] + description[indices, 3] - 1L
  data.table(scaffold = scaffCont[, 1], pos)
}

GRangesFromBed <- function(bed) {
  bed <- fixRanges(bed, 2)
  n <- colnames(bed)
  makeGRangesFromDataFrame(bed, ignore.strand = T, seqnames.field = n[1], start.field = n[2], end.field = n[3])
}



mostFrequent <- function(vect) {
  if (is.integer(vect)) {
    corr = max(1L - min(vect), 0L)
    m <- which.max(tabulate(vect+corr)) - corr
  } else {
    tab = table(vect)
    m <- names(tab)[which.max(tab)]
  }
  m
}



sequenceSizes <- function(fasta, Ns = T) {
  sequences <- readDNAStringSet(fasta)
  if (Ns) counts <- stri_length(sequences) else counts <- stri_count_regex(sequences, "[^N|^n]")
  print(fasta)
  setNames(counts, names(sequences))
}


# for each position of x, gives the number of times the corresponding element is
# found in x. x can be a table, in which case elements correspond to rows
nTimes <- function(x) {
  dt <- as.data.table(x)
  dt[, pos_ := 1:.N]
  counts <- dt[, .(.N, pos_), by = setdiff(colnames(dt), "pos_")]
  res <- integer(nrow(dt))
  res[counts$pos_] <- counts$N
  res
}


# returns the occurrence of each element of a vector or table (e.g. 1 if it's the
# first time it appears, 2 if the 2nd time, etc). If vect is a table, it gets
# the occurrence of each combination of elements at each row
occurences <- function(table) {
  dt <- data.table(table)
  dt[, pos_ := 1:.N]
  counts <- dt[, .(occ = 1:.N, pos_), by = setdiff(colnames(dt), "pos_")]
  res <- integer(nrow(dt))
  res[counts$pos_] <- counts$occ
  res
}



# for all positions in a vector, tells the number of different elements from the
# beginning up to that position. Works also for data tables
nUnique <- function(vect) {
  cumsum(!duplicated(vect))
}

replaceInString <- function(string, start, end, replacement) {
  debs <- stri_sub(string, 1, start - 1)
  ends <- stri_sub(string, end + 1, nchar(string))
  string[] <- stri_join(debs, replacement, ends)
  string
}


# fixes ranges so that start positions are never higher than end positions.
# Ranges is a dataframe or matrix with 2 numeric columns, "starts" is the number
# of the column of start positions
fixRanges <- function(ranges, starts = 1, ends = starts + 1) {
  f <- ranges[, starts] > ranges[, ends]
  ranges[f, c(starts, ends)] <- ranges[f, c(ends, starts)]
  ranges
}



overlappingRegionsGroup <- function(bed) {
  setnames(bed, 1:3, c("scaffold","start","end"))
  depths = genomeCov(bed)
  region = assignToRegion(depths, bed[,.(scaffold, start)])
  depths[region, cov]
}


# assigns positions (like SNPs) to genome regions defined in a data frame in
# bed-style format (scaffold, start, end) and pos is a data frame specifying
# scaffold and position. If a position is in the area of overlap between 2
# regions, it will be assigned to the second one
assignToRegion <- function(bed, pos, olv = T) {
  bed <- as.data.table(bed)
  pos <- as.data.table(pos)
  setnames(bed, 1:3, c("scaffold", "start", "end"))
  setnames(pos, 1:2, c("scaffold", "pos"))

  # we will identify regions by their rows in the original bed
  bed[, id := 1:.N]
  bed <- bed[bed$scaffold %chin% pos$scaffold, ]

  # we convert coordinates into global genome coordinates to call .bincode. This
  # requires estimating the length of scaffold, which is the max coordinate per
  # scaffold for the bed and pos tables combined. "+1" gives a bit or room
  # (probably not needed).
  scaffLengths <- rbind(bed[, .(scaffold, pos = end)], pos)[, max(pos) + 1L, by = scaffold]
  starts <- scaffToGenomeCoord(bed$scaffold, bed$start, scaffLengths$V1, scaffLengths$scaffold)
  pos <- scaffToGenomeCoord(pos$scaffold, pos$pos, scaffLengths$V1, scaffLengths$scaffold)

  # creates a data table with start and ends in the same "boundary" column.
  # These will be intervals for binning. We use "id" to keep track of the rows
  # (regions id) of the original bed
  dt <- data.table(boundary = c(starts, starts + bed[, end - start + 1L]), id = rep(bed$id, 2))
  setorder(dt, boundary)

  # "region" will be the region (id) covering an interval (only for the start of
  # this interval). The idea is that when a region ends, the interval should get
  # the region that main contain the one that ended or be NA (if between 2
  # regions)
  dt[, c("row", "region") := .(1:.N, as.integer(NA))]
  rows <- dt$row[-1]

  # first rows regions that do not overlap with others (whose ids appear
  # successively in the table)
  short <- dt[, which(id[rows - 1L] == id[rows])]

  # all the rows of these regions (first and last)
  shorts <- c(short, short + 1L)

  # if there are overlapping regions
  if (length(shorts) < nrow(dt)) {
    if (length(shorts) == 0) {
      rows <- dt[, .(row = row[1L]:(row[.N] - 1L), size = row[.N] - row[1L]), by = id]

      # for each of these overlapping regions, we create a vector that represent
      # all the rows of dt that are covered by that region
    } else {
      rows <- dt[-shorts, .(row = row[1L]:(row[.N] - 1L), size = row[.N] - row[1L]), by = id]
    }

    # puts the longest regions on top
    setorder(rows, -size)

    # so we can add a region id for each row in dt. If there are duplicates in
    # rows$row, ids of shorter regions will overwrite longer ones. This takes
    # advantage from the fact that r operates in vector order
    dt[rows$row, region := rows$id]
  }

  # we now copy the ids for short regions.
  dt[short, region := id]
  dt[.bincode(pos, boundary, F, F), region]
}

assignToRegion2 = function(scaffold, start, end, scaff2, pos) {

  id = which(scaffold %chin% scaff2)
  dt = data.table(
    scaffold = c(scaffold[id], scaff2, scaffold[id]),
    posi = c(start[id], pos, end[id]),
    bedID = c(id, rep(NA, length(pos)), id),
    posID = c(rep(NA, length(id)), 1:length(pos), rep(NA, length(id)))
    )
  setorder(dt, scaffold, posi)
  
  bed = data.table(scaffold = scaffold[id], start = start[id], end = end[id], id)
  bed[,c("first", "last") := .(match(id, dt$bedID), matchLast(id, dt$bedID))]
  noPos = dt[!is.na(bedID), bedID]
  subdt = bed[last - first > matchLast(id, noPos) - match(id, noPos)]
  
  if(nrow(subdt) >0L) {
    subdt[,l := end-start]
    setorder(subdt, -l)
    rows = subdt[, .(r = first:last), by = id]
    
    dt[rows$r, bedID := rows$id]
  }
  dt[match(1:length(pos), posID), bedID]
}


# generates IRanges from a bed file. The scaffold-based positions are converted
# to absolute positions, and the results are splitted so that positions do not
# exceed 2^31, which is the limit of IRanges
IRangesFromBed <- function(bed) {
  if (!is.data.table(bed)) bed <- data.table(bed)
  setnames(bed, c("scaffold", "start", "end"))

  # we convert coordinates into global genome coordinates so that we can call
  # rangeCov once. "+1" is need to avoid ranges overlapping between scaffold,
  # when some sart at position 1
  scaffLengths <- bed[, max(end) + 1, by = scaffold]
  starts <- scaffToGenomeCoord(bed$scaffold, bed$start, scaffLengths$V1, scaffLengths$scaffold)
  ends <- starts + bed$end - bed$start
}



# splits a string vector into a character vector. May not be compatible with all
# string types. The cumulated length of the string should not exceed the maximum
# length of a vector (2^31-1, typically)
splitToChar <- function(string, charTypes) {
  if (length(string) > 1L) string <- stri_flatten(string)
  raw <- charToRaw(string)
  if (missing(charTypes)) {
    u <- unique(raw)
    charTypes <- unlist(strsplit(rawToChar(u), ""))
  } else {
    u <- charToRaw(charTypes)
    charTypes <- unlist(strsplit(charTypes, ""))
  }
  types <- as.integer(u)
  conv <- rep("", max(types))
  conv[types] <- charTypes
  conv[as.integer(raw)]
}


splitToBases <- function(seqs) {
  splitToChar(seqs, charTypes = "acgtmrwsykvhdbACGTMRWSYKVHDBNn-")
}


# reads a fastq as a data.table that has columns for read name, sequences and quality. Sequence and quality must occupy a single line per read in the fastq file.
fastqToDT <- function(file, desc = F, ...) {
  fastq <- freadLines(file, ...)
  readInd <- seq(1, length(fastq), 4)
  readName <- NULL
  if (desc) {
    readName <- gsub("@", "", fastq[readInd], fixed = T)
  } else {
    readName <- gsub("@", "", testSplit(fastq[readInd], " ", 1), fixed = T)
  }
  data.table(readName, seq = fastq[readInd + 1], qual = fastq[readInd + 3])
}


writeFastq <- function(names, seqs, scores, file, fix = F) {
  fastq <- as.vector(rbind(stri_c("@", names), seqs, "+", scores))
  writeLines(fastq, con = file)
  if (fix) {
    system(stri_c("grep -v -e '^$' ", file, " > ", file, "noBlanckLines"))
    file.rename(stri_c(file, "noBlanckLines"), file)
  }
  invisible("")
}


squareMatrixFromVector <- function(vect, names = "") {
  x <- (1 + sqrt(1 + 8 * length(vect))) / 2
  if (x %% 1 != 0) {
    print("vector is not of adequate length")
  } else {
    mat <- matrix(0, ncol = x, nrow = x)
    if (names != "" & length(names) == x) {
      colnames(mat) <- names
      rownames(mat) <- names
    }
    pos <- 1
    for (i in 1:(x - 1)) {
      mat[(i + 1):x, i] <- vect[pos:(pos + x - i - 1)]
      mat[i, (i + 1):x] <- vect[pos:(pos + x - i - 1)]
      pos <- pos + (x - i)
    }
    mat
  }
}

intersection <- function(start, end, start2, end2, range = F, negative = T, fix = F, olv = 1L) {
  minEnd <- pmin(pmax(start, end), pmax(start2, end2))
  maxStart <- pmax(pmin(start, end), pmin(start2, end2))
  if (range) {
    r <- cbind(maxStart, minEnd)
    r[maxStart > minEnd, ] <- NA
    return(r)
  }
  inter <- minEnd - maxStart + olv
  if (!negative) inter[inter < 0] <- 0L
  inter
}


# like match, but returns the last position of each element of x in vect
matchLast <- function(x, table, ...) {
  if (is.character(x)) {
    return(length(table) + 1L - chmatch(x, rev(table), ...))
  } else {
    return(length(table) + 1L - match(x, rev(table), ...))
  }
}




# gives the distribution of numeric variable "var" according to the numeric
# variable "splitter", with a given number of classes. If regular is true, the
# classes have similar number of observations
distrib <- function(x, splitter, splits, FUN, regular = TRUE) {
  quant <- quantile(splitter, seq(0, 1, 1 / splits), na.rm = T)
  if (!regular) quant <- c(seq(min(splitter), max(splitter), (max(splitter) - mean(splitter)) / (splits - 2)), max(splitter))
  ids <- .bincode(splitter, sort(quant), F)
  res <- tapply(x, ids, FUN = FUN, na.rm = T, simplify = T)
  means <- tapply(splitter, ids, FUN = mean, na.rm = T, simplify = T)
  data.table(mean = as.vector(means), var = as.vector(res), start = quant[-length(quant)], end = quant[-1])
}



# returns the number of positions (an integer) for which each range of "ranges"
# intersect any range of ranges2. Ranges are data.frames in bed-like format with
# closed intervals
intersectsRanges <- function(ranges, ranges2, returnLength = F) {
  ranges <- as.data.table(ranges)
  setnames(ranges, 1:3, c("scaffold", "start", "end"))

  # to keep track of the ranges by their rows
  nums <- 1:nrow(ranges)
  ranges2 <- as.data.table(ranges2)
  setnames(ranges2, 1:3, c("scaffold", "start", "end"))

  # "melts" the ranges and put all starts and ends in the same colums. intervals
  # of ranges1 are given integer numbers >=1  "c" will be useful to compute
  # coverage (+1 when a range starts, -1 when it ends)
  allRanges <- data.table(
      scaffold = rep(c(ranges$scaffold, ranges2$scaffold), 2), pos = c(ranges$start, ranges2$start, ranges$end, ranges2$end), 
      n = rep(c(nums, rep(0L, nrow(ranges2))), 2), c = rep(c(1L, -1L), each = nrow(ranges) + nrow(ranges2))
      )
  setorder(allRanges, scaffold, pos, -c)

  # computes coverage at a given position (start or end of range). Coverage
  # increases by 1 at earch start and decreases by 1 and each end. Note that
  # this considers all ranges (not just ranges2)
  allRanges[, cov := cumsum(c)]

  # to compute coverage while ignoring ranges2
  allRanges[n == 0L, c := 0L]

  # cov2 ignores ranges2
  allRanges[, cov2 := cumsum(c)]

  # so we can flag positions of overlaps between ranges and ranges2 when the two coverages differ and when cov2 is positive
  allRanges[, overlap := cov > cov2 & cov2 > 0L]

  # rows where ranges start and end
  starts <- match(nums, allRanges$n)
  ends <- matchLast(nums, allRanges$n)

  # same, but excluding ranges2
  starts2 <- match(nums, allRanges[n > 0L, n])
  ends2 <- matchLast(nums, allRanges[n > 0L, n])

  # so we can obtain ranges that overlap with ranges2
  f <- (ends - starts > ends2 - starts2) | allRanges[starts, overlap]
  if (!returnLength) {
    return(f)
  }

  # to compute length of overlaps, we remove non-intersecting ranges from the table
  allRanges <- allRanges[n %in% c(0L, nums[f])]

  # and we re-obtain first and last rows for each intersecting range
  starts <- match(nums[f], allRanges$n)
  ends <- matchLast(nums[f], allRanges$n)

  # so we can generate all the rows that each range covers in the allRanges table
  v <- unlist(Map(":", starts, ends))

  # this allows creating a table that list all the rows covered by each range, repeating rows if needed
  dt <- allRanges[v, .(overlap, pos, n2 = n, n = rep(nums[f], ends - starts + 1L))]


  # for a given range, we're only intested in positions related to starts/ends of ranges2 or the end/start or this range (where n2 == n).
  dt <- dt[n2 == 0L | n2 == n, ]

  # at the end of each range, we consider that the overlap is ended (effectively, it is, since the range ends there)
  dt[!duplicated(n, fromLast = T), overlap := F]
  row <- 2:nrow(dt)

  # we remove successive rows where the overlap value doesn't change (positions that don't create/end overlaps). In the end, each row starting and overlap will be followed by one ending the overlap
  redundant <- dt[, c(F, overlap[row] == overlap[row - 1L] & n[row] == n[row - 1L])]
  dt <- dt[!(redundant | (!duplicated(n) & !overlap))]

  # rows where overlaps start
  starts <- which(dt$overlap)

  # so we can generaate overlap blocks for each range
  blocks <- dt[, .(n = n[starts], start = pos[starts], end = pos[starts + 1L])]

  # and use these blocks to compute the total length of overlaps for each range
  perR <- blocks[, sum(end - start + 1L), by = n]

  # intersect will be this length (an interger vector)
  intersect <- integer(nrow(ranges))
  intersect[perR$n] <- perR$V1
  intersect
}


# returns a logical vector telling whether each range (bed-like) contains any position in "pos"
containsPositions <- function(ranges, pos) {
  ranges <- as.data.table(ranges)
  setnames(ranges, 1:3, c("scaffold", "start", "end"))

  # to keep track of the ranges by their rows
  nums <- 1:nrow(ranges)
  pos <- as.data.table(pos)
  setnames(pos, 1:2, c("scaffold", "posi"))

  # "melts" the ranges and put all starts and ends in the same colums.
  allRanges <- data.table(scaffold = c(rep(ranges$scaffold, 2), pos$scaffold), pos = c(ranges$start - 1L, ranges$end, pos$posi), n = c(rep(nums, 2), rep(0L, nrow(pos))))
  setorder(allRanges, scaffold, pos, n)

  # rows where ranges start and end
  starts <- match(nums, allRanges$n)
  ends <- matchLast(nums, allRanges$n)

  # same, but excluding ranges2
  starts2 <- match(nums, allRanges[n > 0L, n])
  ends2 <- matchLast(nums, allRanges[n > 0L, n])

  # so we can obtain ranges that overlap with pos
  ends - starts > ends2 - starts2
}



# returns the intervals that represent the intersections between ranges and ranges2, all in bed-like, data.table format.
intersectionsOfRanges <- function(ranges, ranges2) {
  ranges <- as.data.table(ranges)
  setnames(ranges, 1:3, c("scaffold", "start", "end"))
  nums <- 1:nrow(ranges)
  ranges2 <- as.data.table(ranges2)
  setnames(ranges2, 1:3, c("scaffold", "start", "end"))

  # "melts" the ranges and put all starts and ends in the same colums. intervals of ranges1 are given integer numbers >=1
  allRanges <- data.table(
      scaffold = rep(c(ranges$scaffold, ranges2$scaffold), 2), pos = c(ranges$start, ranges2$start, ranges$end, ranges2$end), 
      n = rep(c(nums, rep(0L, nrow(ranges2))), 2), c = rep(c(1L, -1L), each = nrow(ranges) + nrow(ranges2))
      )
  setorder(allRanges, scaffold, pos, -c)
  allRanges[, cov := cumsum(c)]
  allRanges[n == 0L, c := 0L]
  allRanges[, cov2 := cumsum(c)]
  allRanges[, overlap := cov > cov2 & cov2 > 0L]
  nr <- 2:nrow(allRanges)
  starts <- allRanges[, c(overlap[1], overlap[nr] & !overlap[nr - 1L])]
  ends <- c(F, allRanges[, !overlap[nr] & overlap[nr - 1L]])
  allRanges[, .(scaffold = scaffold[starts], start = pos[starts], end = pos[ends])]
}



# for a set of ranges specified by start and end positions (inclusive), if
# allPos is TRUE, returns the number of time each position is covered by a
# range. If false, returns a set of contiguous ranges
indivCov <- function(starts, ends, len = max(ends), allPos = T) {
  res <- cumsum(tabulate(starts, len + 1L) - tabulate(ends + 1L, len + 1L))[1:len]
  if (allPos) {
    return(res)
  }
  changes <- which(diff(res) != 0L) + 1L
  res <- data.table(start = c(1L, changes), end = c(changes - 1L, len), depth = res[c(1L, changes)])
}



mutate <- function(seq, pos, bases) {
  spl <- unlist(strsplit(stri_flatten(seq), ""))
  spl[pos] <- bases
  stri_flatten(spl)
}


genomeCov <- function(bed, minCov = 1L, seqLength = NULL, seqNames = NULL, combine = F, successive = T) {
  # computes sequencing depth of a genome by aligned reads, given alignment
  # starts and ends of reads in a data table in bed-like format. Outputs
  # successive interval for a given depth. Ignores intervals covered by less than
  # minCov reads. seqLength specify the scaffold/contig length in case one wants
  # to record intervals with 0 depth (required for last intervals of contigs
  # that may have 0 detph). "combine" combines successive intervals covered by
  # at least 1 read (i.e. regions that were sequenced)

  # remembers column names of the bed
  nams <- names(bed)[1:3]
  setnames(bed, 1:3, c("scaffold", "start", "end"))

  # if one wants to record intervals with zero depth, we apply a trick where we
  # add one artificial read totally covering each contig. This adds 1X coverage
  # to all positions, but afterwards we substract this.
  if (minCov == 0L) {
    if(is.null(seqLength)) {
      maxes = bed[,max(end), by = scaffold]
      seqLength = maxes$V1
      seqNames = maxes$scaffold
    }
    bed <- rbind(bed, data.table(scaffold = seqNames, start = 1L, end = seqLength))
  }

  # "melts" the bed to obtain 2 rows per read: alignment start and alignment end
  # + an index (n) of 1 for starts and -1 for ends. We also converts scaffold
  # names to integers to save memory. We also add 1 to end positions, as this
  # will facilitate merging intervals
  dt <- bed[, data.table(scaffID = rep(toInteger(scaffold, unique = F), 2), 
                         pos = c(start, end + 1L), n = rep(c(1L, -1L), each = .N))]

  # orders the table by scaffold and position (regardless of start and end)
  setorder(dt, scaffID, pos)

  # this is the workhorse function to compute coverage at each position since
  # reads starts add 1 to the coverage and ends remove 1
  dt[, cov := cumsum(n)]

  # reads starting or ending at the same position creates duplicate positions.
  # We only retain the last row among each set of duplicates as this is the one
  # with accurate depth
  dt <- dt[!duplicated(data.table(scaffID, pos), fromLast = T)]

  # to combine regions with at least 1X depth, we convert all depth > 1 into 1
  if (combine) dt[cov > 1L, cov := 1L]

  if (successive) {

    # we combine successive rows in the same scaffold with identical coverage
    # (may happen if a read ends just before another read starts). This
    # effectively combines region with ≥ 1X if combine == T
    dt <- dt[c(T, diff(scaffID) != 0L | diff(cov) != 0L)]
  }

  # convert the table back to interval format
  nr <- 2:nrow(dt) - 1L
  res <- dt[, data.table(scaffID = scaffID[nr], start = pos[nr], end = pos[nr + 1L] - 1L, cov = cov[nr])]

  # if we want intervals with 0 depth
  if (minCov == 0L) {
      res[, cov := cov - 1L]
  }
  
  res <- res[cov >= minCov]
  res[, scaffID := bed[scaffID, scaffold]]
  setnames(res, 1:3, nams)
  setnames(bed, 1:3, nams)
  res
}


# returns the "height" of alignments on regions, useful to draw plot were alignments are shown as overlapping horizontal segments
alignmentY <- function(bed) {

  # remembers column names of the bed
  nams <- names(bed)[1:3]
  setnames(bed, 1:3, c("scaffold", "start", "end"))

  # "melts" the bed to obtain 2 rows per read: alignment start and alignment end
  # + an index (n) of 1 for starts and -1 for ends. We also converts scaffold
  # names to integers to save memory. We also add 1 to end positions, as this
  # will facilitate merging intervals
  dt <- bed[, data.table(scaffID = rep(toInteger(scaffold, unique = F), 2), pos = c(start, end), n = rep(c(1L, -1L), each = .N), id = rep(1:.N, 2))]

  # orders the table by scaffold and position (regardless of start and end)
  setorder(dt, scaffID, pos, -n)
  Y <- integer(nrow(bed))
  id <- dt$id
  start <- dt$n == 1L
  taken <- logical(nrow(bed))
  current <- 1L
  for (i in 1:length(start)) {
    if (start[i]) {
      while (taken[current]) {
        current <- current + 1L
      }
      Y[id[i]] <- current
      taken[current] <- T
    } else {
      current <- 1L
      taken[Y[id[i]]] <- F
    }
  }
  Y

  # dt[,cov := cumsum(n)]				#this is the workhorse function to compute coverage at each position since reads starts add 1 to the coverage and ends remove 1
  # starts = dt[,which(!duplicated(id))[-1]]
  # prev = dt[starts-1L, .(n, cov, id)]
  # Y = dt[starts, cov]
  # f = prev[,n == -1L]
  # prevStartCov = dt[match(prev[f, id], id), cov]
  # prevCov = prev[f, cov]
  # Y[f] = pmin(prevStartCov, prevCov)
  # res = integer(nrow(bed))
  # res[dt[c(1L, starts), id]] = c(1L,Y)
  # setnames(bed,1:3, nams)
  # res
}



# combines regions that are distant by a certain 'distance'. By default it will
# aggregate contigous or overlapping regions, but distance could be > 0 (regions
# separated by distance pb or less) or negative (requires a certain amount of
# overlap)
combineRegions <- function(bed, distance = 0L) {
  bed <- as.data.table(bed)
  setnames(bed, c("scaffold", "start", "end"))
  bed[end < start, c("start", "end") := .(end, start)]
  bed[, end := end + distance]
  res <- genomeCov(bed, combine = T)
  res[, end := end - distance]
  res[, -4, with = F]
}



# converts scaffold positions into contig coordinates (note, "contigs" may be any non-overlapping regions of the scaffolds,like exons)
scaffoldToContigPosition <- function(pos, description, paste = F, sep = "_", unicontig = T) {
  # "pos" specifies the scaffold and position to convert and "description"
  # describes the contigs/regions: scaffold, contig identifier, begin and end,
  # in that order. Returns scaffold, contig and new position. If paste == T,
  # scaffold and contig names are pasted together (with separator specified in
  # "sep") and "unicontig" tells if the contig name should be specified for
  # scaffolds that just have one contig in the description


  # see above function
  indices <- assignToRegion(description[, c(1, 3, 4)], pos[, 1:2])

  # converts positions into contig coordinates
  posi <- pos[, 2] - description[indices, 3] + 1
  contig <- description[indices, 2]
  if (paste) {
    contig <- paste(pos[, 1], contig, sep = sep)
    if (!unicontig) {

      # list the number of occurence of a scaffold name in the description
      t <- table(description[, 1])

      # scaffolds that just have one contig
      unicontigs <- unique(names(t)[t == 1])
      filter <- pos[, 1] %in% unicontigs
      pos[, 1] <- as.vector(pos[, 1])
      contig[filter] <- pos[filter, 1]
    }
    return(data.frame(contig, pos = posi))
  }
  data.frame(scaffold = pos[, 1], contig, pos = posi)
}


# reads a multi-column file and makes one colone (the second, by default) the
# values of vector, and the first column, the names of the vector. Other columns
# are ignored. nameCol is the number of the column containing names
readNamedVector <- function(file, nameCol = 1, valueCol = nameCol + 1, ...) {
  dt <- fread(file, ..., data.table = F)
  vect <- dt[, valueCol]
  setNames(vect, dt[, nameCol])
}


# returns a distance matrix based on identity between rows or columns of a
# source matrix (or data frame) mat. dim = 1 means that rows are compared, else
# dim should be specified as 2.
distMat <- function(mat, dim = 1, asDataFrame = F) {
  if (dim == 2) mat <- t(mat)
  n <- nrow(mat)
  ncol <- ncol(mat)
  res <- matrix(NA, n, n)
  rownames(res) <- rownames(mat)
  colnames(res) <- rownames(res)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      res[i, j] <- 1 - sum(mat[i, ] == mat[j, ], na.rm = T) / ncol
    }
  }
  if (asDataFrame == T) {
    res <- as.data.frame(as.table(res), stringsAsFactors = F)
    res <- res[!is.na(res[, 3]), ]
  }
  res
}

compl <- function(bases) {
  chartr("acgtmrwsykvhdbACGTMRWSYKVHDB", "tgcakywsrmbdhvTGCAKYWSRMBDHV", bases)
}


# reverse complements DNA sequences. Ambiguities are not treated
revCom <- function(seqs) {
  com <- compl(seqs)
  seqs[] <- stri_reverse(com)
  seqs
}


Split <- function(x, f, drop = FALSE, sep = ".", recursive = F, ...) {
  # splits a vector/table recursively by each factor of "f" if f is a list
  # and if "recursive" is TRUE

  ls <- split(x, f, drop = drop, sep = sep, ...)

  if (recursive & is.list(f)) {
    for (i in 2:length(f)) {
      fields <- splitToColumns(names(ls), sep)
      names(ls) <- fields[, ncol(fields)]
      ls <- split(
        x = ls,
        f = data.frame(fields[, 1:(ncol(fields) - 1L)]),
        drop = drop,
        sep = sep,
        ...
      )
    }
  }

  ls
}

writeT <- function(data, path, quote = F, row.names = F, na = "NA", sep = "\t", ...) {
  # fwrite with different default values
  fwrite(data, path, quote = quote, row.names = row.names, na = na, sep = sep, ...)
}


# returns the sequence positions corresponding to an alignment position, that is, not counting deletions in the sequence. Start is the position of the first nucleotide of the sequence (in sequence coordinate) that aligns. Reverse must be True if the sequence is in reverse direction in the alignment
alnToSeqCoords <- function(seq, start = 1, reverse = F) {
  f <- seq != "-"
  f[is.na(f)] <- T
  seqCoord <- rep(NA, length(seq))
  close <- -1
  if (!reverse) {
    seqCoord[f] <- start:(start + sum(f) - 1)
  } else {
    close <- 1
    seqCoord[f] <- start:(start - sum(f) + 1)
  }
  w <- which(is.na(seqCoord))
  if (length(w) == 0) {
    return(seqCoord)
  }
  while (length(w) > 0) {
    seqCoord[w] <- seqCoord[w + close]
    w <- which(is.na(seqCoord))
  }
  seqCoord
}

numericClass <- function(x, size) {
  as.integer(ceiling(x / size)) * size
}




# for each numeric position of a vector (pos), returns the window it belongs to (this will be the midpoint of the window).  The step is the distance between successive windows. If windows are overlapping, returns a matrix in which elements of the same row correspond to windows covering the same position
slidingWindows <- function(pos, size, step = size) {

  # number of windows covering a position
  n <- as.integer(ceiling(size / step))
  vapply(1:n - 1L, function(x) numericClass(pos - x * step, size) + x * step - as.integer(size / 2), integer(length(pos)))
}


splitToColumns <- function(vect, split, columns, empty = "", mode = "character") {
  vect <- as.character(vect)
  nc <- stri_length(split)
  out <- NULL
  if (!missing(columns)) {
    cols <- as.integer(columns)
    cols <- cols[cols >= 1]
    out <- matrix(empty, length(vect), length(cols))
  } else {
    cols <- 1:1000
  }
  col <- 1
  repeat {
    temp <- vect
    pos <- as.vector(regexpr(split, vect, fixed = T)) - 1L
    f <- pos >= 0L
    temp[f] <- stri_sub(vect[f], 0L, pos[f])
    if (!missing(columns)) {
      if (col %in% cols) out[, match(col, cols)] <- temp
    } else {
      out <- cbind(out, temp)
    }
    col <- col + 1L
    if (col > max(cols) | !any(f)) break
    vect[f] <- stri_sub(vect[f], pos[f] + nc + 1, stri_length(vect[f]))
    vect[!f] <- empty
  }
  if (ncol(out) == 1L) out <- as.vector(out)
  storage.mode(out) <- mode
  out
}



# same as duplicated (but only for vectors) but marks all elements present more than once as duplicated (even the first occurrence)
allDuplicated <- function(vect) {
  return(duplicated(vect) | duplicated(vect, fromLast = T))
}


# returns true for the first occurrence of an element that is duplicated in the vector
isDuplicated <- function(vect) {
  return(!duplicated(vect) & duplicated(vect, fromLast = T))
}


# returns an object of class fasta (defined in muscle package) from a vector containing sequences and names
makeFasta <- function(seqs, n = names(seqs)) {
  seqs <- data.frame(V1 = n, V2 = as.character(seqs))
  fasta <- list(seqs = seqs, num = nrow(seqs))
  class(fasta) <- "fasta"
  fasta
}


# returns only a subset of data (a data frame or matrix) without repetition of values in the colmumn repeated, based on the value of filter. "filter" and "repeated" are column numbers or column names (with quotes) of the data object
retainBest <- function(data, repeated, filter, decreasing = F, random = F) {
  if (random) data <- data[sample(nrow(data)), ]
  data <- data[order(data[, filter], decreasing = decreasing), ]
  data[!duplicated(data[, repeated]), ]
}


# returns reciprocal best hits from blast outputs, tabular format (and where only the best hists were returned)
rbh <- function(forward, reverse) {

  # select only the best hit for each query. Assumes that the last column is the blast score
  forward <- retainBest(forward, 1, ncol(forward))
  reverse <- retainBest(reverse, 1, ncol(reverse))

  # identifier of hits in the forward blast
  fHits <- paste(forward[, 1], forward[, 2], sep = "___")

  # identifier of hits in the reverse direction, but reversed, so that they appear in the same direction as the forward blast
  revRevHits <- paste(reverse[, 2], reverse[, 1], sep = "___")
  shared <- unique(intersect(fHits, revRevHits))
  splitToColumns(shared, "___")
}


# converts scaffold coordinates to gene coordinates given a "pos" data frame that contains scaffold names and positions, and an "exons" dataframe that contains exon coordinates, including the gene strand at the 4th column (either "+" or "-") and the gene ID as 5th column
scaffToGeneCoord <- function(pos, exons) {
  # we firt determine the start position of each exon in the coordinates of its gene. For this we sort exons in order of appearance in the gene, which depends on the gene strand
  if (!is.data.table(exons)) exons <- as.data.table(exons)
  if (!is.data.table(pos)) pos <- as.data.table(pos)
  setnames(exons, 1:5, c("scaffold", "start", "end", "strand", "gene"))
  setnames(pos, 1:2, c("scaffold", "pos"))

  # for genes in minus strand, the exon that comes last in the scaffold should come first
  exons[, start2 := ifelse(strand == "+", start, -start)]
  setorder(exons, gene, start2)
  exons[, size := end - start + 1L]
  exons[, fromStart := c(0L, size[-length(size)])]
  exons[!duplicated(gene), fromStart := 0L]
  exons$fromStart <- exons[, cumsum(fromStart), by = gene]$V1

  # we determine the exon to which the position belongs (this will be just one exon, even if several overlap)
  exon <- assignToRegion(exons, pos)
  posInGene <- exons[exon, ifelse(strand == "+", pos$pos - start + 1L + fromStart, end - pos$pos + 1L + fromStart)]
  data.table(gene = exons[exon, gene], pos = posInGene)
}



# converts gene coordinates to scaffold coordinates given a "pos" data frame that contains gene names and positions, and an "exons" dataframe that contains exon coordinates, including the gene strand at the 4th column (either "+" or "-") and the gene name as 5th column
geneToScaffCoord <- function(pos, exons) {
  # we firt determine the start and end positions of each exon in the gene coordinate. For this we sort exons in order of appearance in the gene, which depends on the gene strand
  pos <- as.data.table(pos)
  exons <- as.data.table(exons)
  setnames(exons, 1:5, c("scaffold", "start", "end", "strand", "gene"))
  setnames(pos, 1:2, c("gene", "posInGene"))

  # for genes in minus strand, the exon that comes last in the scaffold should come first
  exons[, start2 := ifelse(strand == "-", -end, start)]
  setorder(exons, gene, start2)
  exons[, endInGene := exons[, cumsum(end - start + 1L), by = gene]$V1]
  exons[, startInGene := endInGene - (end - start)]

  # we determine the exon to which each position belongs (this will be just one exon, even if several overlap)
  exon <- assignToRegion(exons[, .(gene, startInGene, endInGene)], pos[, .(gene, posInGene)])
  exons[exon, .(scaffold, posInScaff = abs(start2) + sign(start2) * (pos$posInGene - startInGene))]
}



geneToScaffRanges <- function(geneRanges, exons) {

  # converted in scaffold coordinates
  scaffRanges <- data.frame(geneToScaffCoord(geneRanges[, 1:2], exons), end = geneToScaffCoord(geneRanges[, c(1, 3)], exons)[, 2])

  # ensure that start positions < end positions (due to genes on the reverse strand)
  scaffRanges <- fixRanges(scaffRanges, 2)

  # we replace scaffold names by original gene names for the following instruction
  scaffRanges[, 1] <- geneRanges[, 1]

  # this is to exclude introns, as specified in the exons. We needed to use gene names as certain genes may completely encompass others, so some exons might be retained in full although they shouldn't.
  scaffRanges <- intersectionsOfRanges(scaffRanges, exons[, c(5, 2, 3)])

  # restablishes scaffold names (allowed because a gene encompasses only one scaffold)
  scaffRanges$scaffold <- exons[match(scaffRanges$scaffold, exons[, 5]), 1]
  scaffRanges
}

rangeStarts <- function(lengths) {
  cumsum(as.numeric(c(1, lengths[-length(lengths)])))
}


# converts scaffold/contig positions into unique genome positions assuming contiguous scaffolds whose lengths and name are provided in order
scaffToGenomeCoord <- function(scaffolds, pos, scaffLengths, scaffNames) {
  names(scaffLengths) <- NULL
  scaffStart <- rangeStarts(scaffLengths)
  if (!is.numeric(scaffNames)) {
    return(pos + scaffStart[chmatch(scaffolds, scaffNames)] - 1L)
  } else {
    return(pos + scaffStart[match(scaffolds, scaffNames)] - 1L)
  }
}


# converts scaffold/contig positions into unique genome positions assuming contiguous scaffolds whose lengths and name are provided in order
scaffToGenomeBed <- function(scaffolds, starts, ends, scaffLengths, scaffNames) {
  genomeStarts <- scaffToGenomeCoord(scaffolds, starts, scaffLengths, scaffNames)
  data.table("genome", start = genomeStarts, end = genomeStarts + ends - starts)
}

genomeToScaffCoord <- function(pos, scaffLengths, scaffNames) {
  scaffStart <- cumsum(c(1L, scaffLengths))
  scaff <- .bincode(pos, scaffStart, F, )
  data.frame(scaffold = scaffNames[scaff], pos = pos - scaffStart[scaff] + 1L)
}

patternLeadingGap <- function(x) {
  stopifnot(is(x, "PairwiseAlignments"), type(x) == "global")
  start(subject(x)) - 1L
}

subjectLeadingGap <- function(x) {
  stopifnot(is(x, "PairwiseAlignments"), type(x) == "global")
  start(pattern(x)) - 1L
}


alignWithEndGaps <- function(seq1, seq2, ...) {
  aln <- pairwiseAlignment(seq1, seq2, ...)
  aln <- data.table(pattern = as.character(pattern(aln)), subject = as.character(subject(aln)), pgap = patternLeadingGap(aln), sgap = subjectLeadingGap(aln))
  f <- aln$sgap > 0
  missing <- stri_sub(seq1[f], 1, aln$sgap[f])
  aln[f, subject := stri_join(gsub("[^_]", "-", missing), subject)]
  aln[f, pattern := stri_join(missing, pattern)]
  f <- aln$pgap > 0
  missing <- stri_sub(seq2[f], 1, aln$pgap[f])
  aln[f, pattern := stri_join(gsub("[^_]", "-", missing), pattern)]
  aln[f, subject := stri_join(missing, subject)]

  subjectTrailingGap <- nchar(seq1) - nchar(gsub("-", "", aln$pattern))
  patternTrailingGap <- nchar(seq2) - nchar(gsub("-", "", aln$subject))

  f <- subjectTrailingGap > 0
  missing <- stri_sub(seq1[f], nchar(seq1[f]) - subjectTrailingGap[f] + 1, nchar(seq1[f]))
  aln[f, subject := stri_join(subject, gsub("[^_]", "-", missing))]
  aln[f, pattern := stri_join(pattern, missing)]
  f <- patternTrailingGap > 0
  missing <- stri_sub(seq2[f], nchar(seq2[f]) - patternTrailingGap[f] + 1, nchar(seq2[f]))
  aln[f, pattern := stri_join(pattern, gsub("[^_]", "-", missing))]
  aln[f, subject := stri_join(subject, missing)]
  return(aln[, .(pattern, subject)])
}

mcAlign <- function(seqs1, seqs2, ncores = 4, ...) {
  sizes <- round(length(seqs1) / (ncores * 2))
  starts <- seq(1, length(seqs1), round(length(seqs1) / ncores))
  ends <- c(starts[2:length(starts)] - 1, length(seqs1))
  l <- vector("list", length(starts))
  for (i in seq_along(starts)) l[[i]] <- cbind(seqs1[starts[i]:ends[i]], seqs2[starts[i]:ends[i]])
  aln <- mclapply(l, alignApply, mc.cores = ncores, mc.preschedule = T)
  rbindlist(aln)
}

mcAlign2 <- function(seq1, seq2, ncores = 4, ...) {
  ls <- nchar(seq1)
  parts <- sum(ls) / ncores
  splitter <- as.integer(cumsum(ls) / parts) + 1
  splitter[splitter > ncores] <- ncores
  aln <- mcMap(function(seq1, seq2) alignWithEndGaps(seq1, seq2, ...), 
               split(seq1, splitter), split(seq2, slitter), mc.cores = ncores, mc.preschedule = F)
  rbindlist(aln)
}



# splits a string vector into a matrix with on character per cell. Words are put in rows. The number of columns is the length of the longest word
stringToMatrix <- function(vector) {
  vector <- as.character(vector)
  nchars <- stri_length(vector)
  m <- max(nchars, na.rm = T)
  mat <- matrix(NA, nrow = length(vector), ncol = m, byrow = T)
  pb <- txtProgressBar(char = "*", width = 100, style = 3)
  for (i in 1:m) {
    mat[, i] <- stri_sub(vector, i, i)
    setTxtProgressBar(pb, i / m)
  }
  mat
}


# converts ranges that have specific values (e.g. sequencing depth on a region) to point coordinates so this value can be represented by a line
rangesToPoints <- function(start, end, value) {
  x <- rep(start, each = 2)
  x <- x[-1]
  y <- rep(value, each = 2)
  y <- y[-length(y)]
  cbind(x, y)
}


# Returns a data frame of abba, baba and bbaa indices for each SNP. genotypes is a matrix of 0/1 alleles (columns = individual alleles) and the four taxa are specified as column numbers of the matrix. A taxon can correspond to several individuals, the indices are frequency-based
abbaIndices <- function(genotypes, allo, symp1, symp2, outgroup) {

  freq = vapply(list(allo, symp1, symp2, outgroup), 
                function(col) rowMeans(cbind(genotypes[, col]), na.rm = T)/2, numeric(nrow(genotypes)))
  
  abbaF <- (1 - freq[, 1]) * freq[, 2] * freq[, 3] * (1 - freq[, 4]) + freq[, 1] * (1 - freq[, 2]) * (1 - freq[, 3]) * freq[, 4]

  # site supporting the phylogeny (for reference)
  bbaaF <- (1 - freq[, 1]) * (1 - freq[, 2]) * freq[, 3] * freq[, 4] + freq[, 1] * freq[, 2] * (1 - freq[, 3]) * (1 - freq[, 4])
  babaF <- freq[, 1] * (1 - freq[, 2]) * freq[, 3] * (1 - freq[, 4]) + (1 - freq[, 1]) * freq[, 2] * (1 - freq[, 3]) * freq[, 4]
  res <- data.table(abbaF, babaF, bbaaF)
}

regionSets <- function(regions) {
  pb <- txtProgressBar(char = "*", width = 100, style = 3)
  nrows <- nrow(regions)
  rows <- 2:nrows
  newScaff <- c(FALSE, regions[rows, 1] != regions[rows - 1, 1], T)
  ends <- c(regions[, 3], 0)

  # gives a number ("set" identifier) to all regions that are included in the same genomic region (which should come first in the table)
  set <- rep(NA, nrows)
  end <- ends[1]
  n <- 1
  for (i in 2:(nrows + 1)) {
    if (ends[i] > end | newScaff[i]) {
      set[n:(i - 1)] <- n
      n <- i
      end <- ends[i]
      setTxtProgressBar(pb, i / nrows)
    }
  }
  set
}


# returns gene sequences from a sequence vector having scaffolds/contigs as names, using a dataframe of exon coordinates. Pastes the exons together and reverses complements genes as needed. The exons data.frame contains scaffold, exon starts and ends, the gene strand at the 4th column (either "+" or "-") and the gene name at the 5th column
geneSeqFromScaffolds <- function(scaffolds, exons) {
  exons <- as.data.table(exons)
  setnames(exons, 1:5, c("scaffold", "start", "end", "strand", "gene"))
  setorder(exons, gene, start)

  # extracts exon sequences
  exons[, seq := as.character(subseq(scaffolds[scaffold], start, end))]

  # pastes exons per gene together
  geneSeqs <- exons[, stri_flatten(seq), by = .(gene, strand)]
  geneSeqs[strand == "-", V1 := revCom(V1)]
  setNames(DNAStringSet(geneSeqs$V1), geneSeqs$gene)
}

aas__ <- c("S", "S", "S", "S", "F", "F", "L", "L", "Y", "Y", "*", "*", "C", "C", "*", "W", "L", "L", "L", "L", "P", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M", "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "P", "P", "P", "H", "V", "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G", "G", "G", "-", NA)

names(aas__) <- c("TCA", "TCC", "TCG", "TCT", "TTC", "TTT", "TTA", "TTG", "TAC", "TAT", "TAA", "TAG", "TGC", "TGT", "TGA", "TGG", "CTA", "CTC", "CTG", "CTT", "CCA", "CAT", "CAA", "CAG", "CGA", "CGC", "CGG", "CGT", "ATA", "ATC", "ATT", "ATG", "ACA", "ACC", "ACG", "ACT", "AAC", "AAT", "AAA", "AAG", "AGC", "AGT", "AGA", "AGG", "CCC", "CCG", "CCT", "CAC", "GTA", "GTC", "GTG", "GTT", "GCA", "GCC", "GCG", "GCT", "GAC", "GAT", "GAA", "GAG", "GGA", "GGC", "GGG", "GGT", "---", NA)


# translates codons with the universal genetic code
transCodons <- function(codons, RNA = F) {
  codons <- toupper(codons)
  if (RNA) codons <- chartr("U", "T", codons)
  as.vector(aas__[codons])
}

splitStringInParts <- function(string, size) {
  part <- paste0(".{1,", size, "}")
  setNames(stri_extract_all_regex(string, part), names(string))
}


translateDNA <- function(seqs, frame, join = T) {
  if (!missing(frame)) {
    seqs[frame < 0] <- revCom(seqs[frame < 0])
    frame <- abs(frame)
    newFrame <- frame > 1
    seqs[newFrame] <- stri_sub(seqs[newFrame], frame[newFrame], stri_length(seqs[newFrame]))
  }
  codons <- splitStringInParts(seqs, 3)
  AAs <- lapply(codons, transCodons)
  if (join) {
    AAs <- lapply(AAs, function(x) {
      x[is.na(x)] <- "?"
      x
    })
    AAs <- (sapply(AAs, stri_flatten))
  }
  AAs
}


# gets a nucleotide alignment based on a protein alignment and the corresponding cds
aaToCDS <- function(aas, cds, collapse = T) {
  if (any(stri_length(gsub("-", "", aas, fixed = T)) != stri_length(cds) / 3)) stop("CDS length must be exactly 3 times protein length")
  codons <- splitStringInParts(cds, 3)
  aas <- strsplit(aas, "", fixed = T)
  fillCodons <- function(aas, codons) {
    res <- character(length(aas))
    deleted <- aas == "-"
    res[!deleted] <- codons
    res[deleted] <- "---"
    res
  }
  res <- Map(fillCodons, aas, codons)
  if (collapse) res <- sapply(res, stri_flatten)
  res
}


# builds upon Biostrings' subseq, to automatically reverse complement the sequences (by default) if end coordinates are lower than starts
subSeq <- function(string, start, end, reverse = T) {
  f <- start > end
  temp <- end
  end[f] <- start[f]
  start[f] <- temp[f]
  subs <- subseq(string, start, end)
  if (reverse & any(f)) {
    subs[f] <- reverseComplement(subs[f])
  }
  subs
}


# returns the gene, codons, and amino acids involved in diallelic snps located
# in scaff, pos for the 2 alternative bases (character vectors). scaffSeqs is a
# vector (or DNAStringSet) of scaffold sequences with names as scaffolds. The
# gff data table that contains exon coordinates (3 columns), + the gene strand
# at the 4th column (either "+" or "-") and the gene name/ID as 5th column. the
# "ref" argument enables checking for errors, in case "base" is the base of the
# reference genome (if for example the pileup used a reference genome)
codingChanges <- function(scaff, pos, gff, scaffSeqs, base, base2, ref = F) {
  setnames(gff, 1:5, c("scaffold", "start", "end", "strand", "gene"))
  snps <- data.table(scaffold = scaff, pos, base, base2)
  snps[, c("gene", "genePos") := scaffToGeneCoord(.(scaff, pos), gff)]
  sub <- snps[!is.na(gene)]
  gene <- sub$gene
  posInGene <- sub$genePos
  base <- sub$base
  base2 <- sub$base2
  codon <- ceiling(posInGene / 3L)
  posInCodon <- posInGene %% 3L
  posInCodon[posInCodon == 0L] <- 3L
  gff <- gff[gene %in% sub$gene, ]
  geneSeqs <- geneSeqFromScaffolds(scaffSeqs, gff)
  refCodon <- stri_sub(geneSeqs[gene], codon * 3 - 2, codon * 3)
  toCompl <- gff[match(sub$gene, gene), strand == "-"]
  base[toCompl] <- compl(base[toCompl])
  base2[toCompl] <- compl(base2[toCompl])
  if (ref) {
    refBase <- stri_sub(geneSeqs[gene], posInGene, length = 1L)
    f <- base != refBase
    if (any(f)) warning(paste(sum(f), "mistmatche(s) out of", length(f), "between data and scaffold sequences"))
  }

  stri_sub(refCodon, posInCodon, posInCodon) <- base
  refAA <- transCodons(refCodon)
  altCodon <- refCodon
  stri_sub(altCodon, posInCodon, posInCodon) <- base2
  altAA <- transCodons(altCodon)
  data.table(sub[, .(scaffold, pos)], gene, posInGene, posInCodon, refBase = base, altBase = base2, refCodon, altCodon, refAA, altAA)
}


CDSallelesFromSNPs <- function(scaff, pos, gff, scaffSeqs, base, base2 = NA) {
  # scaff = scaffold names, pos = SNP positions (integer vector in scaffold
  # coordinates), gff = 5-column data.table of CDS (scaffold, start, end,
  # strand, gene name with one line per exon) using gff conventions, hence strand is either "+" or
  # "-". scaffSeqs is a DNAstringSet of scaffold sequences (names must match
  # those of the scaff argument), "base" and "base2" are 1-letter character
  # vectors indicating the mutant base at each position of "pos". base2 can
  # optionally specify another allele. If it is NA (default), it is taken as the base
  # in the scaffold sequence (reference base). scaff, pos, base and base2 must
  # have the same length.
  setnames(gff, 1:5, c("scaffold", "start", "end", "strand", "gene"))
  snps <- data.table(scaff, pos, base, base2)
  snps[, c("gene", "genePos") := scaffToGeneCoord(.(scaff, pos), gff)]
  snps <- snps[!is.na(gene)]
  snps[, strand := gff[match(snps$gene, gene), strand]]
  snps[strand == "-", c("base", "base2") := .(compl(base), compl(base2))]
  geneSeqs <- geneSeqFromScaffolds(scaffSeqs, gff[gene %in% snps$gene, ])

  nucl <- strsplit(as.character(geneSeqs), "")
  mut <- function(gene, pos, base) {
    nuc <- nucl[[gene]]
    if (!any(is.na(base))) {
      nuc[pos] <- base
    }
    stri_flatten(nuc)
  }
  alleles <- snps[, .(n = .N, all1 = mut(gene[1], genePos, base), all2 = mut(gene[1], genePos, base2)), by = gene]
}


pairWisedNdS <- function(alleles) {
  # computes dN and dS from  a data.table generated by the previous function (1
  # column for gene name and 2 other columns containing the sequences)
  require(seqinr)
  alns <- apply(alleles[, cbind(all1, all2)], 1, seqinrAlignment)
  res <- mclapply(alns, kaks, mc.cores = 4)
  res <- lapply(res, unlist)
  data.table(alleles[, .(gene, n)], do.call(rbind, res))
}


# takes a numeric vector and assigns elements to groups based on the maximum distance between them. Returns the groups as integers
clusterize <- function(vect, limit, group = "n") {
  dt <- data.table(g = group, v = vect)
  dt$n <- 1:nrow(dt)
  dt <- dt[order(g, v), ]
  vect <- dt$v
  group <- dt$g
  pos <- 2:length(vect)
  w <- c(1, which(vect[pos] - vect[pos - 1] > limit | group[pos] != group[pos - 1]) + 1, length(vect) + 1)
  pos <- 2:length(w)
  groupSize <- w[pos] - w[pos - 1]
  group <- rep(1:length(groupSize), groupSize)
  group[order(dt$n)]
}


groupsOfMinSize = function(v, minSize,  f = NA) {
  # returns the group to which each element of v belongs, ensuring that 
  # the sum of v in each group is at least minSize. v must be numeric
  # factor is an optional grouping factor f, ensuring that a group does not span several levels
  # f must be ordered!
  l = length(v)
  g = integer(l)
  if(all(is.na(f))) {
    f = g
  } 
  cs = 0
  current = 0L
  newLevel = c(F, f[2:l] != f[2:l-1L])
  for(i in 1:l) {
    if(cs >= minSize | newLevel[i]) {
      current = current+1L
      cs = 0
    }
    cs = cs + v[i]
    g[i] = current
  }
  dt = data.table(v, g)
  n = dt[,sum(v), by = g]
  tooLow = n[V1 < minSize, g]
  i = g %in% tooLow
  g[i] = g[i] - 1L
  g
}


mutateSequences = function(contig, pos, allele, contigSeqs) {
  starts = split(pos, f = contig)
  widths = split(nchar(allele), f = contig)
  bases = split(allele, f = contig)
  
  iRanges = Map(IRanges, start = starts, width = widths)
  iRangeList = do.call(IRangesList, iRanges)
  
  contigNames = names(starts)
  extracted = as.character(unlist(extractAt(contigSeqs[contigNames], at = iRangeList), use.names = F))
  mutated = replaceAt(contigSeqs[contigNames], at = iRangeList, value = bases)  
  contigSeqs[contigNames] = mutated
  contigSeqs
}


# generates a seqinr alignment object from a character vector
seqinrAlignment <- function(seqs, names = NULL) {
  if (is.null(names)) {
    if (!is.null(names(seqs))) {
      names <- names(seqs)
    } else {
      names <- paste(rep("seq", length(seqs)), 1:length(seqs), sep = "_")
    }
  }
  if (var(nchar(seqs)) > 0) stop("sequence lengths differ")
  aln <- list(nb = length(seqs), nam = names, seq = as.character(seqs), com = NA)
  class(aln) <- "alignment"
  aln
}



reorderBlastFasta <- function(fasta, input) {
  input <- data.frame(input)
  gi <- splitToColumns(names(fasta), "|", 2)
  coords <- splitToColumns(names(fasta), ":", 2)
  noCoord <- gi[is.na(coords)]
  inputCoords <- as.numeric(splitToColumns(input[, 2], "-"))
  input$length <- inputCoords[, 2] - inputCoords[, 1]
  longInput <- retainBest(input, 1, 3)
  coords[is.na(coords)] <- longInput[match(noCoord, longInPut[, 1]), 2]
  names(fasta) <- paste(gi, coords)
  fasta[paste(input[, 1], input[, 2])]
}



# returns the first position of the n highest elements of a vector
nHighest <- function(vect, n) {
  sorted <- sort(vect, T)[1:n]
  match(sorted, vect)
}


TwoDHeatMap <- function(x, y, nbins, xlim = c(min(x), max(x)), ylim = c(min(y), max(y))) {
  freq <- as.matrix(table(cut(x, nbins), cut(y, nbins)))
  freq[freq > 500] <- 500
  image((1:nbins) / nbins, (1:nbins) / nbins, freq, col = topo.colors(256), xlab = "SNPs", ylab = "score BLASR")
}


# returns the first position of a value in each row of a matrix
rowMatches <- function(values, mat, nomatch = NA) {
  matches <- rep(nomatch, nrow(mat))
  w <- which(mat == values, arr.ind = T)
  matches[w[, 1]] <- w[, 2]
  matches
}


splitAlignment <- function(aln, coords, fillGaps = F) {
  setnames(coords, c("scaff1", "scaff2", "start", "end", "start2", "end2"))
  nc <- nchar(aln$pattern)
  nuc <- data.table(base1 = unlist(strsplit(stri_flatten(aln$pattern), "")), base2 = unlist(strsplit(stri_flatten(aln$subject), "")), scaff1 = rep(coords$scaff1, nc), scaff2 = rep(coords$scaff2, nc))
  nuc$pos1 <- 0L
  nuc$pos2 <- 0L
  nuc[base1 != "-", pos1 := unlist(Map(":", coords$start, coords$end))]
  nuc[base2 != "-", pos2 := unlist(Map(":", coords$start2, coords$end2))]
  nuc$aln <- rep(1:nrow(coords), nc)
  if (fillGaps) {
    n <- 2:nrow(nuc)
    newAln <- nuc[, c(T, aln[n] != aln[n - 1L], T)]
    fillG <- function(pos) {
      inGap <- pos == 0L
      startGap <- which(inGap[n] & (!inGap[n - 1L] | newAln[n])) + 1L
      endGap <- which(inGap[n - 1L] & (!inGap[n] | newAln[n]))
      if (tail(inGap, 1)) endGap <- c(endGap, max(n))
      previousPos <- pos[startGap - 1L]
      if (startGap[1] == 1L) previouPos <- c(0L, previousPos)
      nextPos <- pos[endGap + 1L]
      if (is.na(tail(nextPos, 1))) nextPos[length(nextPos)] <- 0L
      toTake <- ifelse((previousPos < nextPos & !newAln[startGap]) | newAln[endGap + 1L], previousPos, nextPos)
      pos[inGap] <- rep(toTake, endGap - startGap + 1L)
      pos
    }
    nuc[, c("pos1", "pos2") := .(fillG(pos1), fillG(pos2))]
  }
  nuc
}


# assigns elements to clusters via single-linkage clustering, based on pairs (e.g., query - subject). Pairs are elements of V1 and V2 (of equal length) at the same position
clusterFromPairs <- function(V1, V2, reciprocal = F, int = F) {
  # If an element is found in several pairs (and clusters, iteratively), the elements of these pairs (or clusters) will be assigned to the same cluster
  if (length(V1) != length(V2)) stop("provide vectors of equal length")
  m <- 0
  if (!int) {

    # we assign unique cluster ids (integers) to elements (faster than union())
    uniq <- unique(c(V1, V2))
    if (is.character(V1)) V1 <- chmatch(V1, uniq) else V1 <- match(V1, uniq)
    if (is.character(V2)) V2 <- chmatch(V2, uniq) else V2 <- match(V2, uniq)

    # will give the correspondence between current cluster ids (the indices of m) and new cluter ids (values of m), for replacement. Initially, they are the same
    m <- 1:length(uniq)
  } else {
    m <- 1:max(V1, V2)
  }

  # if there is no reciprocity in pairs
  if (!reciprocal) {

    # we now have reciprocal pairs of clusters
    pairs <- data.table(q = c(V1, V2), s = c(V2, V1))
  } else {
    pairs <- data.table(q = V1, s = V2)
  }
  f <- pairs[, q > s]

  # clustering proceeds until the 2 elements of every pair are assigned to the same cluster
  while (any(f)) {

    # for a given query cluster we need lowest cluster id in all pairs where it is found, including itself. So there's no need to use rows for which the cluster is lower than the query
    mins <- pairs[f, min(s), by = q]

    # and we replace each cluster id by this lowest id in both query and subject clusters (these 2 lines are much faster than a match())
    m[mins$q] <- mins$V1
    pairs[, c("q", "s") := .(m[q], m[s])]
    f <- pairs[, q > s]
    cat("*")
  }

  # we get cluster ids corresponding to V1 elements since they are the same for V2
  res <- pairs$q[1:length(V1)]

  # this removes possible gaps in cluster numbering
  match(res, unique(res))
}


pastePair <- function(string1, string2, sep = "_") {
  f <- string1 < string2
  if (!is.character(string1)) {
    string1 <- as.character(string1)
    string2 <- as.character(string2)
  }
  ifelse(f, stri_join(string1, string2, sep = sep), stri_join(string2, string1, sep = sep))
}


# returns all possible pairs for elements of a vector, or two vectors (2dn vector specified in v2). Much faster than combn(v, 2). reciprocal = T means reciprocal pairs are returned. sort = T means pairs are sorted alphabetically or numerically. same tells if pairs comprising the same element twice should be returned (only valid if v2 = NA)
allPairs <- function(v, sort = T, noDups = T, reciprocal = F, same = F, v2 = NA) {
  if (noDups) {
    v <- unique(v)
    v2 <- unique(v2)
  }
  if (sort) {
    v <- sort(v)
    v2 <- sort(v2)
  }
  
  if (all(is.na(v2))) {
    i = as.integer(!same)
    lv <- length(v)
    V2 <- unlist(lapply((i + 1):lv, function(x) v[x:lv]), use.names = F)
    res = data.table(V1 = rep(v, lv:1 - i), V2)
  }
  else {
    res <- setNames(data.table(expand.grid(v, v2, stringsAsFactors = F)), c("V1","V2"))
  }
  if (reciprocal) res <- rbind(res, res[, .(V1 = V2, V2 = V1)])
  res
}


intMatch <- function(int1, int1b, replacement, nomatch = 0L) {
  corres <- vector(mode(replacement), length = max(int1, int1b))
  # corres = integer(length = max(int1,int1b))
  corres[int1b] <- replacement
  corres[int1]
}

data.tableFromCommunities <- function(comm) {
  cNumbers <- unlist(Map(rep, 1:length(comm), sapply(comm[1:length(comm)], length)))
  members <- unlist(comm[1:length(comm)])
  data.table(member = as.integer(members), community = cNumbers)
}

integerPair <- function(v1, v2) {
  i1 <- chmatch(v1, unique(v1))
  i2 <- chmatch(v2, unique(v2))
  pair <- i1 * (max(i2) + 1L) + i2
  match(pair, unique(pair))
}

HSPgroup <- function(blastResults, maxDist = 50, maxOverlap = 30, maxDiff = 50, blastX = F) {
  blast <- blastResults[, .(query, subject, qStart, qEnd, sStart, sEnd, n = 1:.N)]
  if (blastX) blast[, c("sStart", "sEnd") := .(sStart * 3, sEnd * 3)]
  blast[, pair := integerPair(query, subject)]
  blast[sEnd < sStart, c("sStart", "sEnd") := .(-sStart, -sEnd)]
  blast[qEnd < qStart, c("qStart", "qEnd") := .(-qStart, -qEnd)]
  setorder(blast, pair, qStart, sStart)
  rows <- 2:nrow(blast)
  samePair <- c(F, blast[, pair[rows] == pair[rows - 1L]])
  distQuery <- c(0L, blast[, qStart[rows] - qEnd[rows - 1L]])
  distQuery[!samePair] <- NA
  distSubject <- c(0L, blast[, sStart[rows] - sEnd[rows - 1L]])
  distSubject[!samePair] <- NA
  w <- which(distSubject <= maxDist & distQuery <= maxDist & distSubject >= -maxOverlap & distQuery >= -maxOverlap & abs(distSubject - distQuery) <= maxDiff)
  c <- clusterize(w, 1L)
  group <- integer(nrow(blast))
  group[w] <- c
  group[group == 0L] <- max(group) + 1:sum(group == 0L)
  group[order(blast$n)]
}


combineHits <- function(blast, maxDist = 50, maxOverlap = 30, maxDiff = 50, blastX = F) {
  blast <- blast[, 1:12, with = F]
  setnames(blast, c("query", "subject", "identity", "length", "mismatches", "indels", "qStart", "qEnd", "sStart", "sEnd", "evalue", "score"))
  blast[evalue > 1, evalue := 1]
  if (blastX) blast[, c("sStart", "sEnd") := .(sStart * 3, sEnd * 3)]
  blast[, pair := stri_c(query, subject)]
  blast[sEnd < sStart, c("sStart", "sEnd") := .(-sStart, -sEnd)]
  blast[qEnd < qStart, c("qStart", "qEnd") := .(-qStart, -qEnd)]
  setorder(blast, pair, qStart, sStart)
  rows <- 2:nrow(blast)
  samePair <- c(F, blast[, pair[rows] == pair[rows-1L]])
  distQuery <- c(0L, blast[, qStart[rows] - qEnd[rows-1L]])
  distQuery[!samePair] <- NA
  distSubject <- c(0L, blast[, sStart[rows] - sEnd[rows-1L]])
  distSubject[!samePair] <- NA
  w <- which(distSubject <= maxDist & distQuery <= maxDist & distSubject >= -maxOverlap & distQuery >= -maxOverlap & abs(distSubject - distQuery) <= maxDiff)
  if (length(w) == 0) {
    print("no HSP can be combined")
    return(blast[, -13, with = F])
  }
  c <- clusterize(w, 1)
  blast[w, group := c]
  blast[setdiff(w - 1, w), group := (unique(c))]
  blast[, c("qStart", "qEnd", "sStart", "sEnd") := .(abs(qStart), abs(qEnd), abs(sStart), abs(sEnd))]

  queryRegions <- combineRegions(blast[!is.na(group), data.frame(group, qStart, qEnd)], distance = maxDist * 2)
  subjectRegions <- combineRegions(blast[!is.na(group), data.frame(group, sStart, sEnd)], distance = maxDist * 2)
  regionStats <- blast[!is.na(group), .(query = query[1], subject = subject[1], identity = mean(identity), length = sum(length), mismatches = sum(mismatches), indels = sum(indels) + .N - 1, evalue = prod(evalue), score = sum(score) - .N * 10), by = group]

  combined <- data.table(regionStats[, -1, with = F], queryRegions[, -1, with = F], subjectRegions[, -1, with = F])
  setcolorder(combined, c(1:6, 9:12, 7:8))
  setnames(combined, names(blast)[1:12])
  all <- rbind(combined, blast[is.na(group), -(13:14), with = F])
  if (blastX) all[, c("sStart", "sEnd") := .(sStart / 3, sEnd / 3)]
  all
}



hitBlocks <- function(hits, maxDist = 50, maxOverlap = 30, maxDiff = 50) {
  setnames(hits, 1:6, c("query", "subject", "qStart", "qEnd", "sStart", "sEnd"))
  hits[, c("pair", "n") := .(stri_c(query, subject), 1:.N)]
  hits[sEnd < sStart, c("sStart", "sEnd") := .(-sStart, -sEnd)]
  hits[qEnd < qStart, c("qStart", "qEnd") := .(-qStart, -qEnd)]
  setorder(hits, pair, qStart, sStart)
  rows <- 2:nrow(hits)
  samePair <- c(F, hits[, pair[rows] == pair[rows-1L]])
  distQuery <- c(0L, hits[, qStart[rows] - qEnd[rows-1L]])
  distQuery[!samePair] <- NA
  distSubject <- c(0L, hits[, sStart[rows] - sEnd[rows-1L]])
  distSubject[!samePair] <- NA
  w <- which(distSubject <= maxDist & distQuery <= maxDist & distSubject >= -maxOverlap & distQuery >= -maxOverlap & abs(distSubject - distQuery) <= maxDiff)
  if (length(w) == 0) {
    print("no hit can be combined")
    return(1:nrow(hits))
  }
  c <- clusterize(w, 1)
  hits[w, group := c]
  hits[setdiff(w - 1, w), group := (unique(c))]
  m = max(hits$group, na.rm = T)
  hits[is.na(group), group := 1:.N + nrow(hits)]
  setorder(hits, n)
  toInteger(hits$group)
}


# removes reciprocal and self hits (if removeSelf = T) from blast results, format 6. Score is expected in column 12.
removeReciprocal <- function(blast, removeSelf = T) {
  blastb <- copy(blast)
  names <- colnames(blast)
  setnames(blastb, stri_c("V", 1:ncol(blastb)))
  blastb[, pair := pastePair(V1, V2)]
  setorder(blastb, -V12)
  f <- T
  if (removeSelf) f <- blastb[, V1 != V2]
  blastb <- blastb[!duplicated(pair) & f]
  blastb[, pair := NULL]
  setnames(blastb, names[1:ncol(blastb)])
  blastb
}


cLinkFromPairs <- function(V1, V2, linked) {
  u <- unique(c(V1, V2))
  if (length(V1) != length(V2)) stop("provide vectors of equal length")

  # conversion to integer number
  if (is.integer(V1)) i1 <- match(V1, u) else i1 <- chmatch(V1, u)
  if (is.integer(V2)) i2 <- match(V2, u) else i2 <- chmatch(V2, u)

  # for each element number (index of clique), the element(s) belonging to its clique (values of gr). Originally, an element is just grouped with iself, so this integer vector is suited. But it will become an integer list
  clique <- 1:length(u)

  # all pairs of adjacent elements, with reciprocity.
  allReciprocalPairs <- cbind(c(i1[linked], i2[linked]), c(i2[linked], i1[linked]))

  # we add to these pairs all elements paired with themselves
  allReciprocalPairs <- rbind(allReciprocalPairs, cbind(clique, clique))

  # so the splitting yields, for each element, a "set" (=value of the list) of all the elements that are linked to it (and are potentially allowed to be grouped with it), including itself. Because all element numbers were used to split that list, these sets are ordered by element number (set[[n]] is the set for element number n) so we do not need to fetch a set by name matching => :-)
  set <- split(allReciprocalPairs[, 1], allReciprocalPairs[, 2])
  set <- sapply(set, unique)
  pb <- txtProgressBar(char = "*", width = 100, style = 3)
  n <- 1
  w <- which(linked)

  # for each pair of adjacent elements,
  for (i in w) {

    # gets left element number
    t1 <- i1[i]

    # right element number
    t2 <- i2[i]

    # if the two are not already grouped
    if (!t2 %in% clique[[t1]]) {

      # and if elements of their clique (which includes themselves) are in the set of allowed elements of the other	element
      if (all(clique[[t2]] %in% set[[t1]]) & all(clique[[t1]] %in% set[[t2]])) {

        # we can merge the t1 and t2 cliques into a new clique
        newGroup <- union(clique[[t1]], clique[[t2]])

        # and this new clique becomes the clique of all elements composing it (it's like a square)
        clique[newGroup] <- list(newGroup)

        # then each element of the new clique has its set become the intersection of the t1 and t2 sets, i.e, it now contains only elements that are adjacent to BOTH t1 and t2. Other elements cannot be admitted in the clique as it would break the rule. Note that all sets for elements of the same clique become the same.
        set[newGroup] <- list(intersect(set[[t2]], set[[t1]]))
      }
    }
    setTxtProgressBar(pb, n / length(w))
    n <- n + 1
  }
  cliqueElements <- unique(lapply(clique, sort))

  # makes the correspondance between a element number (vector position) and a clique number (vector value)
  cliqueNumber <- rep(1:length(cliqueElements), sapply(cliqueElements, length))
  names(cliqueNumber) <- u[unlist(cliqueElements)]
  cliqueNumber
}



# for DAVOS, converts numbers to base 10, assuming digits are in separate columns of a matrix, so that rows represent numbers
toBase10 <- function(mat, base = 3) {

  # in case mat is a vector, we convert it to a 1-row matrix (otherwise we would have dimension issues)
  if (!is.matrix(mat)) mat <- rbind(mat)

  # result numbers (a vector)
  res <- 0
  ncol <- ncol(mat)
  for (col in ncol:1) {
    res <- res + mat[, col] * base^(ncol - col)
  }
  res
}


redundantHit = function(query, subject, qStart, qEnd) {
  # flags hits nested into others on a sequence (query)
  dt = data.table(query, subject, qStart, qEnd)
  
}


bestOverlap <- function(qseqid, qStart, qEnd, sStart, sEnd, score, overlapLength = 0L) {
  pb <- txtProgressBar(char = "*", width = 100, style = 3)
  nrow <- length(qStart)
  rows <- 2:nrow
  newScaff <- c(T, qseqid[rows] != qseqid[rows - 1])
  retained <- rep(F, nrow)
  lastRetained <- 1
  for (i in 1:nrow) {
    setTxtProgressBar(pb, i / nrow)
    if (newScaff[i]) {
      lastRetained <- i
      retained[i] <- T
    } else {
      if ((qStart[i] < qEnd[lastRetained] + overlapLength) | (sStart[i] < sEnd[lastRetained] + overlapLength)) {
        if (score[i] > score[lastRetained]) {
          retained[i] <- T
          retained[lastRetained] <- F
          lastRetained <- i
        }
      }
      else {
        retained[i] <- T
        lastRetained <- i
      }
    }
  }
  retained
}


distMode <- function(x, nClasses) {
  h <- hist(x, breaks = nClasses, plot = F)
}


Ysegments <- function(start, end, scaff = 1L) {
  l <- length(start)
  scaff <- rep(scaff, length.out = l)
  dt <- data.table(scaff = rep(scaff, 2), pos = c(start, end), n = rep(1:l, 2), i = c(rep(1L, l), rep(-1L, l)))
  setorder(dt, scaff, pos)
  dt[, cov := cumsum(i)]
  dt <- dt[i == 1L, ]
  dt[match(1:l, n), cov]
}

clusterDiffs <- function(v) {
  n <- 2:length(v)
  w <- c(0L, which(v[n] != v[n - 1L]), length(v))
  n <- 2:length(w)
  rep(1:(length(w) - 1L), w[n] - w[n - 1L])
}



# fades source colors into destination colors, by a certain amount. alpha is not managed and any transparancy is removed
fadeTo <- function(source, dest, amount) {
  dest <- rep(dest, length.out = length(source))
  sourceLevels <- col2rgb(source) / 255
  destLevels <- col2rgb(dest) / 255
  delta <- destLevels - sourceLevels
  levels <- t(sourceLevels) + t(delta) * amount
  rgb(levels[, 1], levels[, 2], levels[, 3])
}


saturate <- function(col, amount) {
  amount <- rep(amount, length.out = length(col))
  amount[amount < -1] <- -1
  amount[amount > 1] <- 1
  f <- amount <= 0
  amount[f] <- amount[f] + 1
  hueSat <- rgb2hsv(col2rgb(col))
  diff <- 1 - hueSat[2, ]
  hueSat[2, f] <- hueSat[2, f] * amount[f]
  hueSat[2, !f] <- hueSat[2, !f] + amount[!f] * diff[!f]
  hueSat[hueSat > 1] <- 1
  hsv(hueSat[1, ], hueSat[2, ], hueSat[3, ])
}


luminosity <- function(col, amount) {
  amount <- rep(amount, length.out = length(col))
  amount[amount < -1] <- -1
  amount[amount > 1] <- 1
  f <- amount <= 0
  amount[f] <- amount[f] + 1
  hueSat <- rgb2hsv(col2rgb(col))
  diff <- 1 - hueSat[2, ]
  hueSat[3, f] <- hueSat[3, f] * amount[f]
  hueSat[3, !f] <- hueSat[3, !f] + amount[!f] * diff[!f]
  hueSat[hueSat > 1] <- 1
  hsv(hueSat[1, ], hueSat[3, ], hueSat[3, ])
}

xPerInch <- function() {
  abs(diff(par("usr")))[1] / par("pin")[1]
}


yPerInch <- function() {
  abs(diff(par("usr")))[3] / par("pin")[2]
}



# ratio of X unit to Y unit on the current plot (with respect to plot dimension in inches) {
xyRatio <- function() {
  xPerInch() / yPerInch()
}


freq_hist = function(x, n, breaks = "Sturges", include.lowest=T, right=T, col = "red", col2 = "blue", border = par("fg"), lwd = par("lwd"),
                     xlab = paste("Proportion of", xname), zeroLine = T, symmetric = T,
                     ylab = "Observations", negative = F, horizontal = negative, ...) {

  breaks = hist(x/n, breaks = breaks, include.lowest = include.lowest, right = right, plot = F)$breaks
  
  class = .bincode(x/n, breaks, include.lowest = include.lowest, right = right)

  dt = data.table(x, n, class)
  dt = dt[!is.na(class), .(n1 = sum(x), n2 = sum(n-x)), by = class]
  dt[,c("left","right") := .(breaks[class], breaks[class+1])]
  
  if(!horizontal) {
    dt[,c("xleft","xright", "xleft2", "xright2") := .(left, right, left, right)]
    if(negative) {
      dt[,c("ybottom","ytop", "ybottom2", "ytop2") := .(-n2, 0, 0, n1)]
    } else {
      dt[,c("ybottom","ytop", "ybottom2", "ytop2") := .(0, n1, n1, n1+n2)]
    }
    
  } else {
    if(negative) {
      dt[,c("xleft","xright", "xleft2", "xright2") := .(-n2, 0, 0, n1)]
    } else {
      dt[,c("xleft","xright", "xleft2", "xright2") := .(0, n1, n1, n1+n2)]
    }
    dt[,c("ybottom","ytop", "ybottom2", "ytop2") := .(left, right, left, right)]
    
  }
  
  rangeX = dt[,range(xleft, xright2)]
  rangeY = dt[,range(ybottom, ytop2)]
  if(negative & symmetric) {
    if(horizontal) {
      rangeX = range(rangeX, -rangeX)
    } else {
      rangeY = range(rangeY, -rangeY)
    }
  }
    
  with(dt, plot(rangeX, rangeY, type = "n", bty = "n", xlab = xlab, ylab = ylab,...))
  if(zeroLine) {
    if(horizontal) {
      abline(v = 0, lwd = 0.5)
    } else {
      abline(h = 0, lwd = 0.5)
    }
  }
  
  with(dt, rect(xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, border = NA, col = col2))
  with(dt, rect(xleft = xleft2, ybottom = ybottom2, xright = xright2, ytop = ytop2, border = NA, col = col))
  
  if(!is.na(border)) {
    with(dt, rect(xleft = xleft, ybottom = ybottom, xright = xright2, ytop = ytop2, border = border, col = NA, lwd = lwd))
  }
}
 

# links sectors between bars of a stacked bar plot
linkedBarPlot <- function(height, space = 1, width = 1, col, junCol = fadeTo(col, "white", 0.5), ...) {
  require(matrixStats)
  b <- barplot(height, space = space, width = width, col = col, plot = T, ...)

  mat <- rbind(0L, height)
  nr <- nrow(mat)
  nc <- ncol(mat)

  x <- as.vector(rbind(b - width / 2, b + width / 2))
  x <- x[-c(1, length(x))]
  x <- matrix(rep(x, nr), ncol = length(x), byrow = T)

  y <- colCumsums(mat)
  if (nc > 2) y <- y[, c(1, rep(2:(nc - 1), each = 2), nc)]

  starts <- seq(1, ncol(x), 2)
  ends <- starts + 1L
  x0 <- as.vector(x[, starts])
  x1 <- as.vector(x[, ends])
  y0 <- as.vector(y[, starts])
  y1 <- as.vector(y[, ends])
  x <- lapply(2:length(x0), function(i) c(x0[i - 1], x0[i], x1[i - 1], x1[i]))
  y <- lapply(2:length(x0), function(i) c(y0[i - 1], y0[i], y1[i], y1[i - 1]))
  if (nc > 2) remove <- seq(nr, length(x), nr) else remove <- length(x) + 1L
  m <- Map(polygon, x[-remove], y[-remove], col = junCol, border = NA)
  args <- list(...)
  if (hasArg("border")) linksCol <- args[["border"]] else linksCol <- par("col")
  if (hasArg("lwd")) linksLWD <- args[["lwd"]] else linksLWD <- par("lwd")
  if (hasArg("lty")) linksLTY <- args[["lty"]] else linksLTY <- par("lty")
  segments(x0, y0, x1, y1, col = linksCol, lwd = linksLWD, lty = linksLTY)
}

closedPolygon <- function(x, y, closeAtY = 0, ...) {
  x <- c(x[1], x, x[length(x)])
  y <- c(closeAtY, y, closeAtY)
  polygon(x, y, ...)
}

readBlastTabular <- function(file, drop = NA) {
  sel <- setdiff(1:12, drop)
  names <- c("query", "subject", "pID", "length", "mismatches", "gapOpen", "qStart", "qEnd", "sStart", "sEnd", "evalue", "score")[sel]
  fread(file, header = F, sep = "\t", select = sel, col.names = names)
}

primerBlast = function(file, maxLen = 2000L) {
  # analyses a blast file of primer pairs against a genome to find possible amplicons
  # maxLen is the maximum length of amplicons to report
  blast = readBlastTabular(file)
  blast = blast[,-c("evalue","score")]
  primerLen = blast[,max(qEnd), by = .(primer = query)] # removes hits that don't include the 3' end of the primers
  blast = blast[qEnd >= primerLen[match(query, primer), V1] - 2L]
  
  blast[,pair := gsub("F$|R$", "", query)]
  blast[,dir := lastChars(query, 1)]
  blast[,rev := sStart > sEnd]
  spl = split(blast, f = list(blast$pair, blast$subject))
  pairs = lapply(spl, findPrimerPair, maxLen = maxLen)
  rbindlist(pairs)  
}


findPrimerPair = function(blast, maxLen = 2000L) {
  if(length(unique(blast$dir)) == 1) {
    return(NULL)
  }
  pairs = blast[,allPairs(v = which(dir == "F"), v2 = which(dir == "R"))]
  pairs = cbind(blast[pairs$V1, -c("query", "dir")], blast[pairs$V2, -c("query", "dir","pair")])
  w = which(duplicated(names(pairs)))
  setnames(pairs, w, stri_c(names(pairs)[w], ".R"))
  pairs = pairs[rev != rev.R]
  pairs[,len := ifelse(rev, sStart - sStart.R, sStart.R - sStart)]
  pairs[len < maxLen & len > 0, -c("subject.R","rev","rev.R")]
}


seqtk <- function(fas, bed, out, ex = "seqtk subseq", formated = F) {
  if (missing(out)) {
    #seqs <- system(paste(ex, fas, bed), intern = T)
    seqs <- fread(cmd = paste(ex, fas, bed), sep = "\n", header = F)$V1
    if (formated) {
      f <- stri_sub(seqs, 1, 1) == ">"
      dt <- data.table(content = seqs[!f], id = cumsum(f)[!f])
      concat <- DNAStringSet(dt[, stri_flatten(content), by = id]$V1)
      names(concat) <- stri_sub(seqs[f], 2L, nchar(seqs[f]))
      return(concat)
    } else {
      return(seqs)
    }
  } else {
    system(paste(ex, fas, bed, ">", out))
  }
}

pileupParser = function(pile) {
  pile = as.data.frame(pile)
  cigars = pile[,seq(5, ncol(pile) -1, 3)]
  bases = setNames(c("A|a","C|c","G|g","T|t"),c("A","C","G","T"))
  counts = lapply(cigars, function(cigar) vapply(bases, function(base) stri_count(cigar, regex = base), integer(length(cigar))))
  res = list(dt = data.table(contig = pile[,1], pos = as.integer(pile[,2])), counts = counts)
  
}


diallelicCounts = function(counts, bases = c("A","C","G","T")) {
  # for a n-colmumn matrix of integer read counts for n bases (in principle, n = 4),
  # returns a table listing the base and the count of the two most frequent alleles
  # in case of equality in counts, the first base in the bases vector is picked
  maxes = rowMaxs(counts)
  wMaxs = rowMatches(maxes, counts)
  counts[cbind(1:nrow(counts), wMaxs)] = 0L
  secs = rowMaxs(counts)
  wSecs = rowMatches(secs, counts)
  data.table(A1 = bases[wMaxs],A2 = bases[wSecs], c1 = maxes, c2 = secs)
}


# used to convert count data from pileup parser to integer matrices
stringToIntegerMatrix <- function(string, split = ",", NA.value = 0L) {
  counts <- stri_split(string, fixed = split, simplify = T)
  storage.mode(counts) <- "integer"
  counts[is.na(counts)] <- NA.value
  counts
}


plotHitsBetween2 <- function(qStart, qEnd, sStart, sEnd, add = T) {
  qR <- range(qStart, qEnd)
  sR <- range(sStart, sEnd)
  s <- sum(sEnd - sStart + 1L)
  if (s < 0L) {
    sStart <- sR[2] - sStart + 1L
    sEnd <- sR[2] - sEnd + 1L
  }
  factor <- (qR[2] - qR[1]) / (sR[2] - sR[1])
  if (factor > 1) factor <- 1
  plot(qR, c(0, 1), type = "n", yaxt = "n", xlab = "pos", ylab = "")
  x1 <- par("usr")[1]
  x2 <- par("usr")[2]
  abline(h = 0:1)
  for (i in 1:length(qStart)) {
    polygon(c(qStart[i], (sStart[i] - sR[1]) * factor + x1, (sEnd[i] - sR[1]) * factor + x1, qEnd[i]), c(0, 1, 1, 0), col = grey(0, 0.2), border = NA)
  }
  ticks <- (x2 - x1) * 0:4 / 4 + x1
  axis(3, at = ticks, labels = round((ticks - x1) / factor + sR[1]))
  invisible(NULL)
}

closedPolygon <- function(x, y, closeAt = 0L, ...) {
  x <- c(x[1], x, tail(x, 1))
  y <- c(closeAt, y, closeAt)
  polygon(x, y, ...)
}


# stolen from poolfstat
poolFST <- function(refCounts, altCounts, poolsize) {
  npop <- ncol(refCounts)
  nsnp <- nrow(refCounts)
  YY <- as.matrix(refCounts)
  NN <- as.matrix(refCounts) + as.matrix(altCounts)
  mtrx.n_i <- matrix(poolsize, nrow = nsnp, ncol = npop, byrow = TRUE)
  C_1 <- rowSums(NN)
  C_2 <- rowSums(NN^2)
  D_2 <- rowSums(NN / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i, na.rm = TRUE)
  D_2.star <- rowSums(NN * (NN / mtrx.n_i + (mtrx.n_i - 1) / mtrx.n_i), na.rm = TRUE) / C_1
  n_c <- (C_1 - C_2 / C_1) / (D_2 - D_2.star)
  YY_alt <- NN - YY
  SSI <- rowSums(YY - YY^2 / NN + YY_alt - YY_alt^2 / NN, na.rm = TRUE)
  SSP <- rowSums(NN * ((YY / NN) - (rowSums(YY) / C_1))^2 + NN * ((YY_alt / NN) - (rowSums(YY_alt) / C_1))^2, na.rm = TRUE)
  MSI <- SSI / (C_1 - D_2)
  MSP <- SSP / (D_2 - D_2.star)
  FST <- (MSP - MSI) / (MSP + (n_c - 1) * MSI)
  #   keep <- !is.na(F_ST)
  #   F_ST_multi <- sum(MSP[keep] - MSI[keep])/sum(MSP[keep] + (n_c[keep] - 1) * MSI[keep])
  #   Q_1 <- 1 - MSI
  #   Q_2 <- 1 - MSI - (MSP - MSI)/n_c
  rslt <- data.table(FST, MSP, MSI, n_c)
}



reflowComments <- function(file, out = file) {
  rFile <- data.table(string = readLines(file))
  rFile[, commentStart := stri_locate_first(string, fixed = "#")[, 1]]
  lengths <- nchar(rFile$string)
  rFile[commentStart > 5L, extractedComment := stri_sub(string, commentStart, nchar(string))]
  rFile[commentStart > 5L, string := stri_sub(string, 1L, commentStart - 1L)]
  rFile[!grepl("[^ ]", string) & !is.na(extractedComment), c("string", "extractedComment") := .(extractedComment, NA)]
  rFile[!is.na(extractedComment), extractedComment := stri_c("\n", extractedComment)]
  concat <- rFile[, as.vector(rbind(extractedComment, string))]
  concat <- concat[!is.na(concat)]
  gsub("#", '# ', concat, fixed = T)
  gsub("#  ", '# ', concat, fixed = T)
  writeLines(concat, out)
  invisible(NULL)
}


densityPolygon <- function(x, Yscale = 1, col = NA, border = NULL, lty = par("lty"), ...) {
  d <- density(x, ...)
  closedPolygon(x = d$x, y = d$y * Yscale, col = col, border = border, lty = lty)
}


ggstyle <- function(...) {
  par(bty = "n", lend = 1, cex.axis = 0.8, tcl = -0.2, mgp = c(1.5, 0.3, 0), ...)
}


ggBackground <- function(x, y = NULL, main = "", axes = T, xaxt = par("xaxt"), yaxt = par("yaxt"), xgrid = T, ygrid = T, las = 1, bg.col = grey(0.92), grid.lwd = 1.5, subgrid.lwd = 0.5, ...) {
  ggstyle()
  if (!is.vector(x)) y <- as.matrix(x)[, 2]
  x <- as.matrix(x)[, 1]
  plot(x, y, axes = F, main = main, type = "n", las = las, ...)
  do.call(rect, c(as.list(par("usr")[c(1, 3, 2, 4)]), col = bg.col, border = NA))
  # grid(lty =1, col = "white", lwd = grid.lwd)

  ticks <- par("xaxp")
  step <- (ticks[2] - ticks[1]) / ticks[3]
  mainXTicks <- seq(ticks[1], ticks[2], length.out = ticks[3] + 1)
  if (xgrid) {
    abline(v = mainXTicks, col = "white", lwd = grid.lwd)
    secondaryTicks <- setdiff(seq(ticks[1] - step, ticks[2] + step, step / 2), mainXTicks)
    secondaryTicks <- secondaryTicks[secondaryTicks > par("usr")[1] & secondaryTicks < par("usr")[2]]
    abline(v = secondaryTicks, col = "white", lwd = subgrid.lwd)
  }
  ticks <- par("yaxp")
  step <- (ticks[2] - ticks[1]) / ticks[3]
  mainYTicks <- seq(ticks[1], ticks[2], length.out = ticks[3] + 1)
  if (ygrid) {
    abline(h = mainYTicks, col = "white", lwd = grid.lwd)
    secondaryTicks <- setdiff(seq(ticks[1] - step, ticks[2] + step, step / 2), mainXTicks)
    secondaryTicks <- secondaryTicks[secondaryTicks > par("usr")[3] & secondaryTicks < par("usr")[4]]
    abline(h = secondaryTicks, col = "white", lwd = subgrid.lwd)
  }

  if (axes) {
    if (xaxt != "n") axis(1, mainXTicks, col = NA, col.ticks = 1, col.axis = grey(0.3), las = las, lwd = grid.lwd)
    if (yaxt != "n") axis(2, mainYTicks, col = NA, col.ticks = 1, col.axis = grey(0.3), las = las, lwd = grid.lwd)
  }
}

barplotText = function(mat, text = mat, xPos, horiz = F, ...) {
  # draws labels at the middle of sectors of a stacked barplot. xPos represent values returned by barplot()
  mids = colCumsums(mat) -mat/2
  xPos = rep(xPos, each = nrow(mat))
  x = xPos
  y = mids
  if(horiz) {
    x = mids
    y = xPos
  } 
  text(x, y, text,...)
}

noTrailing0 = function(x) {
  x = as.character(x)
  f = x != "0" & grepl(".", x, fixed = T)
  x[f] = sub("0+$", "", x[f])
  x
}

angleForPoint = function(x, y, centerX, centerY) {
  dX = x - centerX
  dY = y - centerY
  atan(dY/dX) + pi * (dX < 0) 
}


rotatePoints = function(x, y, centerX, centerY, angle = 0) {
  d = xyDistance(x, y, centerX, centerY)
  alpha = angleForPoint(x, y, centerX, centerY) + angle
  nX = centerX + d * cos(alpha)
  nY = centerY + d * sin(alpha)
  data.frame(x = nX, y = nY)
}


rectangles = function(x, y, width, height = width, angle = 0, circ = F, ...) {
  # draws rectangles centered at x and y with a rotation angle in degrees
  # all arguments can be vectors
  # width and height are in inches
  # ... are arguments that can be passed to polygon() for colors, borders, line types, etc.
  ma = max(length(x), length(y), length(angle))
  x = rep(x, length.out = ma)
  y = rep(y, length.out = ma)
  angle = rep(angle, length.out = ma)
  
  if(exists("circ__") & circ) {
    angle = angleOfPoint(x) - circ__[4] + angle
    coords = convertCoords(x, y)
    x = coords$x
    y = coords$y
  }
  x = coordsToInches(x, 1)
  y = coordsToInches(y, 2)
  w = width/2
  h = height/2
  
  polyX = as.vector(rbind(x-w, x-w, x+w, x+w, NA))
  polyY = as.vector(rbind(y-h, y+h, y+h, y-h, NA))

  rotated = rotatePoints(polyX, polyY, rep(x, each = 5), rep(y, each = 5), angle = rep(angle/180*pi, each = 5))
  
  polygon(x = inchesToCoords(rotated$x, 1), y = inchesToCoords(rotated$y, 2), ...)
}

# Functions related to trees (phylo-class objects) ----------------------------------------

nestedClades <- function(tree, symmetrical = T) {
  # from a tree ("phylo" class), returns a matrix were the value is TRUE if the
  # clade at the column (the "child") is the same as, or nested in, the clade at
  # the row (the "parent"). If symmetrical is TRUE, the matrix is also TRUE if
  # the clade at the row is nested within the clade at the column. Uses the fact
  # that nodes of a phylo object are integer numbers
  
  nodes <- unique(as.vector(tree$edge))
  
  childs <- lapply(nodes, function(x) {
    subNodes(tree, x)
  })
  
  childOf <- cbind(
    parent = rep(nodes, sapply(childs, length)),
    child = unlist(childs)
  )
  
  mat <- matrix(F, max(childOf), max(childOf))
  diag(mat) <- T
  mat[childOf] <- T
  
  if (symmetrical) {
    mat[childOf[, 2:1]] <- T
  }
  
  mat
}



MRCA <- function(tree, tips) {
  # the same as getMRCA() from ape, except it returns 
  # the tip if only one tip is specified, instead of NULL
  
  if (length(tips) == 1) {
    return(match(tips, tree$tip.label))
  }
  getMRCA(tree, tips)
}


nodeDepth <- function(tree, nodes = sort(unique(as.vector(tree$edge)))) {
  # returns the age, or depth, of nodes of a timetree (nodes specified as integers)
  
  len <- node.depth.edgelength(tree)
  max(len) - len[nodes]
}



tipsForNode <- function(tree,
                        node,
                        names = F,
                        itself = F) {
  # returns all tips for a node from tree, or NULL if node is a tip (except if itself is TRUE). 
  # Returns tip names if names is TRUE, else tip numbers
  
  tips <- NULL
  n <- node
  
  while (length(n) > 0) {
    descendants <- tree$edge[tree$edge[, 1] %in% n, 2]
    tips <- c(tips, descendants[descendants <= length(tree$tip)])
    n <- setdiff(descendants, tips)
  }
  
  if (itself & length(tips) == 0) {
    tips <- node
  }
  
  if (names) {
    return(tree$tip.label[tips])
  }
  
  tips
}



tipsForNodes <- function(tree,
                         nodes,
                         names = F,
                         itself = T) {
  # returns all tips of several tree nodes, based on tipsForNode()
  # returns a data.table whose first column is the node number and second the tips
  
  tips <- lapply(nodes, function(x) {
    tipsForNode(tree, x, names = names, itself = itself)
  })
  
  nSP <- sapply(tips, length)
  dt <- data.table(node = rep(nodes, nSP), tip = unlist(tips))
  
  dt[order(-rep(nSP, nSP))]
}



subNodes <- function(tree, node) {
  # returns all subnodes and tips a tree node
  
  descendants <- node
  subNodes <- NULL
  
  repeat {
    descendants <- tree$edge[tree$edge[, 1] %in% descendants, 2]
    
    if (length(descendants) == 0) {
      break
    }
    
    subNodes <- c(subNodes, descendants)
  }
  
  subNodes
}


cladesOfAge <- function(tree,
                        age,
                        withTips = F,
                        names = T,
                        singletons = T) {
  # returns the oldest clades of age ≤ age, from a tree. None of these clades are nested within another
  
  ages <- nodeDepth(tree)
  
  # all nodes younger than the requested age
  clades <- which(ages <= age)
  
  # to obtain children of nodes
  edges <- as.data.table(tree$edge)
  
  # removes clades that are children of nodes younger than the age. This leaves only the oldest clade of age ≤ minAge
  clades <- setdiff(clades, edges[V1 %in% clades, ]$V2)
  
  if (!singletons) {
    clades <- setdiff(clades, 1:length(tree$tip.label))
  }
  
  if (!withTips) {
    return(clades)
  }
  
  tipsForNodes(tree, clades, names = names)
}


cladesForTips = function(tree, tips) {
  # return the oldest clades (nod numbers) leading to tips of a tree, 
  # ensuring that each clade contains only these tips (and no other)
  if(length(tips) == 0) {
    return(integer(0))
  }
  if(is.character(tips)) {
    tips = match(tips, tree$tip.label)
  }
  edges = data.table(tree$edge)
  test = T
  nodes = NULL
  currentNodes = tips
  while(test) {
    anc = edges[V2 %in% currentNodes]
    descs = tipsForNodes(tree, anc$V1)
    good = descs[,all(tip %in% tips), by = node]
    f = good$V1
    test = any(f)
    if(test) {
      nodes = c(nodes, anc[V1 %in% good[!f, node], V2])
      currentNodes = good[f, node]
    }
  }
  nodes
}


linesFromTree <- function(tree, angular = T) {
  # returns a list of lines (x, y coordinates of points) that can be used to draw a tree on a plot
  
  if(angular) {
    edges = data.table(tree$edge, len = tree$edge.length,  x = node.depth.edgelength(tree)[tree$edge[,1]], y= node.height(tree)[tree$edge[,2]])
    setorder(edges, V1, y)
    outter = edges[allDuplicated(V1),]
    outter = outter[!duplicated(V1) | !duplicated(V1, fromLast = T)]
    cups = split(outter, f = outter$V1)
   
     makeCup = function(cup) {
      l = cup$len
      x = cup$x
      y = cup$y
      if(l[1]> 0) {
        x = c(x[1]+l[1], x)
        y = c(y[1], y)
      }
      if(l[2]> 0) {
        n = length(x)
        x = c(x, x[n]+l[2])
        y = c(y, y[n])
      }
      data.table(x = c(x,NA), y = c(y, NA))
    }
    
    cupLines = rbindlist(lapply(cups, makeCup))
    
    inner = edges[duplicated(V1) & duplicated(V1, fromLast = T) & len > 0, .(len, x, y)]
    id = 1:nrow(inner)
    lines = rbind(inner[,.(x, y, id, n = 1)], 
                  inner[,.(x = x + len, y, id, n = 2)],
                  data.table(x = NA, y = NA, id, n = 3))
    setorder(lines, id, n)
    
    return(rbind(lines[,.(x, y)], cupLines))
    
  }
  
  # the path to each tip
  paths <- nodepath(tree)
  
  # the X coordinates of nodes/tips on the plot
  X <- node.depth.edgelength(tree)
  
  # the Y coordinates
  Y <- node.height(tree)
  
  # we scan paths to remove branches already covered by previous paths. We start with the root node
  temp <- getMRCA(tree, tree$tip.label)
  
  for (i in 1:length(paths)) {
    
    # position of the last node already covered by previous paths
    m <- max(match(temp, paths[[i]]), na.rm = T)
    
    # node connections that are specific to this path (i.e. between this node and the tip)
    nodes <- paths[[i]][m:length(paths[[i]])]
    temp <- union(temp, nodes)
    x <- X[nodes]
    
    # we get their coordinates
    y <- Y[nodes]
    if (angular) {
      # if this is an angular tree (not V-shaped), we add points for corners
      x <- rep(x, each = 2)
      x <- x[-length(x)]
      y <- rep(y, each = 2)
      y <- y[-1]
    }
    
    # we use the path list to store our coordinates
    dt = unique(data.table(x, y))
    paths[[i]] <- dt
  }
  paths[sapply(paths, nrow)>1]
  
}


groupOfTip = function(tree) {
  # returns the "group" of each tip, where a "group" comprise tips that are at same
  # location on the tree (same character states). Typically, a group would represent
  # a molecular sequence (e.g., a haplotype). The group is named after the ancestor node
  edges = data.table(tree$edge, len = tree$edge.length)
  edges = edges[V2 <= length(tree$tip.label)]
  edges[,parent := paste(V1, 1:.N)]
  edges[len == 0, parent := as.character(V1)]
  setorder(edges, V2)
  edges$parent
}



N50 = function(l, frac = 0.5, genomeSize) {
  if(!is.numeric(l)) l = nchar(l)
  l = sort(l, decreasing = T)
  cs = cumsum(as.numeric(l))
  if(missing(genomeSize)) ref = max(cs)*frac else ref = genomeSize*frac
  n = sapply(ref, function(ref) l[which.min(abs(cs-ref))])
  names(n) = frac
  n
}

meltPCRplate = function(plate, rowNames = LETTERS[1:8], colNames = 1:12) {
  nc = ncol(plate)
  setnames(plate, 1:nc, as.character(colNames[1:nc]))
  plate[,row := rowNames[1:.N]]
  melted = melt(plate, id.vars = "row", measure.vars = colNames[1:nc], variable.name = "col", value.name = "sample")
  melted[,.(well = paste(row, col, sep = ""), sample)]
}
