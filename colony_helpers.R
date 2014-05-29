# Downloaded from:
# https://github.com/randomexecutable
# MIT license
#

# Searches for column names starting with C and followed by numbers
# with dot, which is again followed by numbers
# e.g. column name C3.2 would be found
#
ColonyCoefficientNames <- function(names) {
  return(names[grep("C\\d+[.]\\d+", names)])
}

# Returns the number of technical replicates per dilution level
# in a vector with dilution levels as names
# e.g.
# C4 C5  <- dilution levels (10^4, 10^5)
# 4  4   <- number of technical replicates in the data 
# file (C4.1-C4.4, C5.1-C5.4)
#
ColonyCoefficientTechrVector <- function(coefs) {
  res <- numeric()
  for(i in seq_along(coefs)) {
    split.mark <- grep("[.]", strsplit(coefs,"")[[i]])
    dilution.level <- substr(coefs[[i]], 1, (split.mark-1))
    if (is.null(names(res)) || !((dilution.level %in% names(res)))) {
      res <- append(res, 1)
      names(res)[length(res)] <- dilution.level
    } else {
      res[dilution.level] <- res[dilution.level] + 1
    }
  }
  return(res)
}

# Calculates the average from technical replicates 
# and selects the best candidate for the further
# calculations
#
ColonySelectBestAvg <- function(data, coefs.vect) {
  data['selected'] <- NA
  for(i in seq_along(names(coefs.vect))) {
    treplicates <- character(0)
    dilution.label <- names(coefs.vect[i])
    number.treplicates <- as.numeric(coefs.vect[dilution.label])
    for(j in seq_along(1:number.treplicates)) {
      current <- paste(names(coefs.vect)[i], as.character(j), sep=".")
      treplicates <- append(treplicates, current) # collects all labels for one dilution
    }
    data <- cbind(data, rowMeans(data[treplicates] * data$Coefficient, na.rm=TRUE))
    names(data)[length(data)] <- dilution.label   
    index <- is.na(data[['selected']])
    data[['selected']][index] <- dilution.label
    tempcol <- which(colnames(data) == dilution.label)
    index <- (data[[tempcol]] > 0)
    data[['selected']][index] <- dilution.label
  }  
  return(data)
}

# Calculates the absolute number of colonies 
# per sample based on the dilution level.
# Will return the calculations inside a list
# with 'Samples' and 'Controls' labels.
#
ColonyCalcAbsolute <- function(data) {
  dilution.levels <- data[['selected']]
  controls <- numeric(0)
  samples <- numeric(0)
  for(i in seq_along(1:length(dilution.levels))) {
    ifelse(grepl("[C]", data[i, 'Sample']), controls <- append(controls, data[i, dilution.levels[i]] * 10^as.numeric(substring(dilution.levels[i], 2))), samples <- append(samples, data[i, dilution.levels[i]] * 10^as.numeric(substring(dilution.levels[i], 2))))
  }
  result <- list(controls, samples)
  names(result) <- c("Controls", "Samples")
  return(result)
}

# Wraps up the other functions
#
ColonyCalculateBReplicates <- function(data) {
  coef.names <- ColonyCoefficientNames(names(data))
  technical.replicates <- ColonyCoefficientTechrVector(coef.names)
  avgs <- ColonySelectBestAvg(data, technical.replicates)  # selection works but could be better
  return(ColonyCalcAbsolute(avgs))
}

