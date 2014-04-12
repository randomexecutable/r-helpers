# https://github.com/randomexecutable
# MIT license
#

colony.coefficient.names <- function(names) {
	colony.coefs <- names[grep("C\\d+[.]\\d+", names)]
	return(colony.coefs)
}

# Returns the number of technical replicates per dilution level
# in a vector with dilution levels as names
# e.g.
# C4 C5  <- dilution levels (10^4, 10^5)
# 4  4   <- number of technical replicates in the data file (C4.1-C4.4, C5.1-C5.4)
#
colony.coefficient.techr.vector <- function(coefs) {
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


