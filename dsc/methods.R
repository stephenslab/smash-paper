source_dir("methods")

parse.period.delimited.name <- function(name) {
  return(strsplit(name, ".", fixed = TRUE)[[1]])
}

build.period.delimited.name <- function(words) {
  return(paste(words, collapse = "."))
}

suffix <- function(name) {
  words <- parse.period.delimited.name(name)
  return(words[length(words)])
}

remove.suffix <- function(name) {
  words <- parse.period.delimited.name(name)
  return(build.period.delimited.name(words[1:(length(words) - 1)]))
}

fxn.wrapper <- function(name, named.method) {
    if ((suffix(name) == "s8") || (suffix(name) == "haar")) {
        if (named.method) {
            stem <- remove.suffix(remove.suffix(name))
        } else {
            stem <- remove.suffix(name)
        }
    } else {
        stem <- name
    }
    return(parse(text = paste(stem, "wrapper", sep = ".")))
}

arg.list <- function(name, named.method, blank.arglist) {
    args <- list()

    if (!blank.arglist) {
        if (named.method) {
            args$method <- suffix(remove.suffix(name))
        }

        if (suffix(name) == "haar") {
            args$filter.number <- 1
            args$family <- "DaubExPhase"
        } else {
            
            # This covers a suffix of "s8" but also the case
            # "smash.jash".
            args$filter.number <- 8
            args$family <- "DaubLeAsymm"
        }
    }

    return(args)
}

add.method.by.flags <- function(row) {
  method.name <- as.character(row$name)
  add_method(dsc_smash,
             name = method.name,
             fn   = eval(fxn.wrapper(method.name,row$named.method)),
             args = arg.list(method.name,row$named.method,row$blank.arglist))
}

# Define methods and flags used to name the function wrappers and
# argument lists. The first seven of this list use a blank argument
# list. The tithresh.rmad and tithresh.smash methods put the
# penultimate word in the argument list instead of using separate
# wrappers.
methods <- data.frame(list(name = c("ebayesthresh",
                                    "bams.homo.s8",
                                    "blockjs.homo.s8",
                                    "neighblock.homo.s8",
                                    "sure.homo.s8",
                                    "postmean.homo.s8",
                                    "tithresh.homo.s8",
                                    "smash.homo.haar",
                                    "smash.homo.s8",
                                    "smash.haar",
                                    "smash.s8",
                                    "smash.jash",
                                    "tithresh.rmad.haar",
                                    "tithresh.rmad.s8",
                                    "tithresh.smash.haar",
                                    "tithresh.smash.s8",
                                    "smash.true.haar",
                                    "smash.true.s8",
                                    "tithresh.true.haar",
                                    "tithresh.true.s8")))

methods$named.method <- FALSE
methods$named.method[13:16] <- TRUE

methods$blank.arglist <- FALSE
methods$blank.arglist[1:7] <- TRUE

# Add the methods to the DSC.
by(methods,1:nrow(methods),add.method.by.flags)
