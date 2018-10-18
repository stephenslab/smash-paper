# Load the MACS peaks data from a text file, returning all peaks
# within the selected range of base-pair positions.
read.macs.peaks <- function (file, min.pos, max.pos) {

  # Read the MACS data from the text file.
  peaks <- read.table(file)

  # Order the peaks by the starting base-pair position.
  i           <- order(peaks[,2])
  peaks       <- peaks[i,]

  # Output the start and end positions of the MACS peaks, within the
  # desired range of base-pair positions.
  i <- which(peaks[,2] >= min.pos & peaks[,2] <= max.pos)
  return(list(start = peaks[i,2],
              end    = peaks[i,3]))
}

# Create a plot showing. This is a very *ad hoc* implementation that
# will only work for the specific data set that was analyzed in the
# "chipseq.Rmd" example.
create.chipseq.plot <- function (pos, counts, smash.est, peaks, nbreaks) {

  # Sum the read counts in equally sized "bins" (small chromosome
  # intervals).
  breaks  <- seq(min(pos),max(pos),length.out = nbreaks)
  x       <- cut(pos,breaks)
  n       <- table(x,factor(counts))

  # Create a data frame containing the mean base-pair position and
  # total count for each bin.
  pdat <- data.frame(position = rep(breaks[-nbreaks],times = 8) +
                                  diff(breaks[1:2])/2,
                     count    = rep(0:7,each = nbreaks - 1),
                     n        = sqrt(as.vector(n)),
                     is.zero  = as.vector(n) == 0)
  return(ggplot(pdat,aes_string(x = "position",y = "count",size = "n",
                              fill = "is.zero")) +
         geom_point(color = "white",shape = 21) +
         scale_fill_manual(values = c("deepskyblue","white"),guide = FALSE) +
         scale_size(range = c(0.5,3.5)) +
         geom_line(data = data.frame(position = pos,y = 4*smash.est),
                   mapping = aes_string(x = "position",y = "y"),
                   color = "darkorange",inherit.aes = FALSE) +
         geom_point(data = data.frame(position = peaks,y = -1),
                    mapping = aes_string(x = "position",y = "y"),color = "red",
                    inherit.aes = FALSE,shape = 2,size = 2) +
         labs(x = "base-pair position (Mb)",
              y = "read count"))
}
