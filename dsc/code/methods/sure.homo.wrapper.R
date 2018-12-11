# I'm assuming that this code is being run from the "dsc" directory
# inside this git repository; that is, when this code is called, the
# working directory is the "dsc" directory.
sure.homo.wrapper = function (input, args) {
  in.file  <- "results/temp/ml_in.txt"
  out.file <- "results/temp/ml_out.csv"
  cat("     - Writing data to file.\n")
  write(input$x,in.file,ncolumns = length(input$x))
  cat("     - Running recbams from WavLab in MATLAB.\n")
  system(paste(matlab.exec,"-r \"run('code/methods/sure.m')\""),
         ignore.stdout = TRUE)
  while (!file.exists(out.file))
    Sys.sleep(5)
  cat("     - Loading recbams estimate from file.\n")
  mu.est <- as.vector(read.csv(out.file,header = FALSE))
  file.remove(c(in.file,out.file))
  return(mu.est)
} 
