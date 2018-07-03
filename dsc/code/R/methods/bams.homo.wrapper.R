bams.homo.wrapper = function(input, args) {
  write(input$x, "input/ml_in.txt", ncolumns = length(input$x))
  system(paste(matlab.exec,"-r \"run('methods/bams_matlab.m')\""),
         ignore.stdout = TRUE)
  if (Sys.info()["sysname"] == "Windows") {
    while (!file.exists("input/ml_out.csv")) {
      Sys.sleep(5)
    }
  }
  mu.est = as.vector(read.csv("input/ml_out.csv", header = FALSE))
  system("rm input/ml_in.txt input/ml_out.csv")
  return(mu.est)
} 
