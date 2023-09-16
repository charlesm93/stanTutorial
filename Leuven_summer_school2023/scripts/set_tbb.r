
set_tbb <- function(cmdstan_path, model_path) {
  command <- 
    paste0("cp ", 
           cmdstan_path, "/stan/lib/stan_math/lib/tbb/tbb.dll",
           " ", model_path, "/tbb.dll")
  
  system(command)
}
