.First.lib <- function(lib, pkg) {
  if(as.numeric(R.Version()$minor) < 7) {
    autoload("confint", "MASS")
  }
}
