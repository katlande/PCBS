lin <- function(x, l){
  m <- unname(l$coefficients[2])
  b <- unname(l$coefficients[1])
  return((m*x)+b)
}
