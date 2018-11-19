checkloss <- function(y, tau){
  obj <- y*tau
  obj[y<0] <- y[y<0]*(tau-1)
  return(obj)
}
