f_raw_data <- function(file = "",
                       rangeL = c(minL,maxL),
                       rangeJ = c(minJ,maxJ)){
  # library(RandomFields)
  library(INLA)
  library(raster)
  library(RANN)
  
  x <- read.csv(file)
  # x$a <- 1
  #Rename column to original names
  names(x) <- c("l","j","y","s","a","diff")
  x <- x %>%
    mutate(l = ifelse(a == 1, l + round(diff/7,0)*2, l)) %>% #growth fish above dam
    na.omit() %>% #get rid of NAs
    group_by(l,j,y,a) %>%  #get rid of diff column
    summarise(ns = sum(s), 
              nt = length(s)) %>% #survivors and sample size
    filter(!(a == 0 & (l < rangeL[1] | l > rangeL[2] ))) %>% #remove bad dam fish
    filter(!(a == 1 & (l < rangeL[1] | l > rangeL[2] ))) %>% #remove bad fish above dam
    filter(!(j < rangeJ[1] | j > rangeJ[2] )) %>%#remove bad days
    mutate(sl = scale(l))
  return(x)
}