f_management_projection <-function(df = df,
                                   m_proj = m_proj,
                                   rangeL = rangeL,
                                   rangeJ = rangeL){
  grid <- expand.grid('l' = m_proj$l,
                             'j' = sort(m_proj$j))
  
  # print(nrow(df))
  tmp <- array(NA , dim = c(nrow(df),nrow(grid),ncol(grid)+1))
  # print(dim(tmp))
  dimnames(tmp)[[3]] <- c("l", "j", "s_i")
  
  for(i in 1:nrow(grid)){
    for(j in c('l','j')){
      if(j=='l'){
        tmp[,i,j] <- as.data.frame(df[,j])[,1] + grid[i,j] #pain in the ass
        tmp[tmp[,i,j]>max(rangeL),i,j] <- max(rangeL)
        tmp[tmp[,i,j]<min(rangeL),i,j] <- min(rangeL)
      }
      if(j=='j'){
        tmp[,i,j] <- as.data.frame(df[,j])[,1] + grid[i,j] #pain in the ass
        tmp[tmp[,i,j]>max(rangeJ),i,j] <- max(rangeJ)
        tmp[tmp[,i,j]<min(rangeJ),i,j] <- min(rangeJ)
      }
    }
  }
  return(list(proj = tmp,
         grid = grid))
}
