#Super annoying. This keeps TMB from recompile if there's been changes to the cpp
f_cpp <- function(version){
  cppFile <- paste0("spde_aniso_wt_",version,".cpp")
  cppLog <- paste0("lenWt_",version,".log")
  cppVersion <- paste0("spde_aniso_wt_",version)
  if(length(list.files()[grep(cppLog,list.files())])==0){
    write.table(file.info(cppFile)$mtime,cppLog, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  cppLoc <- list.files()[grep(cppFile,list.files())]
  cppTime <- file.info(cppLoc)$mtime
  mTime <- apply(read.table(cppLog,header = FALSE),1,paste,collapse=" ")
  write.table(file.info(cppFile)$mtime,cppLog, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  if(as.character(cppTime)!=as.character(mTime)){
    try(dyn.unload(cppVersion))
    compile(paste0(cppVersion,".cpp"), CPPFLAGS="-Wno-ignored-attributes")
  }
  dyn.load(cppVersion)
}

