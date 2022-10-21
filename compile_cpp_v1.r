#Super annoying. This keeps TMB from recompile if there's been changes to the cpp
cppLog <- paste0("lenWt_",version,".log")
if(length(list.files()[grep(cppLog,list.files())])==0){
  write.table(file.info(cppFile)$mtime,cppLog, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
cppLoc <- list.files()[grep(cppFile,list.files())]
cppTime <- file.info(cppLoc)$mtime
mTime <- apply(read.table(cppLog,header = FALSE),1,paste,collapse=" ")
write.table(file.info(cppFile)$mtime,cppLog, quote = FALSE, row.names = FALSE, col.names = FALSE)

if(as.character(cppTime)!=as.character(mTime)){
  source("compile_spde_aniso_wt_v2.r")
}
dyn.load(paste0("spde_aniso_wt_",version))

