args<-commandArgs(TRUE)
a<-normalizePath(args[1],winslash='/',mustWork=TRUE)
out<-normalizePath(args[2],winslash='/',mustWork=TRUE)
d<-readRDS(file.path(a,'phase_a_data_candidate.rds'))
z<-list()
for(fr in c('main','common'))for(sx in c('Overall','Female','Male')){
 keep<-if(fr=='main')d$domain_main else d$domain_common
 if(sx!='Overall')keep<-keep & d$sex %in% sx
 x<-d[keep,]
 z[[length(z)+1]]<-data.frame(framework=fr,sex=sx,N=nrow(x),events=sum(x$stroke==1))
 if(sx!='Overall')writeLines(as.character(x$SEQN),file.path(out,paste0(sx,'_',if(fr=='main')'main'else'MCNP','_SEQN.txt')))
}
write.csv(do.call(rbind,z),file.path(out,'core_counts.csv'),row.names=FALSE)
for(div in c(8,7)){
 keep<-d$base_ok & (div==8 | d$cycle!='2003_2004')
 writeLines(sprintf('%.0f',d$SEQN[keep]),file.path(out,paste0('base',div,'_SEQN.txt')))
}
for(ex in c('MCNP','MCOP'))for(sx in c('Female','Male')){
 keep<-d[[paste0('domain_',ex)]] & d$sex %in% sx
 writeLines(sprintf('%.0f',d$SEQN[keep]),file.path(out,paste0(sx,'_',ex,'_SEQN.txt')))
}
cat('CLEAN_CORE_COUNTS_EXPORTED\n')
