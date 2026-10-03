root0<-commandArgs(TRUE)[1];source(file.path(root0,'01_corrected_helpers.R'),encoding='UTF-8')
safe_run({
 rows<-list()
 for(sex in c('Female','Male'))for(ex in met)for(mod in names(cv)){
  key<-paste('primary',sex,ex,mod,sep='_');f<-reformulate(c(paste0(ex,'_ln'),cv[[mod]]),response='stroke')
  fit<-fit_checked(f,getdomain(ex,sex),key)
  rows[[length(rows)+1]]<-single(fit,paste0(ex,'_ln'),'primary',sex,ex,mod)
  wc(do.call(rbind,rows),'results/primary_partial.csv')
 }
 ans<-do.call(rbind,rows);ii<-ans$model=='Model3';ans[ii,]<-bh(ans[ii,])
 ck(nrow(ans)==60,'Primary completeness');wc(ans,'results/primary_all_models.csv')
 writeLines(c('CORRECTED_EXPOSURE_STROKE_MODELS_RUN = YES','CORRECTED_OR_CI_P_Q_GENERATED = YES','PHASE_B_EXECUTED = YES'),file.path(root,'execution_flags.txt'))
 writeLines('PRIMARY_EXECUTION = COMPLETE',file.path(root,'primary_completion.txt'))
 cat('PRIMARY_EXECUTION = COMPLETE; FDR_SURVIVORS=',sum(ans$FDR[ii]),'/20\n',sep='')
})
