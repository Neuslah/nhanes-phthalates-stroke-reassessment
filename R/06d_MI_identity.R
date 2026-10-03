args<-commandArgs(TRUE)
out<-normalizePath(args[1],winslash='/',mustWork=TRUE)
s<-normalizePath(args[2],winslash='/',mustWork=TRUE)
inputs<-readRDS(file.path(out,'MI_preflight_inputs.rds'))
records<-list()
for(k in names(inputs)){
 fp<-file.path(out,'completed',k,paste0(k,'_mids_candidate.rds'))
 fit<-readRDS(fp);z<-inputs[[k]]
 stopifnot(identical(fit$data,z$data),identical(fit$method,z$method),identical(fit$predictorMatrix,z$predictorMatrix),fit$m==50,fit$iteration==20,fit$seed==20260925,is.null(fit$loggedEvents))
 records[[k]]<-data.frame(framework=k,mids_path=fp,mids_sha256=toupper(digest::digest(file=fp,algo='sha256')),m=fit$m,maxit=fit$iteration,seed=fit$seed,N=nrow(fit$data),fresh_imputation=TRUE)
}
write.csv(do.call(rbind,records),file.path(s,'mi_identity_recheck.csv'),row.names=FALSE)
cat('FRESH_MI_IDENTITY_VERIFIED\n')
