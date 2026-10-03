# Authorized FCS only. No exposure-outcome model; no fitted coefficients exported.
out<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
suppressPackageStartupMessages(library(mice))
stopifnot(readLines(file.path(out,"initialization_status.txt"))[1]=="INITIALIZATION_PASS")
initial<-readRDS(file.path(out,"mice_initialization_only.rds"))
parent_out<-out
key_requested<-commandArgs(TRUE)[2]; stopifnot(key_requested %in% c("Female_eight_2003_2018","Male_eight_2003_2018","Female_ten_2005_2018","Male_ten_2005_2018")); initial<-initial[key_requested]
out<-file.path(parent_out,"completed",key_requested); stopifnot(!dir.exists(out)); dir.create(out)
dir.create(out,showWarnings=FALSE)
active_target<-""
active_framework<-"";calls_visited<-0L;fit_checks<-0L
wc<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na="")
event_boundary<-function(event){
 wc(data.frame(framework=active_framework,target=active_target,phase="FCS_ITERATION",conditional_calls_entered=calls_visited,event=paste(event,collapse="; "),status="STOP_BEFORE_ASSOCIATION"),"blocking_iteration_event.csv")
 stop(structure(list(message="MI_TECHNICAL_DECISION_BOUNDARY",call=NULL),class=c("mi_preflight_boundary","error","condition")))
}
tick<-function(target){
 active_target<<-target
 if(file.exists(file.path(out,"STOP_REQUEST")))event_boundary("OPERATOR_TECHNICAL_STOP")
 calls_visited<<-calls_visited+1L
 if(calls_visited==1L||calls_visited%%45L==0L)wc(data.frame(framework=active_framework,target=active_target,conditional_calls_entered=calls_visited,fit_convergence_checks=fit_checks,time=as.character(Sys.time())),"runtime_progress.csv")
}
check_fit_convergence<-function(ok,kind){fit_checks<<-fit_checks+1L;if(!isTRUE(ok))event_boundary(paste(kind,"CONVERGENCE_FAILURE"))}
# Guards examine only variable names, boolean convergence and allowed technical logs.
invisible(trace("updateLog",where=asNamespace("mice"),tracer=quote(get("event_boundary",.GlobalEnv)(out)),print=FALSE))
invisible(trace("remove.lindep",where=asNamespace("mice"),exit=quote({k<-returnValue();if(length(k)&&any(!k))get("event_boundary",.GlobalEnv)(paste("UNEXPECTED_PREDICTOR_REMOVAL",paste(colnames(x)[!k],collapse=",")))}),print=FALSE))
invisible(trace("sampler.univ",where=asNamespace("mice"),tracer=quote(get("tick",.GlobalEnv)(yname)),print=FALSE))
invisible(trace("glm.fit",where=parent.env(asNamespace("mice")),exit=quote(get("check_fit_convergence",.GlobalEnv)(returnValue()$converged,"INTERNAL_GLM")),print=FALSE))
invisible(trace("multinom",where=asNamespace("nnet"),exit=quote(get("check_fit_convergence",.GlobalEnv)(returnValue()$convergence==0,"INTERNAL_MULTINOM")),print=FALSE))
completed<-list()
result<-tryCatch({
 for(key in names(initial)){
  active_framework<-key;calls_visited<-0L;fit_checks<-0L
  cat("START_FRAMEWORK",key,"m=50 maxit=20\n")
  stopifnot(initial[[key]]$iteration==0, initial[[key]]$m==50, initial[[key]]$seed==20260925, as.character(packageVersion("mice"))=="3.19.0")
  saveRDS(initial[[key]]$lastSeedValue,file.path(out,"initial_rng_state.rds"))
  fit<-mice.mids(initial[[key]],maxit=20,printFlag=FALSE,nnet.maxit=1000)
  stopifnot(is.null(fit$loggedEvents),identical(fit$predictorMatrix,initial[[key]]$predictorMatrix),identical(fit$method,initial[[key]]$method))
  # mids contains imputed values and marginal chain summaries, not fitted conditional coefficients.
  saveRDS(fit,file.path(out,paste0(key,"_mids_candidate.rds")))
  completed[[key]]<-data.frame(framework=key,m=fit$m,iterations=fit$iteration,conditional_calls=calls_visited,fit_convergence_checks=fit_checks,logged_events=0,status="COMPUTED_PENDING_DIAGNOSTIC_REVIEW")
  wc(do.call(rbind,completed),"completed_frameworks.csv")
  cat("COMPLETED_FRAMEWORK",key,"PENDING_DIAGNOSTIC_REVIEW\n")
 }
 "ALL_FCS_COMPUTED_PENDING_DIAGNOSTIC_REVIEW"
},mi_preflight_boundary=function(e)"STOPPED_AT_MI_TECHNICAL_DECISION_BOUNDARY",error=function(e){
 wc(data.frame(framework=active_framework,phase="FCS_ITERATION",conditional_calls=calls_visited,error_class=paste(class(e),collapse=";"),message="UNSCREENED_ERROR_TEXT_NOT_EXPORTED"),"guarded_error.csv")
 "STOPPED_AT_TECHNICAL_ERROR"
})
for(n in c("updateLog","remove.lindep","sampler.univ"))invisible(untrace(n,where=asNamespace("mice")))
invisible(untrace("glm.fit",where=parent.env(asNamespace("mice"))));invisible(untrace("multinom",where=asNamespace("nnet")))
writeLines(result,file.path(out,"guarded_run_status.txt"))
capture.output(sessionInfo(),file=file.path(out,"guarded_run_sessionInfo.txt"))
cat(result,"\n")
