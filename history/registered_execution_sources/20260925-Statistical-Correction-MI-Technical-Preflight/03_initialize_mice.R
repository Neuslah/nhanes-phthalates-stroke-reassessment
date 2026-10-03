# Frozen prediction specification and immediate stop on automatic variable changes.
out<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
.libPaths(c(.Library,"C:/Users/Hal Suen/AppData/Local/R/win-library/4.5"))
suppressPackageStartupMessages(library(mice))
stopifnot(as.character(packageVersion("mice"))=="3.19.0")
inputs<-readRDS(file.path(out,"MI_preflight_inputs.rds"))
wc<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na="")
active_framework<-"";preflight_phase<-"INITIALIZATION"
# Only the logging interface is instrumented; no coefficient is extracted.
event_boundary<-function(event){
 text<-paste(as.character(event),collapse="; ")
 wc(data.frame(framework=active_framework,phase=preflight_phase,event=text,status="STOP_UNEXPECTED_AUTOMATIC_MODEL_CHANGE"),"blocking_mice_event.csv")
 stop(structure(list(message="AUTOMATIC_MODEL_CHANGE_BOUNDARY",call=NULL),class=c("mi_preflight_boundary","error","condition")))
}
invisible(trace("updateLog",where=asNamespace("mice"),tracer=quote(get("event_boundary",envir=.GlobalEnv)(out)),print=FALSE))
on.exit<-NULL
initial<-list();inventory<-list()
result<-tryCatch({
 for(key in names(inputs)){
  active_framework<-key;z<-inputs[[key]]
  # Default initialization diagnostics retained; any automatic change stops before use.
  init<-mice(z$data,m=50,maxit=0,method=z$method,predictorMatrix=z$predictorMatrix,seed=20260925,printFlag=FALSE)
  stopifnot(identical(init$predictorMatrix,z$predictorMatrix),identical(init$method,z$method),is.null(init$loggedEvents))
  initial[[key]]<-init
  inventory[[key]]<-data.frame(framework=key,initialization="PASS",m=init$m,iterations=init$iteration,predictor_matrix_unchanged=TRUE,methods_unchanged=TRUE,logged_events=0)
 }
 saveRDS(initial,file.path(out,"mice_initialization_only.rds"))
 wc(do.call(rbind,inventory),"initialization_QC.csv")
 "INITIALIZATION_PASS"
},mi_preflight_boundary=function(e)"STOPPED_AT_AUTOMATIC_MODEL_CHANGE",error=function(e){
 # Do not emit unscreened package errors that could contain fitted quantities.
 wc(data.frame(framework=active_framework,phase=preflight_phase,error_class=paste(class(e),collapse=";"),status="UNSCREENED_ERROR_TEXT_NOT_EXPORTED"),"initialization_error.csv")
 "STOPPED_AT_TECHNICAL_ERROR"
})
invisible(untrace("updateLog",where=asNamespace("mice")))
loaded<-loadedNamespaces();wc(data.frame(package=loaded,version=vapply(loaded,function(p)as.character(packageVersion(p)),character(1)),path=vapply(loaded,find.package,character(1))),"loaded_package_versions.csv")
capture.output(sessionInfo(),file=file.path(out,"mice_sessionInfo.txt"))
writeLines(deparse(get("remove.lindep",asNamespace("mice"))),file.path(out,"mice_remove_lindep_source.txt"))
writeLines(deparse(get("mice.edit.setup",asNamespace("mice"))),file.path(out,"mice_setup_source.txt"))
writeLines(result,file.path(out,"initialization_status.txt"))
cat(result,"\n")
