# Locked unweighted exploratory mixtures; inference intentionally retains non-design implementation.
root0<-commandArgs(TRUE)[1];source(file.path(root0,'01_corrected_helpers.R'),encoding='UTF-8')
suppressPackageStartupMessages(library(gWQS))
future::plan('sequential')
progress<-list();issues<-list();effects<-list();weights<-list();diagnostics<-list()
persist_mix<-function(){if(length(progress))wc(do.call(rbind,progress),'QC/mixture_runtime_progress_retry.csv');if(length(issues))wc(do.call(rbind,issues),'QC/mixture_issues.csv');if(length(diagnostics))wc(do.call(rbind,diagnostics),'QC/mixture_model_diagnostics.csv')}
mark<-function(key,status){progress[[length(progress)+1]]<<-data.frame(time=as.character(Sys.time()),key=key,status=status);persist_mix();cat(format(Sys.time()),' ',key,' ',status,'\n',sep='');flush.console()}
capture_mix<-function(expr,key){
 mark(key,'STARTED')
 value<-withCallingHandlers(force(expr),warning=function(w){msg<-conditionMessage(w);issues[[length(issues)+1]]<<-data.frame(key=key,type='WARNING',message=msg,time=as.character(Sys.time()));persist_mix();if(grepl('converg|singular|numer|overflow|fitted probabilities|iteration|infinite|rank.deficien',msg,ignore.case=TRUE))stop(paste('Mixture numerical warning',key,msg));invokeRestart('muffleWarning')})
 saveRDS(value,file.path(root,'models',paste0(key,'.rds')));value
}
effect<-function(b,se,p,key,sex,method,direction,N,events){
 ck(is.finite(b)&&is.finite(se)&&se>0&&is.finite(p),'Mixture invalid target inference')
 data.frame(key=key,sex=sex,method=method,direction=direction,N=N,events=events,beta=b,SE=se,OR=exp(b),CI_lower=exp(b-1.96*se),CI_upper=exp(b+1.96*se),P_full=p,analysis_role='UNWEIGHTED EXPLORATORY / NON_DESIGN_BASED / CURRENT INFERENCE RETAINED')
}
glm_qc<-function(fit,N,key){
 ck(isTRUE(fit$converged),paste('Mixture final GLM nonconvergence',key));ck(nobs(fit)==N,paste('Mixture unexpected N',key));ck(all(is.finite(coef(fit)))&&all(is.finite(vcov(fit)))&&fit$rank==ncol(model.matrix(fit)),paste('Mixture final fit singularity/variance',key))
}
safe_run({
 expvars<-paste0(met,'_ln');covars<-cv$Model3
 sample<-getdomain('HMW')$variables;sample<-sample[c('SEQN','cycle','sex','stroke',expvars,covars)]
 ck(nrow(sample)==9366&&sum(sample$sex=='Female')==4642&&sum(sample$sex=='Male')==4724&&!anyNA(sample),'Mixture sample identity')
 for(sex in c('Female','Male'))for(pos in c(TRUE,FALSE)){
  dir<-if(pos)'Positive'else'Negative';key<-paste('WQS',sex,dir,sep='_');ds<-sample[sample$sex==sex,];f<-reformulate(c('wqs',covars),response='stroke');set.seed(2026)
  fit<-capture_mix(gWQS::gwqs(formula=f,mix_name=expvars,data=ds,q=4,validation=0,b=1000,b1_pos=pos,family=binomial(link='logit'),seed=2026,signal='t2',plan_strategy='sequential'),key)
  glm_qc(fit$fit,nrow(ds),key);nonconv<-sum(fit$conv!=0,na.rm=TRUE);w<-fit$final_weights
  ck(setequal(w$mix_name,expvars)&&all(is.finite(w$mean_weight))&&all(w$mean_weight>=0)&&abs(sum(w$mean_weight)-1)<=1e-8,'WQS component weight identity')
  diagnostics[[length(diagnostics)+1]]<-data.frame(key=key,N=nrow(ds),events=sum(ds$stroke),final_converged=TRUE,bootstrap_nonconverged=nonconv,requested_bootstraps=1000,returned_bootstraps=length(fit$conv));persist_mix()
  ck(length(fit$conv)==1000&&nonconv==0,paste('WQS bootstrap convergence boundary',key))
  ct<-summary(fit)$coefficients;effects[[key]]<-effect(ct['wqs','Estimate'],ct['wqs','Std. Error'],ct['wqs','Pr(>|z|)'],key,sex,'WQS',dir,nrow(ds),sum(ds$stroke))
  weights[[key]]<-data.frame(key=key,sex=sex,direction=dir,component=w$mix_name,weight=w$mean_weight,rank=rank(-w$mean_weight,ties.method='first'))
  wc(do.call(rbind,effects),'results/WQS_effects.csv');wc(do.call(rbind,weights),'results/WQS_weights.csv');mark(key,'COMPLETE')
 }
 ck(length(effects)==4,'WQS4model completeness');writeLines('WQS = COMPLETE',file.path(root,'WQS_completion.txt'))
 qeff<-list();qweights<-list()
 for(sex in c('Female','Male'))for(boot in c(FALSE,TRUE)){
  ds<-sample[sample$sex==sex,];f<-reformulate(c(expvars,covars),response='stroke');key<-paste('QGCOMP',sex,if(boot)'Boot'else'NoBoot',sep='_');set.seed(2026)
  fit<-if(boot)capture_mix(qgcomp::qgcomp.boot(f=f,data=ds,expnms=expvars,family=binomial(),q=4,B=1000,seed=2026),key)else capture_mix(qgcomp::qgcomp.noboot(f=f,data=ds,expnms=expvars,family=binomial(),q=4),key)
  actual<-if(inherits(fit$fit,'glm'))fit$fit else fit$fit$fit;glm_qc(actual,nrow(ds),key)
  b<-as.numeric(fit$psi[1]);se<-sqrt(as.numeric(fit$var.psi[1]));qeff[[key]]<-effect(b,se,2*pnorm(-abs(b/se)),key,sex,if(boot)'qgcomp.boot'else'qgcomp.noboot','Both',nrow(ds),sum(ds$stroke))
  if(!boot){
   expand<-function(v){out<-setNames(rep(0,10),expvars);ck(all(names(v)%in%expvars),'qgcomp component identity');out[names(v)]<-v;out}
   pos<-expand(fit$pos.weights);neg<-expand(fit$neg.weights);ck(all(is.finite(c(pos,neg)))&&all(c(pos,neg)>=0)&&(sum(pos)==0||abs(sum(pos)-1)<=1e-8)&&(sum(neg)==0||abs(sum(neg)-1)<=1e-8),'qgcomp weight QC')
   qweights[[sex]]<-data.frame(sex=sex,component=expvars,positive_weight=pos,negative_weight=neg,net_weight=pos-neg);wc(do.call(rbind,qweights),'results/qgcomp_component_weights.csv')
  }
  wc(do.call(rbind,qeff),'results/qgcomp_effects.csv');mark(key,'COMPLETE')
 }
 ck(length(qeff)==4,'qgcomp4model completeness');writeLines('QGCOMP = COMPLETE',file.path(root,'qgcomp_completion.txt'))
 writeLines('CURRENT_MIXTURE_EXECUTION = COMPLETE; RH EXECUTION DISABLED',file.path(root,'mixture_completion.txt'));cat('CURRENT_MIXTURE_EXECUTION = COMPLETE\n')
})
