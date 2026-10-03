# Locked unweighted exploratory mixtures; inference intentionally retains non-design implementation.
root0<-commandArgs(TRUE)[1];source(file.path(root0,'01_corrected_helpers.R'),encoding='UTF-8')
suppressPackageStartupMessages(library(gWQS))
future::plan('sequential')
progress<-list();issues<-list();effects<-list();weights<-list();diagnostics<-list()
persist_mix<-function(){if(length(progress))wc(do.call(rbind,progress),'QC/mixture_runtime_progress_retry.csv');if(length(issues))wc(do.call(rbind,issues),'QC/mixture_issues.csv');if(length(diagnostics))wc(do.call(rbind,diagnostics),'QC/mixture_model_diagnostics.csv')}
mark<-function(key,status){progress[[length(progress)+1]]<<-data.frame(time=as.character(Sys.time()),key=key,status=status);persist_mix();cat(format(Sys.time()),' ',key,' ',status,'\n',sep='');flush.console()}
capture_mix<-function(expr,key){
 ck(!file.exists(file.path(root,'models',paste0(key,'.rds'))),'Refuse overwriting completed model');mark(key,'STARTED')
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
# Inherit only the three completed, independently audited WQS fits; no refitting.
completed_keys<-c('WQS_Female_Positive','WQS_Female_Negative','WQS_Male_Positive')
old_effects<-readRDS(file.path(previous_root,'results/WQS_effects.rds'))
old_weights<-readRDS(file.path(previous_root,'results/WQS_weights.rds'))
old_diagnostics<-readRDS(file.path(previous_root,'QC/mixture_model_diagnostics.rds'))
ck(identical(as.character(old_effects$key),completed_keys),'Completed WQS keys changed')
effects<-setNames(lapply(completed_keys,function(k)old_effects[old_effects$key==k,,drop=FALSE]),completed_keys)
weights<-setNames(lapply(completed_keys,function(k)old_weights[old_weights$key==k,,drop=FALSE]),completed_keys)
diagnostics<-lapply(seq_len(nrow(old_diagnostics)),function(i)old_diagnostics[i,,drop=FALSE])
preservation<-read.csv(file.path(root,'QC/resume_preservation_preflight.csv'))
ck(all(vapply(seq_len(nrow(preservation)),function(i)hash(preservation$path[i])==preservation$expected_sha256[i],logical(1))),'Resume preservation hash gate failed')
 # Reuse the pre-fit frozen sample and partitions; no regeneration or overwrite.
 sample<-readRDS(file.path(previous_root,'mixtures/mixture_analysis_sample_2005_2018.rds'))
 expvars<-paste0(met,'_ln');covars<-cv$Model3
 expected<-getdomain('HMW')$variables;expected<-expected[c('SEQN','cycle','sex','stroke',expvars,covars)]
 ck(identical(sample,expected),'Previously frozen mixture sample changed')
 ck(nrow(sample)==9366&&sum(sample$sex=='Female')==4642&&sum(sample$sex=='Male')==4724&&!anyNA(sample),'Mixture sample identity')
 partition<-lapply(setNames(c('Female','Male'),c('Female','Male')),function(sex)readRDS(file.path(previous_root,'mixtures',paste0(sex,'_validation_rows.rds'))))
 ph<-read.csv(file.path(previous_root,'QC/holdout_partition_identity.csv'))
 ck(all(vapply(seq_len(nrow(ph)),function(i)hash(file.path(previous_root,'mixtures',paste0(ph$sex[i],'_validation_rows.rds')))==ph$sha256[i],logical(1))),'Holdout frozen hash changed')
 ck(!any(file.exists(file.path(root,'models',paste0(c('WQS_Female_Positive','WQS_Female_Negative','WQS_Male_Positive','WQS_Male_Negative'),'.rds')))),'Existing WQS fit must be preserved')
 for(sex in c('Female','Male'))for(pos in c(TRUE,FALSE)){
  dir<-if(pos)'Positive'else'Negative';key<-paste('WQS',sex,dir,sep='_');if(key %in% completed_keys)next;ds<-sample[sample$sex==sex,];f<-reformulate(c('wqs',covars),response='stroke');set.seed(2026)
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
 rheff<-list();rhweights<-list();rhsummary<-list()
 for(sex in c('Female','Male'))for(pos in c(TRUE,FALSE)){
  dir<-if(pos)'Positive'else'Negative';key<-paste('WQS_RH',sex,dir,sep='_');ds<-sample[sample$sex==sex,];f<-reformulate(c('wqs',covars),response='stroke');set.seed(2026)
  fit<-capture_mix(gWQS::gwqs(formula=f,data=ds,mix_name=expvars,rh=100,b=100,b1_pos=pos,b_constr=FALSE,q=4,validation=.6,validation_rows=partition[[sex]],family=binomial(link='logit'),signal='t2',rs=FALSE,seed=2026,plan_strategy='sequential'),key)
  ck(length(fit$gwqslist)==100,paste('Incomplete repeated holdouts',key))
  map<-vapply(fit$gwqslist,function(ch){ii<-which(vapply(partition[[sex]],function(v)identical(as.logical(ch$validation_rows),v),logical(1)));if(length(ii)==1)ii else NA_integer_},integer(1));ck(!anyNA(map)&&!anyDuplicated(map)&&setequal(map,seq_len(100)),'Holdout identity mapping')
  for(j in seq_len(100)){
   ch<-fit$gwqslist[[match(j,map)]];ct<-coef(summary(ch$fit));expectedN<-sum(partition[[sex]][[j]]);glm_qc(ch$fit,expectedN,paste(key,j));nonconv<-sum(ch$conv!=0,na.rm=TRUE)
   ck(length(ch$conv)==100&&nonconv==0,paste('Repeated-holdout bootstrap convergence boundary',key,j))
   row<-effect(ct['wqs','Estimate'],ct['wqs','Std. Error'],ct['wqs',grep('Pr',colnames(ct),value=TRUE)[1]],key,sex,'WQS_RH',dir,expectedN,sum(ds$stroke[partition[[sex]][[j]]]))
   row$holdout<-j;rheff[[length(rheff)+1]]<-row;w<-ch$final_weights;ck(setequal(w$mix_name,expvars)&&all(is.finite(w$mean_weight))&&all(w$mean_weight>=0)&&abs(sum(w$mean_weight)-1)<=1e-8,'RH component QC')
   rhweights[[length(rhweights)+1]]<-data.frame(key=key,sex=sex,direction=dir,holdout=j,component=w$mix_name,weight=w$mean_weight,rank=rank(-w$mean_weight,ties.method='first'))
  }
  wc(do.call(rbind,rheff),'results/WQS_RH_effects_long.csv');wc(do.call(rbind,rhweights),'results/WQS_RH_weights_long.csv')
  e<-do.call(rbind,rheff);e<-e[e$key==key,];rhsummary[[key]]<-data.frame(key=key,sex=sex,direction=dir,successful_holdouts=nrow(e),OR_median=median(e$OR),OR_p025=quantile(e$OR,.025,type=7,names=FALSE),OR_p975=quantile(e$OR,.975,type=7,names=FALSE),beta_median=median(e$beta),beta_mean=mean(e$beta),beta_SD=sd(e$beta),proportion_OR_lt1=mean(e$OR<1),proportion_OR_gt1=mean(e$OR>1),proportion_P_lt005_descriptive=mean(e$P_full<.05))
  wc(do.call(rbind,rhsummary),'results/WQS_RH_effect_summary.csv');mark(key,'COMPLETE')
 }
 w<-do.call(rbind,rhweights);stability<-do.call(rbind,lapply(split(w,paste(w$key,w$component)),function(z)data.frame(key=z$key[1],sex=z$sex[1],direction=z$direction[1],component=z$component[1],N=nrow(z),weight_mean=mean(z$weight),weight_median=median(z$weight),weight_SD=sd(z$weight),weight_p025=quantile(z$weight,.025,type=7,names=FALSE),weight_p975=quantile(z$weight,.975,type=7,names=FALSE),above_010_count=sum(z$weight>.1),top1_count=sum(z$rank==1),top3_count=sum(z$rank<=3),rank_median=median(z$rank),rank_min=min(z$rank),rank_max=max(z$rank))))
 wc(stability,'results/WQS_RH_component_stability.csv');ck(length(rheff)==400&&nrow(stability)==40,'RH completeness')
 writeLines('MIXTURE_EXECUTION = COMPLETE',file.path(root,'mixture_completion.txt'));cat('MIXTURE_EXECUTION = COMPLETE\n')
})
