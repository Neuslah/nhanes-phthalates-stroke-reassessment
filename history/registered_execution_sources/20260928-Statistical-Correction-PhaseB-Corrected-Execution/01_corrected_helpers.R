options(stringsAsFactors=FALSE,digits=17)
.libPaths(c('C:/Users/Hal Suen/AppData/Local/R/win-library/4.5',.libPaths()))
root<-normalizePath(commandArgs(TRUE)[1],winslash='/',mustWork=TRUE)
ev<-dirname(root);a<-file.path(ev,'20260925-Statistical-Correction-PhaseA/resume_cholesterol_proxy')
s<-file.path(ev,'20260928-Statistical-Correction-PhaseB-Step0-Identity-Freeze')
m<-file.path(ev,'20260925-Statistical-Correction-MI-Technical-Preflight')
for(pkg in c('survey','digest','mice','mitools','rms','gWQS','qgcomp','future')){
 if(!requireNamespace(pkg,quietly=TRUE))stop('Missing locked runtime package: ',pkg)
}
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust',survey.adjust.domain.lonely=TRUE)
for(k in c('models','results','QC','logs','mixtures'))dir.create(file.path(root,k),showWarnings=FALSE)
hash<-function(p)toupper(digest::digest(file=p,algo='sha256'))
wc<-function(x,p){write.csv(x,file.path(root,p),row.names=FALSE,na='',fileEncoding='UTF-8');saveRDS(x,file.path(root,sub('\\.csv$','.rds',p)))}
ck<-function(b,msg){if(!isTRUE(b))stop(msg)}
ck(hash(file.path(root,'PHASE_B_AUTHORIZED_ANALYSIS_INVENTORY.csv'))=='856A7CEC9A9A99F7F43B194FA7438826969313A84A79F986DD6D29D7576B3A2D','Frozen inventory changed')
protected<-read.csv(file.path(root,'input_integrity_before.csv'))
ck(all(vapply(seq_len(nrow(protected)),function(i)file.exists(protected$path[i])&&hash(protected$path[i])==protected$expected_sha256[i],logical(1))),'Protected input identity changed')
versions<-do.call(rbind,lapply(c('survey','mice','mitools','rms','gWQS','qgcomp','future','digest','Matrix','nnet'),function(p)data.frame(package=p,version=as.character(packageVersion(p)),path=find.package(p))))
wc(versions,'QC/software_versions.csv');capture.output({print(R.version.string);print(RNGkind());print(.libPaths());sessionInfo()},file=file.path(root,'QC/sessionInfo.txt'))
capture.output(print(getS3method('MIcombine','default',envir=asNamespace('mitools'))),file=file.path(root,'QC/MIcombine_actual_implementation.txt'))
d<-readRDS(file.path(a,'phase_a_data_candidate.rds'))
rawmap<-c(MBP='URXMBP',MBzP='URXMZP',MECPP='URXECP',MEHHP='URXMHH',MEOHP='URXMOH',MCPP='URXMC1',MEP='URXMEP',MiBP='URXMIB',MCNP='URXCNP',MCOP='URXCOP')
met<-names(rawmap)
for(ex in met)d[[paste0(ex,'_ln')]]<-log(d[[rawmap[[ex]]]]/d$URXUCR*1000)
cv<-list(Model1=c('RIDAGEYR','race_eth_f'),Model2=c('RIDAGEYR','race_eth_f','edu_f','marital_f','pir'))
cv$Model3<-c(cv$Model2,'bmi','smoking_f','alcohol_harmonized_f','hypertension_f','diabetes_f','hyperlipidemia_f')
d$cycle_f<-factor(d$cycle);d$era<-as.integer(d$cycle %in% c('2011_2012','2013_2014','2015_2016','2017_2018'))
d$age_group_f<-factor(ifelse(d$RIDAGEYR<60,'<60','>=60'),levels=c('<60','>=60'))
d$bmi_group_f<-factor(ifelse(d$bmi<25,'<25',ifelse(d$bmi<30,'25-30','>=30')),levels=c('<25','25-30','>=30'))
d$race_subgrp<-factor(ifelse(d$race_eth_f=='Non-Hispanic White','NH White',ifelse(d$race_eth_f=='Non-Hispanic Black','NH Black',ifelse(d$race_eth_f %in% c('Mexican American','Other Hispanic'),'Hispanic','Other'))),levels=c('NH White','NH Black','Hispanic','Other'))
MW<-c(MEP=194.18,MBP=222.24,MiBP=222.24,MBzP=312.36,MEHHP=294.35,MEOHP=292.33,MECPP=308.33,MCPP=252.22,MCNP=336.38,MCOP=322.35)
for(z in list(LMW=c('MEP','MBP','MiBP'),HMW=c('MBzP','MEHHP','MEOHP','MECPP','MCNP','MCOP','MCPP'),DEHP=c('MEHHP','MEOHP','MECPP'))){
 nm<-if(identical(z,c('MEP','MBP','MiBP')))'LMW'else if(length(z)==7)'HMW'else'DEHP'
 d[[paste0('ln_sigma_',nm)]]<-log(rowSums(sapply(z,function(ex)exp(d[[paste0(ex,'_ln')]])/MW[[ex]])))
}
rownames(d)<-as.character(d$SEQN)
bases<-list()
for(div in c(8,7)){
 x<-d[d$base_ok & (div==8 | d$cycle!='2003_2004'),,drop=FALSE];x$analysis_weight<-x$ph_weight/div
 des<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~analysis_weight,data=x,nest=TRUE)
 bases[[as.character(div)]]<-des
}
design_info<-function(des){
 x<-des$variables;u<-unique(x[c('cycle','SDMVSTRA','SDMVPSU')]);st<-unique(x[c('cycle','SDMVSTRA')]);ct<-table(interaction(u$cycle,u$SDMVSTRA,drop=TRUE))
 list(N=nrow(x),events=sum(x$stroke==1),PSU=nrow(u),strata=nrow(st),df=nrow(u)-nrow(st),lonely=sum(ct==1))
}
for(div in c(8,7)){z<-design_info(bases[[as.character(div)]]);ck(z$lonely==0 && z$df==if(div==8)124 else 109,'Unexpected base design')}
getdomain<-function(ex,sex=NULL,common=FALSE){
 div<-if(common||ex %in% c('MCNP','MCOP','HMW'))7 else 8
 base<-bases[[as.character(div)]];x<-base$variables
 keep<-if(common||ex=='HMW')x$domain_common else if(ex %in% c('MCNP','MCOP'))x[[paste0('domain_',ex)]]else x$domain_main
 if(!is.null(sex))keep<-keep & x$sex==sex
 ck(!anyNA(keep),'NA domain membership');des<-base[keep,]
 if(!is.null(sex)){
  fr<-if(div==8)'main'else'MCNP';f<-file.path(s,paste0(sex,'_',fr,'_SEQN.txt'));ids<-scan(f,quiet=TRUE)
  ck(setequal(ids,des$variables$SEQN),'Step0 frozen domain membership mismatch')
 }
 z<-design_info(des);ck(z$df==if(div==8)124 else 109,'Primary domain df changed');des
}
fit_qc<-list();warn_qc<-list();row_qc<-list();result_store<-list()
persist_qc<-function(){if(length(fit_qc))wc(do.call(rbind,fit_qc),'QC/model_identity.csv');if(length(warn_qc))wc(do.call(rbind,warn_qc),'QC/model_warnings.csv');if(length(row_qc))wc(do.call(rbind,row_qc),'QC/model_rows.csv')}
fit_checked<-function(f,des,key){
 ck(!file.exists(file.path(root,'models',paste0(key,'.rds'))),'Refuse model overwrite')
 x<-des$variables;mf<-model.frame(f,data=x,na.action=na.pass);mm<-model.matrix(f,mf)
 ck(nrow(mf)==nrow(x)&&!anyNA(mf)&&all(is.finite(mm)),'Unexpected model-frame missingness/non-finite matrix')
 ck(qr(mm)$rank==ncol(mm),paste('Pre-fit singular model matrix',key))
 z<-design_info(des);ck(z$df>0,'Nonpositive domain df')
 writeLines(c('CORRECTED_EXPOSURE_STROKE_MODELS_RUN = YES','CORRECTED_OR_CI_P_Q_GENERATED = PENDING','PHASE_B_EXECUTED = YES'),file.path(root,'execution_flags.txt'))
 observed<-character();cat(format(Sys.time()),' START ',key,' N=',z$N,' df=',z$df,'\n',sep='');flush.console()
 fit<-withCallingHandlers(svyglm(f,design=des,family=quasibinomial(link='logit'),na.action=na.fail),warning=function(w){observed<<-c(observed,conditionMessage(w));invokeRestart('muffleWarning')})
 saveRDS(fit,file.path(root,'models',paste0(key,'.rds')))
 if(length(observed))warn_qc[[length(warn_qc)+1]]<<-data.frame(key=key,warning=unique(observed))
 b<-coef(fit);v<-vcov(fit)
 fit_qc[[length(fit_qc)+1]]<<-data.frame(key=key,formula=paste(deparse(f),collapse=' '),N=z$N,events=z$events,PSU=z$PSU,strata=z$strata,DOMAIN_DF=z$df,degf=degf(des),lonely=z$lonely,nobs=nobs(fit),default_residual_df=fit$df.residual,rank=fit$rank,matrix_columns=ncol(mm),converged=isTRUE(fit$converged),finite=all(is.finite(b))&&all(is.finite(v)),warning_count=length(observed))
 row_qc[[length(row_qc)+1]]<<-data.frame(key=key,SEQN=x$SEQN)
 persist_qc()
 ck(isTRUE(fit$converged),'Model failed convergence');ck(fit$rank==ncol(mm)&&length(b)==ncol(mm)&&all(is.finite(b))&&all(is.finite(v))&&all(diag(v)>0),'Singularity/non-finite variance')
 ck(nobs(fit)==nrow(x)&&identical(rownames(fit$model),rownames(x)),'Unexpected model row deletion/order')
 ck(!any(grepl('converg|singular|numer|overflow|fitted probabilities|iteration|infinite|rank.deficien',observed,ignore.case=TRUE)),paste('Numerical warning',key))
 cat(format(Sys.time()),' COMPLETE ',key,'\n',sep='');flush.console()
 list(fit=fit,info=z)
}
single<-function(obj,term,analysis,sex,ex,model){
 b<-unname(coef(obj$fit)[term]);se<-sqrt(unname(vcov(obj$fit)[term,term]));df<-obj$info$df
 ck(is.finite(b)&&is.finite(se)&&se>0,'Invalid target coefficient');p<-2*pt(-abs(b/se),df=df);crit<-qt(.975,df=df)
 data.frame(analysis=analysis,sex=sex,exposure=ex,model=model,N=obj$info$N,events=obj$info$events,DOMAIN_DF=df,beta=b,SE=se,OR=exp(b),CI_lower=exp(b-crit*se),CI_upper=exp(b+crit*se),P_full=p,P_report=if(p<.001)'<0.001'else sprintf('%.3f',p),direction=if(b>0)'POSITIVE'else if(b<0)'NEGATIVE'else'ZERO',nominal=p<.05,q=NA_real_,FDR=NA)
}
bh<-function(x){ck(nrow(x)==20&&!anyNA(x$P_full)&&!anyDuplicated(paste(x$sex,x$exposure)),'BH family membership failure');x$q<-p.adjust(x$P_full,method='BH');x$FDR<-x$q<.05;ord<-order(x$P_full);q2<-pmin(1,rev(cummin(rev(x$P_full[ord]*20/seq_len(20)))));ck(max(abs(x$q[ord]-q2))<1e-14,'Independent BH check failure');x}
safe_run<-function(expr){tryCatch(force(expr),error=function(e){persist_qc();writeLines(c('PHASE_B_EXECUTION = BLOCKED',paste0('BLOCKER = ',conditionMessage(e)),paste0('time = ',Sys.time()),'No model/specification repair performed; subsequent analyses stopped.'),file.path(root,'decision_boundary.txt'));stop(e)})}
