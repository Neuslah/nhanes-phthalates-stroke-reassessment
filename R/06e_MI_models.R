root0<-commandArgs(TRUE)[1];source(file.path(root0,'01_corrected_helpers.R'),encoding='UTF-8')
persist_qc<-function(){if(length(fit_qc))wc(do.call(rbind,fit_qc),'QC/MI_model_identity.csv');if(length(warn_qc))wc(do.call(rbind,warn_qc),'QC/MI_model_warnings.csv')}
safe_run({
 identity<-read.csv(file.path(s,'mi_identity_recheck.csv'));inputs<-readRDS(file.path(m,'MI_preflight_inputs.rds'))
 targets<-c('edu_f','marital_f','pir','bmi','smoking_f','alcohol_harmonized_f','hypertension_f','diabetes_f','hyperlipidemia_f')
 pooled_rows<-list();framework_rows<-list();all_fits<-list()
 for(k in names(inputs)){
  z<-inputs[[k]];id<-identity[identity$framework==k,];fp<-id$mids_path;ck(length(fp)==1&&hash(fp)==id$mids_sha256,'Completed mids hash changed')
  fitmi<-readRDS(fp);ck(identical(fitmi$data,z$data)&&fitmi$m==50&&fitmi$iteration==20&&is.null(fitmi$loggedEvents),'Completed MI object identity')
  sex<-strsplit(k,'_',fixed=TRUE)[[1]][1];div<-z$weight_divisor;exs<-if(div==8)met[1:8]else met[9:10];base<-bases[[as.character(div)]]
  ids<-z$SEQNs;ii<-match(ids,base$variables$SEQN);ck(!anyNA(ii)&&!anyDuplicated(ids),'MI participant join identity')
  estimates<-setNames(lapply(exs,function(e)vector('list',50)),exs);variances<-estimates;perfit<-list()
  for(j in seq_len(50)){
   completed<-mice::complete(fitmi,action=j);ck(nrow(completed)==length(ids),'Completed MI N mismatch'); attr(completed,'row.names')<-attr(z$data,'row.names')
   fixed<-setdiff(names(completed),targets);ck(identical(completed[fixed],z$data[fixed]),'Fixed MI exposure/outcome/design changed')
   for(tg in targets){obs<-!is.na(z$data[[tg]]);ck(identical(completed[[tg]][obs],z$data[[tg]][obs]),paste('Observed/structural proxy overwritten',tg))}
   des<-base;x<-des$variables
   for(tg in targets)x[[tg]][ii]<-completed[[tg]]
   des$variables<-x;keep<-x$SEQN %in% ids;des<-des[keep,];x<-des$variables
   ck(setequal(x$SEQN,ids)&&all(x$sex==sex),'MI framework membership/sex changed')
   info<-design_info(des);ck(info$df==if(div==8)124 else 109,'MI design df mismatch')
   if(j==1)framework_rows[[k]]<-data.frame(framework=k,N=info$N,PSU=info$PSU,strata=info$strata,df_complete=info$df,lonely=info$lonely,mids_sha256=hash(fp),m=50,maxit=20,seed=fitmi$seed,imputation_rerun=TRUE)
   for(ex in exs){
    key<-paste('MI',k,ex,j,sep='_');e<-paste0(ex,'_ln');f<-reformulate(c(e,cv$Model3),response='stroke');mf<-model.frame(f,x,na.action=na.pass);mm<-model.matrix(f,mf)
    ck(!anyNA(mf)&&all(is.finite(mm))&&qr(mm)$rank==ncol(mm),'MI analysis-frame missingness/singularity')
    warnings<-character();fitted<-withCallingHandlers(svyglm(f,design=des,family=quasibinomial(link='logit'),na.action=na.fail),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    b<-coef(fitted);v<-vcov(fitted);rows_same<-identical(rownames(fitted$model),rownames(x))
    record<-list(key=key,formula=f,beta=b,vcov=v,N=info$N,events=info$events,DOMAIN_DF=info$df,member_SEQN_sha256=digest::digest(sort(ids),algo='sha256'),converged=fitted$converged,rank=fitted$rank,warnings=warnings,model_rows_identical=rows_same)
    perfit[[key]]<-record
    fit_qc[[length(fit_qc)+1]]<-data.frame(key=key,framework=k,exposure=ex,imputation=j,N=info$N,events=info$events,PSU=info$PSU,strata=info$strata,DOMAIN_DF=info$df,lonely=info$lonely,nobs=nobs(fitted),rank=fitted$rank,matrix_columns=ncol(mm),converged=isTRUE(fitted$converged),finite=all(is.finite(b))&&all(is.finite(v)),rows_identical=rows_same,warning_count=length(warnings))
    if(length(warnings))warn_qc[[length(warn_qc)+1]]<-data.frame(key=key,warning=unique(warnings))
    saveRDS(perfit,file.path(root,'models',paste0('MI_',k,'_per_imputation_evidence.rds')))
    ck(isTRUE(fitted$converged)&&fitted$rank==ncol(mm)&&all(is.finite(b))&&all(is.finite(v))&&all(diag(v)>0),'MI association convergence/singularity/variance failure')
    ck(nobs(fitted)==nrow(x)&&rows_same,'MI unexpected model row deletion')
    ck(!any(grepl('converg|singular|numer|overflow|fitted probabilities|iteration|infinite|rank.deficien',warnings,ignore.case=TRUE)),'MI association numerical warning')
    estimates[[ex]][[j]]<-setNames(unname(b[e]),e);variances[[ex]][[j]]<-matrix(v[e,e],1,1,dimnames=list(e,e))
   }
   persist_qc();cat(format(Sys.time()),' MI ',k,' completed imputation ',j,'/50\n',sep='');flush.console()
  }
  for(ex in exs){
   cmb<-mitools::MIcombine(estimates[[ex]],variances[[ex]],df.complete=info$df)
   betas<-vapply(estimates[[ex]],function(x)x[[1]],numeric(1));vv<-vapply(variances[[ex]],function(x)x[1,1],numeric(1));bar<-mean(betas);within<-mean(vv);between<-var(betas);total<-within+(1+1/50)*between;r<-(1+1/50)*between/within
   olddf<-(50-1)*(1+1/r)^2;obsdf<-((info$df+1)/(info$df+3))*info$df*within/(within+between);newdf<-1/(1/obsdf+1/olddf)
   ck(isTRUE(all.equal(as.numeric(cmb$coefficients),bar,tolerance=1e-12))&&isTRUE(all.equal(as.numeric(cmb$variance),total,tolerance=1e-12))&&isTRUE(all.equal(as.numeric(cmb$df),newdf,tolerance=1e-12)),'Rubin/finite MIcombine independent check failure')
   ck(is.finite(total)&&total>0&&is.finite(newdf)&&newdf>0,'MI pooled inference failure')
   se<-sqrt(total);p<-2*pt(-abs(bar/se),df=newdf);crit<-qt(.975,df=newdf)
   pooled_rows[[length(pooled_rows)+1]]<-data.frame(analysis='MI_Model3',sex=sex,exposure=ex,model='Model3',framework=k,N=info$N,events=info$events,df_complete=info$df,df_pooled=newdf,m=50,beta=bar,within_variance=within,between_variance=between,total_variance=total,SE=se,OR=exp(bar),CI_lower=exp(bar-crit*se),CI_upper=exp(bar+crit*se),P_full=p,P_report=if(p<.001)'<0.001'else sprintf('%.3f',p),fraction_missing_information=as.numeric(cmb$missinfo),direction=if(bar>0)'POSITIVE'else if(bar<0)'NEGATIVE'else'ZERO',nominal=p<.05,q=NA_real_,FDR=NA)
   saveRDS(cmb,file.path(root,'models',paste0('MI_pooled_',sex,'_',ex,'.rds')))
  }
  wc(do.call(rbind,pooled_rows),'results/MI_Model3_partial.csv');wc(do.call(rbind,framework_rows),'QC/MI_framework_identity.csv')
 }
 out<-bh(do.call(rbind,pooled_rows));wc(out,'results/MI_Model3.csv');ck(length(fit_qc)==1000,'MI1000 association fits completeness')
 cc<-readRDS(file.path(root,'results/primary_all_models.rds'));cc<-cc[cc$model=='Model3',];cmp<-merge(cc,out,by=c('sex','exposure'),suffixes=c('_CC','_MI'))
 cmp$direction_match<-cmp$direction_CC==cmp$direction_MI;cmp$nominal_match<-cmp$nominal_CC==cmp$nominal_MI;cmp$FDR_match<-cmp$FDR_CC==cmp$FDR_MI
 wc(cmp,'results/CC_MI_mechanical_comparison.csv');writeLines('MI_ASSOCIATION_EXECUTION = COMPLETE',file.path(root,'MI_completion.txt'));cat('MI_ASSOCIATION_EXECUTION = COMPLETE\n')
})
