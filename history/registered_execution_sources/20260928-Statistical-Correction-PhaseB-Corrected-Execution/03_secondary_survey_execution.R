root0<-commandArgs(TRUE)[1];source(file.path(root0,'01_corrected_helpers.R'),encoding='UTF-8')
persist_qc<-function(){if(length(fit_qc))wc(do.call(rbind,fit_qc),'QC/secondary_model_identity.csv');if(length(warn_qc))wc(do.call(rbind,warn_qc),'QC/secondary_model_warnings.csv');if(length(row_qc))wc(do.call(rbind,row_qc),'QC/secondary_model_rows.csv')}
safe_run({
 for(analysis in c('common_cycle','cycle_adjusted')){
  rows<-list()
  for(sex in c('Female','Male'))for(ex in met){
   des<-getdomain(ex,sex,common=analysis=='common_cycle');f<-reformulate(c(paste0(ex,'_ln'),cv$Model3,if(analysis=='cycle_adjusted')'cycle_f'),response='stroke')
   key<-paste(analysis,sex,ex,sep='_');obj<-fit_checked(f,des,key)
   rows[[length(rows)+1]]<-single(obj,paste0(ex,'_ln'),analysis,sex,ex,'Model3')
   wc(do.call(rbind,rows),paste0('results/',analysis,'_partial.csv'))
  }
  wc(bh(do.call(rbind,rows)),paste0('results/',analysis,'.csv'))
 }
 rows<-list()
 for(ex in c('MBzP','MCNP','MCOP')){
  e<-paste0(ex,'_ln');f<-as.formula(paste('stroke ~',e,'* sex +',paste(cv$Model3,collapse=' + ')))
  obj<-fit_checked(f,getdomain(ex),paste0('sex_interaction_',ex));terms<-names(coef(obj$fit));trm<-terms[grepl(':',terms,fixed=TRUE)&grepl(e,terms,fixed=TRUE)&grepl('sex',terms,fixed=TRUE)]
  ck(length(trm)==1,'Sex interaction term identity');rows[[length(rows)+1]]<-single(obj,trm,'sex_interaction','Overall',ex,'Model3')
 }
 wc(do.call(rbind,rows),'results/sex_interaction.csv')
 rows<-list()
 for(idx in c('LMW','HMW','DEHP'))for(sx in c('Overall','Male','Female')){
  des<-getdomain(idx,if(sx=='Overall')NULL else sx);ex<-paste0('ln_sigma_',idx);f<-reformulate(c(ex,cv$Model3),response='stroke')
  obj<-fit_checked(f,des,paste('grouped',idx,sx,sep='_'));rows[[length(rows)+1]]<-single(obj,ex,'grouped',sx,idx,'Model3')
 }
 wc(do.call(rbind,rows),'results/grouped.csv');wc(data.frame(metabolite=names(MW),MW=MW),'QC/grouped_molecular_weights.csv')
 rows<-list();knots<-list()
 for(sex in c('Female','Male'))for(ex in met){
  des<-getdomain(ex,sex);v<-des$variables[[paste0(ex,'_ln')]];knot<-as.numeric(quantile(v,c(.1,.5,.9),type=7,names=FALSE));ck(all(diff(knot)>0),'Invalid RCS knots')
  b<-rms::rcs(v,parms=knot);ck(ncol(b)==2&&identical(as.logical(attr(b,'nonlinear')),c(FALSE,TRUE))&&identical(as.numeric(attr(b,'parms')),knot),'RCS basis identity')
  des$variables$RCS_linear<-b[,1];des$variables$RCS_nonlinear<-b[,2];f<-reformulate(c('RCS_linear','RCS_nonlinear',cv$Model3),response='stroke')
  obj<-fit_checked(f,des,paste('RCS',sex,ex,sep='_'))
  ov<-survey::regTermTest(obj$fit,~RCS_linear+RCS_nonlinear,method='Wald',df=obj$info$df);nl<-survey::regTermTest(obj$fit,~RCS_nonlinear,method='Wald',df=obj$info$df)
  ck(ov$df==2&&nl$df==1&&ov$ddf==obj$info$df&&nl$ddf==obj$info$df,'RCS explicit F df identity')
  ck(is.finite(ov$p)&&is.finite(nl$p),'RCS invalid global P')
  rows[[length(rows)+1]]<-data.frame(analysis='RCS',sex=sex,exposure=ex,N=obj$info$N,events=obj$info$events,DOMAIN_DF=obj$info$df,overall_df=ov$df,overall_F=as.numeric(ov$Ftest),P_overall_full=ov$p,nonlinear_df=nl$df,nonlinear_F=as.numeric(nl$Ftest),P_nonlinear_full=nl$p)
  knots[[length(knots)+1]]<-data.frame(sex=sex,exposure=ex,N=obj$info$N,knot10=knot[1],knot50=knot[2],knot90=knot[3],rule='ordinary unweighted type7; 3 explicit knots')
  wc(do.call(rbind,rows),'results/RCS_partial.csv');wc(do.call(rbind,knots),'QC/RCS_knots.csv')
 }
 wc(do.call(rbind,rows),'results/RCS.csv')
 rows<-list();specs<-list(Age=list(var='age_group_f',cov=cv$Model3,df=1),BMI=list(var='bmi_group_f',cov=setdiff(cv$Model3,'bmi'),df=2),Race=list(var='race_subgrp',cov=setdiff(cv$Model3,'race_eth_f'),df=3))
 for(ex in c('MCNP','MCOP'))for(mod in names(specs)){
  z<-specs[[mod]];e<-paste0(ex,'_ln');f<-as.formula(paste('stroke ~',e,'*',z$var,'+',paste(z$cov,collapse=' + ')))
  obj<-fit_checked(f,getdomain(ex,'Male'),paste('male_subgroup_global',ex,mod,sep='_'))
  test<-regTermTest(obj$fit,as.formula(paste('~',e,':',z$var)),method='Wald',df=obj$info$df)
  ck(test$df==z$df&&test$ddf==obj$info$df&&is.finite(test$p),'Subgroup global-test identity')
  rows[[length(rows)+1]]<-data.frame(analysis='male_subgroup_global',sex='Male',exposure=ex,modifier=mod,N=obj$info$N,events=obj$info$events,DOMAIN_DF=obj$info$df,numerator_df=test$df,F=as.numeric(test$Ftest),P_full=test$p)
 }
 wc(do.call(rbind,rows),'results/male_subgroup_global.csv')
 rows<-list();dropped<-list()
 for(sex in c('Female','Male'))for(ex in met)for(era in 0:1){
  div<-if(ex %in% c('MCNP','MCOP'))7 else 8;x<-bases[[as.character(div)]]$variables;x$period_weight<-x$ph_weight/4
  base<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~period_weight,data=x,nest=TRUE)
  keep<-x$sex==sex&x$era==era & if(div==8)x$domain_main else x[[paste0('domain_',ex)]]
  des<-base[keep,];x<-des$variables
  covkeep<-cv$Model3[vapply(cv$Model3,function(v){val<-x[[v]];if(is.factor(val)||is.character(val))length(unique(na.omit(val)))>=2 else sd(val,na.rm=TRUE)>0},logical(1))]
  dropped[[length(dropped)+1]]<-data.frame(sex=sex,exposure=ex,era=era,dropped=paste(setdiff(cv$Model3,covkeep),collapse='|'),rule='original deterministic constant-variable rule')
  f<-reformulate(c(paste0(ex,'_ln'),covkeep),response='stroke');obj<-fit_checked(f,des,paste('era_stratified',sex,ex,era,sep='_'))
  row<-single(obj,paste0(ex,'_ln'),'era_stratified',sex,ex,'Model3');row$era<-if(era==0)'2003-2010'else'2011-2018';rows[[length(rows)+1]]<-row
  wc(do.call(rbind,rows),'results/era_stratified_partial.csv')
 }
 wc(do.call(rbind,rows),'results/era_stratified.csv');wc(do.call(rbind,dropped),'QC/era_constant_covariate_rule.csv')
 rows<-list()
 for(ex in c('MCNP','MCOP')){
  e<-paste0(ex,'_ln');f<-as.formula(paste('stroke ~',e,'* era +',paste(cv$Model3,collapse=' + ')))
  obj<-fit_checked(f,getdomain(ex,'Male'),paste0('era_interaction_',ex));trm<-paste0(e,':era')
  rows[[length(rows)+1]]<-single(obj,trm,'era_interaction','Male',ex,'Model3')
 }
 wc(do.call(rbind,rows),'results/era_interaction.csv')
 ck(nrow(do.call(rbind,fit_qc))==120,'Secondary survey model completeness')
 writeLines('SECONDARY_SURVEY_EXECUTION = COMPLETE',file.path(root,'secondary_survey_completion.txt'));cat('SECONDARY_SURVEY_EXECUTION = COMPLETE\n')
})
