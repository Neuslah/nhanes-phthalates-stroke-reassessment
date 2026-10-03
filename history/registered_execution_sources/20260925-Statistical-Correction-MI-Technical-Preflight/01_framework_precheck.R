# Coefficient-blind preparation only: membership, design information, donor support.
out<-normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
phase<-file.path(dirname(out),"20260925-Statistical-Correction-PhaseA")
src<-file.path(phase,"resume_cholesterol_proxy/phase_a_data_candidate.rds")
stopifnot(toupper(digest::digest(file=src,algo="sha256"))=="43826EA86EAA72F598D35BE1B36FC3BEA282B1BB0C255A54E8A0F19025C3FE16")
raw<-readRDS(src)
wc<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na="")
covs<-c("RIDAGEYR","race_eth_f","edu_f","marital_f","pir","bmi","smoking_f","alcohol_harmonized_f","hypertension_f","diabetes_f","hyperlipidemia_f")
met<-c(URXMBP="MBP_ln",URXMZP="MBzP_ln",URXECP="MECPP_ln",URXMHH="MEHHP_ln",URXMOH="MEOHP_ln",URXMC1="MCPP_ln",URXMEP="MEP_ln",URXMIB="MiBP_ln",URXCNP="MCNP_ln",URXCOP="MCOP_ln")
borderline<-tolower(as.character(raw$DIQ010))%in%"borderline"
membership<-list();summary<-list();support<-list();designmap<-list();inputs<-list();issues<-list()
for(window in c("eight_2003_2018","ten_2005_2018"))for(sex in c("Female","Male")){
 key<-paste(sex,window,sep="_")
 eligible<-raw$precc & raw$sex%in%sex & if(window=="eight_2003_2018")TRUE else raw$cycle!="2003_2004" & raw$all10_ok
 d<-raw[eligible & !borderline,,drop=FALSE]
 nmet<-if(window=="eight_2003_2018")8 else 10;divisor<-if(nmet==8)8 else 7
 dat<-d[,c("stroke",covs),drop=FALSE]
 for(v in names(met)[seq_len(nmet)])dat[[met[[v]]]]<-log(d[[v]]/d$URXUCR*1000)
 stopifnot(all(vapply(dat[,unname(met[seq_len(nmet)]),drop=FALSE],function(z)all(is.finite(z)),logical(1))))
 psu<-interaction(d$cycle,d$SDMVSTRA,d$SDMVPSU,drop=TRUE,lex.order=TRUE)
 basis<-model.matrix(~psu);design<-basis[,-1,drop=FALSE]
 colnames(design)<-sprintf("PSU_FE_%03d",seq_len(ncol(design)))
 dat$log_analysis_weight<-log(d$ph_weight/divisor)
 dat<-cbind(dat,as.data.frame(design))
 # Membership/nesting verification does not compute any exposure/outcome association.
 unit<-unique(data.frame(level=as.character(psu),cycle=d$cycle,stratum=d$SDMVSTRA,PSU=d$SDMVPSU))
 stopifnot(!anyDuplicated(unit$level))
 unit$framework<-key;unit$dummy_column<-ifelse(match(unit$level,levels(psu))==1,"REFERENCE",sprintf("PSU_FE_%03d",match(unit$level,levels(psu))-1))
 designmap[[key]]<-unit
 cycle_strata<-model.matrix(~factor(d$cycle)+factor(d$SDMVSTRA))
 # A full PSU dummy basis exactly represents every coarser cycle/stratum indicator.
 psu_index<-as.integer(psu)
 nesting_equal<-all(vapply(seq_len(ncol(cycle_strata)),function(j){z<-cycle_strata[,j];all(vapply(split(z,psu_index),function(x)length(unique(x))==1,logical(1)))},logical(1)))
 stopifnot(nesting_equal)
 methods<-vapply(dat,function(x)if(!anyNA(x))"" else if(is.factor(x)){if(nlevels(x)==2)"logreg"else"polyreg"}else"pmm",character(1))
 targets<-names(methods)[methods!=""]
 stopifnot(all(targets%in%covs),!anyNA(dat$stroke),!anyNA(dat$log_analysis_weight))
 pred<-matrix(0,ncol(dat),ncol(dat),dimnames=list(names(dat),names(dat)))
 pred[targets,]<-1;diag(pred)<-0
 for(v in targets){
  obs<-!is.na(dat[[v]]);nobs<-table(factor(psu[obs],levels=levels(psu)));nmis<-table(factor(psu[!obs],levels=levels(psu)))
  absent<-which(nobs==0 & nmis>0)
  support[[length(support)+1]]<-data.frame(framework=key,target=v,method=methods[v],n=nrow(dat),observed=sum(obs),missing=sum(!obs),observed_categories=if(is.factor(dat[[v]]))length(unique(dat[[v]][obs]))else NA_integer_,PSUs_missing_target_without_observed_target=length(absent),design_only_rank=ncol(basis)-sum(nobs==0),design_only_columns=ncol(basis))
  if(length(absent))issues[[length(issues)+1]]<-data.frame(framework=key,target=v,issue="PSU_TARGET_NO_OBSERVED_INFORMATION",PSU_level=levels(psu)[absent])
 }
 nstr<-nrow(unique(unit[,c("cycle","stratum")]))
 lonely<-sum(table(interaction(unit$cycle,unit$stratum,drop=TRUE))==1)
 summary[[key]]<-data.frame(framework=key,preCC_before_borderline=sum(eligible),borderline_excluded=sum(eligible&borderline),MI_eligible=nrow(d),CC_in_MI_domain=sum(d$cov_complete),PSUs=nrow(unit),strata=nstr,manual_design_df=nrow(unit)-nstr,lonely_strata=lonely,weight_divisor=divisor,n_exposures=nmet,imputation_targets=length(targets),PSU_reference=levels(psu)[1],design_basis_preserves_cycle_strata=nesting_equal)
 membership[[key]]<-data.frame(framework=key,SEQN=d$SEQN,cycle=d$cycle,SDMVSTRA=d$SDMVSTRA,SDMVPSU=d$SDMVPSU)
 inputs[[key]]<-list(data=dat,method=methods,predictorMatrix=pred,SEQNs=d$SEQN,weight_divisor=divisor)
 wc(data.frame(variable=names(methods),method=unname(methods)),paste0(key,"_methods.csv"))
 write.csv(pred,file.path(out,paste0(key,"_predictor_matrix.csv")))
}
wc(do.call(rbind,summary),"MI_framework_counts.csv")
wc(do.call(rbind,support),"target_information_support.csv")
wc(do.call(rbind,designmap),"design_encoding_map.csv")
wc(do.call(rbind,membership),"framework_membership.csv")
saveRDS(inputs,file.path(out,"MI_preflight_inputs.rds"))
ip<-installed.packages();wc(data.frame(package=ip[,"Package"],version=ip[,"Version"],library=ip[,"LibPath"]),"initial_package_versions.csv")
capture.output(sessionInfo(),file=file.path(out,"framework_sessionInfo.txt"))
if(length(issues)){wc(do.call(rbind,issues),"precheck_blocking_issues.csv");cat("STOP / TARGET_INFORMATION_SUPPORT_BOUNDARY\n")}else cat("FRAMEWORK_SUPPORT_PRECHECK_PASS\n")
print(do.call(rbind,summary),row.names=FALSE)
