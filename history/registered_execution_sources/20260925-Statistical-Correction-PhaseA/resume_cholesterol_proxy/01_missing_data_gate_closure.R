# Bounded Phase A author-adjudicated cholesterol correction and descriptive audit.
out <- normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
ev <- dirname(out); prior <- file.path(ev,"resume_released_values")
suppressPackageStartupMessages(library(survey))
options(stringsAsFactors=FALSE)
wc <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na="")
bind <- function(x) do.call(rbind,x)
qc <- list()
check <- function(name,ok) {qc[[length(qc)+1]] <<- data.frame(check=name,pass=isTRUE(ok)); if(!isTRUE(ok))stop(name)}
# Verify immutable inputs before computing any new descriptive results.
verify <- function(m) {
 m$actual_bytes <- vapply(m$path,function(p)file.info(p)$size,numeric(1))
 m$actual_sha256 <- vapply(m$path,function(p)toupper(digest::digest(file=p,algo="sha256")),character(1))
 m$unchanged <- m$bytes==m$actual_bytes & toupper(m$sha256)==m$actual_sha256
 m
}
protected <- verify(read.csv(file.path(ev,"input_and_protected_manifest.csv")))
inputs <- verify(read.csv(file.path(prior,"raw_input_manifest.csv")))
wc(protected,"protected_preflight.csv");wc(inputs,"raw_preflight.csv")
check("protected_preflight",all(protected$unchanged));check("raw_preflight",all(inputs$unchanged))
orig <- readRDS(file.path(prior,"phase_a_data_candidate.rds"));raw <- orig
check("unique_person",!anyDuplicated(raw$SEQN))
early <- raw$cycle %in% c("2003_2004","2005_2006","2007_2008","2009_2010")
b60 <- tolower(as.character(raw$BPQ060));b80 <- tolower(as.character(raw$BPQ080))
proxy <- early & b60 %in% "no" & is.na(raw$BPQ080) & raw$BPQ_component_present
upstream_missing <- early & b60 %in% c("refused","don't know")
# Observed answers on a skip path would require review, never silently overridden.
check("no_early_observed_answer_on_no_or_unknown_screen_path",
 !any(early & b60 %in% c("no","refused","don't know") & !is.na(raw$BPQ080)))
h <- ifelse(b80 %in% "yes","Yes",ifelse(b80 %in% "no","No",NA_character_))
h[proxy] <- "No"; h[upstream_missing] <- NA_character_
raw$hyperlipidemia_f <- factor(h,levels=levels(orig$hyperlipidemia_f))
reason <- rep("ITEM_NONRESPONSE_OR_UNEXPLAINED_MISSING",nrow(raw))
reason[b80 %in% "yes"] <- "OBSERVED_BPQ080_YES"
reason[b80 %in% "no"] <- "OBSERVED_BPQ080_NO"
reason[upstream_missing] <- "UPSTREAM_REFUSED_DONT_KNOW"
reason[b80 %in% c("refused","don't know")] <- "BPQ080_REFUSED_DONT_KNOW"
reason[!raw$BPQ_component_present] <- "COMPONENT_ABSENCE"
reason[proxy] <- "STRUCTURAL_PROXY_NO_NEVER_CHECKED"
raw$hyperlipidemia_reason <- reason
covs <- c("sex","RIDAGEYR","race_eth_f","edu_f","marital_f","pir","bmi","smoking_f","alcohol_harmonized_f","hypertension_f","diabetes_f","hyperlipidemia_f")
raw$cov_complete <- complete.cases(raw[,covs])
raw$domain_main <- raw$precc & raw$cov_complete
raw$domain_MCNP <- raw$domain_main & raw$cycle!="2003_2004" & raw$MCNP_ok
raw$domain_MCOP <- raw$domain_main & raw$cycle!="2003_2004" & raw$MCOP_ok
raw$domain_common <- raw$domain_main & raw$cycle!="2003_2004" & raw$all10_ok
allowed <- c("hyperlipidemia_f","cov_complete","domain_main","domain_MCNP","domain_MCOP","domain_common")
fieldqc <- data.frame(variable=names(orig),unchanged=vapply(names(orig),function(v)identical(orig[[v]],raw[[v]]),logical(1)))
fieldqc$authorized_change <- fieldqc$variable %in% allowed
wc(fieldqc,"candidate_field_integrity.csv")
check("only_authorized_fields_changed",all(fieldqc$unchanged | fieldqc$authorized_change))
check("nonproxy_cholesterol_unchanged",identical(orig$hyperlipidemia_f[!proxy],raw$hyperlipidemia_f[!proxy]))
check("proxy_all_no",all(as.character(raw$hyperlipidemia_f[proxy])=="No"))
check("no_loss_of_prior_complete_cases",!any(orig$domain_main & !raw$domain_main))
label <- function(x)ifelse(is.na(x),"<NA>",as.character(x))
paths <- aggregate(rep(1,nrow(raw)),list(cycle=raw$cycle,BPQ060=label(raw$BPQ060),BPQ080=label(raw$BPQ080),reason=reason),sum)
names(paths)[5]<-"n";wc(paths,"cholesterol_mapping_paths.csv")
cycles <- sort(unique(raw$cycle));groups <- c("OVERALL",cycles,"Female","Male")
gm <- function(d,g)if(g=="OVERALL")rep(TRUE,nrow(d))else if(g %in% c("Female","Male"))d$sex %in% g else d$cycle==g
flow <- list();counts <- list()
stages <- list(pooled=rep(TRUE,nrow(raw)),adult=raw$age_ok,nonpregnant=raw$age_ok & !raw$pregnant,valid_stroke=raw$clinical,
 positive_weight_design=raw$clinical & raw$base_ok,creatinine_available=raw$clinical & raw$base_ok & raw$ucr_ok,
 eight_analytes_available_preCC=raw$precc,covariates_complete=raw$domain_main,common_cycle_complete=raw$domain_common)
for(g in groups) {
 for(s in names(stages))flow[[length(flow)+1]]<-data.frame(group=g,stage=s,n=sum(gm(raw,g)&stages[[s]]))
 for(f in c("main","common","MCNP","MCOP")) {
  ix <- gm(raw,g)&raw[[paste0("domain_",f)]]
  counts[[length(counts)+1]]<-data.frame(framework=f,group=g,n=sum(ix),stroke_events=sum(raw$stroke[ix]))
 }
}
wc(bind(flow),"participant_flow_counts.csv");wc(bind(counts),"candidate_domain_counts.csv")
classes <- list()
sourcevar<-c(sex="RIAGENDR",RIDAGEYR="RIDAGEYR",race_eth_f="RIDRETH1",edu_f="DMDEDUC2",marital_f="DMDMARTL",pir="INDFMPIR",bmi="BMXBMI",smoking_f="SMQ020",hypertension_f="BPQ020",diabetes_f="DIQ010",hyperlipidemia_f="BPQ080")
component<-c(sex="DEMO",RIDAGEYR="DEMO",race_eth_f="DEMO",edu_f="DEMO",marital_f="DEMO",pir="DEMO",bmi="BMX",smoking_f="SMQ",hypertension_f="BPQ",diabetes_f="DIQ",hyperlipidemia_f="BPQ")
for(v in covs) {
 miss<-is.na(raw[[v]]);cl<-ifelse(miss,"TRUE_ITEM_OR_COMPONENT_MISSING","OBSERVED")
 if(v=="alcohol_harmonized_f") {
  cl[!miss & raw$alcohol_reason %in% "STRUCTURAL_PROXY_LT12_LIFETIME"]<-"STRUCTURAL_PROXY_LT12_LIFETIME"
  cl[!miss & raw$alcohol_reason %in% "LIFETIME_NEVER_STRUCTURAL"]<-"LIFETIME_NEVER_STRUCTURAL"
 } else {
  sv<-sourcevar[[v]]
  derived<-miss & !is.na(raw[[sv]]) & !(tolower(as.character(raw[[sv]]))%in%c("refused","don't know"))
  cl[derived]<-"DERIVED_VARIABLE_MISSING"
  if(v=="hyperlipidemia_f")cl[proxy]<-"STRUCTURAL_PROXY_NO_NEVER_CHECKED"
 }
 classes[[v]]<-cl
}
missrows<-list();reasonrows<-list()
refs<-list(POSITIVE_WEIGHT_CLINICAL_DOMAIN=raw$base_ok&raw$clinical,
 MAIN_PRE_CC=raw$precc,COMMON_PRE_CC=raw$precc&raw$cycle!="2003_2004"&raw$all10_ok)
for(ref in names(refs))for(g in groups)for(v in covs) {
 ix<-refs[[ref]]&gm(raw,g);den<-sum(ix);if(!den)next
 cl<-classes[[v]];m<-is.na(raw[[v]])
 missrows[[length(missrows)+1]]<-data.frame(reference=ref,group=g,variable=v,denominator=den,
 missing_n=sum(ix&m),missing_proportion=sum(ix&m)/den,
 true_item_component_missing_n=sum(ix&cl=="TRUE_ITEM_OR_COMPONENT_MISSING"),
 true_item_component_missing_proportion=sum(ix&cl=="TRUE_ITEM_OR_COMPONENT_MISSING")/den,
 derived_missing_n=sum(ix&cl=="DERIVED_VARIABLE_MISSING"),
 structural_proxy_n=sum(ix&cl%in%c("STRUCTURAL_PROXY_NO_NEVER_CHECKED","STRUCTURAL_PROXY_LT12_LIFETIME")),
 lifetime_never_structural_n=sum(ix&cl=="LIFETIME_NEVER_STRUCTURAL"))
 tt<-as.data.frame(table(cl[ix]));names(tt)<-c("classification","n")
 tt$reference<-ref;tt$group<-g;tt$variable<-v;reasonrows[[length(reasonrows)+1]]<-tt
}
mi<-bind(missrows);wc(mi,"covariate_missingness.csv");wc(bind(reasonrows),"missingness_classification_counts.csv")
# Full positive-weight base designs precede domain selection; same fixed pooled weights.
dat<-raw[raw$base_ok,];dat$pooled_weight_8<-dat$ph_weight/8;dat$pooled_weight_7<-dat$ph_weight/7
base8<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~pooled_weight_8,data=dat,nest=TRUE)
base7<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~pooled_weight_7,data=dat[dat$cycle!="2003_2004",],nest=TRUE)
descvars<-c("stroke","RIDAGEYR","sex","race_eth_f","edu_f","marital_f","pir","bmi","smoking_f","hypertension_f","diabetes_f","hyperlipidemia_f")
check("selection_variable_allowlist_no_concentrations",!any(grepl("^URX|ph_weight|alcohol",descvars)))
response<-list();sel<-list();weighted_qc<-list()
# Only point descriptive estimates retained. No tests, intervals or regression fits.
survey_point<-function(des,z,tag) {
 ok<-!is.na(z);if(!any(ok))return(NA_real_)
 des$variables$.audit_value<-z;dd<-des[ok,]
 ans<-unname(coef(svymean(~.audit_value,dd,na.rm=TRUE))[1])
 w<-as.numeric(weights(dd));direct<-sum(w*z[ok])/sum(w)
 weighted_qc[[length(weighted_qc)+1]]<<-data.frame(estimate=tag,absolute_error=abs(ans-direct),pass=isTRUE(abs(ans-direct)<1e-10))
 ans
}
for(fr in c("MAIN_8_CYCLES","COMMON_7_CYCLES")) {
 base<-if(fr=="MAIN_8_CYCLES")base8 else base7
 d<-base$variables;elig<-d$precc & (if(fr=="COMMON_7_CYCLES")d$all10_ok else TRUE)
 for(g in groups) {
  ix<-elig&gm(d,g);if(!any(ix))next
  ds<-base[ix,];z<-ds$variables;inc<-z$cov_complete;w<-as.numeric(weights(ds))
  for(v in c("COMPLETE_CASE",covs)) {
   observed<-if(v=="COMPLETE_CASE")inc else !is.na(z[[v]])
   p<-survey_point(ds,as.numeric(observed),paste(fr,g,v,"response"))
   response[[length(response)+1]]<-data.frame(framework=fr,group=g,variable=v,denominator=nrow(z),response_n=sum(observed),unweighted_response_proportion=mean(observed),weight_sum=sum(w),respondent_weight_sum=sum(w[observed]),weighted_response_proportion=p,weighted_missing_proportion=1-p)
  }
  for(v in descvars) {
   x<-z[[v]];levs<-if(is.factor(x))levels(x)else if(v=="stroke")"1"else NA_character_
   for(lv in levs) {
    a<-if(is.na(lv))as.numeric(x)else ifelse(is.na(x),NA_real_,as.numeric(as.character(x)==lv))
    vals<-list()
    for(method in c("UNWEIGHTED","SURVEY_WEIGHTED")) {
     estimates<-lapply(c(TRUE,FALSE),function(which) {
      keep<-inc==which & !is.na(a);n<-sum(keep)
      if(!n)return(c(n=0,wsum=0,mean=NA,variance=NA))
      ww<-if(method=="UNWEIGHTED")rep(1,n)else w[keep]
      mu<-if(method=="UNWEIGHTED")mean(a[keep])else survey_point(ds[keep,],a[keep],paste(fr,g,v,lv,which))
      vv<-if(is.na(lv)) {if(method=="UNWEIGHTED")var(a[keep])else sum(ww*(a[keep]-mu)^2)/sum(ww)}else mu*(1-mu)
      c(n=n,wsum=sum(ww),mean=mu,variance=vv)
     })
     aa<-estimates[[1]];bb<-estimates[[2]];sdpool<-sqrt((aa['variance']+bb['variance'])/2)
     smd<-if(is.finite(sdpool)&&sdpool>0)abs(aa['mean']-bb['mean'])/sdpool else if(isTRUE(aa['mean']==bb['mean']))0 else NA_real_
     sel[[length(sel)+1]]<-data.frame(framework=fr,group=g,variable=v,level=if(is.na(lv))"CONTINUOUS"else lv,weighting=method,
      included_total=sum(inc),excluded_total=sum(!inc),included_nonmissing=unname(aa['n']),excluded_nonmissing=unname(bb['n']),
      included_valid_weight_sum=if(method=="SURVEY_WEIGHTED")unname(aa['wsum'])else NA_real_,excluded_valid_weight_sum=if(method=="SURVEY_WEIGHTED")unname(bb['wsum'])else NA_real_,
      included_mean_or_proportion=unname(aa['mean']),excluded_mean_or_proportion=unname(bb['mean']),
      included_variance=unname(aa['variance']),excluded_variance=unname(bb['variance']),absolute_SMD=unname(smd))
    }
   }
  }
 }
}
resp<-bind(response);selection<-bind(sel);wqc<-bind(weighted_qc)
wc(resp,"weighted_unweighted_response_proportions.csv");wc(selection,"included_excluded_selection.csv");wc(wqc,"weighted_estimate_QC.csv")
check("survey_points_match_direct_weighted_means",all(wqc$pass))
check("flow_partition_main",sum(raw$precc)==sum(raw$domain_main)+sum(raw$precc&!raw$cov_complete))
check("missingness_classification_partition",all(mi$missing_n==mi$true_item_component_missing_n+mi$derived_missing_n))
check("response_proportions_valid",all(resp$weighted_response_proportion>=0 & resp$weighted_response_proportion<=1))
gate_missing<-mi[mi$true_item_component_missing_proportion>0.10,]
gate_smd<-selection[is.finite(selection$absolute_SMD)&selection$absolute_SMD>=0.10,]
wc(gate_missing,"gate_true_missing_over10pct.csv");wc(gate_smd,"gate_selection_SMD_ge010.csv")
saveRDS(raw,file.path(out,"phase_a_data_candidate.rds"))
saveRDS(list(base8=base8,base7=base7),file.path(out,"phase_a_base_designs_candidate.rds"))
protected_after<-verify(read.csv(file.path(ev,"input_and_protected_manifest.csv")))
raw_after<-verify(read.csv(file.path(prior,"raw_input_manifest.csv")))
wc(protected_after,"protected_after.csv");wc(raw_after,"raw_after.csv")
check("protected_after",all(protected_after$unchanged));check("raw_after",all(raw_after$unchanged))
wc(bind(qc),"closure_QC.csv")
jsonlite::write_json(list(status="STOPPED_AT_DECISION_BOUNDARY",MISSING_DATA_DECISION_GATE="MISSING_DATA_DECISION_REQUIRED",
 author_cholesterol_decision="APPLIED_TO_PHASE_A_CANDIDATE_ONLY",main_preCC=sum(raw$precc),main_complete_case=sum(raw$domain_main),main_excluded=sum(raw$precc&!raw$cov_complete),
 common_preCC=sum(raw$precc&raw$cycle!="2003_2004"&raw$all10_ok),common_complete_case=sum(raw$domain_common),
 main_added_complete_cases=sum(raw$domain_main&!orig$domain_main),cholesterol_structural_proxy_preCC=sum(proxy&raw$precc),
 true_missing_over10pct_rows=nrow(gate_missing),selection_SMD_ge010_rows=nrow(gate_smd),
 missing_data_route="NOT_SELECTED_REQUIRES_AUTHOR_DECISION",CORRECTED_EXPOSURE_STROKE_MODELS_RUN="NO",CORRECTED_OR_CI_P_Q_GENERATED="NO",PHASE_B_EXECUTED="NO",
 protected_files_unchanged=nrow(protected_after),raw_files_unchanged=nrow(raw_after),qc_checks=nrow(bind(qc)),weighted_point_checks=nrow(wqc)),
 file.path(out,"phase_a_checkpoint.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat(readLines(file.path(out,"phase_a_checkpoint.json")),sep="\n")
