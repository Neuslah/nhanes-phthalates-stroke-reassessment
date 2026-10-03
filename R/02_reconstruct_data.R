# Phase A data/design/missingness audit. No exposure-outcome model or association.
out <- normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
inputs<-normalizePath(commandArgs(TRUE)[2],winslash="/",mustWork=TRUE); docs<-file.path(inputs,"official_docs"); package_root<-normalizePath(commandArgs(TRUE)[3],winslash="/",mustWork=TRUE)
rawroot<-file.path(inputs,"translated")
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(survey))
wtcsv<-function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na="")
cyc<-data.frame(suffix=LETTERS[3:10],year=seq(2003,2017,2));cyc$cycle<-paste(cyc$year,cyc$year+1,sep="_")
cyc$ph<-c("L24PH_C",paste0("PHTHTE_",LETTERS[4:10]))
met<-c(URXMBP="MBP",URXMZP="MBzP",URXECP="MECPP",URXMHH="MEHHP",URXMOH="MEOHP",URXMC1="MCPP",URXMEP="MEP",URXMIB="MiBP",URXCNP="MCNP",URXCOP="MCOP")
basic<-names(met)[1:8]
flags<-c("URDMBPLC","URDMZPLC","URDECPLC","URDMHHLC","URDMOHLC","URDMC1LC","URDMEPLC","URDMIBLC","URDCNPLC","URDCOPLC")
alcohol<-readRDS(file.path(out,"alcohol_harmonized_candidate.rds"))
parts<-list();pregdocs<-list()
rd<-function(cy,comp,cols) {
 z<-readRDS(file.path(rawroot,cy,paste0(comp,".rds")))
 stopifnot(!anyDuplicated(z$SEQN)); z[,intersect(c("SEQN",cols),names(z)),drop=FALSE]
}
for(i in 1:8) {
 cy<-cyc[i,];s<-cy$suffix
 d<-rd(cy$cycle,paste0("DEMO_",s),c("RIAGENDR","RIDAGEYR","RIDRETH1","DMDEDUC2","DMDMARTL","INDFMPIR","SDMVPSU","SDMVSTRA","WTMEC2YR","RIDEXPRG"))
 d$cycle<-cy$cycle;d$demographic_mec_weight<-d$WTMEC2YR;d$WTMEC2YR<-NULL
 comps<-list(MCQ="MCQ160F",BPQ=c("BPQ020","BPQ040A","BPQ060","BPQ080"),DIQ="DIQ010",SMQ="SMQ020",BMX="BMXBMI")
 for(comp in names(comps)) {
  z<-rd(cy$cycle,paste0(comp,"_",s),comps[[comp]])
  z[[paste0(comp,"_component_present")]]<-TRUE
  d<-left_join(d,z,by="SEQN")
 }
 ph<-haven::read_xpt(file.path(inputs,"official_xpt",paste0(cy$ph,".xpt")))
 wt<-if(i==5)"WTSA2YR" else "WTSB2YR"
 ph$ph_weight<-as.numeric(ph[[wt]]);ph$ph_component_present<-TRUE
 ph<-as.data.frame(ph[,intersect(c("SEQN","ph_weight","ph_component_present",names(met),flags),names(ph))])
 for(v in setdiff(c(names(met),flags),names(ph))) ph[[v]]<-NA_real_
 for(v in intersect(c(names(met),flags),names(ph))) ph[[v]]<-as.numeric(ph[[v]])
 d<-left_join(d,ph,by="SEQN")
 if(i>=7) ucr<-rd(cy$cycle,paste0("ALB_CR_",s),"URXUCR") else ucr<-rd(cy$cycle,cy$ph,"URXUCR")
 ucr$ucr_component_present<-TRUE
 d<-left_join(d,ucr,by="SEQN")
 a<-alcohol[alcohol$cycle==cy$cycle,];a$cycle<-NULL;a$ALQ_component_present<-TRUE
 d<-left_join(d,a,by="SEQN")
 # Historical alcohol implementation, for participant-flow provenance ONLY.
 olda<-rd(cy$cycle,paste0("ALQ_",s),c("ALQ101","ALQ111","ALQ120Q"))
 if(!"ALQ101"%in%names(olda)) olda$ALQ101<-olda$ALQ111
 if(!"ALQ120Q"%in%names(olda)) olda$ALQ120Q<-NA_real_
 aq<-as.character(olda$ALQ101);fq<-as.numeric(olda$ALQ120Q)
 olda$alcohol_legacy<-ifelse(aq=="No","No",ifelse(aq=="Yes"&!is.na(fq)&fq==0,"No",ifelse(aq=="Yes","Yes",NA_character_)))
 d<-left_join(d,olda[,c("SEQN","alcohol_legacy")],by="SEQN")
 d$weight_source<-wt;d$ucr_source<-if(i>=7)paste0("ALB_CR_",s) else cy$ph
 parts[[i]]<-d
 doc<-xml2::read_html(file.path(docs,paste0("DEMO_",s,".html")))
 sec<-xml2::xml_find_first(doc,"//*[@id='RIDEXPRG']/parent::*")
 stopifnot(!inherits(sec,"xml_missing"))
 pregdocs[[i]]<-data.frame(cycle=cy$cycle,variable="RIDEXPRG",official_section=rvest::html_text2(sec))
}
raw<-bind_rows(parts);stopifnot(!anyDuplicated(raw$SEQN))
for(v in grep("_component_present$",names(raw),value=TRUE)) raw[[v]][is.na(raw[[v]])]<-FALSE
raw$alcohol_reason[!raw$ALQ_component_present]<-"ALQ_COMPONENT_ABSENCE"
raw$age_ok<-!is.na(raw$RIDAGEYR)&raw$RIDAGEYR>=20
raw$pregnant<-as.character(raw$RIDEXPRG) %in% "Yes, positive lab pregnancy test or self-reported pregnant at exam"
raw$stroke_ok<-as.character(raw$MCQ160F)%in%c("Yes","No")
raw$clinical<-raw$age_ok&!raw$pregnant&raw$stroke_ok
raw$stroke<-ifelse(raw$stroke_ok,as.integer(as.character(raw$MCQ160F)=="Yes"),NA_integer_)
raw$weight_ok<-!is.na(raw$ph_weight)&is.finite(raw$ph_weight)&raw$ph_weight>0
raw$design_vars_ok<-!is.na(raw$SDMVPSU)&!is.na(raw$SDMVSTRA)
raw$base_ok<-raw$ph_component_present&raw$weight_ok&raw$design_vars_ok
raw$ucr_ok<-!is.na(raw$URXUCR)&is.finite(raw$URXUCR)&raw$URXUCR>0
avail<-sapply(names(met),function(v)!is.na(raw[[v]])&is.finite(raw[[v]])&raw[[v]]>0)
raw$basic_ok<-rowSums(avail[,basic,drop=FALSE])==8
raw$all10_ok<-rowSums(avail)==10
raw$MCNP_ok<-avail[,"URXCNP"];raw$MCOP_ok<-avail[,"URXCOP"]
# Upstream logic unchanged; enforce before downstream auditing.
anchor<-data.frame(metric=c("initial_pooled","age_excluded","pregnancy_excluded_after_age","invalid_stroke_after_age_preg"),
 historical=c(80312,35522,941,68),
 observed=c(nrow(raw),sum(!raw$age_ok),sum(raw$age_ok&raw$pregnant),sum(raw$age_ok&!raw$pregnant&!raw$stroke_ok)),
 class="HARD_UPSTREAM_ANCHOR")
anchor$match<-anchor$historical==anchor$observed
wtcsv(anchor,"upstream_anchors.csv")
if(!all(anchor$match))stop("STOP / upstream hard anchor mismatch")
# Covariates preserve the registered builder definitions, except authorized alcohol mapping.
f<-function(x,lev)factor(x,levels=lev)
raw$sex<-f(ifelse(as.character(raw$RIAGENDR)%in%c("Female","Male"),as.character(raw$RIAGENDR),NA),c("Female","Male"))
race<-as.character(raw$RIDRETH1)
raw$race_eth_f<-f(ifelse(race%in%c("Mexican American","Other Hispanic","Non-Hispanic White","Non-Hispanic Black"),race,
 ifelse(race%in%c("Non-Hispanic Asian","Other Race - Including Multi-Racial","Other/Multi-racial"),"Other/Multiracial",NA)),
 c("Non-Hispanic White","Mexican American","Other Hispanic","Non-Hispanic Black","Other/Multiracial"))
edu<-as.character(raw$DMDEDUC2)
raw$edu_f<-f(ifelse(edu%in%c("Less than 9th grade","Less Than 9th Grade","9-11th grade (Includes 12th grade with no diploma)","9-11th Grade (Includes 12th grade with no diploma)"),"Less than high school",
 ifelse(edu%in%c("High school graduate/GED or equivalent","High School Grad/GED or Equivalent"),"High school graduate",
 ifelse(edu%in%c("Some college or AA degree","Some College or AA degree","College graduate or above","College Graduate or above"),"College or above",NA))),
 c("College or above","High school graduate","Less than high school"))
mar<-as.character(raw$DMDMARTL)
raw$marital_f<-f(ifelse(mar%in%c("Married","Living with partner"),"Married/partnered",ifelse(mar=="Never married","Never married",
 ifelse(mar%in%c("Widowed","Divorced","Separated"),"Separated/divorced/widowed",NA))),c("Married/partnered","Never married","Separated/divorced/widowed"))
raw$pir<-raw$INDFMPIR;raw$bmi<-raw$BMXBMI
sm<-as.character(raw$SMQ020);raw$smoking_f<-f(ifelse(sm=="Yes","Ever",ifelse(sm=="No","Never",NA)),c("Never","Ever"))
raw$alcohol_harmonized_f<-f(raw$alcohol_harmonized_f,c("No","Yes"))
bp<-as.character(raw$BPQ020);bp4<-as.character(raw$BPQ040A)
raw$hypertension_f<-f(ifelse(bp=="Yes","Yes",ifelse(!is.na(bp4)&bp4=="Yes","Yes",ifelse(bp=="No","No",NA))),c("No","Yes"))
di<-as.character(raw$DIQ010);raw$diabetes_f<-f(ifelse(di=="Yes","Yes",ifelse(di=="No","No",NA)),c("No","Yes"))
hp<-as.character(raw$BPQ080);raw$hyperlipidemia_f<-f(ifelse(hp=="Yes","Yes",ifelse(hp=="No","No",NA)),c("No","Yes"))
covs<-c("sex","RIDAGEYR","race_eth_f","edu_f","marital_f","pir","bmi","smoking_f","alcohol_harmonized_f","hypertension_f","diabetes_f","hyperlipidemia_f")
othercovs<-setdiff(covs,"alcohol_harmonized_f")
raw$cov_complete<-complete.cases(raw[,covs])
raw$legacy_cov_complete<-complete.cases(raw[,othercovs])&!is.na(raw$alcohol_legacy)
raw$precc<-raw$clinical&raw$base_ok&raw$ucr_ok&raw$basic_ok
raw$domain_main<-raw$precc&raw$cov_complete
raw$legacy_main<-raw$clinical&raw$ucr_ok&raw$basic_ok&raw$legacy_cov_complete&raw$base_ok
raw$domain_MCNP<-raw$domain_main&raw$cycle!="2003_2004"&raw$MCNP_ok
raw$domain_MCOP<-raw$domain_main&raw$cycle!="2003_2004"&raw$MCOP_ok
raw$domain_common<-raw$domain_main&raw$cycle!="2003_2004"&raw$all10_ok
# Pregnancy counts retain denominator/stage explicitly.
preg<-raw%>%group_by(cycle)%>%summarise(pooled_n=n(),pregnant_all_ages=sum(pregnant),adults=sum(age_ok),excluded_pregnant_after_age=sum(age_ok&pregnant),RIDEXPRG_missing_adults=sum(age_ok&is.na(RIDEXPRG)),.groups="drop")
wtcsv(preg,"pregnancy_counts.csv");wtcsv(bind_rows(pregdocs),"pregnancy_official_definitions.csv")
stages<-list(pooled=rep(TRUE,nrow(raw)),age=raw$age_ok,nonpregnant=raw$age_ok&!raw$pregnant,valid_stroke=raw$clinical,
 urinary_creatinine=raw$clinical&raw$ucr_ok,basic_metabolites=raw$clinical&raw$ucr_ok&raw$basic_ok,
 legacy_complete_case=raw$legacy_main,corrected_complete_case=raw$domain_main)
flow<-list()
for(gr in c("OVERALL",cyc$cycle)) {
 g<-if(gr=="OVERALL")rep(TRUE,nrow(raw)) else raw$cycle==gr
 for(nm in names(stages)) flow[[length(flow)+1]]<-data.frame(group=gr,stage=nm,n=sum(g&stages[[nm]]))
}
wtcsv(bind_rows(flow),"participant_flow_counts.csv")
legacy_steps<-list(excluded_23403=raw$clinical&!raw$ucr_ok,
 excluded_7049=raw$clinical&raw$ucr_ok&!raw$basic_ok,
 excluded_4200=raw$clinical&raw$ucr_ok&raw$basic_ok&!raw$legacy_main)
raw$sampling_status<-ifelse(!raw$ph_component_present,"OUTSIDE_RELEASED_PHT_COMPONENT",
 ifelse(!raw$weight_ok,"PHT_COMPONENT_NONPOSITIVE_OR_MISSING_WEIGHT",ifelse(!raw$design_vars_ok,"INVALID_DESIGN_VARIABLES","PHT_COMPONENT_POSITIVE_WEIGHT")))
raw$exposure_status<-ifelse(!raw$ph_component_present,"COMPONENT_ABSENCE",
 ifelse(raw$basic_ok,"BASIC_ANALYTES_AVAILABLE","COMPONENT_PRESENT_ANALYTE_UNAVAILABLE"))
raw$ucr_status<-ifelse(raw$ucr_ok,"AVAILABLE",ifelse(!raw$ucr_component_present,"SOURCE_COMPONENT_ABSENCE","SOURCE_PRESENT_UCR_MISSING_OR_INVALID"))
decomp<-list()
for(nm in names(legacy_steps)) {
 z<-raw[legacy_steps[[nm]],]
 tt<-z%>%group_by(cycle,sampling_status,exposure_status,ucr_status,legacy_cov_complete)%>%summarise(n=n(),.groups="drop")
 tt$legacy_exclusion_stage<-nm;decomp[[nm]]<-tt
}
wtcsv(bind_rows(decomp),"historical_exclusion_decomposition.csv")
sampling<-raw%>%group_by(cycle,sampling_status)%>%summarise(pooled_n=n(),clinical_n=sum(clinical),.groups="drop")
wtcsv(sampling,"sampling_framework_counts.csv")
# Corrected flow through a positive-weight design framework, no imputation of unavailable concentrations.
masks<-list(base=raw$base_ok,adult=raw$base_ok&raw$age_ok,nonpregnant=raw$base_ok&raw$age_ok&!raw$pregnant,
 stroke_valid=raw$base_ok&raw$clinical,ucr_valid=raw$base_ok&raw$clinical&raw$ucr_ok,
 eight_metabolites_available=raw$precc,covariates_complete=raw$domain_main,common_cycle=raw$domain_common)
cf<-list()
for(gr in c("OVERALL",cyc$cycle))for(nm in names(masks))cf[[length(cf)+1]]<-data.frame(group=gr,stage=nm,n=sum(masks[[nm]]&(if(gr=="OVERALL")TRUE else raw$cycle==gr)))
wtcsv(bind_rows(cf),"corrected_design_domain_flow.csv")
# Historical downstream references are comparisons, never pass conditions.
h<-data.frame(metric=c("creatinine_excluded","basic_metabolites_excluded","covariate_or_design_excluded","main_n","female_n","male_n","stroke_events","common_n","common_events"),
 historical=c(23403,7049,4200,9129,4611,4518,375,8311,339),
 recreated_legacy=c(sapply(legacy_steps,sum),sum(raw$legacy_main),sum(raw$legacy_main&raw$sex=="Female",na.rm=TRUE),sum(raw$legacy_main&raw$sex=="Male",na.rm=TRUE),sum(raw$stroke[raw$legacy_main]),sum(raw$legacy_main&raw$cycle!="2003_2004"&raw$all10_ok),sum(raw$stroke[raw$legacy_main&raw$cycle!="2003_2004"&raw$all10_ok])),
 corrected=c(sum(legacy_steps[[1]]),sum(legacy_steps[[2]]),sum(raw$clinical&raw$ucr_ok&raw$basic_ok&!raw$domain_main),sum(raw$domain_main),sum(raw$domain_main&raw$sex=="Female",na.rm=TRUE),sum(raw$domain_main&raw$sex=="Male",na.rm=TRUE),sum(raw$stroke[raw$domain_main]),sum(raw$domain_common),sum(raw$stroke[raw$domain_common])),
 class="HISTORICAL_V1.0.1_REFERENCE")
wtcsv(h,"historical_anchor_comparison.csv")
# Categorize missingness without treating phthalate subsampling as covariate nonresponse.
sourcevar<-c(sex="RIAGENDR",RIDAGEYR="RIDAGEYR",race_eth_f="RIDRETH1",edu_f="DMDEDUC2",marital_f="DMDMARTL",pir="INDFMPIR",bmi="BMXBMI",smoking_f="SMQ020",hypertension_f="BPQ020",diabetes_f="DIQ010",hyperlipidemia_f="BPQ080")
component<-c(sex="DEMO",RIDAGEYR="DEMO",race_eth_f="DEMO",edu_f="DEMO",marital_f="DEMO",pir="DEMO",bmi="BMX",smoking_f="SMQ",hypertension_f="BPQ",diabetes_f="DIQ",hyperlipidemia_f="BPQ")
missing_rows<-list(); reason_counts<-list()
groups<-c("OVERALL",cyc$cycle,"Female","Male")
groupmask<-function(gr)if(gr=="OVERALL")rep(TRUE,nrow(raw))else if(gr%in%c("Female","Male"))raw$sex%in%gr else raw$cycle==gr
for(ref in c("POSITIVE_WEIGHT_CLINICAL_DOMAIN","EXPOSURE_ELIGIBLE_BEFORE_CC")) {
 mask<-if(ref=="POSITIVE_WEIGHT_CLINICAL_DOMAIN")raw$base_ok&raw$clinical else raw$precc
 for(v in covs) {
  miss<-is.na(raw[[v]])
  cls<-rep("OBSERVED",nrow(raw));cls[miss]<-"ITEM_NONRESPONSE"
  if(v=="alcohol_harmonized_f") {
   cls[miss]<-raw$alcohol_reason[miss]
   cls[!miss&raw$alcohol_reason=="STRUCTURAL_PROXY_LT12_LIFETIME"]<-"OBSERVED_VIA_STRUCTURAL_PROXY"
   cls[!miss&raw$alcohol_reason=="LIFETIME_NEVER_STRUCTURAL"]<-"OBSERVED_VIA_LIFETIME_NEVER_SKIP"
  } else {
   sv<-sourcevar[[v]]
   derived<-miss&!is.na(raw[[sv]])&!(tolower(as.character(raw[[sv]]))%in%c("refused","don't know"))
   cls[derived]<-"DERIVED_VARIABLE_MISSING"
   comp<-component[[v]]
   if(comp!="DEMO")cls[miss&!raw[[paste0(comp,"_component_present")]]]<-"COMPONENT_ABSENCE"
   if(v=="hyperlipidemia_f") {
    early<-raw$cycle%in%cyc$cycle[1:4]
    no_test<-tolower(as.character(raw$BPQ060))%in%"no"
    refused_screen<-tolower(as.character(raw$BPQ060))%in%c("refused","don't know")
    cls[miss&is.na(raw$BPQ080)&early&no_test]<-"STRUCTURAL_SKIP_NO_CHOLESTEROL_TEST"
    cls[miss&is.na(raw$BPQ080)&early&refused_screen]<-"UPSTREAM_ITEM_NONRESPONSE_CAUSING_SKIP"
   }
  }
  for(gr in groups) {
   g<-mask&groupmask(gr);den<-sum(g)
   missing_rows[[length(missing_rows)+1]]<-data.frame(reference=ref,group=gr,variable=v,denominator=den,missing_n=sum(g&miss),missing_proportion=if(den)sum(g&miss)/den else NA_real_,
    true_item_or_component_missing_n=sum(g&miss&!(cls%in%c("DERIVED_VARIABLE_MISSING","STRUCTURAL_SKIP_NO_CHOLESTEROL_TEST"))),
    derived_missing_n=sum(g&cls=="DERIVED_VARIABLE_MISSING"),
    structural_missing_n=sum(g&miss&cls=="STRUCTURAL_SKIP_NO_CHOLESTEROL_TEST"),
    structural_proxy_observed_n=sum(g&cls=="OBSERVED_VIA_STRUCTURAL_PROXY"))
   t<-as.data.frame(table(cls[g]));names(t)<-c("reason","n");t<-t[t$n>0,]
   t$reference<-ref;t$group<-gr;t$variable<-v;reason_counts[[length(reason_counts)+1]]<-t
  }
 }
}
mi<-bind_rows(missing_rows);mi$true_item_or_component_proportion<-mi$true_item_or_component_missing_n/mi$denominator
wtcsv(mi,"covariate_missingness.csv");wtcsv(bind_rows(reason_counts),"missingness_reason_counts.csv")
# Design-variable availability is reported before base-design selection as well.
designmiss<-list()
for(gr in groups)for(v in c("SDMVPSU","SDMVSTRA","ph_weight")) {
 g<-raw$clinical&groupmask(gr)
 designmiss[[length(designmiss)+1]]<-data.frame(group=gr,variable=v,denominator=sum(g),missing_n=sum(g&is.na(raw[[v]])),
  sampling_component_absence_n=sum(g&!raw$ph_component_present&is.na(raw[[v]])),nonpositive_n=sum(g&!is.na(raw[[v]])&raw[[v]]<=0))
}
wtcsv(bind_rows(designmiss),"design_variable_availability.csv")
# Base designs created BEFORE all clinical, exposure and covariate domain restrictions.
dat<-raw[raw$base_ok,]
dat$pooled_weight_8<-dat$ph_weight/8
dat$pooled_weight_7<-dat$ph_weight/7
base8<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~pooled_weight_8,data=dat,nest=TRUE)
dat7<-dat[dat$cycle!="2003_2004",]
base7<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~pooled_weight_7,data=dat7,nest=TRUE)
# Verify that original stratum IDs are unique across cycles, before using them pooled.
cross<-raw%>%filter(base_ok)%>%distinct(cycle,SDMVSTRA)%>%count(SDMVSTRA)
if(any(cross$n>1))stop("STOP: stratum identifiers repeat across cycles; design identity review required")
design_diag<-list();units<-list()
forms<-~RIDAGEYR+race_eth_f+edu_f+marital_f+pir+bmi+smoking_f+alcohol_harmonized_f+hypertension_f+diabetes_f+hyperlipidemia_f
for(fr in c("basic_2003_2018","MCNP_2005_2018","MCOP_2005_2018","common_2005_2018"))for(sex in c("Female","Male")) {
 base<-if(fr=="basic_2003_2018")base8 else base7
 v<-base$variables
 dm<-switch(fr,basic_2003_2018="domain_main",MCNP_2005_2018="domain_MCNP",MCOP_2005_2018="domain_MCOP",common_2005_2018="domain_common")
 keep<-v[[dm]]&as.character(v$sex)==sex
 ss<-base[keep,]
 v<-ss$variables
 unit<-v%>%distinct(cycle,SDMVSTRA,SDMVPSU)
 strat<-unit%>%count(cycle,SDMVSTRA,name="represented_psu")
 manual<-nrow(unit)-nrow(strat)
 nuisance_rank<-qr(model.matrix(forms,data=v))$rank
 design_diag[[length(design_diag)+1]]<-data.frame(framework=fr,sex=sex,n=nrow(v),stroke_events=sum(v$stroke),
  base_n=nrow(base$variables),represented_strata=nrow(strat),represented_PSU=nrow(unit),lonely_PSU_strata=sum(strat$represented_psu==1),
  manual_domain_df=manual,degf_base=degf(base),degf_subset=degf(ss),nuisance_design_matrix_rank=nuisance_rank,
  hypothetical_default_residual_df_if_exposure_adds_one_rank=degf(ss)-nuisance_rank,
  actual_svyglm_residual_df="NOT_FITTED",primary_df_rule="MANUALLY_VERIFIED_DOMAIN_DESIGN_DF")
 unit$framework<-fr;unit$sex<-sex;units[[length(units)+1]]<-unit
 stopifnot(manual==degf(ss))
}
dd<-bind_rows(design_diag);wtcsv(dd,"survey_design_domain_diagnostics.csv");wtcsv(bind_rows(units),"domain_represented_design_units.csv")
body<-deparse(getS3method("svyglm","survey.design",envir=asNamespace("survey")))
writeLines(c(paste("survey version:",packageVersion("survey")),"No model was fitted.",
 "Default rule (unless special degf option): degf(design)+1-number of non-aliased coefficients.",
 "Hypothetical diagnostic assumes exposure contributes one additional independent column; its rank is deliberately not checked.",
 grep("df.residual|degf",body,value=TRUE)),file.path(out,"svyglm_df_rule_source.txt"))
# Included/excluded comparisons only for pre-CC eligible domain and non-exposure variables.
# Unweighted descriptive statistics and one-level-vs-rest SMD; no significance tests.
selection<-list()
descvars<-c("stroke","RIDAGEYR","sex","race_eth_f","edu_f","marital_f","pir","bmi","smoking_f","hypertension_f","diabetes_f","hyperlipidemia_f")
for(gr in groups) {
 g<-raw$precc&groupmask(gr)
 inc<-g&raw$cov_complete;exc<-g&!raw$cov_complete
 for(v in descvars) {
  x<-raw[[v]]
  levs<-if(is.factor(x))levels(x)else if(v=="stroke")"1"else NA_character_
  for(lv in levs) {
   z<-if(is.na(lv))as.numeric(x)else ifelse(is.na(x),NA_real_,as.numeric(as.character(x)==lv))
   a<-z[inc&!is.na(z)];b<-z[exc&!is.na(z)]
   ma<-if(length(a))mean(a)else NA_real_;mb<-if(length(b))mean(b)else NA_real_
   pooled<-if(is.na(lv))sqrt((var(a)+var(b))/2)else sqrt((ma*(1-ma)+mb*(1-mb))/2)
   smd<-if(is.finite(pooled)&&pooled>0)abs(ma-mb)/pooled else if(isTRUE(ma==mb))0 else NA_real_
   selection[[length(selection)+1]]<-data.frame(group=gr,variable=v,level=if(is.na(lv))"CONTINUOUS"else lv,
    included_total=sum(inc),excluded_total=sum(exc),included_nonmissing=length(a),excluded_nonmissing=length(b),
    included_mean_or_proportion=ma,excluded_mean_or_proportion=mb,
    included_SD=if(is.na(lv))sd(a)else NA_real_,excluded_SD=if(is.na(lv))sd(b)else NA_real_,absolute_SMD=smd,weighting="UNWEIGHTED_DESCRIPTIVE")
  }
 }
}
sel<-bind_rows(selection);wtcsv(sel,"included_excluded_selection.csv")
# LOD QC is source/flag validation only; concentration fields remain unchanged.
lod<-read.csv(file.path(package_root,"config/lod_exact_comparison.csv"),check.names=FALSE)
lod$adjudicated_status<-ifelse(lod$exact_mismatch_n>0,"RELEASED_VALUE_RETAINED_RECONSTRUCTION_DIFFERENCE_NONBLOCKING","NO_RECONSTRUCTION_DIFFERENCE_OR_NO_BELOW_LOD_RECORDS")
lod$secondary_reimputation<-"NOT_PERFORMED"
det<-list()
for(i in 1:8) {
 ph<-haven::read_xpt(file.path(inputs,"official_xpt",paste0(cyc$ph[i],".xpt")))
 for(j in seq_along(met))if(names(met)[j]%in%names(ph)) {
  val<-as.numeric(ph[[names(met)[j]]]);flag<-as.numeric(ph[[flags[j]]])
  det[[length(det)+1]]<-data.frame(cycle=cyc$cycle[i],variable=names(met)[j],component_n=nrow(ph),
   nonmissing_value_n=sum(!is.na(val)),missing_value_n=sum(is.na(val)),nonmissing_flag_n=sum(!is.na(flag)),
   below_LOD_n=sum(flag==1,na.rm=TRUE),above_LOD_n=sum(flag==0,na.rm=TRUE),
   detection_proportion=sum(flag==0,na.rm=TRUE)/sum(!is.na(flag)),
   flag1_value_missing_n=sum(flag==1&is.na(val),na.rm=TRUE))
 }
}
wtcsv(lod,"lod_adjudicated_qc.csv");wtcsv(bind_rows(det),"detection_and_analyte_availability.csv")
# Data foundation is a candidate only. No logarithms, regressions or target cross-tabs generated.
attr(raw,"IDENTITY")<-"PHASE_A_DATA_DESIGN_CANDIDATE / NOT_STATISTICS_AUTHORITY / PHASE_B_NOT_AUTHORIZED"
attr(raw,"LOD_RULE")<-"OFFICIAL_RELEASED_ANALYTE_VALUE=SOURCE_OF_TRUTH; SECONDARY_LOD_SQRT2_REIMPUTATION=NOT_PERFORMED"
saveRDS(raw,file.path(out,"phase_a_data_candidate.rds"))
saveRDS(list(base8=base8,base7=base7),file.path(out,"phase_a_base_designs_candidate.rds"))
# Missing-data decision is evaluated after all permitted descriptive diagnostics are written.
gate_rows<-mi[mi$true_item_or_component_proportion>0.10,]
gate_smd<-sel[!is.na(sel$absolute_SMD)&sel$absolute_SMD>=0.10,]
wtcsv(gate_rows,"missing_data_gate_over10pct.csv");wtcsv(gate_smd,"missing_data_gate_SMD_ge010.csv")
triggered<-nrow(gate_rows)>0||nrow(gate_smd)>0
jsonlite::write_json(list(status=if(triggered)"STOPPED_AT_DECISION_BOUNDARY"else"PHASE_A_COMPLETE / RETURN_FOR_ADJUDICATION",
 missing_data_gate=if(triggered)"MISSING_DATA_DECISION_REQUIRED"else"REQUIRES_AUTHOR_REVIEW_OF_NONNUMERIC_CRITERIA",
 over10pct_rows=nrow(gate_rows),SMD_ge010_rows=nrow(gate_smd),
 candidate_main_n=sum(raw$domain_main),candidate_main_events=sum(raw$stroke[raw$domain_main]),
 candidate_common_n=sum(raw$domain_common),candidate_common_events=sum(raw$stroke[raw$domain_common]),
 corrected_exposure_stroke_models_run="NO",corrected_OR_CI_P_q_generated="NO",phase_B_executed="NO",
 explanatory_note="No MI/IPW selection. Important-variable SMD, cycle clustering and stroke prevalence are reported for adjudication; no scientific threshold invented for prevalence."),
 file.path(out,"phase_a_machine_summary.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,"data_design_sessionInfo.txt"))
cat("Upstream anchors:\n");print(anchor)
cat("Historical / corrected counts:\n");print(h)
cat("Design diagnostics:\n");print(dd)
cat("Missing-data gate:",if(triggered)"MISSING_DATA_DECISION_REQUIRED"else"RETURN_FOR_REVIEW","\n")
