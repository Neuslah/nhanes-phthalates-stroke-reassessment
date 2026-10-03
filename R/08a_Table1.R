# Task B-R1: only author-locked Table 1 descriptives and two-sample tests.
options(stringsAsFactors=FALSE, digits=17)
suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu='adjust', survey.adjust.domain.lonely=TRUE)
w <- normalizePath(commandArgs(TRUE)[1],winslash='/',mustWork=TRUE)
run <- dirname(w); package_root<-normalizePath(commandArgs(TRUE)[2],winslash='/',mustWork=TRUE)
a <- file.path(run,'resume_cholesterol_proxy');s <- file.path(run,'identity')
out <- file.path(w,'03_OUTPUT','R1_FINAL'); dir.create(out,showWarnings=FALSE); qcdir <- file.path(w,'04_QC')
stopifnot(!file.exists(file.path(out,'table1_canonical.csv')))
hash <- function(p) toupper(digest::digest(file=p,algo='sha256'))
write_csv <- function(z,p) write.csv(z,p,row.names=FALSE,na='',fileEncoding='UTF-8')
checks <- list()
check <- function(name,pass,detail='') {
 checks[[length(checks)+1L]] <<- data.frame(check=name,pass=isTRUE(pass),detail=detail)
 if(!isTRUE(pass)) {write_csv(do.call(rbind,checks),file.path(qcdir,'table1_R1_execution_QC.csv'));stop(name)}
}
warnings <- character()
withCallingHandlers({
d <- readRDS(file.path(a,'phase_a_data_candidate.rds'))
covs <- c('sex','RIDAGEYR','race_eth_f','edu_f','marital_f','pir','bmi','smoking_f','alcohol_harmonized_f','hypertension_f','diabetes_f','hyperlipidemia_f')
check('frozen_CC_indicator',identical(d$domain_main,d$precc & complete.cases(d[,covs])))
check('base_domain_membership',!anyNA(d$base_ok) && !anyNA(d$domain_main) && all(d$base_ok[d$domain_main]))
x <- d[d$base_ok,,drop=FALSE]
check('positive_weight_base',all(is.finite(x$ph_weight) & x$ph_weight>0) && !anyNA(x[,c('SDMVPSU','SDMVSTRA')]))
x$analysis_weight <- x$ph_weight/8
base <- svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~analysis_weight,data=x,nest=TRUE)
saved <- readRDS(file.path(a,'phase_a_base_designs_candidate.rds'))$base8
check('saved_base_membership_order',identical(base$variables$SEQN,saved$variables$SEQN))
check('saved_base_weights',isTRUE(all.equal(as.numeric(weights(base)),as.numeric(weights(saved)),tolerance=1e-14)))
check('saved_base_PSU_strata',identical(base$cluster,saved$cluster) && identical(base$strata,saved$strata))
des <- subset(base,domain_main)
dat <- des$variables
check('domain_counts',nrow(dat)==10528 && sum(dat$stroke==1)==385 && sum(dat$stroke==0)==10143)
check('Female_counts',sum(dat$sex=='Female')==5222 && sum(dat$sex=='Female' & dat$stroke==1)==198)
check('Male_counts',sum(dat$sex=='Male')==5306 && sum(dat$sex=='Male' & dat$stroke==1)==187)
check('unique_ids',!anyDuplicated(dat$SEQN))
for(sex in c('Female','Male')) {
 ids <- scan(file.path(s,paste0(sex,'_main_SEQN.txt')),quiet=TRUE)
 check(paste0(sex,'_frozen_ID_identity'),setequal(ids,dat$SEQN[dat$sex==sex]))
}
check('no_missing_Table1_variables',!anyNA(dat[,covs]) && !anyNA(dat$stroke))
check('no_borderline_in_CC',all(as.character(dat$DIQ010) %in% c('Yes','No')) && identical(as.character(dat$DIQ010),as.character(dat$diabetes_f)))
des <- update(des,stroke_f=factor(stroke,levels=c(0,1),labels=c('No stroke','Stroke')))
design_info <- function(z,name) {
 v <- z$variables; u <- unique(v[c('cycle','SDMVSTRA','SDMVPSU')]); st <- unique(v[c('cycle','SDMVSTRA')])
 counts <- table(interaction(u$cycle,u$SDMVSTRA,drop=TRUE))
 data.frame(domain=name,N=nrow(v),events=sum(v$stroke==1,na.rm=TRUE),PSU=nrow(u),strata=nrow(st),represented_df=nrow(u)-nrow(st),survey_degf=degf(z),lonely=sum(counts==1),weight_sum=sum(weights(z)),weight_variable='ph_weight / 8',nest=TRUE)
}
design_rows <- rbind(design_info(base,'positive_weight_base'),design_info(des,'primary_CC'),design_info(subset(des,stroke==1),'stroke'),design_info(subset(des,stroke==0),'no_stroke'))
check('primary_design_identity',all(design_rows[2,c('PSU','strata','represented_df','survey_degf','lonely')]==c(244,120,124,124,0)))
write_csv(design_rows,file.path(out,'table1_design_identity.csv'))
cv <- c(RIDAGEYR='Age (years)',bmi='BMI (kg/m²)',pir='Poverty-income ratio')
cats <- list(
 sex=list(label='Male (%)',levels=c('Female','Male'),display=c('Male'),labels=c('Female','Male'),ref='Female'),
 race_eth_f=list(label='Race/ethnicity',levels=c('Non-Hispanic White','Non-Hispanic Black','Mexican American','Other Hispanic','Other/Multiracial'),display=c('Non-Hispanic White','Non-Hispanic Black','Mexican American','Other Hispanic','Other/Multiracial'),labels=c('Non-Hispanic White','Non-Hispanic Black','Mexican American','Other Hispanic','Other/Multiracial'),ref='Non-Hispanic White'),
 edu_f=list(label='Education',levels=c('College or above','High school graduate','Less than high school'),display=c('College or above','High school graduate','Less than high school'),labels=c('Some college or above','High school graduate','Less than high school'),ref='College or above'),
 marital_f=list(label='Marital status',levels=c('Married/partnered','Never married','Separated/divorced/widowed'),display=c('Married/partnered','Never married','Separated/divorced/widowed'),labels=c('Married/partnered','Never married','Separated/divorced/widowed'),ref='Married/partnered'),
 smoking_f=list(label='Ever smoking, %',levels=c('Never','Ever'),display='Ever',labels=c('Never','Ever'),ref='Never'),
 alcohol_harmonized_f=list(label='Alcohol use, yes (%)',levels=c('No','Yes'),display='Yes',labels=c('No','Yes'),ref='No'),
 hypertension_f=list(label='Hypertension, yes (%)',levels=c('No','Yes'),display='Yes',labels=c('No','Yes'),ref='No'),
 diabetes_f=list(label='Diabetes, yes (%)',levels=c('No','Yes'),display='Yes',labels=c('No','Yes'),ref='No'),
 hyperlipidemia_f=list(label='Hyperlipidemia, yes (%)',levels=c('No','Yes'),display='Yes',labels=c('No','Yes'),ref='No')
)
coding <- list()
for(v in names(cats)) {
 cfg <- cats[[v]]
 check(paste0(v,'_level_identity'),setequal(levels(dat[[v]]),cfg$levels) && all(as.character(dat[[v]]) %in% cfg$levels))
 check(paste0(v,'_reference'),levels(dat[[v]])[1]==cfg$ref)
 coding[[v]] <- data.frame(variable=v,level=cfg$levels,display_label=cfg$labels,is_reference=cfg$levels==cfg$ref,shown_in_Table1=cfg$levels %in% cfg$display)
}
write_csv(do.call(rbind,coding),file.path(out,'table1_category_definition.csv'))
raw <- list(); testrows <- list(); tests <- list()
for(v in c(names(cv),names(cats))) {
 if(v %in% names(cv)) tt <- svyttest(as.formula(paste(v,'~ stroke_f')),des)
 else tt <- svychisq(as.formula(paste('~',v,'+ stroke_f')),des,statistic='F')
 tests[[v]] <- tt
 params <- as.numeric(tt$parameter)
 testrows[[v]] <- data.frame(variable=v,method=if(v %in% names(cv))'survey::svyttest()' else 'survey::svychisq(statistic="F")',statistic=unname(tt$statistic),statistic_name=names(tt$statistic),numerator_df=if(v %in% names(cv))NA_real_ else params[1],denominator_df=tail(params,1),full_precision_P=tt$p.value)
 check(paste0(v,'_test_finite'),all(is.finite(c(tt$statistic,params,tt$p.value))) && tt$p.value>=0 && tt$p.value<=1)
 for(gr in c('Stroke','No stroke')) {
  dz <- if(gr=='Stroke')subset(des,stroke==1) else subset(des,stroke==0)
  xx <- dz$variables[[v]]; nn <- length(xx)
  if(v %in% names(cv)) {
   ff <- as.formula(paste('~',v));mu <- unname(coef(svymean(ff,dz,na.rm=TRUE))[1]);sdv <- sqrt(unname(coef(svyvar(ff,dz,na.rm=TRUE))[1]))
   raw[[length(raw)+1]] <- data.frame(variable=v,level='',group=gr,type='CONTINUOUS',N=nn,count=NA_integer_,mean=mu,SD=sdv,weighted_percent=NA_real_)
   check(paste0(v,'_',gr,'_finite'),is.finite(mu) && is.finite(sdv) && sdv>=0)
  } else {
   pp <- numeric()
   for(lv in cats[[v]]$levels) {
    iz <- as.numeric(xx==lv);dd <- update(dz,.indicator=iz)
    pct <- unname(coef(svymean(~.indicator,dd,na.rm=TRUE))[1])*100;pp <- c(pp,pct)
    raw[[length(raw)+1]] <- data.frame(variable=v,level=lv,group=gr,type='CATEGORICAL',N=nn,count=sum(xx==lv),mean=NA_real_,SD=NA_real_,weighted_percent=pct)
   }
   check(paste0(v,'_',gr,'_counts'),sum(vapply(cats[[v]]$levels,function(lv)sum(xx==lv),numeric(1)))==nn)
   check(paste0(v,'_',gr,'_percentage_closure'),all(is.finite(pp) & pp>=0 & pp<=100) && abs(sum(pp)-100)<1e-10)
  }
 }
}
raw <- do.call(rbind,raw);testrows <- do.call(rbind,testrows);rownames(testrows)<-NULL
check('tests_allowlist',nrow(testrows)==12 && sum(grepl('svychisq',testrows$method))==9 && sum(grepl('svyttest',testrows$method))==3)
format_p <- function(p) ifelse(p<0.001,'<0.001',sprintf('%.3f',p))
testrows$displayed_P <- format_p(testrows$full_precision_P)
hist <- read.csv(file.path(package_root,'config/table1-layout.csv'),check.names=FALSE)
can <- list(); displays <- list()
for(i in 2:nrow(hist)) {
 h <- hist[i,];label <- h$characteristic;v <- '';lv <- '';type <- '';hasp <- FALSE;ref <- FALSE
 if(label=='N') type <- 'N'
 else if(label %in% unname(cv)) {v <- names(cv)[match(label,cv)];type <- 'CONTINUOUS';hasp<-TRUE}
 else {
  parent <- names(cats)[vapply(cats,function(z)z$label==label,logical(1))]
  if(length(parent)) {v <- parent;hasp<-TRUE;if(length(cats[[v]]$display)>1)type<-'CATEGORY_HEADER' else {type<-'CATEGORICAL';lv<-cats[[v]]$display}}
  else for(z in names(cats)) if(label %in% cats[[z]]$labels) {v<-z;lv<-cats[[z]]$levels[match(label,cats[[z]]$labels)];type<-'CATEGORICAL';break}
 }
 check(paste0('row_',i,'_mapping'),nzchar(type))
 if(nzchar(lv))ref<-lv==cats[[v]]$ref
 rr <- data.frame(row_order=h$row_order,characteristic=label,variable=v,level=lv,row_type=type,is_reference=ref,stroke_N=NA_integer_,stroke_count=NA_integer_,stroke_mean=NA_real_,stroke_SD=NA_real_,stroke_weighted_percent=NA_real_,no_stroke_N=NA_integer_,no_stroke_count=NA_integer_,no_stroke_mean=NA_real_,no_stroke_SD=NA_real_,no_stroke_weighted_percent=NA_real_,test_statistic=NA_real_,numerator_df=NA_real_,denominator_df=NA_real_,full_precision_P=NA_real_)
 ds <- data.frame(row_order=h$row_order,Characteristic=label,Stroke='',No_stroke='',P_value='—',is_reference=ref)
 if(type=='N') {rr$stroke_N<-385;rr$no_stroke_N<-10143;ds$Stroke<-'385';ds$No_stroke<-'10,143'}
 if(type %in% c('CONTINUOUS','CATEGORICAL')) for(gr in c('Stroke','No stroke')) {
  zz<-raw[raw$variable==v & raw$level==lv & raw$group==gr,,drop=FALSE];check(paste0('row_',i,'_',gr,'_source'),nrow(zz)==1)
  prefix<-if(gr=='Stroke')'stroke_' else 'no_stroke_'
  rr[[paste0(prefix,'N')]]<-zz$N;rr[[paste0(prefix,'count')]]<-zz$count;rr[[paste0(prefix,'mean')]]<-zz$mean;rr[[paste0(prefix,'SD')]]<-zz$SD;rr[[paste0(prefix,'weighted_percent')]]<-zz$weighted_percent
  value<-if(type=='CONTINUOUS')sprintf('%.1f (%.1f)',zz$mean,zz$SD) else sprintf('%.1f',zz$weighted_percent)
  ds[[if(gr=='Stroke')'Stroke' else 'No_stroke']]<-value
 }
 if(hasp) {tt<-testrows[testrows$variable==v,];rr$test_statistic<-tt$statistic;rr$numerator_df<-tt$numerator_df;rr$denominator_df<-tt$denominator_df;rr$full_precision_P<-tt$full_precision_P;ds$P_value<-tt$displayed_P}
 can[[length(can)+1]]<-rr;displays[[length(displays)+1]]<-ds
}
can<-do.call(rbind,can);displays<-do.call(rbind,displays)
check('historical_display_order',identical(can$characteristic,hist$characteristic[-1]) && identical(can$row_order,hist$row_order[-1]))
check('one_P_per_variable',sum(!is.na(can$full_precision_P))==12 && !anyDuplicated(can$variable[!is.na(can$full_precision_P)]))
write_csv(raw,file.path(out,'table1_raw_descriptives.csv'))
write_csv(testrows,file.path(out,'table1_test_details.csv'))
write_csv(can,file.path(out,'table1_canonical.csv'))
write_csv(displays,file.path(out,'table1_display.csv'))
saveRDS(list(raw=raw,tests=tests,test_details=testrows,canonical=can,display=displays),file.path(out,'table1_raw_objects.rds'))
check('only_locked_domain_lonely_warnings',all(grepl('^Stratum \\(.*\\) has only one PSU at stage 1$',warnings)),paste('Recorded warnings:',length(warnings)))
write_csv(data.frame(warning_id=seq_along(warnings),message=warnings),file.path(out,'table1_warnings.csv'))
write_csv(do.call(rbind,checks),file.path(qcdir,'table1_R1_execution_QC.csv'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
capture.output({print(getS3method('svyvar','survey.design'));print(getS3method('svyttest','default'));print(RNGkind());print(.libPaths())},file=file.path(out,'installed_descriptive_implementation.txt'))
write_csv(data.frame(package=c('R','survey','digest'),version=c(as.character(getRversion()),as.character(packageVersion('survey')),as.character(packageVersion('digest')))),file.path(out,'software_versions.csv'))
writeLines('TABLE1_EXECUTION = RETURNED_PENDING_INDEPENDENT_QC',file.path(qcdir,'table1_execution_status.txt'))
cat('Table 1 only: 10528 / 385; 3 continuous tests; 9 uniform categorical F tests; pending independent QC.\n')
},warning=function(z){
 msg <- conditionMessage(z); warnings<<-c(warnings,msg)
 known <- grepl('^Stratum \\(.*\\) has only one PSU at stage 1$',msg) && identical(getOption('survey.lonely.psu'),'adjust') && isTRUE(getOption('survey.adjust.domain.lonely'))
 if(!known)stop(msg)
})
