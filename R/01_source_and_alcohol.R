# Source and fixed alcohol mapping verification; no outcome/exposure associations.
out <- normalizePath(commandArgs(TRUE)[1],winslash="/",mustWork=TRUE)
run<-dirname(out); inputs<-normalizePath(commandArgs(TRUE)[2],winslash="/",mustWork=TRUE); docs<-file.path(inputs,"official_docs")
rawroot<-file.path(inputs,"translated")
options(stringsAsFactors=FALSE)
write_tab <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na="")
cycles <- data.frame(suffix=LETTERS[3:10],year=seq(2003,2017,2))
cycles$cycle <- paste(cycles$year,cycles$year+1,sep="_")
cycles$ph <- c("L24PH_C",paste0("PHTHTE_",LETTERS[4:10]))
met <- c(URXMBP="MBP",URXMZP="MBzP",URXECP="MECPP",URXMHH="MEHHP",URXMOH="MEOHP",URXMC1="MCPP",URXMEP="MEP",URXMIB="MiBP",URXCNP="MCNP",URXCOP="MCOP")
flags <- c(URXMBP="URDMBPLC",URXMZP="URDMZPLC",URXECP="URDECPLC",URXMHH="URDMHHLC",URXMOH="URDMOHLC",URXMC1="URDMC1LC",URXMEP="URDMEPLC",URXMIB="URDMIBLC",URXCNP="URDCNPLC",URXCOP="URDCOPLC")
manifest <- list()
for(i in 1:8) {
 cy<-cycles[i,]
 components<-c(paste0(c("DEMO","MCQ","BPQ","DIQ","SMQ","ALQ","BMX"),"_",cy$suffix),cy$ph)
 if(i>=7) components<-c(components,paste0("ALB_CR_",cy$suffix))
 for(comp in components) {
  p<-file.path(rawroot,cy$cycle,paste0(comp,".rds"))
  stopifnot(file.exists(p))
  manifest[[length(manifest)+1]]<-data.frame(cycle=cy$cycle,component=comp,path=p,bytes=file.info(p)$size,sha256=toupper(digest::digest(file=p,algo="sha256")))
 }
}
write_tab(do.call(rbind,manifest),"raw_input_manifest.csv")
check <- list(); mapping <- list(); alco <- list(); docs_check <- list()
record <- function(cy,comp,var,local,official) {
 same <- (is.na(local)&is.na(official)) | (!is.na(local)&!is.na(official)&local==official)
 mismatch <- sum(!same)
 check[[length(check)+1]] <<- data.frame(cycle=cy,component=comp,variable=var,n=length(local),mismatch_n=mismatch)
 if(mismatch>0) {
  write_tab(do.call(rbind,check),"source_value_comparison.csv")
  stop("STOP: local raw vs official XPT mismatch: ",cy," ",comp," ",var)
 }
}
for(i in 1:8) {
 cy<-cycles[i,]
 ph<-readRDS(file.path(rawroot,cy$cycle,paste0(cy$ph,".rds")))
 x<-haven::read_xpt(file.path(inputs,"official_xpt",paste0(cy$ph,".xpt")))
 stopifnot(!anyDuplicated(ph$SEQN),!anyDuplicated(x$SEQN),setequal(ph$SEQN,x$SEQN))
 ph<-ph[match(x$SEQN,ph$SEQN),]
 wt<-if(i==5)"WTSA2YR" else "WTSB2YR"
 for(v in intersect(c(names(met),"URXUCR",wt),names(x))) record(cy$cycle,cy$ph,v,as.numeric(ph[[v]]),as.numeric(x[[v]]))
 for(v in intersect(unname(flags),names(x))) {
  z<-as.character(ph[[v]])
  stopifnot(all(is.na(z)|z %in% c("At or above the detection limit","Below lower detection limit")))
  zz<-ifelse(z=="Below lower detection limit",1,ifelse(z=="At or above the detection limit",0,NA_real_))
  record(cy$cycle,cy$ph,v,zz,as.numeric(x[[v]]))
 }
 # Fresh official alcohol XPT is compared to the translated local cache.
 comp<-paste0("ALQ_",cy$suffix)
 a<-haven::read_xpt(file.path(inputs,"official_xpt",paste0(comp,".xpt")))
 al<-readRDS(file.path(rawroot,cy$cycle,paste0(comp,".rds")))
 stopifnot(!anyDuplicated(a$SEQN),!anyDuplicated(al$SEQN),setequal(a$SEQN,al$SEQN))
 al<-al[match(a$SEQN,al$SEQN),]
 doc<-xml2::read_html(file.path(docs,paste0(comp,".html")))
 vars<-if(i<8)c("ALQ101","ALQ110","ALQ120Q","ALQ120U")else c("ALQ111","ALQ121")
 for(v in vars) {
  stopifnot(v %in% names(a),v %in% names(al))
  sec<-xml2::xml_find_first(doc,paste0("//*[@id='",v,"']/parent::*"))
  tab<-rvest::html_table(rvest::html_elements(sec,"table"),fill=TRUE)[[1]]
  code<-suppressWarnings(as.numeric(tab[[1]]))
  label<-trimws(as.character(tab[[2]]))
  docs_check[[length(docs_check)+1]]<-data.frame(cycle=cy$cycle,variable=v,codes=paste(tab[[1]],collapse="|"),labels=paste(label,collapse="|"),skip=paste(tab[[ncol(tab)]],collapse="|"))
  if(is.factor(al[[v]])) {
   z<-as.character(al[[v]]); converted<-code[match(z,label)]
   if(any(!is.na(z)&is.na(converted))) stop("Unmatched official alcohol label ",cy$cycle," ",v)
  } else converted<-as.numeric(al[[v]])
  record(cy$cycle,comp,v,converted,as.numeric(a[[v]]))
 }
 aa<-function(v) as.numeric(a[[v]])
 n<-nrow(a); value<-rep(NA_character_,n); reason<-rep("ITEM_OR_UNEXPLAINED_MISSING",n)
 if(i<8) {
  # Route verified against each codebook: ALQ101=1 -> frequency;
  # otherwise ALQ110=1 -> frequency; ALQ110=2 -> end of section.
  z1<-aa("ALQ101"); z2<-aa("ALQ110"); q<-aa("ALQ120Q"); u<-aa("ALQ120U")
  sec<-xml2::xml_find_first(doc,"//*[@id='ALQ110']/parent::*")
  tb<-rvest::html_table(rvest::html_elements(sec,"table"),fill=TRUE)[[1]]
  stopifnot(grepl("End of Section",as.character(tb[[ncol(tb)]][as.character(tb[[1]])=="2"]),fixed=TRUE))
  route <- z1 %in% 1 | (!(z1 %in% 1)&z2 %in% 1)
  proxy <- !(z1 %in% 1)&z2 %in% 2&is.na(q)
  zero <- route&!is.na(q)&q==0
  positive <- route&!is.na(q)&q>0&q<=366
  value[proxy|zero]<-"No"; value[positive]<-"Yes"
  reason[proxy]<-"STRUCTURAL_PROXY_LT12_LIFETIME";reason[zero]<-"OBSERVED_PAST_YEAR_ZERO";reason[positive]<-"OBSERVED_POSITIVE_FREQUENCY"
  refused <- q %in% c(777,999) | (!route&z2 %in% c(7,9))
  reason[refused]<-"REFUSED_OR_DONT_KNOW";value[refused]<-NA_character_
  impossible <- (z1 %in% 1 & !is.na(z2)) | (!route & !is.na(q)) | (route&!is.na(q)&!(q %in% c(777,999))&(q<0|q>366))
  # Do not force a zero for a non-structural missing response.
  value[impossible]<-NA_character_;reason[impossible]<-"UNEXPLAINED_ROUTE"
  src1<-z1;src2<-z2;freq<-q;unit<-u
 } else {
  z1<-aa("ALQ111"); q<-aa("ALQ121")
  never<-z1 %in% 2&is.na(q);zero<-z1 %in% 1&q %in% 0;positive<-z1 %in% 1&q %in% 1:10
  value[never|zero]<-"No";value[positive]<-"Yes"
  reason[never]<-"LIFETIME_NEVER_STRUCTURAL";reason[zero]<-"OBSERVED_PAST_YEAR_ZERO";reason[positive]<-"OBSERVED_POSITIVE_FREQUENCY"
  refused<-z1 %in% c(7,9)|q %in% c(77,99)
  value[refused]<-NA_character_;reason[refused]<-"REFUSED_OR_DONT_KNOW"
  impossible<-(!(z1 %in% 1)&!is.na(q))|(z1 %in% 1&!is.na(q)&!(q %in% c(0:10,77,99)))
  value[impossible]<-NA_character_;reason[impossible]<-"UNEXPLAINED_ROUTE"
  src1<-z1;src2<-rep(NA_real_,n);freq<-q;unit<-rep(NA_real_,n)
 }
 row<-data.frame(SEQN=as.numeric(a$SEQN),cycle=cy$cycle,alcohol_harmonized_f=value,alcohol_reason=reason)
 alco[[i]]<-row
 path<-data.frame(cycle=cy$cycle,screen_variable=vars[1],screen_value=src1,lifetime_variable=if(i<8)"ALQ110" else "ALQ111",lifetime_value=if(i<8)src2 else src1,
 frequency_variable=if(i<8)"ALQ120Q" else "ALQ121",frequency_value=freq,frequency_unit=unit,harmonized_value=ifelse(is.na(value),"NA",value),reason=reason)
 path[is.na(path)]<-"<MISSING>"
 path$n<-1L
 mapping[[i]]<-aggregate(n~.,data=path,sum)
}
write_tab(do.call(rbind,check),"source_value_comparison.csv")
write_tab(do.call(rbind,docs_check),"alcohol_codebook_validation.csv")
write_tab(do.call(rbind,mapping),"alcohol_mapping_paths.csv")
a<-do.call(rbind,alco)
write_tab(a,"alcohol_harmonized_candidate.csv")
saveRDS(a,file.path(out,"alcohol_harmonized_candidate.rds"))
cat("SOURCE_VALUE_COMPARISON PASS: ",length(check)," fields; no outcome data read.\n",sep="")
cat("Alcohol mapping path counts:\n")
print(with(a,table(cycle,alcohol_reason)))
capture.output(sessionInfo(),file=file.path(out,"source_mapping_sessionInfo.txt"))
