# Serialization only: no survey calculations or fitting.
options(stringsAsFactors=FALSE,digits=17)
w<-normalizePath(commandArgs(TRUE)[1],winslash='/',mustWork=TRUE)
o<-file.path(w,'03_OUTPUT/R1_FINAL');q<-file.path(w,'04_QC')
z<-readRDS(file.path(o,'table1_raw_objects.rds'))
sets<-list(table1_raw_descriptives_full_precision=z$raw,table1_test_details_full_precision=z$test_details,table1_canonical_full_precision=z$canonical)
checks<-list()
for(nm in names(sets)) {
 x<-sets[[nm]];xx<-x
 for(v in names(x))if(is.numeric(x[[v]]))xx[[v]]<-vapply(x[[v]],function(a)if(is.na(a))'' else sprintf('%.17g',a),character(1))
 p<-file.path(o,paste0(nm,'.csv'));stopifnot(!file.exists(p))
 write.csv(xx,p,row.names=FALSE,na='',fileEncoding='UTF-8')
 y<-read.csv(p,check.names=FALSE,na.strings='')
 for(v in names(x))if(is.numeric(x[[v]])){
  ok<-identical(as.numeric(x[[v]]),as.numeric(y[[v]]));stopifnot(ok)
  checks[[length(checks)+1]]<-data.frame(file=nm,column=v,exact_numeric_roundtrip=ok)
 }
}
write.csv(do.call(rbind,checks),file.path(q,'full_precision_serialization_QC.csv'),row.names=FALSE)
cat('Full precision CSVs: exact double-precision numeric roundtrip PASS. No recalculation.\n')
