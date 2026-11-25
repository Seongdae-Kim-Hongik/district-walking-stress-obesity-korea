# +1 offset sensitivity of the green-zone covariate: the main model refitted with only that covariate changed (five fits, same seed).
# Sources 81fh_obes.R without running its steps. Usage: PC_RES=<folder with green_variants.csv> ST_SHP=... ST_GRAPH=... ST_THREADS=1 PC_ONLY=<list> Rscript 84_refit_green_plus1.R OBES
PCRES<-Sys.getenv("PC_RES");stopifnot(nzchar(PCRES))
Sys.setenv(ST_FH_ONLY="NONE",ST_RES=PCRES)
HERE<-dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))))
source(file.path(HERE,"81fh_obes.R"))
RES<-PCRES
GV<-read.csv(file.path(RES,"green_variants.csv"),stringsAsFactors=FALSE);GV$fbd_region<-enc2utf8(GV$fbd_region)
ix<-match(paste(DD$year,DD$fbd_region),paste(GV$year,GV$fbd_region));stopifnot(!anyNA(ix));GV<-GV[ix,]
GK<-"eloc__dt_1yl20421e__1___1";stopifnot(max(abs(GV$green_p1-DD[[GK]]))<1e-9)
FITS<-list(
  orig_rerun=list(col="green_p1",drop=FALSE),
  A_halfmin=list(col="green_halfmin",drop=FALSE),
  B_zeroNA_fill=list(col="green_zeroNA",drop=FALSE),
  B2_zeroNA_drop=list(col="green_zeroNA",drop=TRUE),
  C_plus1_drop16=list(col="green_p1_raw",drop=TRUE))
SEL<-Sys.getenv("PC_ONLY");if(nzchar(SEL))FITS<-FITS[strsplit(SEL,",")[[1]]]
for(nm in names(FITS)){fo<-file.path(RES,sprintf("fit_%s_OR.csv",nm));if(file.exists(fo))next
  d<-DD;d[[GK]]<-GV[[FITS[[nm]]$col]]
  if(FITS[[nm]]$drop){k<-GV$green_zero_flag==0;d<-d[k,];x<-d[[GK]];d[[GK]]<-(x-mean(x))/sqrt(mean((x-mean(x))^2));d$id<-seq_len(nrow(d))}
  set.seed(20261110)
  m<-timed(sprintf("PC %s (n=%d)",nm,nrow(d)),safe(MAINV,data=d));stopifnot(!is.null(m))
  ir<-tab(m,MAINV,list(fit=nm,n_obs=nrow(d)));write.csv(ir,fo,row.names=FALSE)
  write.csv(cbind(fit=nm,n_obs=nrow(d),fstat(m),mlik=round(m$mlik[1,1],2)),file.path(RES,sprintf("fit_%s_stat.csv",nm)),row.names=FALSE)
  rm(m);gc()}
writeLines(capture.output(sessionInfo()),file.path(RES,"sessionInfo.txt"));LOG("plus1 all done")
