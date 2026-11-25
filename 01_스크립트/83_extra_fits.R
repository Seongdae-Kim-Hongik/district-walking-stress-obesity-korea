# Additional analyses A-F (ICAR and non-spatial models, two-way decomposition, facilities per km2, boundary-change exclusions, imputed-cell exclusions and complete cases, multiple imputation).
# Sources 81fh_obes.R without running its steps. Usage: EX_RES=<result folder with extra_inputs.csv> ST_SHP=... ST_GRAPH=... ST_THREADS=1 EX_ONLY=<list> Rscript 83_extra_fits.R OBES
EXRES<-Sys.getenv("EX_RES");stopifnot(nzchar(EXRES))
Sys.setenv(ST_FH_ONLY="NONE",ST_RES=EXRES)
HERE<-dirname(normalizePath(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))))
source(file.path(HERE,"81fh_obes.R"))
RES<-EXRES;WD<-file.path(RES,"waic_local");dir.create(WD,showWarnings=FALSE)
EX<-read.csv(file.path(RES,"extra_inputs.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);EX$fbd_region<-enc2utf8(EX$fbd_region)
ix<-match(paste(DD$year,DD$fbd_region),paste(EX$year,EX$fbd_region));stopifnot(!anyNA(ix));EX<-EX[ix,]
for(kk in setdiff(names(EX),c("year","fbd_region")))DD[[kk]]<-EX[[kk]]
PUB<-"public_sports_facilities_per_100k";VIL<-"village_sports_facilities_per_100k";MV<-"moderate_pa_rate";TE<-"tertiary_educated_share";BL<-"basic_livelihood_rate"
LAB["pub_km2_z"]<-"Public sports facilities per km2 (log)";LAB["vil_km2_z"]<-"Village sports facilities per km2 (log)"
ICAR<-"f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=PRI_IID)"
reid<-function(d){d$id<-seq_len(nrow(d));d}
bm<-function(m,v){r<-m$summary.fixed[v,,drop=FALSE];data.frame(beta_mean=r$mean)}
FITS<-list(
  orig_rerun=list(vars=MAINV),
  A1_icar=list(vars=MAINV,sp=ICAR),
  A2_nospatial=list(vars=MAINV,sp="NONE"),
  B_twoway=list(twoway=TRUE),
  C1_km2=list(vars=c(setdiff(HYPV,c(PUB,VIL)),"pub_km2_z","vil_km2_z",FORCED)),
  C2_both=list(vars=c(HYPV,"pub_km2_z","vil_km2_z",FORCED)),
  D1_bnd_terr=list(vars=MAINV,keep=function(d)d$bnd_terr==0),
  D2_bnd_any=list(vars=MAINV,keep=function(d)d$bnd_any==0),
  E1_drop_mvpa=list(vars=setdiff(MAINV,MV)),
  E1cc_mvpa=list(vars=MAINV,keep=function(d)d$imp_mvpa==0),
  E2_drop_fac=list(vars=setdiff(MAINV,c(PUB,VIL))),
  E2cc_fac=list(vars=MAINV,keep=function(d)d$imp_fac==0&d$imp_vil==0),
  E3_drop_tert=list(vars=setdiff(MAINV,TE)),
  E3cc_tert=list(vars=MAINV,keep=function(d)d$imp_tert==0),
  E4_drop_bl=list(vars=setdiff(MAINV,BL)))
MIM<-as.integer(Sys.getenv("MI_M",unset="20"))
for(mm in seq_len(MIM))FITS[[sprintf("F_MI_%02d",mm)]]<-list(vars=MAINV,mi=mm)
SEL<-Sys.getenv("EX_ONLY");if(nzchar(SEL)){p<-strsplit(SEL,",")[[1]];FITS<-FITS[names(FITS)%in%p|sapply(names(FITS),function(z)any(startsWith(z,sub("\\*$","",p[endsWith(p,"*")]))))]}
MIMAP<-list(imp_mvpa=MV,imp_fac=PUB,imp_vil=VIL,imp_tert=TE,imp_bl=BL)
for(nm in names(FITS)){fo<-file.path(RES,"fits",sprintf("fit_%s_OR.csv",nm));if(file.exists(fo))next;dir.create(dirname(fo),showWarnings=FALSE)
  Fi<-FITS[[nm]];d<-DD;extra<-character(0);vars<-Fi$vars;spx<-if(is.null(Fi$sp))SPf() else if(Fi$sp=="NONE")NULL else Fi$sp
  if(!is.null(Fi$keep))d<-reid(d[Fi$keep(d),])
  if(isTRUE(Fi$twoway)){ex<-character(0)
    for(cvar in HYPV){b<-ave(d[[cvar]],d$area);tt<-ave(d[[cvar]],d$year);d[[paste0(cvar,"__B")]]<-b;d[[paste0(cvar,"__W2")]]<-d[[cvar]]-b-tt+mean(d[[cvar]]);ex<-c(ex,bt(paste0(cvar,"__B")),bt(paste0(cvar,"__W2")))}
    extra<-ex;vars<-FORCED}
  if(!is.null(Fi$mi)){set.seed(20261125+Fi$mi);for(fl in names(MIMAP)){v<-MIMAP[[fl]];k<-d[[fl]]==1;d[[v]][k]<-d[[v]][k]+rnorm(sum(k),0,d[[paste0("sd_",v)]][k])}}
  set.seed(20261110)
  m<-timed(sprintf("EX %s (n=%d)",nm,nrow(d)),safe(vars,data=d,extra=extra,sp=spx));stopifnot(!is.null(m))
  if(nrow(d)==nrow(DD))SAVEW(m,paste0("EX_",nm),d)
  if(isTRUE(Fi$twoway)){ir<-NULL;for(cvar in HYPV)for(kk in c("B","W2")){vv<-paste0(cvar,"__",kk);r<-cbind(orr(m,vv),bm(m,vv));r$cov<-cvar;r$component<-ifelse(kk=="B","between","within_twoway");r$sd_component<-sd(d[[vv]]);ir<-rbind(ir,r)}
    ir$label_en<-LAB[ir$cov];ir$fit<-nm;ir$n_obs<-nrow(d)
    fr<-cbind(tab(m,FORCED),bm(m,FORCED));fr$component<-"forced";fr$sd_component<-NA;fr$fit<-nm;fr$n_obs<-nrow(d);ir<-rbind(ir,fr[,names(ir)])
  } else {ir<-cbind(tab(m,vars,list(fit=nm,n_obs=nrow(d))),bm(m,vars))}
  write.csv(ir,fo,row.names=FALSE)
  write.csv(cbind(fit=nm,n_obs=nrow(d),fstat(m),mlik=round(m$mlik[1,1],2)),file.path(RES,"fits",sprintf("fit_%s_stat.csv",nm)),row.names=FALSE)
  if(nm%in%c("A1_icar","orig_rerun")){ur<-m$summary.random$area;if(nrow(ur)>=n){so<-data.frame(area=1:n,fbd_region=G$region,OR=round(exp(ur$mean[1:n]),4));write.csv(so,file.path(RES,"fits",sprintf("spatial_%s.csv",nm)),row.names=FALSE)}}
  rm(m);gc()}
writeLines(capture.output(sessionInfo()),file.path(RES,"code","sessionInfo.txt"));LOG("extra all done")
