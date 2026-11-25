# SUPERSEDED in v1.2.0: replaced by 01_스크립트/84_refit_green_plus1.R, which sources 81fh_obes.R instead of copying its model code.
# Sensitivity of the main model to the +1 offset in the log green-zone covariate: the main model refitted with only that covariate changed. Model code is the same as 81fh_obes.R.
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);TAG<-args[1];stopifnot(TAG=="OBES")
ROOT<-Sys.getenv("X_ROOT")
SHP<-Sys.getenv("ST_SHP")
SUF<-"_std3";OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,"results_plus1_check")
NT<-as.integer(Sys.getenv("ST_THREADS",unset="1"));inla.setOption(num.threads=sprintf("%d:1",NT))
ONLY<-Sys.getenv("ST_FH_ONLY",unset="");RUN<-function(s)!nzchar(ONLY)||s%in%strsplit(ONLY,",")[[1]]
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE);RU<-S$rules
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
META<-read.csv(file.path(OUTD,sprintf("panel_META_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
LAB<-setNames(META$label_en,META$variable);for(k in names(c("eloc__dt_1yl20421e__1___1"=1,"ext__restaurant_per10k"=1,"ext__instant_food_per10k"=1)))LAB[k]<-c("eloc__dt_1yl20421e__1___1"="Green zone area per capita (log)","ext__restaurant_per10k"="General restaurants per 10,000 (log)","ext__instant_food_per10k"="Instant-food manufacturing-retail outlets per 10,000 (log)")[k]
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
FO<-read.csv(file.path(ROOT,"02_패널","OBES_fh_outcomes.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);FO$fbd_region<-enc2utf8(FO$fbd_region)
DD<-merge(DD,FO[,c("year","fbd_region","SR","SR_SE","N","CR","CR_SE","y_logit","v_logit","y_logit_cr","v_logit_cr","n_eff","y_eff","meas_SR","meas_SE","y_logit_meas","v_logit_meas","report_gap")],by=c("year","fbd_region"),all.x=TRUE,sort=FALSE)
stopifnot(nrow(DD)==229*13,!any(DD$year==2019),all(!is.na(DD$y_logit)),max(abs(DD$SR-DD$cases))<1e-9)
FORCED<-S$forced;NEG<-S$negative_controls;SETS<-lapply(S$sets_after_vif,as.character);SETS<-SETS[setdiff(names(SETS),"S1b_obesity_sensitivity")];stopifnot(!("obesity_rate"%in%unlist(SETS)));HYPV<-SETS[[RU$main_set]];MAINV<-unique(c(HYPV,FORCED))
CAND<-setdiff(unique(unlist(SETS)),FORCED);NMV<-RU$multiverse_n;if(nzchar(Sys.getenv("ST_MV_N")))NMV<-as.integer(Sys.getenv("ST_MV_N"))
WIN<-INFO$windows;POSTY<-RU$covid_post_from
LOG<-function(...)cat(sprintf("[81fh %s %s] ",TAG,format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")
stopifnot(all(c(MAINV,NEG)%in%names(DD)))

G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);stopifnot(all(unique(DD$fbd_region)%in%G$region))
nbq<-poly2nb(G,queen=TRUE);isl<-which(card(nbq)==0)
if(length(isl)){cen<-st_coordinates(st_centroid(st_geometry(G)));kn<-knearneigh(cen,k=1)$nn
  for(i in isl){j<-kn[i,1];nbq[[i]]<-sort(unique(c(nbq[[i]][nbq[[i]]>0],j)));nbq[[j]]<-sort(unique(c(nbq[[j]][nbq[[j]]>0],i)))}}
nb<-nbq;GP<-Sys.getenv("ST_GRAPH");stopifnot(nzchar(GP))
sg<-inla.read.graph(GP);stopifnot(sg$n==nrow(G));for(i in seq_len(sg$n))nb[[i]]<-as.integer(sort(sg$nbs[[i]]));attr(nb,"sym")<-TRUE
GF<-file.path(tempdir(),sprintf("nb_%s_%d.graph",TAG,Sys.getpid()));nb2INLA(GF,nb);g<-inla.read.graph(GF);lw<-nb2listw(nb,style="W",zero.policy=TRUE)
GFQ<-file.path(tempdir(),sprintf("nbqueen_%s_%d.graph",TAG,Sys.getpid()));nb2INLA(GFQ,nbq);gq<-inla.read.graph(GFQ)
n<-nrow(G);mI<-function(x,alt="greater"){mm<-moran.mc(x,lw,999,zero.policy=TRUE,alternative=alt);c(I=unname(mm$statistic),p=mm$p.value,EI=-1/(length(x)-1))}

DD$area<-match(DD$fbd_region,G$region);DD<-DD[order(DD$year,DD$area),];Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time;DD$id<-seq_len(nrow(DD))
DD$w<-1/DD$v_logit;DD$w_cr<-1/DD$v_logit_cr;DD$p<-DD$SR/100
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
bt<-function(v)sprintf("`%s`",v)
SPf<-function(gr="g",hp="pc")sprintf("f(area,model='bym2',graph=%s,scale.model=TRUE,constr=TRUE,hyper=%s)",gr,hp)
TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(yiid,model='iid',hyper=pcp)";CELL<-"f(id,model='iid',hyper=pcp)"
FIXP<-list(hyper=list(prec=list(initial=0,fixed=TRUE)))
fitfh<-function(vars,data=DD,lik="fh",extra=character(0),re=TRUE,cell=TRUE,cpo=TRUE,sp=SPf(),robust=FALSE){
  yv<-switch(lik,fh="y_logit",gfree="y_logit",cr="y_logit_cr",meas="y_logit_meas",binom="y_eff",beta="p")
  usecell<-cell&&lik%in%c("fh","cr","meas","binom")
  TIv<-sprintf("f(time,model='rw1',scale.model=TRUE,hyper=pcp,values=1:%d)",max(data$time))
  rhs<-c("1",bt(vars),extra,if(re)c(sp,TIv,IT),if(usecell)CELL);form<-as.formula(paste(yv,"~",paste(rhs,collapse="+")))
  ci<-if(robust)list(control.vb=list(emergency=50),int.strategy="eb",strategy="gaussian") else list(control.vb=list(emergency=50))
  cc<-list(dic=TRUE,waic=TRUE,cpo=cpo,config=FALSE)
  if(lik%in%c("fh","cr","meas")){sc<-switch(lik,fh=data$w,cr=data$w_cr,meas=1/data$v_logit_meas)
    return(inla(form,family="gaussian",data=data,scale=sc,control.family=FIXP,control.fixed=CF,control.compute=cc,control.inla=ci,verbose=FALSE))}
  if(lik=="gfree")return(inla(form,family="gaussian",data=data,control.fixed=CF,control.compute=cc,control.inla=ci,verbose=FALSE))
  if(lik=="binom")return(inla(form,family="binomial",Ntrials=data$n_eff,data=data,control.fixed=CF,control.compute=cc,control.inla=ci,verbose=FALSE))
  if(lik=="beta")return(inla(form,family="beta",data=data,control.fixed=CF,control.compute=cc,control.inla=ci,verbose=FALSE))
  stop(lik)}
safe<-function(...)tryCatch(fitfh(...),error=function(e){message("  ✗ fit failed: ",conditionMessage(e)," → retry robust");tryCatch(fitfh(...,robust=TRUE),error=function(e2){message("  ✗ retry failed: ",conditionMessage(e2));NULL})})
orr<-function(m,v){r<-m$summary.fixed[v,,drop=FALSE];data.frame(cov=v,OR=round(exp(r$mean),4),lo=round(exp(r$`0.025quant`),4),hi=round(exp(r$`0.975quant`),4),
  P_gt1=round(1-INLA::inla.pmarginal(0,m$marginals.fixed[[v]]),4),cred=(r$`0.025quant`>0|r$`0.975quant`<0),stringsAsFactors=FALSE)}
fstat<-function(m){cp<-m$cpo$cpo;hp<-m$summary.hyperpar;data.frame(WAIC=round(m$waic$waic,1),DIC=round(m$dic$dic,1),LCPO=if(is.null(cp))NA else round(-mean(log(cp[is.finite(cp)&cp>0])),4),
  phi=round(if("Phi for area"%in%rownames(hp))hp["Phi for area","mean"] else NA,3),
  sd_area=round(if("Precision for area"%in%rownames(hp))1/sqrt(hp["Precision for area","mean"]) else NA,4),sd_cell=round(if("Precision for id"%in%rownames(hp))1/sqrt(hp["Precision for id","mean"]) else NA,4))}
tab<-function(m,vars,extra=list()){r<-do.call(rbind,lapply(vars,function(v)orr(m,v)));r$label_en<-LAB[r$cov];r$role<-ifelse(r$cov%in%FORCED,"forced","hypothesis");for(k in names(extra))r[[k]]<-extra[[k]];r}
timed<-function(label,expr){t0<-Sys.time();r<-expr;LOG("%s · %.1f min",label,as.numeric(difftime(Sys.time(),t0,units="mins")));r}
clusters<-function(nodes){nodes<-as.integer(nodes);seen<-integer(0);best<-0;for(s in nodes){if(s%in%seen)next;st<-s;comp<-integer(0)
  while(length(st)){v<-st[1];st<-st[-1];if(v%in%seen)next;seen<-c(seen,v);comp<-c(comp,v);st<-c(st,setdiff(intersect(nb[[v]],nodes),seen))};best<-max(best,length(comp))};best}
sidoOf<-sapply(G$region,function(r)sub("^(서울시|부산시|대구시|인천시|광주시|대전시|울산시|세종시|경기도|강원도|충청북도|충청남도|전라북도|전라남도|경상북도|경상남도|제주도).*$","\\1",r))
spatialOR<-function(m){ur<-m$summary.random$area[1:n,];pep<-sapply(m$marginals.random$area[1:n],function(mg)1-inla.pmarginal(0,mg))
  rr<-data.frame(area=1:n,fbd_region=G$region,sido=sidoOf,OR=round(exp(ur$mean),4),lo=round(exp(ur$`0.025quant`),4),hi=round(exp(ur$`0.975quant`),4),P_gt1=round(pep,4))
  rr$class<-ifelse(rr$P_gt1>0.95,"hotspot",ifelse(rr$P_gt1<0.05,"coldspot","neutral"));rr}

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
writeLines(capture.output(sessionInfo()),file.path(RES,"sessionInfo.txt"));LOG("all done")
