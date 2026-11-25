# Round-2 additional fits (stages J1-J7): determinism rerun, binomial check, posterior intervals and predictive checks, unadjusted model, Moran's I under two nulls, persistence null,
# graph comparison, residual-Moran specifications, access and urbanicity, imputed covariates, instant-food years, regional slopes, calendar lags, null hot spots, split interactions, gate thresholds.
# Sources 01_스크립트/81fh_obes.R without running its steps (same data, functions, priors and fit ladder).
# Usage: E2_RES=<results_extra2> E2_MAINRES=<results_std3, read only> E2_RUNSRC=<01_스크립트/81fh_obes.R> E2_STAGE=J1..J7 ST_SHP=... ST_GRAPH=... ST_THREADS=1 ST_OLD65_COL=pop65_share_raw_fbd5 Rscript e2_fits.R OBES
E2<-Sys.getenv("E2_RES");stopifnot(nzchar(E2));STG<-Sys.getenv("E2_STAGE");stopifnot(nzchar(STG))
MAINRES<-Sys.getenv("E2_MAINRES");stopifnot(dir.exists(MAINRES))
RUNSRC<-Sys.getenv("E2_RUNSRC");stopifnot(file.exists(RUNSRC))
Sys.setenv(ST_FH_ONLY="NONE",ST_RES=file.path(E2,"_side"))
commandArgs<-function(trailingOnly=FALSE){a<-base::commandArgs(trailingOnly);if(trailingOnly)a else c(a[!grepl("^--file=",a)],paste0("--file=",RUNSRC))}
source(RUNSRC);rm(commandArgs)
RES<-E2;FD<-file.path(E2,"fits");dir.create(FD,showWarnings=FALSE);WD<-file.path(E2,"waic_local");dir.create(WD,showWarnings=FALSE)
EX<-read.csv(file.path(E2,"inputs_extra2.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);EX$fbd_region<-enc2utf8(EX$fbd_region)
ix<-match(paste(DD$year,DD$fbd_region),paste(EX$year,EX$fbd_region));stopifnot(!anyNA(ix));EX<-EX[ix,];for(kk in setdiff(names(EX),c("year","fbd_region")))DD[[kk]]<-EX[[kk]]
LAB[c("road_pavement_rate","outside_care_share","care_institutions_per_100k","hospitals_per_100k")]<-c("Road pavement rate (%)","Outside-district care share (%)","Healthcare institutions per 100k","Hospitals (general + hospital) per 100k")
ROB6<-c("walking_practice_rate","tertiary_educated_share","perceived_stress_rate","high_risk_drinking_rate_drinkers","ext__restaurant_per10k","current_smoking_rate")
reid<-function(d){d$id<-seq_len(nrow(d));d}
retime<-function(d){d<-d[order(d$year,d$area),];d$time<-d$year-min(d$year)+1;d$yiid<-d$time;d$id<-seq_len(nrow(d));d}
full<-function(d)nrow(d)==nrow(DD)
RF<-function(nm,vars,data=DD,extra=character(0),sp=SPf(),lik="fh",ti=NULL,cpo=FALSE,keepm=FALSE,note=""){
  fo<-file.path(FD,sprintf("fit_%s_OR.csv",nm));if(file.exists(fo)&&!keepm)return(invisible(NULL))
  set.seed(20261110);m<-timed(sprintf("E2 %s (n=%d)",nm,nrow(data)),safe(vars,data=data,extra=extra,sp=sp,lik=lik,ti=ti,cpo=cpo))
  if(is.null(m)){write.csv(data.frame(fit=nm,note="fit failed after ladder"),fo,row.names=FALSE);return(invisible(NULL))}
  ir<-if(length(vars))tab(m,vars,list(fit=nm,n_obs=nrow(data))) else data.frame(fit=nm,n_obs=nrow(data),note="no covariates");write.csv(ir,fo,row.names=FALSE)
  rm_<-if(lik=="fh")residMoran(m,data) else c(I=NA,p=NA)
  rr<-if(!is.null(sp)&&!is.null(m$summary.random$area)&&nrow(m$summary.random$area)>=n)spatialOR(m) else NULL
  st<-cbind(fit=nm,n_obs=nrow(data),fstat(m),resid_moran_I=round(unname(rm_["I"]),4),resid_moran_p=unname(rm_["p"]),n_hot=if(is.null(rr))NA else sum(rr$class=="hotspot"),n_cold=if(is.null(rr))NA else sum(rr$class=="coldspot"),note=note)
  write.csv(st,file.path(FD,sprintf("fit_%s_stat.csv",nm)),row.names=FALSE)
  if(!is.null(rr))write.csv(rr,file.path(FD,sprintf("spatial_%s.csv",nm)),row.names=FALSE)
  if(full(data))SAVEW(m,paste0("E2_",nm),data)
  if(keepm)return(m);rm(m);gc();invisible(NULL)}
LOG("stage %s start; MAINV %d (HYP %d, forced %d)",STG,length(MAINV),length(HYPV),length(FORCED))

if(STG=="J1"){
  m<-timed("determinism main rerun",safe(MAINV));stopifnot(!is.null(m))
  a<-tab(m,MAINV);b<-read.csv(file.path(MAINRES,"S3_main_OR.csv"),stringsAsFactors=FALSE);b<-b[match(a$cov,b$cov),]
  hp<-m$summary.hyperpar;hb<-read.csv(file.path(MAINRES,"S3_main_hyper.csv"),stringsAsFactors=FALSE)
  so<-spatialOR(m);sb<-read.csv(file.path(MAINRES,"S3_spatial_OR.csv"),stringsAsFactors=FALSE);sb<-sb[match(so$area,sb$area),]
  fs<-fstat(m);fb<-read.csv(file.path(MAINRES,"S3_main_fit.csv"))
  det<-data.frame(quantity=c("fixed beta (18)","fixed beta_sd","OR","OR lo","OR hi","P_gt1","hyperpar mean","hyperpar 2.5%","hyperpar 97.5%","spatial OR (229)","spatial P_gt1","hot/cold class differ","WAIC","DIC","LCPO","phi","sd_area","sd_cell"),
    max_abs_diff=c(max(abs(a$beta-b$beta)),max(abs(a$beta_sd-b$beta_sd)),max(abs(a$OR-b$OR)),max(abs(a$lo-b$lo)),max(abs(a$hi-b$hi)),max(abs(a$P_gt1-b$P_gt1)),
      max(abs(round(hp$mean,4)-hb$mean)),max(abs(round(hp$`0.025quant`,4)-hb$X0.025quant)),max(abs(round(hp$`0.975quant`,4)-hb$X0.975quant)),
      max(abs(so$OR-sb$OR)),max(abs(so$P_gt1-sb$P_gt1)),sum(so$class!=sb$class),abs(fs$WAIC-fb$WAIC),abs(fs$DIC-fb$DIC),abs(fs$LCPO-fb$LCPO),abs(fs$phi-fb$phi),abs(fs$sd_area-fb$sd_area),abs(fs$sd_cell-fb$sd_cell)))
  det$identical<-det$max_abs_diff==0;det$seed_note<-"INLA num.threads 1:1, 같은 입력·같은 호출(S3와 동일 순서); 재적합 1회 더(results_extra orig_rerun)와 합쳐 3회 일치";write.csv(det,file.path(E2,"determinism_determinism_S3_rerun.csv"),row.names=FALSE)
  LOG("determinism max diff overall %.3g",max(det$max_abs_diff))
  eta<-m$summary.fitted.values$mean[seq_len(nrow(DD))];cm<-m$summary.random$id$mean[DD$id]
  write.csv(data.frame(year=DD$year,fbd_region=DD$fbd_region,sigungu_cd=DD$sigungu_cd,y_logit=DD$y_logit,fit_logit_nocell=eta-cm,fit_logit_withcell=eta,resid_logit=DD$y_logit-(eta-cm),SR=DD$SR,fit_SR=round(100*plogis(eta-cm),3)),file.path(E2,"obsfit_obs_fit_resid_cells.csv"),row.names=FALSE)
  rm(m);gc()
  pp<-DD$y_eff/DD$n_eff;stopifnot(all(is.finite(pp)),all(pp>0&pp<1),all(DD$n_eff>=10),all(DD$n_eff==round(DD$n_eff)),all(DD$y_eff==round(DD$y_eff)))
  b12<-read.csv(file.path(MAINRES,"S12_likelihood.csv"),stringsAsFactors=FALSE);b12<-b12[b12$likelihood=="binom",]
  mb<-timed("binomchk binom with assert",safe(MAINV,lik="binom",cpo=FALSE));ab<-tab(mb,MAINV,list(likelihood="binom"));write.csv(ab,file.path(FD,"fit_binomchk_binom_assert_OR.csv"),row.names=FALSE)
  write.csv(data.frame(check=c("pp min","pp max","n_eff min","assert passed","max |beta - S12 binom beta|","max |OR - S12|","cred differ"),value=c(min(pp),max(pp),min(DD$n_eff),TRUE,max(abs(ab$beta-b12$beta[match(ab$cov,b12$cov)])),max(abs(ab$OR-b12$OR[match(ab$cov,b12$cov)])),sum(ab$cred!=b12$cred[match(ab$cov,b12$cov)]))),file.path(E2,"binomchk_binom_assert.csv"),row.names=FALSE);rm(mb)
  cfgfit<-function(vars){rhs<-c("1",bt(vars),SPf(),sprintf("f(time,model='rw1',scale.model=TRUE,hyper=PRI_IID,values=1:%d)",max(DD$time)),IT,CELL)
    inla(as.formula(paste("y_logit~",paste(rhs,collapse="+"))),family="gaussian",data=DD,scale=DD$w,control.family=FIXP,control.fixed=PRI_FIX,control.compute=list(config=TRUE,waic=TRUE,dic=TRUE),control.inla=list(control.vb=list(emergency=50)),verbose=FALSE)}
  out<-NULL;PPC<-NULL
  for(mn in c("main","forced_only")){vv<-if(mn=="main")MAINV else FORCED;mc<-timed(sprintf("intervals %s config",mn),cfgfit(vv));hp<-mc$summary.hyperpar
    set.seed(20261125);ps<-inla.posterior.sample(1000,mc);cn<-rownames(ps[[1]]$latent);ia<-which(cn%in%sprintf("area:%d",1:n));ip<-grep("^Predictor:",cn);ii<-grep("^id:",cn)
    sdA<-sapply(ps,function(s)sd(s$latent[ia,1]));pm<-mc$summary.random$area$mean[1:n]
    out<-rbind(out,data.frame(model=mn,phi_mean=hp["Phi for area","mean"],phi_lo=hp["Phi for area","0.025quant"],phi_hi=hp["Phi for area","0.975quant"],
      bym2_sd_mean=1/sqrt(hp["Precision for area","mean"]),bym2_sd_lo=1/sqrt(hp["Precision for area","0.975quant"]),bym2_sd_hi=1/sqrt(hp["Precision for area","0.025quant"]),
      sd_of_posterior_mean_logOR=sd(pm),sd_district_effects_post_median=median(sdA),sd_district_effects_lo=quantile(sdA,.025),sd_district_effects_hi=quantile(sdA,.975),n_samples=1000,seed=20261125))
    if(mn=="main"){
      ih<-sapply(ps,function(s)s$hyperpar["Precision for id"]);set.seed(20261126)
      Tst<-function(y){da<-tapply(y,DD$area,mean);c(moran_pooled=unname(moran.test(as.numeric(da),lw,zero.policy=TRUE)$estimate[1]),sd_district_mean=sd(da),sd_all=sd(y),p90_p10_SR=unname(quantile(plogis(da),.9)/quantile(plogis(da),.1)),sd_year_mean=sd(tapply(y,DD$year,mean)))}
      T0<-Tst(DD$y_logit);TR<-t(sapply(seq_along(ps),function(k){s<-ps[[k]];mu<-s$latent[ip,1]-s$latent[ii,1][DD$id];yr<-mu+rnorm(nrow(DD),0,1/sqrt(ih[k]))+rnorm(nrow(DD),0,sqrt(DD$v_logit));Tst(yr)}))
      PPC<-data.frame(statistic=names(T0),observed=round(T0,4),rep_mean=round(colMeans(TR),4),rep_lo=round(apply(TR,2,quantile,.025),4),rep_hi=round(apply(TR,2,quantile,.975),4),ppp=round(colMeans(TR>=matrix(T0,nrow(TR),length(T0),byrow=TRUE)),3),n_rep=nrow(TR))
      write.csv(PPC,file.path(E2,"ppc_ppc_main.csv"),row.names=FALSE)}
    rm(mc,ps);gc()}
  write.csv(out,file.path(E2,"intervals_intervals_phi_spatialSD.csv"),row.names=FALSE)
  m0<-RF("nocov_nocov",character(0),keepm=TRUE,note="no covariates, same random effects");r0<-spatialOR(m0);rm(m0)
  rM<-read.csv(file.path(MAINRES,"S3_spatial_OR.csv"));rF<-read.csv(file.path(MAINRES,"S0_forced_only_spatialOR.csv"));raw<-tapply(DD$y_logit,DD$area,mean)
  write.csv(data.frame(pair=c("nocov smoothed vs main","nocov smoothed vs forced-only","nocov smoothed vs raw district mean logit","main vs raw district mean logit","forced-only vs main (manuscript .78)"),
    spearman=round(c(cor(log(r0$OR),log(rM$OR[match(r0$area,rM$area)]),method="spearman"),cor(log(r0$OR),log(rF$OR[match(r0$area,rF$area)]),method="spearman"),cor(log(r0$OR),raw[as.character(r0$area)],method="spearman"),cor(log(rM$OR),raw[as.character(rM$area)],method="spearman"),cor(log(rF$OR),log(rM$OR[match(rF$area,rM$area)]),method="spearman")),3),
    n_hot_nocov=sum(r0$class=="hotspot"),n_cold_nocov=sum(r0$class=="coldspot")),file.path(E2,"nocov_unadjusted_smoothed_correlation.csv"),row.names=FALSE)
  pooled<-as.numeric(tapply(DD$y_logit,DD$area,mean));set.seed(20261110);mc1<-moran.mc(pooled,lw,999,zero.policy=TRUE)
  tn<-moran.test(pooled,lw,randomisation=FALSE,zero.policy=TRUE);tr<-moran.test(pooled,lw,randomisation=TRUE,zero.policy=TRUE)
  mu_t<-tapply(DD$y_logit*DD$w,DD$year,sum)/tapply(DD$w,DD$year,sum);set.seed(20261127)
  Isim<-replicate(999,{ys<-mu_t[as.character(DD$year)]+rnorm(nrow(DD),0,sqrt(DD$v_logit));as.numeric(moran.test(as.numeric(tapply(ys,DD$area,mean)),lw,zero.policy=TRUE)$estimate[1])})
  yr<-sort(unique(DD$year));yI<-t(sapply(yr,function(y){d<-DD[DD$year==y,];x<-d$y_logit[order(d$area)];a_<-moran.test(x,lw,randomisation=FALSE,zero.policy=TRUE);set.seed(y);s_<-replicate(499,{xs<-mu_t[as.character(y)]+rnorm(n,0,sqrt(d$v_logit[order(d$area)]));as.numeric(moran.test(xs,lw,zero.policy=TRUE)$estimate[1])})
    c(year=y,I=unname(a_$estimate[1]),p_normal=a_$p.value,null_size_mean=mean(s_),null_size_p975=unname(quantile(s_,.975)),p_size=(1+sum(s_>=a_$estimate[1]))/500)}))
  write.csv(rbind(data.frame(target="pooled logit SR (2011-24 mean)",I=unname(mc1$statistic),p_permutation999=mc1$p.value,p_normal=tn$p.value,p_randomisation=tr$p.value,z_normal=unname(tn$statistic),null_size_mean=mean(Isim),null_size_p975=unname(quantile(Isim,.975)),p_size_preserving=(1+sum(Isim>=mc1$statistic))/1000,n_sim=999,seed="20261110/20261127")),file.path(E2,"moran2_moran_two_nulls.csv"),row.names=FALSE)
  write.csv(as.data.frame(yI),file.path(E2,"moran2_moran_by_year_two_nulls.csv"),row.names=FALSE)
  SRm<-sapply(yr,function(y)DD$y_logit[DD$year==y][order(DD$area[DD$year==y])]);Vm<-sapply(yr,function(y)DD$v_logit[DD$year==y][order(DD$area[DD$year==y])])
  above<-rowSums(SRm>matrix(mu_t,n,length(yr),byrow=TRUE));Wf<-function(M){R<-apply(M,2,rank);Rs<-rowSums(R);12*sum((Rs-mean(Rs))^2)/(ncol(R)^2*(nrow(R)^3-nrow(R)))}
  set.seed(20261128);NS<-2000;AB<-matrix(0,n,NS);WS<-numeric(NS)
  for(kk in 1:NS){Ms<-matrix(mu_t,n,length(yr),byrow=TRUE)+matrix(rnorm(n*length(yr)),n)*sqrt(Vm);AB[,kk]<-rowSums(Ms>matrix(mu_t,n,length(yr),byrow=TRUE));WS[kk]<-Wf(Ms)}
  p_hi<-(1+rowSums(AB>=above))/(NS+1);p_lo<-(1+rowSums(AB<=above))/(NS+1);ph<-pmin(1,2*pmin(p_hi,p_lo))
  dist<-data.frame(area=1:n,fbd_region=G$region,years_above_national=above,null_mean=round(rowMeans(AB),2),p_two_sided=round(ph,4),q_BH=round(p.adjust(ph,"BH"),4))
  dist$persistent<-ifelse(dist$q_BH<0.05&dist$years_above_national>dist$null_mean,"persistently above",ifelse(dist$q_BH<0.05,"persistently below","not persistent"))
  write.csv(dist,file.path(E2,"nullhot_persistence_size_null_district.csv"),row.names=FALSE)
  write.csv(data.frame(kendall_W_obs=round(Wf(SRm),4),null_size_W_mean=round(mean(WS),4),null_size_W_p975=round(quantile(WS,.975),4),p=(1+sum(WS>=Wf(SRm)))/(NS+1),n_persist_above=sum(dist$persistent=="persistently above"),n_persist_below=sum(dist$persistent=="persistently below"),n_sim=NS,seed=20261128),file.path(E2,"nullhot_persistence_size_null_summary.csv"),row.names=FALSE)
  isl_link<-function(nbx){isl<-which(card(nbx)==0);if(length(isl)){cen<-st_coordinates(st_centroid(st_geometry(G)));kn<-knearneigh(cen,k=1)$nn;for(i in isl){j<-kn[i,1];nbx[[i]]<-sort(unique(c(nbx[[i]][nbx[[i]]>0],j)));nbx[[j]]<-sort(unique(c(nbx[[j]][nbx[[j]]>0],i)))}};nbx}
  ed<-function(nbx){e<-do.call(rbind,lapply(seq_along(nbx),function(i){j<-nbx[[i]][nbx[[i]]>0];if(length(j))cbind(i,j) else NULL}));e<-e[e[,1]<e[,2],,drop=FALSE];paste(e[,1],e[,2])}
  nbr<-isl_link(poly2nb(G,queen=FALSE));nbr0<-poly2nb(G,queen=FALSE);nbq0<-poly2nb(G,queen=TRUE)
  E_g593<-ed(nb);EQ<-ed(nbq);ER<-ed(nbr);EQ0<-ed(nbq0);ER0<-ed(nbr0)
  sn<-c(0.1,1,10,100,500);snE<-lapply(sn,function(s)ed(isl_link(poly2nb(G,queen=FALSE,snap=s))))
  allE<-unique(c(E_g593,EQ,ER));ij<-do.call(rbind,strsplit(allE," "));ij<-matrix(as.integer(ij),ncol=2)
  cen<-st_coordinates(st_centroid(st_geometry(G)))
  seglen<-function(x){if(length(x)==0)return(0);ty<-as.character(st_geometry_type(x));if(all(ty%in%c("POINT","MULTIPOINT")))return(0)
    if(any(ty=="GEOMETRYCOLLECTION")){x<-tryCatch(st_collection_extract(x,"LINESTRING"),error=function(e)NULL);if(is.null(x)||length(x)==0)return(0)}
    x<-x[as.character(st_geometry_type(x))%in%c("LINESTRING","MULTILINESTRING")];if(length(x)==0)0 else as.numeric(sum(st_length(x)))}
  shared<-sapply(seq_len(nrow(ij)),function(k)seglen(suppressWarnings(st_intersection(st_geometry(G)[ij[k,1]],st_geometry(G)[ij[k,2]]))))
  EDF<-data.frame(i=ij[,1],j=ij[,2],region_i=G$region[ij[,1]],region_j=G$region[ij[,2]],in_593=allE%in%E_g593,in_queen_isl=allE%in%EQ,in_rook_isl=allE%in%ER,in_queen_raw=allE%in%EQ0,in_rook_raw=allE%in%ER0,
    shared_boundary_m=round(shared,1),centroid_km=round(sqrt((cen[ij[,1],1]-cen[ij[,2],1])^2+(cen[ij[,1],2]-cen[ij[,2],2])^2)/1000,1))
  for(kk in seq_along(sn))EDF[[sprintf("in_rook_snap%g",sn[kk])]]<-allE%in%snE[[kk]]
  EDF$diff_type<-with(EDF,ifelse(in_593&!in_queen_isl,ifelse(in_queen_raw,"593 only (queen raw has, island rule differs)","593 only (not touching in shapefile)"),ifelse(!in_593&in_queen_isl,ifelse(!in_queen_raw,"queen+island only (island nearest link)",ifelse(in_rook_raw,"queen only, shared edge in rook","queen only (point contact, not rook)")),"same")))
  write.csv(EDF[EDF$diff_type!="same",],file.path(E2,"graphcmp_graph_edge_differences.csv"),row.names=FALSE)
  write.csv(data.frame(graph=c("593 (analysis)","queen+island (S13)","rook+island",sprintf("rook+island snap %g",sn),"queen raw","rook raw"),n_edges=c(length(E_g593),length(EQ),length(ER),sapply(snE,length),length(EQ0),length(ER0)),
    common_with_593=c(length(E_g593),length(intersect(EQ,E_g593)),length(intersect(ER,E_g593)),sapply(snE,function(e)length(intersect(e,E_g593))),length(intersect(EQ0,E_g593)),length(intersect(ER0,E_g593)))),file.path(E2,"graphcmp_graph_comparison.csv"),row.names=FALSE)
  LOG("J1 done")}

if(STG=="J2"){
  source(file.path(ROOT,"01_스크립트","kh_constraints.R"))
  for(ty in c("II","III","IV")){z<-kh_matrices(ty,g,14);d<-DD;d$st<-(d$time-1L)*n+d$area
    assign("KHQ",INLA::inla.as.sparse(z$Q),envir=.GlobalEnv);assign("KHC",list(A=z$A,e=rep(0,z$rankdef)),envir=.GlobalEnv)
    ex<-sprintf("f(st,model='generic0',Cmatrix=KHQ,rankdef=%d,constr=FALSE,extraconstr=KHC,hyper=PRI_IID,values=1:%d)",z$rankdef,n*14);RF(paste0("residspec_KH_",ty),MAINV,data=d,extra=ex)}
  d<-DD;ex<-character(0);for(cvar in HYPV){b<-ave(d[[cvar]],d$area);d[[paste0(cvar,"__B")]]<-b;d[[paste0(cvar,"__W")]]<-d[[cvar]]-b;ex<-c(ex,bt(paste0(cvar,"__B")),bt(paste0(cvar,"__W")))};RF("residspec_S9_between_within",FORCED,data=d,extra=ex)
  RF("residspec_S13_queen_s577",MAINV,sp=SPf("gq"),note="= sets_OBES_s577 (queen 577 edges + island links)")
  assign("pc2",list(prec=list(prior="pc.prec",param=c(1,0.01)),phi=list(prior="pc",param=c(0.5,0.5))),envir=.GlobalEnv);RF("residspec_S15_pc_alt",MAINV,sp=SPf("g","pc2"))
  Xb<-apply(as.matrix(DD[,MAINV]),2,function(x)tapply(x,DD$area,mean))[as.character(1:n),];Xb<-scale(Xb,scale=FALSE);assign("RSRC",list(A=cbind(t(Xb),matrix(0,ncol(Xb),n)),e=rep(0,ncol(Xb))),envir=.GlobalEnv)
  RF("residspec_S16_RSR",MAINV,sp="f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=PRI_BYM2,extraconstr=RSRC)")
  cen<-st_coordinates(st_centroid(st_geometry(G)));for(kk in 3:5){nk<-make.sym.nb(knn2nb(knearneigh(cen,k=kk)));assign("gk",mkG(nk,sprintf("knn%d",kk)),envir=.GlobalEnv);RF(sprintf("residspec_S20_knn%d",kk),MAINV,sp=SPf("gk"))}
  Am<-matrix(0,n,n);for(ii in seq_len(n))Am[ii,nb[[ii]][nb[[ii]]>0]]<-1;Am<-(Am+t(Am))>0;Lp<-diag(rowSums(Am))-Am*1;ev<-eigen(Lp,symmetric=TRUE);oo<-order(ev$values);nz<-sum(ev$values<1e-8);Ev<-ev$vectors[,oo[(nz+1):(nz+23)]]
  d<-DD;EE<-Ev[d$area,];for(cvar in MAINV){r<-resid(lm(d[[cvar]]~EE));d[[cvar]]<-(r-mean(r))/sd(r)};RF("residspec_S24_spatialplus",MAINV,data=d)
  LOG("J2 done")}

if(STG=="J3"){
  RF("access_plus_road",c(MAINV,"road_pavement_rate"))
  d<-DD;for(kk in 2:5)d[[sprintf("urbQ%d",kk)]]<-as.integer(d$urban_q5==kk);V0<-setdiff(MAINV,"urban_population_share")
  RF("access_urban_q5",c(V0,sprintf("urbQ%d",2:5)),data=d);RF("access_urban_q5_plus_road",c(V0,sprintf("urbQ%d",2:5),"road_pavement_rate"),data=d)
  RF("access_urban_rw2",V0,extra="f(urban_rw2_bin,model='rw2',scale.model=TRUE,hyper=PRI_IID)")
  RF("access_plus_outside_care",c(MAINV,"outside_care_share"))
  RF("care_plus_care_institutions",c(MAINV,"care_institutions_per_100k"));RF("care_plus_hospitals",c(MAINV,"hospitals_per_100k"))
  RF("care_plus_outside_and_hospitals",c(MAINV,"outside_care_share","hospitals_per_100k"))
  IMP9<-c("cars_per_capita","clinics_per_100k","eloc__dt_1yl20421e__1___1","fiscal_independence","moderate_pa_rate","public_sports_facilities_per_100k","sex_ratio","urban_population_share","village_sports_facilities_per_100k")
  RF("imputed_drop_all9_imputed_vars",setdiff(MAINV,IMP9));RF("imputed_drop_all11_incl_census_interp",setdiff(MAINV,c(IMP9,"tertiary_educated_share","basic_livelihood_rate")))
  RF("imputed_cc_no_imputed_cell",MAINV,data=reid(DD[DD$imp_any_main==0,]))
  RF("instant_years_2016_2024",MAINV,data=retime(DD[DD$year>=2016,]),note="instant-food licence stock within 1% of MFDS yearbook only from 2016")
  RF("s2_no_pop80_forced",setdiff(MAINV,"pop80_share"),note="sets_OBES_s2 (post hoc label)")
  LOG("J3 done")}

if(STG=="J4"){
  WR<-NULL;WF<-NULL;WM<-NULL
  for(wn in names(WIN)){yy<-WIN[[wn]];d<-retime(DD[DD$year>=yy[1]&DD$year<=yy[2],]);m<-RF(paste0("W85_",wn),MAINV,data=d,keepm=TRUE);rr<-spatialOR(m);rr$window<-wn;WR<-rbind(WR,rr)
    hp<-m$summary.hyperpar;WF<-rbind(WF,data.frame(window=wn,n_obs=nrow(d),WAIC=round(m$waic$waic,1),phi=hp["Phi for area","mean"],phi_lo=hp["Phi for area","0.025quant"],phi_hi=hp["Phi for area","0.975quant"],n_hot=sum(rr$class=="hotspot"),n_cold=sum(rr$class=="coldspot"),OR_min=min(rr$OR),OR_max=max(rr$OR),OR_p10=quantile(rr$OR,.1),OR_p90=quantile(rr$OR,.9)))
    set.seed(20261110);mm<-moran.mc(log(rr$OR),lw,999,zero.policy=TRUE);pw<-as.numeric(tapply(d$y_logit,d$area,mean)[as.character(1:n)]);mo<-moran.mc(pw,lw,999,zero.policy=TRUE)
    WM<-rbind(WM,data.frame(window=wn,quantity=c("bym2_logOR_window","observed_logitSR_window_mean"),moran_I=round(c(mm$statistic,mo$statistic),4),p_perm=c(mm$p.value,mo$p.value)));rm(m);gc()}
  w<-names(WIN);H<-sapply(w,function(k)WR$class[WR$window==k]=="hotspot");C<-sapply(w,function(k)WR$class[WR$window==k]=="coldspot")
  WF$hot_kept_from_prev<-c(NA,sapply(2:3,function(k)sum(H[,k]&H[,k-1])));WF$hot_new<-c(NA,sapply(2:3,function(k)sum(H[,k]&!H[,k-1])));WF$hot_lost<-c(NA,sapply(2:3,function(k)sum(!H[,k]&H[,k-1])))
  WF$cold_kept_from_prev<-c(NA,sapply(2:3,function(k)sum(C[,k]&C[,k-1])));WF$hot_all_windows<-sum(rowSums(H)==3);WF$cold_all_windows<-sum(rowSums(C)==3)
  WF$spearman_logOR_vs_prev<-c(NA,sapply(2:3,function(k)round(cor(log(WR$OR[WR$window==w[k]]),log(WR$OR[WR$window==w[k-1]]),method="spearman"),3)))
  write.csv(WR,file.path(E2,"W85_window_spatialOR.csv"),row.names=FALSE);write.csv(WF,file.path(E2,"W85_window_fit.csv"),row.names=FALSE);write.csv(WM,file.path(E2,"W85_window_moran.csv"),row.names=FALSE)
  RG<-c("서울시"="Capital","인천시"="Capital","경기도"="Capital","대전시"="Chungcheong","세종시"="Chungcheong","충청북도"="Chungcheong","충청남도"="Chungcheong","광주시"="Honam","전라북도"="Honam","전라남도"="Honam","부산시"="Yeongnam","대구시"="Yeongnam","울산시"="Yeongnam","경상북도"="Yeongnam","경상남도"="Yeongnam","강원도"="Gangwon_Jeju","제주도"="Gangwon_Jeju")
  DD$rg<-unname(RG[sidoOf[DD$area]]);stopifnot(!anyNA(DD$rg));m0<-safe(MAINV,cpo=FALSE);w0<-m0$waic$waic;rm(m0);out<-NULL
  for(cvar in c("moderate_pa_rate","eloc__dt_1yl20421e__1___1","cars_per_capita")){D2<-DD;nv<-character(0);for(rgn in unique(DD$rg)){cn<-paste0(cvar,"__",rgn);D2[[cn]]<-D2[[cvar]]*(D2$rg==rgn);nv<-c(nv,cn)}
    set.seed(20261110);m<-safe(c(setdiff(MAINV,cvar),nv),data=D2,cpo=FALSE);if(is.null(m))next;r<-do.call(rbind,lapply(nv,function(x)orr(m,x)));r$WAIC_region<-round(m$waic$waic,1);r$WAIC_common<-round(w0,1);r$covariate<-cvar;r$region<-sub(".*__","",r$cov);r$dWAIC_vs_common<-round(m$waic$waic-w0,1);out<-rbind(out,r);rm(m)}
  write.csv(out,file.path(E2,"slopes_S17_region_slopes_remaining3.csv"),row.names=FALSE)
  for(L_ in 1:2)for(vr in c("A","B")){d<-DD;for(cvar in HYPV)d[[cvar]]<-d[[sprintf("cl%s%d_%s",vr,L_,cvar)]];k<-stats::complete.cases(d[,HYPV]);d<-retime(d[k,])
    RF(sprintf("callag_cal_lag%d_%s",L_,vr),MAINV,data=d,note=sprintf("calendar lag %d, variant %s",L_,vr))
    if(vr=="A")RF(sprintf("callag_nolag_matched_lag%d_A",L_),MAINV,data=retime(DD[k,]),note="same rows, no lag")}
  d<-DD;d$sexr_bin<-as.integer(cut(rank(d$sex_ratio,ties.method="first"),20,labels=FALSE));d$p80_bin<-as.integer(cut(rank(d$pop80_share,ties.method="first"),20,labels=FALSE))
  RF("G88_smooth_sexratio_pop80",setdiff(MAINV,c("sex_ratio","pop80_share")),data=d,extra=c("f(sexr_bin,model='rw2',scale.model=TRUE,hyper=PRI_IID)","f(p80_bin,model='rw2',scale.model=TRUE,hyper=PRI_IID)"))
  d<-DD;KC<-c("walking_practice_rate","moderate_pa_rate","high_risk_drinking_rate_drinkers","current_smoking_rate","perceived_stress_rate");for(cvar in KC)d[[cvar]]<-d[[paste0("ma3_",cvar)]];RF("G88_kchs_cov_ma3",MAINV,data=d)
  d<-DD;for(cvar in HYPV)d[[cvar]]<-d[[paste0("xt_",cvar)]];RF("G88_extreme_5IQR_missing",MAINV,data=d)
  fo<-file.path(FD,"fit_G88_poisson_lognormal_OR.csv");if(!file.exists(fo)){rhs<-c("1",bt(MAINV),SPf(),sprintf("f(time,model='rw1',scale.model=TRUE,hyper=PRI_IID,values=1:%d)",max(DD$time)),IT,CELL)
    set.seed(20261110);FB$last<-"default(direct inla)";m<-timed("G88 poisson-lognormal",inla(as.formula(paste("y_eff~",paste(rhs,collapse="+"))),family="poisson",E=DD$n_eff,data=DD,control.fixed=PRI_FIX,control.compute=list(waic=TRUE,dic=TRUE),control.inla=list(control.vb=list(emergency=50)),verbose=FALSE))
    write.csv(tab(m,MAINV,list(fit="G88_poisson_lognormal",n_obs=nrow(DD))),fo,row.names=FALSE);rr<-spatialOR(m);write.csv(cbind(fit="G88_poisson_lognormal",n_obs=nrow(DD),fstat(m),resid_moran_I=NA,resid_moran_p=NA,n_hot=sum(rr$class=="hotspot"),n_cold=sum(rr$class=="coldspot"),note="y_eff ~ Poisson(n_eff·θ), log link: effects are prevalence ratios"),file.path(FD,"fit_G88_poisson_lognormal_stat.csv"),row.names=FALSE);rm(m)}
  LOG("J4 done")}

if(STG%in%c("J5","J6")){
  ks<-if(STG=="J5")1:10 else 11:20;mu_t<-tapply(DD$y_logit*DD$w,DD$year,sum)/tapply(DD$w,DD$year,sum);out<-NULL;fo<-file.path(E2,sprintf("nullhot_null_hotspots_%s.csv",STG))
  for(kk in ks){set.seed(20261200+kk);d<-DD;d$y_logit<-mu_t[as.character(d$year)]+rnorm(nrow(d),0,sqrt(d$v_logit))
    m<-timed(sprintf("nullhot null sim %d",kk),safe(MAINV,data=d,cpo=FALSE));if(is.null(m))next;rr<-spatialOR(m);fs<-fstat(m)
    out<-rbind(out,data.frame(sim=kk,seed=20261200+kk,n_hot=sum(rr$class=="hotspot"),n_cold=sum(rr$class=="coldspot"),largest_hot=clusters(rr$area[rr$class=="hotspot"]),largest_cold=clusters(rr$area[rr$class=="coldspot"]),phi=fs$phi,sd_area=fs$sd_area,sd_logOR=sd(log(rr$OR)),OR_min=min(rr$OR),OR_max=max(rr$OR),fit_setting=FB$last));write.csv(out,fo,row.names=FALSE);rm(m);gc()}
  LOG("%s done",STG)}

if(STG=="J7"){
  out<-NULL;for(kk in c(2013:2018,2020:2023)){d<-DD;d$postk<-as.integer(d$year>=kk);ex<-"postk";for(cvar in HYPV){nm<-sprintf("%s_x_postk",cvar);d[[nm]]<-d[[cvar]]*d$postk;ex<-c(ex,bt(nm))}
    set.seed(20261110);m<-timed(sprintf("splitperm split %d",kk),safe(MAINV,data=d,extra=ex,cpo=FALSE));if(is.null(m))next
    for(cvar in HYPV){r<-orr(m,sprintf("%s_x_postk",cvar));r$cov<-cvar;r$split_year<-kk;r$z<-r$beta/r$beta_sd;out<-rbind(out,r)};write.csv(out,file.path(E2,"splitperm_inla_split_interactions.csv"),row.names=FALSE);rm(m);gc()}
  LOG("J7 done")}
if(STG=="J8"){
  AX<-read.csv(file.path(E2,"inputs_extra2_a205.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE,check.names=FALSE);AX$fbd_region<-enc2utf8(AX$fbd_region);ia<-match(paste(DD$year,DD$fbd_region),paste(AX$year,AX$fbd_region));stopifnot(!anyNA(ia));AX<-AX[ia,]
  for(tg in c("a205k25","a205k20")){d<-DD;for(cvar in MAINV)d[[cvar]]<-AX[[paste0(tg,"__",cvar)]];stopifnot(!anyNA(d[,MAINV]));RF(paste0("gatethr_",tg),MAINV,data=d,note=sprintf("scan alternative %s cells set missing and interpolated",tg))}
  LOG("J8 done")}
writeLines(capture.output(sessionInfo()),file.path(E2,"code",sprintf("sessionInfo_%s.txt",STG)))
LOG("stage %s ALL DONE",STG)
