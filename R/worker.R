# Portable adapter around the byte-preserved CRNP v0.2.0 functions.
if(.Platform$OS.type=="windows") {
  utf8_locale<-suppressWarnings(Sys.setlocale("LC_CTYPE",".UTF-8"))
  if(!nzchar(utf8_locale))stop("A UTF-8 capable Windows R runtime is required; use R 4.4.3")
}
suppressPackageStartupMessages({library(jsonlite);library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
job <- fromJSON(args[[1L]], simplifyVector=FALSE)
for (name in c("core.R","fidelity.R","mapping.R","model.R","evaluation.R","simulate.R")) {
  source(file.path(job$vendor_r,name),local=.GlobalEnv,encoding="UTF-8")
}
if (!requireNamespace("WGCNA",quietly=TRUE)) stop("Install WGCNA; see docs/INSTALL.md")
WGCNA::disableWGCNAThreads()
if (isTRUE(job$strict_environment) &&
    (getRversion()!="4.4.3" || as.character(packageVersion("WGCNA"))!="1.73")) {
  stop("Frozen R/WGCNA environment differs; use the documented versions or explicitly allow drift")
}
env <- list(R=as.character(getRversion()),WGCNA=as.character(packageVersion("WGCNA")),LC_CTYPE=Sys.getlocale("LC_CTYPE"),
            crnp="0.2.0",data.table=as.character(packageVersion("data.table")),jsonlite=as.character(packageVersion("jsonlite")))
dir.create(job$output,recursive=TRUE,showWarnings=FALSE)
write_json(env,file.path(job$output,"R_environment.json"),pretty=TRUE,auto_unbox=TRUE)

read_locked <- function(path) {
  con <- if(grepl("\\.gz$",path,ignore.case=TRUE)) gzfile(path,"rt") else file(path,"rt",encoding="UTF-8")
  on.exit(close(con))
  d <- as.data.table(read.csv(con,check.names=FALSE,stringsAsFactors=FALSE,
                             colClasses=c(participant_id="character",record_id="character")))
  nm <- names(d)[1L];nm <- sub("^\ufeff","",nm);setnames(d,1L,nm)
  stopifnot(identical(names(d)[1:2],c("participant_id","record_id")))
  x <- as.matrix(d[,-c("participant_id","record_id")]);storage.mode(x)<-"double"
  list(participant=d$participant_id,record=d$record_id,x=x)
}

if (job$task %in% c("recovery","timing")) {
  weights <- fread(job$weights,check.names=FALSE)
  masks <- fread(job$masks,check.names=FALSE)
  units <- fread(job$units,check.names=FALSE)
  seeds <- fread(job$seeds,check.names=FALSE)
  batches <- lapply(job$batches,read_locked)
  make_case <- function(i) {
    u<-units[i];wd<-weights[module==u$module]
    w<-setNames(as.numeric(wd$weight),as.character(wd$protein))
    observed<-names(w)[names(w)%in%masks[mask_id==u$mask_unit_id,protein]]
    source<-batches[[u$train_batch]]$x[,names(w),drop=FALSE]
    source<-source[complete.cases(source),,drop=FALSE]
    if(nrow(source)<50L) stop("Fewer than 50 complete source rows")
    target<-batches[[u$test_batch]]$x[,observed,drop=FALSE]
    stream_name<-paste0("fold::",u$direction_unit_id)
    seed<-seeds[stream==stream_name & replicate==0L,seed]
    if(length(seed)!=1L)stop("Non-unique frozen seed")
    list(unit=u,source=source,target=target,w=w,observed=observed,seed=as.integer(seed),
         ids=batches[[u$test_batch]])
  }
  fit_case <- function(ca) {
    model<-fit_crnp(ca$source,ca$w,ca$observed,
                    lambdas=c(1e-4,1e-3,1e-2,.1,1,10),calibration_replicates=2L,folds=3L,seed=ca$seed)
    full<-as.numeric(ca$source%*%ca$w)/sum(abs(ca$w))
    list(model=model,mu=mean(full),sigma=sd(full))
  }
  score_case <- function(ca,fit) {
    z<-project_crnp(fit$model,ca$target,center=rep(0,length(ca$observed)),
                    scale=rep(1,length(ca$observed)),standardize_score=FALSE)
    fit$mu+fit$sigma*as.numeric(z)
  }
}

if(job$task=="recovery") {
  scores<-list();audit<-list()
  for(i in seq_len(nrow(units))) {
    ca<-make_case(i);fit<-fit_case(ca);pred<-score_case(ca,fit)
    if(any(!is.finite(pred)))stop("Nonfinite CRNP predictions")
    scores[[i]]<-data.table(direction_unit_id=ca$unit$direction_unit_id,
       participant_id=ca$ids$participant,record_id=ca$ids$record,`CRNP_v0.2.0`=pred)
    audit[[i]]<-data.table(direction_unit_id=ca$unit$direction_unit_id,
       lambda=fit$model$lambda,recovery_fraction=fit$model$recovery_fraction,
       source_score_center=fit$mu,source_score_scale=fit$sigma,output_scale="normalized_frozen_weight_via_source_inverse_z")
  }
  fwrite(rbindlist(scores),file.path(job$output,"crnp_scores.csv"))
  fwrite(rbindlist(audit),file.path(job$output,"crnp_fit_audit.csv"))
} else if(job$task=="timing") {
  if(!requireNamespace("microbenchmark",quietly=TRUE))stop("Install microbenchmark")
  cases<-lapply(seq_len(nrow(units)),make_case)
  nano<-microbenchmark::get_nanotime
  reply<-function(x){cat("RESULT:",toJSON(x,auto_unbox=TRUE,null="null",digits=17),"\n",sep="");flush.console()}
  block<-function(ca,fit,calls){t<-nano();for(j in seq_len(calls))last<-score_case(ca,fit);list(elapsed=(nano()-t)/1e9,last=last)}
  expected<-lapply(cases,function(ca)score_case(ca,fit_case(ca)))
  reply(list(action="ready",cases=length(cases),environment=env))
  input<-file("stdin",open="r")
  repeat {
    line<-readLines(input,n=1L,warn=FALSE);if(!length(line))break
    cmd<-fromJSON(line);if(cmd$action=="stop"){reply(list(action="stopped"));break}
    ca<-cases[[cmd$case_index]];gc();t<-nano();fit<-fit_case(ca);fit_s<-(nano()-t)/1e9
    initial<-score_case(ca,fit)
    difference<-max(abs(initial-expected[[cmd$case_index]]))
    if(!is.finite(difference)||difference>1e-7)stop("Repeated CRNP fitting changed predictions")
    if(cmd$action=="warmup") {
      calls<-16L
      repeat {b<-block(ca,fit,calls);if(b$elapsed>=.025||calls>=131072L)break;calls<-calls*2L}
      calls<-as.integer(min(131072,max(32,2^ceiling(log2(.1/(b$elapsed/calls))))))
      reply(list(action="warmup",prediction_calls=calls,max_abs_score_difference=difference))
    } else {
      gc();b<-block(ca,fit,as.integer(cmd$prediction_calls))
      difference<-max(abs(b$last-expected[[cmd$case_index]]))
      reply(list(action="measure",fit_s=fit_s,prediction_ms=1000*b$elapsed/cmd$prediction_calls,
                 prediction_calls=cmd$prediction_calls,max_abs_score_difference=difference))
    }
  }
  close(input)
} else if(job$task=="emory") {
  weights<-fread(job$weights);mapping<-fread(job$mapping)
  source_dt<-fread(job$source_matrix,check.names=FALSE,colClasses=list(character=1L))
  source_ids<-as.character(source_dt[[1L]]);source_dt[[1L]]<-NULL
  if(anyNA(source_ids)||anyDuplicated(source_ids)||anyDuplicated(names(source_dt)))stop("Source IDs and protein columns must be unique and complete")
  source<-as.matrix(source_dt);storage.mode(source)<-"double";rownames(source)<-source_ids
  if(nrow(source)<50L||any(!is.finite(source)))stop("Source definition matrix needs >=50 finite rows")
  modules<-unique(weights$ModuleColor)
  branch_data<-lapply(job$branches,function(b){d<-fread(b$path,check.names=FALSE,colClasses=list(character=1L));ids<-as.character(d[[1L]]);d[[1L]]<-NULL;if(anyNA(ids)||anyDuplicated(ids)||anyDuplicated(names(d)))stop("Assay IDs and feature columns must be unique and complete");x<-as.matrix(d);storage.mode(x)<-"double";rownames(x)<-ids;list(spec=b,x=x)})
  common<-sort(Reduce(intersect,lapply(branch_data,function(b)rownames(b$x))))
  if(length(common)<5L)stop("Insufficient paired Emory participants")
  if(isTRUE(job$strict_paper_design)&&length(common)!=35L)stop("Paper Emory intersection requires 35 paired participants")
  mats<-list();all_scores<-list();score_pos<-0L;qc<-list()
  for(b in branch_data) {
    x<-b$x[common,,drop=FALSE];sp<-b$spec
    if(sp$input_scale=="raw") {
      if(sp$platform%in%c("SomaScan","TMT-MS")){x[!is.finite(x)|x<=0]<-NA_real_;x<-log2(x)}
    } else if(sp$input_scale!="transformed_features")stop("input_scale must be raw or transformed_features")
    mp<-mapping[input_id==sp$input_id & conservative_eligible==TRUE,.(feature_id,protein_id,mapping_route)]
    ag<-aggregate_features(x,mp,weights$ProteinID,min_finite_fraction=.8,sample_missing="median",aggregate="median")
    pm<-ag$protein_matrix;mats[[sp$input_id]]<-pm
    qc[[sp$input_id]]<-as.data.table(ag$feature_qc)[,input_id:=sp$input_id]
    for(mod in modules) {
      wd<-weights[ModuleColor==mod];w<-setNames(wd$kME_ownModule,wd$ProteinID)
      obs<-names(w)[names(w)%in%colnames(pm)]
      if(length(obs)<2L)next
      if(length(obs)==length(w))stop("Frozen CRNP requires partial coverage; complete-coverage branches need a separately specified analysis")
      model<-fit_crnp(source[,names(w),drop=FALSE],w,obs,seed=93017L)
      sc<-as.data.table(compare_methods(model,pm[,obs,drop=FALSE],methods=c("direct_kME","CRNP"),standardize_scores=TRUE))
      sc[,`:=`(input_id=sp$input_id,platform=sp$platform,biofluid=sp$biofluid,module=mod)]
      score_pos<-score_pos+1L;all_scores[[score_pos]]<-sc
    }
  }
  scores<-rbindlist(all_scores);fwrite(scores,file.path(job$output,"restricted_emory_scores.csv"))
  fwrite(rbindlist(qc,fill=TRUE),file.path(job$output,"feature_qc.csv"))
  wide<-dcast(scores,sample_id+biofluid+module~platform+method,value.var="score")
  pairs<-list("Olink-TMT-MS"=c("Olink","TMT-MS"),"SomaScan-Olink"=c("SomaScan","Olink"),"SomaScan-TMT-MS"=c("SomaScan","TMT-MS"))
  keys<-CJ(biofluid=c("CSF","plasma"),module=modules,platform_pair=names(pairs),sorted=TRUE)
  safe<-function(a,b){k<-is.finite(a)&is.finite(b);if(sum(k)<5L||sd(a[k])==0||sd(b[k])==0)return(NA_real_);suppressWarnings(cor(a[k],b[k],method="spearman"))}
  fisher<-function(x)atanh(pmax(-.999999,pmin(.999999,x)))
  cell_delta<-function(k,idx) {
    ps<-pairs[[k$platform_pair]];d<-wide[biofluid==k$biofluid & module==k$module][order(sample_id)]
    if(nrow(d)!=length(common))return(NA_real_)
    cols<-c(paste0(ps,"_CRNP"),paste0(ps,"_direct_kME"))
    if(!all(cols%in%names(d)))return(NA_real_)
    fisher(safe(d[[cols[1L]]][idx],d[[cols[2L]]][idx]))-fisher(safe(d[[cols[3L]]][idx],d[[cols[4L]]][idx]))
  }
  dv<-vapply(seq_len(nrow(keys)),function(i)cell_delta(keys[i],seq_along(common)),numeric(1))
  cells<-copy(keys);cells[,delta_z:=dv];fwrite(cells,file.path(job$output,"cell_differences.csv"))
  if(!any(is.finite(dv)))stop("No evaluable concordance cells")
  if(isTRUE(job$strict_paper_design)&&(sum(is.finite(dv)&keys$biofluid=="CSF")!=27L||sum(is.finite(dv)&keys$biofluid=="plasma")!=25L))stop("Paper Emory analysis requires 27 CSF and 25 plasma evaluable cells")
  summarize<-function(x)c(all=mean(x,na.rm=TRUE),CSF=mean(x[keys$biofluid=="CSF"],na.rm=TRUE),plasma=mean(x[keys$biofluid=="plasma"],na.rm=TRUE))
  point<-summarize(dv);B<-as.integer(job$n_bootstrap);set.seed(2026081601L)
  boot<-matrix(NA_real_,B,3L)
  for(bb in seq_len(B)) {
    idx<-sample.int(length(common),length(common),replace=TRUE)
    dd<-vapply(seq_len(nrow(keys)),function(i)cell_delta(keys[i],idx),numeric(1))
    boot[bb,]<-summarize(dd)
  }
  ans<-data.table(level=names(point),estimate=as.numeric(point),ci_low=NA_real_,ci_high=NA_real_,cells=c(sum(is.finite(dv)),sum(is.finite(dv)&keys$biofluid=="CSF"),sum(is.finite(dv)&keys$biofluid=="plasma")),participants=length(common))
  for(i in seq_len(3L)) {
    v<-boot[,i];v<-v[is.finite(v)]
    if(length(v)>=.95*B){ans$ci_low[i]<-quantile(v,.025,names=FALSE);ans$ci_high[i]<-quantile(v,.975,names=FALSE)}
  }
  fwrite(ans,file.path(job$output,"emory_concordance.csv"))
  fwrite(as.data.table(boot),file.path(job$output,"bootstrap_aggregates.csv"))
} else if(job$task=="contracts") {
  set.seed(37L);x<-scale(matrix(rnorm(160*10),160,10));colnames(x)<-paste0("SYN_P",seq_len(10))
  w<-setNames(seq(.3,1.2,length.out=10),colnames(x));obs<-names(w)[1:6]
  m<-fit_crnp(x,w,obs,seed=91L)
  y<-x[1:40,,drop=FALSE];p<-project_crnp(m,y,center=rep(0,6),scale=rep(1,6),standardize_score=FALSE)
  manual<-m$calibration_intercept+m$calibration_slope*as.numeric(y[,obs,drop=FALSE]%*%m$effective_weights)
  stopifnot(max(abs(p-manual))<1e-12)
  altered<-y;altered[,setdiff(colnames(y),obs)]<-1e6
  p2<-project_crnp(m,altered,center=rep(0,6),scale=rep(1,6),standardize_score=FALSE)
  stopifnot(identical(p,p2))
  permuted<-y[,rev(colnames(y)),drop=FALSE]
  p3<-project_crnp(m,permuted,center=rep(0,6),scale=rep(1,6),standardize_score=FALSE)
  stopifnot(identical(p,p3))
  expected_fraction<-max(0,min(1,.6*(.6-.2)))
  stopifnot(abs(m$recovery_fraction-expected_fraction)<1e-12)
  mu<-mean(as.numeric(x%*%w)/sum(abs(w)));sigma<-sd(as.numeric(x%*%w)/sum(abs(w)))
  stopifnot(max(abs((mu+sigma*p-mu)/sigma-p))<1e-12)
  checks<-list(status="PASS",data_kind="synthetic",checks=c("projection_linear_identity","unobserved_target_features_do_not_change_predictions","column_name_order_invariance","coverage_attenuation","source_inverse_z"),installed_crnp_parity="not_requested")
  if(requireNamespace("crnp",quietly=TRUE)&&as.character(packageVersion("crnp"))=="0.2.0") {
    installed<-crnp::fit_crnp(x,w,obs,seed=91L)
    ip<-crnp::project_crnp(installed,y,center=rep(0,6),scale=rep(1,6),standardize_score=FALSE)
    stopifnot(max(abs(p-ip))<1e-12);checks$installed_crnp_parity<-"PASS"
  }
  write_json(checks,file.path(job$output,"R_CONTRACTS.json"),auto_unbox=TRUE,pretty=TRUE)
} else stop("Unknown R task")
