#if (exists(".DEF.COMMON")) stop ("_common.R has already been loaded") else .DEF.COMMON=TRUE
library(methods)
library(dplyr)
library(marginalizedRisk)
library(survival)
library(parallel)
library(kyotil)

if(Sys.getenv("TRIAL")=="") stop("Environmental variable TRIAL not defined!!!!!!!!!!!!!!")
TRIAL=Sys.getenv("TRIAL")

set.seed(98109)

if(!exists("verbose")) verbose=0
if (Sys.getenv("VERBOSE") %in% c("T","TRUE")) verbose=1
if (Sys.getenv("VERBOSE") %in% c("1", "2", "3")) verbose=as.integer(Sys.getenv("VERBOSE"))



###################################################################################################
# disable lower level parallelization in favor of higher level of parallelization

library(RhpcBLASctl)
blas_get_num_procs()
blas_set_num_threads(1L)
#stopifnot(blas_get_num_procs() == 1L) # Commented this out as it does not work as expected any more!
omp_set_num_threads(1L)
stopifnot(omp_get_max_threads() == 1L)
    
# for mclapply etc
numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows", 1, as.integer(future::availableCores()/2)))


###################################################################################################
# load config, assay metadata, define common variables and labels

config <- config::get(config = Sys.getenv("TRIAL"))
study_name=config$study_name


# created named lists for assay metadata to easier access, e.g. assay_labels_short["bindSpike"]
assay_metadata = read.csv(config$assay_metadata)
assays=assay_metadata$assay
assay_labels=assay_metadata$assay_label; names(assay_labels)=assays
assay_labels_short=assay_metadata$assay_label_short; names(assay_labels_short)=assays
llox_labels=assay_metadata$llox_label; names(llox_labels)=assays
lloqs=assay_metadata$lloq; names(lloqs)=assays
lods=assay_metadata$lod; names(lods)=assays
lloxs=ifelse(llox_labels=="lloq", lloqs, lods)


form.s = Surv(EventTimePrimary, EventIndPrimary) ~ 1
form.0 = update (form.s, as.formula(config$covariates))
print(form.0)

# read analysis ready data
dat = read.csv(config$data_cleaned)


###############################################################################
# compute tfinal.tpeak

# returns smaller of the two: 1) time of the last case, 2) last time to have 15 at risk
# use case: 1) dat.ph2 is the ph2 dataset from a case control study
get.tfinal.tpeak.1 = function(dat.ph2, event.ind.col="EventIndPrimary", event.time.col="EventTimePrimary") {
  min(
    max (dat.ph2[[event.time.col]] [dat.ph2[[event.ind.col]]==1]),
    sort(dat.ph2[[event.time.col]], decreasing=T)[15]-1
  )
}

if (TRIAL=="moderna_boost") {
  # compute tfinal.tpeak as the minimum of the four quadrants and no larger than 105 days
  tfinal.tpeaks=c(
    get.tfinal.tpeak.1(subset(dat, Trt==1 &  naive & ph2.BD29), event.ind.col="EventIndOmicronBD29", event.time.col="EventTimeOmicronBD29"),
    get.tfinal.tpeak.1(subset(dat, Trt==1 &  naive & ph2.BD29), event.ind.col="EventIndOmicronBD29", event.time.col="EventTimeOmicronBD29"),
    get.tfinal.tpeak.1(subset(dat, Trt==1 & !naive & ph2.BD29), event.ind.col="EventIndOmicronBD29", event.time.col="EventTimeOmicronBD29"),
    get.tfinal.tpeak.1(subset(dat, Trt==1 & !naive & ph2.BD29), event.ind.col="EventIndOmicronBD29", event.time.col="EventTimeOmicronBD29"),
    105)
  myprint(tfinal.tpeaks)
  tfinal.tpeak=min(tfinal.tpeaks)
  myprint(tfinal.tpeak)
}




###################################################################################################
# shared functions: survival analysis

# get marginalized risk to the followup followup.day without marker
get.marginalized.risk.no.marker=function(formula, dat.ph1, followup.day){
  if (!is.list(formula)) {
    # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
    fit.risk = coxph(formula, dat.ph1, model=T) 
    dat.ph1$EventTimePrimary=followup.day
    risks = 1 - exp(-predict(fit.risk, newdata=dat.ph1, type="expected"))
    mean(risks)
  } else {
    # competing risk estimation
    out=pcr2(formula, dat.ph1, followup.day)
    mean(out)
  }
}




###################################################################################################
# shared functions: misc

add.trichotomized.markers=function(dat, markers, ph2.col.name="ph2", wt.col.name="wt") {
  
  if(verbose) print("add.trichotomized.markers ...")
  
  marker.cutpoints <- list()    
  for (a in markers) {
    if (verbose) myprint(a, newline=F)
    tmp.a=dat[[a]]
    
    # if we estimate cutpoints using all non-NA markers, it may have an issue when a lot of subjects outside ph2 have non-NA markers
    flag=dat[[ph2.col.name]]
    
    if(startsWith(a, "Delta")) {
      # fold change
      q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/3, 2/3))
    } else {
      # not fold change
      uloq=assay_metadata$uloq[assay_metadata$assay==marker.name.to.assay(a)]
      
      uppercut=log10(uloq); uppercut=uppercut*ifelse(uppercut>0,.9999,1.0001)
      lowercut=min(tmp.a, na.rm=T)*1.0001; lowercut=lowercut*ifelse(lowercut>0,1.0001,.9999)
      if (mean(tmp.a>uppercut, na.rm=T)>1/3) {
        # if more than 1/3 of vaccine recipients have value > ULOQ, let q.a be (median among those < ULOQ, ULOQ)
        if (verbose) cat("more than 1/3 of vaccine recipients have value > ULOQ\n")
        q.a=c(wtd.quantile(tmp.a[dat[[a]]<=uppercut & flag], weights = dat[[wt.col.name]][tmp.a<=uppercut & flag], probs = c(1/2)),  uppercut)
      } else if (mean(tmp.a<lowercut, na.rm=T)>1/3) {
        # if more than 1/3 of vaccine recipients have value at min, let q.a be (min, median among those > LLOQ)
        if (verbose) cat("more than 1/3 of vaccine recipients have at min\n")
        q.a=c(lowercut, wtd.quantile(tmp.a[dat[[a]]>=lowercut & flag], weights = dat[[wt.col.name]][tmp.a>=lowercut & flag], probs = c(1/2))  )
      } else {
        # this implementation uses all non-NA markers, which include a lot of subjects outside ph2, and that leads to uneven distribution of markers between low/med/high among ph2
        #q.a <- wtd.quantile(tmp.a, weights = dat[[wt.col.name]], probs = c(1/3, 2/3))
        q.a <- wtd.quantile(tmp.a[flag], weights = dat[[wt.col.name]][flag], probs = c(1/3, 2/3))
      }
    }
    tmp=try(factor(cut(tmp.a, breaks = c(-Inf, q.a, Inf))), silent=T)
    
    do.cut=FALSE # if TRUE, use cut function which does not use weights
    # if there is a huge point mass, an error would occur, or it may not break into 3 groups
    if (inherits(tmp, "try-error")) do.cut=TRUE else if(length(table(tmp)) != 3) do.cut=TRUE
    
    if(!do.cut) {
      dat[[a %.% "cat"]] <- tmp
      marker.cutpoints[[a]] <- q.a
    } else {
      cat("\nfirst cut fails, call cut again with breaks=3 \n")
      # cut is more robust but it does not incorporate weights
      tmp=cut(tmp.a, breaks=3)
      stopifnot(length(table(tmp))==3)
      dat[[a %.% "cat"]] = tmp
      # extract cut points from factor level labels
      tmpname = names(table(tmp))[2]
      tmpname = substr(tmpname, 2, nchar(tmpname)-1)
      marker.cutpoints[[a]] <- as.numeric(strsplit(tmpname, ",")[[1]])
    }
    stopifnot(length(table(dat[[a %.% "cat"]])) == 3)
    if(verbose) {
      print(table(dat[[a %.% "cat"]]))
      cat("\n")
    }
  }
  
  attr(dat, "marker.cutpoints")=marker.cutpoints
  dat
  
}



###################################################################################################
# shared functions: make bootstrap samples

# bootstrap for COVE boost
# Within each quadrant (2 Trt * 2 naive status):
#   1. resample the cohort and count the number of cases and controls: n1 and n0
#         if n1 < 32 or n0 < 32, redo
#   2. resample 32 cases and 32 controls from ph2 samples 
#   3. resample n1-32 cases and n0-32 controls from non-ph2 samples
# 4. Collapse strata if needed and compute inverse probability sampling weights
# Thus, the number of cases may vary across bootstrap replicates, but the ph2 sample size remains constant

# dat.ph1 is ph1 data and need to have, in addition to markers and covariates columns:
#   Ptid, Trt, naive, EventIndPrimary, ph2, demo.stratum, CalendarBD1Interval
# return a dataframe with wt column

bootstrap.cove.boost=function(dat.ph1, seed) {
  
  set.seed(seed)
  
  # perform bootstrap within each quadrant (2 Trt * 2 naive status)
  dat.b=NULL
  for (idat in 1:4) {
    if (idat==1) {dat.tmp = subset(dat.ph1, Trt==1 & naive==1)}
    if (idat==2) {dat.tmp = subset(dat.ph1, Trt==0 & naive==1)}
    if (idat==3) {dat.tmp = subset(dat.ph1, Trt==1 & naive==0)}
    if (idat==4) {dat.tmp = subset(dat.ph1, Trt==0 & naive==0)}
    if (nrow(dat.tmp)==0) next
    
    dat.tmp.nph2=subset(dat.tmp, !ph2)
    dat.tmp.ph2=subset(dat.tmp, ph2)
    
    # n1.ph2 and n0.ph2 are expected to be 32 in COVE Boost
    # we make it data-dependent here to be more flexible
    n1.ph2 = sum(dat.tmp.ph2$EventIndPrimary)
    n0.ph2 = sum(1-dat.tmp.ph2$EventIndPrimary)
    
    # 1. 
    dat.2=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),]
    n1 = nrow(subset(dat.2, EventIndPrimary==1))
    n0 = nrow(subset(dat.2, EventIndPrimary==0))
    
    while(n1<n1.ph2 | n0<n0.ph2) {   
      dat.2=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),]
      n1 = nrow(subset(dat.2, EventIndPrimary==1))
      n0 = nrow(subset(dat.2, EventIndPrimary==0))
    }
    
    # 2.
    dat.ph2.cases=subset(dat.tmp.ph2, EventIndPrimary==1)
    dat.ph2.cases.b=dat.ph2.cases[sample.int(nrow(dat.ph2.cases), size=n1.ph2, r=TRUE),]
    
    dat.ph2.ctrls=subset(dat.tmp.ph2, EventIndPrimary==0)
    dat.ph2.ctrls.b=dat.ph2.ctrls[sample.int(nrow(dat.ph2.ctrls), size=n0.ph2, r=TRUE),]
    
    # 3.
    dat.nph2.cases=subset(dat.tmp.nph2, EventIndPrimary==1)
    dat.nph2.cases.b=dat.nph2.cases[sample.int(nrow(dat.nph2.cases), size=n1-n1.ph2, r=TRUE),]
    
    dat.nph2.ctrls=subset(dat.tmp.nph2, EventIndPrimary==0)
    dat.nph2.ctrls.b=dat.nph2.ctrls[sample.int(nrow(dat.nph2.ctrls), size=n0-n0.ph2, r=TRUE),]
    
    dat.b=rbind(dat.b, dat.ph2.cases.b, dat.ph2.ctrls.b, dat.nph2.cases.b, dat.nph2.ctrls.b)
  }
  
  
  # 4. 
  n.demo = length(table(dat.b$demo.stratum))
  assertthat::assert_that(n.demo==6, msg = "n.demo != 6")
  
  # adjust Wstratum
  # this call may throw exceptions
  ret = cove.boost.collapse.strata (dat.b, n.demo)
  
  # compute inverse probability sampling weights
  tmp = with(ret, ph1)
  wts_table <- with(ret[tmp,], table(Wstratum, ph2))
  wts_norm <- rowSums(wts_table) / wts_table[, 2]
  ret[["wt"]] = ifelse(ret$ph1, wts_norm[ret$Wstratum %.% ""], NA)
  
  assertthat::assert_that(
    all(!is.na(subset(ret, tmp & !is.na(Wstratum))[["wt"]])),
    msg = "missing wt.BD for D analyses ph1 subjects")
  
  return (ret)
}




# a second version, simpler, faster, results are close to bootstrap.cove.boost
bootstrap.cove.boost.2=function(dat.ph1, seed) {
  
  set.seed(seed)
  
  # perform bootstrap within each quadrant (2 Trt * 2 naive status)
  dat.b=NULL
  for (idat in 1:4) {
    if (idat==1) {dat.tmp = subset(dat.ph1, Trt==1 & naive==1)}
    if (idat==2) {dat.tmp = subset(dat.ph1, Trt==0 & naive==1)}
    if (idat==3) {dat.tmp = subset(dat.ph1, Trt==1 & naive==0)}
    if (idat==4) {dat.tmp = subset(dat.ph1, Trt==0 & naive==0)}
    if (nrow(dat.tmp)==0) next
    
    dat.b=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),]
  }
  
  # 4. adjust Wstratum
  n.demo = length(table(dat.b$demo.stratum))
  assertthat::assert_that(n.demo==6, msg = "n.demo != 6")
  ret = cove.boost.collapse.strata (dat.b, n.demo)
  
  # compute inverse probability sampling weights
  tmp = with(ret, ph1)
  wts_table <- with(ret[tmp,], table(Wstratum, ph2))
  wts_norm <- rowSums(wts_table) / wts_table[, 2]
  ret[["wt"]] = ifelse(ret$ph1, wts_norm[ret$Wstratum %.% ""], NA)
  
  assertthat::assert_that(
    all(!is.na(subset(ret, tmp & !is.na(Wstratum))[["wt"]])),
    msg = "missing wt.BD for D analyses ph1 subjects")
  
  return (ret)
}




# when all cases are sampled, bootstrap from case control studies is done by resampling cases, ph2 controls, and non-ph2 controls separately. 
# e.g. hvtn705

# Across bootstrap replicates, the number of cases does not stay constant, neither do the numbers of ph2 controls by demographics strata. 
# Specifically,
# 1) sample with replacement to get dat.b. From this dataset, take the cases and count ph2 and non-ph2 controls by strata
# 2) sample with replacement ph2 and non-ph2 controls by strata

bootstrap.case.control.samples=function(dat.ph1, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2", min.cell.size=1) {
  #dat.ph1=dat.tmp; delta.name="EventIndPrimary"; strata.name="tps.stratum"; ph2.name="ph2"; min.cell.size=0
  
  set.seed(seed)
  
  dat.tmp=data.frame(ptid=1:nrow(dat.ph1), delta=dat.ph1[,delta.name], strata=dat.ph1[,strata.name], ph2=dat.ph1[,ph2.name])
  
  nn.ph1=with(dat.tmp, table(strata, delta))
  strat=rownames(nn.ph1); names(strat)=strat
  # ctrl.ptids is a list of lists
  ctrl.ptids = with(subset(dat.tmp, delta==0), lapply(strat, function (i) list(ph2=ptid[strata==i & ph2], nonph2=ptid[strata==i & !ph2])))
  
  # 1. resample dat.ph1 to get dat.b, but only take the cases 
  dat.b=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),]
  
  # re-do resampling if the bootstrap dataset has too few samples in a cell in nn.ctrl.b
  while(TRUE) {   
    nn.ctrl.b=with(subset(dat.b, !delta), table(strata, ph2))
    if (min(nn.ctrl.b)<min.cell.size | ncol(nn.ctrl.b)<2) dat.b=dat.tmp[sample.int(nrow(dat.tmp), r=TRUE),] else break
  }
  
  # take the case ptids
  case.ptids.b = dat.b$ptid[dat.b$delta==1]
  
  # 2. resample controls in dat.ph1 (numbers determined by dat.b) stratified by strata and ph2/nonph2
  # ph2 and non-ph2 controls by strata
  nn.ctrl.b=with(subset(dat.b, !delta), table(strata, ph2))
  # sample the control ptids
  ctrl.ptids.by.stratum.b=lapply(strat, function (i) {
    c(sample(ctrl.ptids[[i]]$ph2, nn.ctrl.b[i,2], r=T),
      sample(ctrl.ptids[[i]]$nonph2, nn.ctrl.b[i,1], r=T))
  })
  ctrl.ptids.b=do.call(c, ctrl.ptids.by.stratum.b)    
  
  # return data frame
  dat.ph1[c(case.ptids.b, ctrl.ptids.b), ]
}

## testing
#dat.b=bootstrap.case.control.samples(dat.vac.seroneg)
#with(dat.vac.seroneg, table(ph2, tps.stratum, EventIndPrimary))
#with(dat.b, table(ph2, tps.stratum, EventIndPrimary))
#> with(dat.vac.seroneg, table(ph2, tps.stratum, EventIndPrimary))
#, , EventIndPrimary = 0
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE 1483  915  759  439 1677 1138  894  591 3018 1973 1559 1051 1111  693  511  329
#  TRUE    57   53   55   57   56   57   57   56   58   55   55   57   57   56   56   56
#
#, , EventIndPrimary = 1
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE    1    0    0    1    0    1    0    0    2    1    2    1    0    0    0    1
#  TRUE     3    7    7   10    8   11    2   13   17   23   15   23    5    6    4    6
#
#> with(dat.b, table(ph2, tps.stratum, EventIndPrimary))
#, , EventIndPrimary = 0
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE 1487  911  750  462 1675 1181  884  570 3058 2023 1499 1034 1094  694  487  329
#  TRUE    47   57   65   62   50   53   50   64   55   61   65   53   64   53   54   60
#
#, , EventIndPrimary = 1
#
#       tps.stratum
#ph2       33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48
#  FALSE    0    0    0    0    0    2    0    0    1    1    3    3    0    0    0    2
#  TRUE     2    6    8    5    9   13    0   11   20   26   10   20    4    3    4    5


# for bootstrap use
get.ptids.by.stratum.for.bootstrap = function(data) {
  strat=sort(unique(data$tps.stratum))
  ptids.by.stratum=lapply(strat, function (i) 
    list(subcohort=subset(data, tps.stratum==i & SubcohortInd==1, Ptid, drop=TRUE), nonsubcohort=subset(data, tps.stratum==i & SubcohortInd==0, Ptid, drop=TRUE))
  )    
  # add a pseudo-stratum for subjects with NA in tps.stratum (not part of Subcohort). 
  # we need this group because it contains some cases with missing tps.stratum
  # if data is ph2 only, then this group is only cases because ph2 = subcohort + cases
  tmp=list(subcohort=subset(data, is.na(tps.stratum), Ptid, drop=TRUE),               nonsubcohort=NULL)
  ptids.by.stratum=append(ptids.by.stratum, list(tmp))    
  ptids.by.stratum
}


# bootstrap case cohort samples
# data is assumed to contain only ph1 ptids
get.bootstrap.data.cor = function(data, ptids.by.stratum, seed) {
  set.seed(seed)    
  
  # For each sampling stratum, bootstrap samples in subcohort and not in subchort separately
  tmp=lapply(ptids.by.stratum, function(x) c(sample(x$subcohort, r=TRUE), sample(x$nonsubcohort, r=TRUE)))
  
  dat.b=data[match(unlist(tmp), data$Ptid),]
  
  # compute weights
  tmp=with(dat.b, table(Wstratum, ph2))
  weights=rowSums(tmp)/tmp[,2]
  dat.b$wt=weights[""%.%dat.b$Wstratum]
  # we assume data only contains ph1 ptids, thus weights is defined for every bootstrapped ptids
  
  dat.b
}




###################################################################################################
# shared functions: plotting

# get plotting range
get.xlim=function(dat, marker) {
  assay=marker.name.to.assay(a)
  
  # the default
  ret=range(dat[[marker]], log10(lloxs[assay]/2), na.rm=T)
  
  # may be customized, e.g. to have the same xlim for different variants in the same type of assay
  # if (TRIAL=="moderna_boost") {
  #   if(assay %in% c("bindSpike", "bindRBD")) {
  #     ret=range(dat[["Day"%.%time%.%"bindSpike"]], 
  #               dat[["Day"%.%time%.%"bindRBD"]], 
  #               log10(lloxs[c("bindSpike","bindRBD")]/2), na.rm=T)
  #     
  #   } 
  # }

  delta=(ret[2]-ret[1])/20     
  c(ret[1]-delta, ret[2]+delta)
}

draw.x.axis.cor=function(xlim, llox, llox.label){
        
    xx=seq(ceiling(xlim[1]), floor(xlim[2]))        
    if (is.na(llox)) {
        for (x in xx) {
            axis(1, at=x, labels=if (x>=3) bquote(10^.(x)) else 10^x )    
        }
    } else if (llox.label=="delta") {
        for (x in xx) {
            axis(1, at=x, labels=if (x>=3 | x<=-3) bquote(10^.(x)) else 10^x )    
        }    
    } else {
        axis(1, at=log10(llox), labels=llox.label)
        for (x in xx[xx>log10(llox*1.8)]) {
            axis(1, at=x, labels= if(x>=3) bquote(10^.(x)) else 10^x)
        }
    }
    
    # add e.g. 30 between 10 and 100
    if (length(xx)<=3 & length(xx)>1) { 
        # a hack for prevent19 ID50 to not draw 3 b/c it is too close to LOD
        tmp=2:length(xx)
        if (study_name=="PREVENT19") tmp=3:length(xx)
        for (i in tmp) {
            x=xx[i-1]
            axis(1, at=x+log10(3), labels=if (x>=3) bquote(3%*%10^.(x)) else 3*10^x )
        }
    }
    
}

##### Copy of draw.x.axis.cor but returns the x-axis ticks and labels
# This is necessary if one works with ggplot as the "axis" function does not work.
get.labels.x.axis.cor=function(xlim, llox){
  xx=seq(floor(xlim[1]), ceiling(xlim[2]))
  if (!is.na(llox)) xx=xx[xx>log10(llox*2)]
  x_ticks <- xx
  if (is.na(llox)) {
      labels <- sapply(xx, function(x) {
        if (x>=3) bquote(10^.(x)) else 10^x
      })
  } else {
      labels <- sapply(xx, function(x) {
        if (log10(llox)==x) config$llox_label else if (x>=3) bquote(10^.(x)) else 10^x
      })
      #if(!any(log10(llox)==x_ticks)){
        x_ticks <- c(log10(llox), x_ticks)
        labels <- c(config$llox_label, labels)
      #}
  }
  return(list(ticks = x_ticks, labels = labels))
}




# get histogram object to add to VE plots etc
get.marker.histogram=function(marker, wt, trial, marker.break=marker) {
  # first call hist to get breaks, then call weighted.hist
  tmp.1=hist(marker.break,breaks=ifelse(trial=="moderna_real",25,15),plot=F)  # 15 is treated as a suggestion and the actual number of breaks is determined by pretty()
  tmp=weighted.hist(marker,wt, breaks=tmp.1$breaks, plot=F)
  attr(tmp,"class")="histogram" 
  tmp
}



###################################################################################################
# shared functions: for tables

# x is the marker values
# assay is one of assays, e.g. pseudoneutid80
report.assay.values=function(x, assay){
    lars.quantiles=seq(0,1,length.out=30) [round(seq.int(1, 30, length.out = 10))]
    sens.quantiles=c(0.15, 0.85)
    # cannot have different lengths for different assays, otherwise downstream code may break
    fixed.values = log10(c("500"=500, "1000"=1000))
    # if we want to add "llox/2"=unname(lloxs[assay]/2))) to fixed.values, we have to get assay right, which will take some thought because marker.name.to.assay is hardcoded
    out=sort(c(quantile(x, c(lars.quantiles,sens.quantiles), na.rm=TRUE), fixed.values[fixed.values<max(x, na.rm=T) & fixed.values>min(x, na.rm=T)]))    
    out
    #out[!duplicated(out)] # unique strips away the names. But don't take out duplicates because 15% may be needed and because we may want the same number of values for each assay
}
#report.assay.values (dat.vac.seroneg[["Day57pseudoneutid80"]], "pseudoneutid80")



# a function to print tables of cases counts with different marker availability
# note that D57 cases and intercurrent cases may add up to more than D29 cases because ph1.D57 requires EarlyendpointD57==0 while ph1.D29 requires EarlyendpointD29==0
make.case.count.marker.availability.table=function(dat) {
    if (study_name=="COVE" | study_name=="MockCOVE" ) {
        idx.trt=1:0
        names(idx.trt)=c("vacc","plac")
        cnts = sapply (idx.trt, simplify="array", function(trt) {
             idx=1:3
             names(idx)=c("Day 29 Cases", "Day 57 Cases", "Intercurrent Cases")
             tab=t(sapply (idx, function(i) {           
                tmp.1 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(BbindSpike)     | is.na(BbindRBD) )
                tmp.2 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(Day29bindSpike) | is.na(Day29bindRBD))
                tmp.3 = with(subset(dat, Trt==trt & Bserostatus==0 & (if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29) &   (if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases)), is.na(Day57bindSpike) | is.na(Day57bindRBD))    
                
                c(sum(tmp.1 & tmp.2 & tmp.3), sum(tmp.1 & tmp.2 & !tmp.3), sum(tmp.1 & !tmp.2 & tmp.3), sum(tmp.1 & !tmp.2 & !tmp.3), 
                  sum(!tmp.1 & tmp.2 & tmp.3), sum(!tmp.1 & tmp.2 & !tmp.3), sum(!tmp.1 & !tmp.2 & tmp.3), sum(!tmp.1 & !tmp.2 & !tmp.3))
            }))
            colnames(tab)=c("---", "--+", "-+-", "-++", "+--", "+-+", "++-", "+++")
            tab
        })
        cnts
    } else if (study_name=="ENSEMBLE" | study_name=="MockENSEMBLE" ) {
        idx.trt=1:0
        names(idx.trt)=c("vacc","plac")
        cnts = sapply (idx.trt, simplify="array", function(trt) {
             idx=1:1
             tab=t(sapply (idx, function(i) {           
                tmp.1 = with(subset(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29 &   if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases), is.na(BbindSpike)     | is.na(BbindRBD) )
                tmp.2 = with(subset(dat, Trt==trt & Bserostatus==0 & if(i==2) EventIndPrimaryD57 else EventIndPrimaryD29 &   if(i==2) ph1.D57 else if(i==1) ph1.D29 else ph1.intercurrent.cases), is.na(Day29bindSpike) | is.na(Day29bindRBD))
                
                c(sum(tmp.1 & tmp.2), sum(!tmp.1 & tmp.2), sum(tmp.1 & !tmp.2), sum(!tmp.1 & !tmp.2))
             }))
             colnames(tab)=c("--", "+-", "-+", "++")
             tab
        })
        t(drop(cnts))
    } else {
        NA
    }
}
#make.case.count.marker.availability.table(dat.mock)




###############################################################################
# theme options

# fixed knitr chunk options
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  out.width = "80%",
  out.extra = "",
  fig.pos = "H",
  fig.show = "hold",
  fig.align = "center",
  fig.width = 6,
  fig.asp = 0.618,
  fig.retina = 0.8,
  dpi = 600,
  echo = FALSE,
  message = FALSE,
  warning = FALSE
)

# global options
options(
  digits = 6,
  #scipen = 999,
  dplyr.print_min = 6,
  dplyr.print_max = 6,
  crayon.enabled = FALSE,
  bookdown.clean_book = TRUE,
  knitr.kable.NA = "NA",
  repos = structure(c(CRAN = "https://cran.rstudio.com/"))
)

# no complaints from installation warnings
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

# overwrite options by output type
if (knitr:::is_html_output()) {
  #options(width = 80)
  
  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}
if (knitr:::is_latex_output()) {
  #knitr::opts_chunk$set(width = 67)
  #options(width = 67)
  options(cli.unicode = TRUE)
  
  # automatically create a bib database for R packages
  knitr::write_bib(c(
    .packages(), "bookdown", "knitr", "rmarkdown"
  ), "packages.bib")
}

# create and set global ggplot theme
# borrowed from https://github.com/tidymodels/TMwR/blob/master/_common.R
theme_transparent <- function(...) {
  # use black-white theme as base
  ret <- ggplot2::theme_bw(...)
  
  # modify with transparencies
  trans_rect <- ggplot2::element_rect(fill = "transparent", colour = NA)
  ret$panel.background  <- trans_rect
  ret$plot.background   <- trans_rect
  ret$legend.background <- trans_rect
  ret$legend.key        <- trans_rect
  
  # always have legend below
  ret$legend.position <- "bottom"
  return(ret)
}

library(ggplot2)
theme_set(theme_transparent())
theme_update(
  text = element_text(size = 25),
  axis.text.x = element_text(colour = "black", size = 30),
  axis.text.y = element_text(colour = "black", size = 30)
)

# custom ggsave function with updated defaults
ggsave_custom <- function(filename = default_name(plot),
                          height= 15, width = 21, ...) {
  ggsave(filename = filename, height = height, width = width, ...)
}



###############################################################################

# race labeling
labels.race <- c(
  "White", 
  "Black or African American",
  "Asian", 
  if ((study_name=="ENSEMBLE" | study_name=="MockENSEMBLE") & startsWith(attr(config, "config"),"janssen_la")) "Indigenous South American" else "American Indian or Alaska Native",
  "Native Hawaiian or Other Pacific Islander", 
  "Multiracial",
  if ((study_name=="COVE" | study_name=="MockCOVE")) "Other", 
  "Not reported and unknown"
)

# ethnicity labeling
labels.ethnicity <- c(
  "Hispanic or Latino", "Not Hispanic or Latino",
  "Not reported and unknown"
)

# baseline stratum labeling
if (study_name=="COVEBoost") {
  demo.stratum.labels <- c(
    "Age >= 65, URM",
    "Age < 65, At risk, URM",
    "Age < 65, Not at risk, URM",
    "Age >= 65, White non-Hisp",
    "Age < 65, At risk, White non-Hisp",
    "Age < 65, Not at risk, White non-Hisp"
  )
} else stop("unknown study_name 1")

