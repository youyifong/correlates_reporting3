renv::activate(project = here::here(".."))     

#Sys.setenv(TRIAL = "moderna_boost"); Sys.setenv(VERBOSE = 1)
source(here::here("..", "_common.R"))
source("code/params.R")

library(survey)
library(plotrix) # weighted.hist
library(forestplot)
library(Hmisc) # wtd.quantile, cut2

time.start=Sys.time()
myprint(TRIAL)
myprint(verbose)
begin=Sys.time()
print(date())

# need this function b/c svycoxh may error due to singularity if, e.g. all cases have the same marker value
run.svycoxph=function(f, design) {
  fit=try(svycoxph(f, design=design), silent=T)
  if (class(fit)[1]=="try-error") NA else fit
}


# path for figures and tables etc
save.results.to.0 = here::here("output");                                if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)
save.results.to.0 = paste0(save.results.to.0, "/", attr(config,"config")); if (!dir.exists(save.results.to.0))  dir.create(save.results.to.0)

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for one analysis
B <-       config$num_boot_replicates 
numPerm <- config$num_perm_replicates # number permutation replicates 1e4
myprint(B, numPerm)


# uloq censoring, done here b/c should not be done for immunogenicity reports
# note that if delta are used, delta needs to be recomputed
for (a in assays) {
  uloq=assay_metadata$uloq[assay_metadata$assay==a]
  for (t in c("BD1", "BD29")  ) {
    dat[[t %.% a]] <- ifelse(dat[[t %.% a]] > log10(uloq), log10(uloq), dat[[t %.% a]])
  }
}    


dat$ph1=dat$ph1.BD29
dat$ph2=dat$ph2.BD29
dat$wt=dat$wt.BD29
dat$EventIndPrimary =dat$EventIndOmicronBD29
dat$EventTimePrimary=dat$EventTimeOmicronBD29
dat$yy=dat$EventIndPrimary

# subset to ph1
dat=subset(dat, ph1.BD29)

# add trichotomized markers. use the same cutpoints for naive and nnaive
obj.assays=c("bindSpike_BA.1", "pseudoneutid50_BA.1", "bindSpike", "pseudoneutid50")  
all.markers = c(paste0("BD29", obj.assays), paste0("DeltaBD29overBD1", obj.assays))
names(all.markers)=all.markers
dat = add.trichotomized.markers (dat, all.markers)
marker.cutpoints=attr(dat, "marker.cutpoints")
  
# define naive and nnaive datasets
dat.naive=subset(dat,   naive & ph1.BD29)
dat.nnaive=subset(dat, !naive & ph1.BD29)

# compute overall risks
prev.naive  = get.marginalized.risk.no.marker(form.0, dat.naive,  tfinal.tpeak)
prev.nnaive = get.marginalized.risk.no.marker(form.0, dat.nnaive, tfinal.tpeak)
overall.risks=list(prev.naive, prev.nnaive)
myprint(prev.naive, prev.nnaive)



###################################################################################################
# Obj 1: To assess BD29 omicron Ab as a correlate of risk (CoR) against omicron COVID-19 
# Obj 2: To assess fold-rise in omicron Ab from BD1/pre-booster to BD29 as a CoR against omicron COVID-19


for (iObj in 1:2) {
# iObj=1
  
  if (iObj==1) {
    all.markers = paste0("BD29", obj.assays)
  } else if (iObj==2) {
    all.markers = paste0("DeltaBD29overBD1", obj.assays)
  }
  names(all.markers)=all.markers

  # loop through naive and nonnaive
  for (idat in 2:2) {
    # idat=2
    
    myprint(idat)
    if (idat==1) {dat.ph1 = dat.naive;  save.results.to = glue("{save.results.to.0}/obj{iObj}_naive/")}
    if (idat==2) {dat.ph1 = dat.nnaive; save.results.to = glue("{save.results.to.0}/obj{iObj}_nnaive/")}

    if (!dir.exists(save.results.to))  dir.create(save.results.to)
    print(paste0("save results to ", save.results.to))
    
    # save info
    write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk"))
    for (a in all.markers) {        
      q.a=marker.cutpoints[[a]]
      to.write = paste0(gsub("_", "\\\\_",a), " [", concatList(round(q.a, 2), ", "), ")%")
      write(to.write, file=paste0(save.results.to, "cutpoints_", a))
    }
    
    # create data objects
    dat.ph2 = subset(dat.ph1, ph2)
    design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
    with(dat.ph1, table(Wstratum, ph2, useNA="ifany"))
    
    # table of ph1 and ph2 cases
    tab=with(dat.ph1, table(ph2, EventIndPrimary))
    names(dimnames(tab))[2]="Event Indicator"
    print(tab)
    mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)

    # prepare for Cox models runs and margialized risks
    tpeak=29
    # the origin of followup days, may be different from tpeak, e.g., D43start48
    tpeak1 = 29
    
    source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))
    
    source(here::here("code", "cor_coxph_ph.R"))
    
    # unit testing
    if (TRIAL == "") {
      tmp.1=c(rv$tab.1[,4], rv$tab.2[,"overall.p.0"])
      tmp.2=c("0.162","0.079","0.006",      "0.498","   ","   ","0.162","   ","   ","0.003","   ","   ")
      assertthat::assert_that(all(tmp.1==tmp.2), msg = "failed cor_coxph unit testing")    
      print("Passed cor_coxph unit testing")    
    } 
    
    # marginalized risk and controlled VE
    comp.risk=F
    COR=idat # only used in table figure labels
    
    source(here::here("code", "cor_coxph_marginalized_risk_bootstrap.R"))
    
    source(here::here("code", "cor_coxph_marginalized_risk_plotting.R"))
    
    # source(here::here("code", "cor_coxph_samplesizeratio.R"))

  } 
}



###################################################################################################
# Obj 3: To assess whether the CoR in 1. or 2. is modified by SARS-CoV-2 naive/non-naive status

for (iObj in 1:2) {
  
  if (iObj==1) {
    all.markers = paste0("BD29", obj.assays)
  } else if (iObj==2) {
    all.markers = paste0("DeltaBD29overBD1", obj.assays)
  }
  names(all.markers)=all.markers

  dat.ph1 = dat

  ilabel="all"
  
  save.results.to = glue("{save.results.to.0}/obj3_{iObj}/")
  if (!dir.exists(save.results.to))  dir.create(save.results.to)
  print(paste0("save results to ", save.results.to))
  
  # save info
  write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk"))
  
  # create data objects
  design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
  with(dat.ph1, table(Wstratum, ph2, useNA="ifany"))
  
  # prepare for Cox models runs and margialized risks
  tpeak=29
  # the origin of followup days, may be different from tpeak, e.g., D43start48
  tpeak1 = 29
  
  # fit the interaction model and save regression results to a table
  for (a in all.markers) {
    f= update(form.0, as.formula(paste0("~.+naive*", a)))
    fits=list(svycoxph(f, design=design.1))
    est=getFormattedSummary(fits, exp=T, robust=T, type=1)
    ci= getFormattedSummary(fits, exp=T, robust=T, type=13)
    est = paste0(est, " ", ci)
    p=  getFormattedSummary(fits, exp=T, robust=T, type=10)
    # # generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
    # var.ind=5:7
    # stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
    # p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
    # put together the table
    tab=cbind(est, p)
    colnames(tab)=c("HR", "P value")
    # tab=rbind(tab, "Generalized Wald Test for Markers"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
    tab
    mytex(tab, file.name=paste0("CoR_itxnnaive_",a), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
    # itxn.pvals=c(itxn.pvals, last(getFixedEf(fit)[,"p.value"]))
  }

  # names(itxn.pvals)=config$interaction
  # itxn.pvals=itxn.pvals[!contain(config$interaction, "ICS4AnyEnv")] # remove the ones with ICS4AnyEnv
  # itx.pvals.adj.fdr=p.adjust(itxn.pvals, method="fdr")
  # itx.pvals.adj.hol=p.adjust(itxn.pvals, method="holm")
  # tab=cbind(itxn.pvals, itx.pvals.adj.hol, itx.pvals.adj.fdr)
  # colnames(tab)=c("interaction P value", "FWER", "FDR")
  # mytex(tab, file.name="CoR_itxn_multitesting", align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
}


###################################################################################################
# Obj 4: To assess whether the CoR in 1. is modified by the BD1 antibody value

for (idat in 1:2) {
  # idat=1
  
  all.markers = paste0("BD29", obj.assays)
  names(all.markers)=all.markers
  
  myprint(idat)
  if (idat==1) {dat.ph1 = dat.naive;  save.results.to = glue("{save.results.to.0}/obj4_naive/")}
  if (idat==2) {dat.ph1 = dat.nnaive; save.results.to = glue("{save.results.to.0}/obj4_nnaive/")}
  
  if (!dir.exists(save.results.to))  dir.create(save.results.to)
  print(paste0("save results to ", save.results.to))
  
  # save info
  write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk"))

  # create data objects
  design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
  with(dat.ph1, table(Wstratum, ph2, useNA="ifany"))
  
  # prepare for Cox models runs and margialized risks
  tpeak=29
  # the origin of followup days, may be different from tpeak, e.g., D43start48
  tpeak1 = 29
  
  # fit the interaction model and save regression results to a table
  for (a in all.markers) {
    f= update(form.0, as.formula(paste0("~.+scale(",sub("BD29","BD1",a),")*scale(", a, ")")))
    fits=list(svycoxph(f, design=design.1))
    est=getFormattedSummary(fits, exp=T, robust=T, type=1)
    ci= getFormattedSummary(fits, exp=T, robust=T, type=13)
    est = paste0(est, " ", ci)
    p=  getFormattedSummary(fits, exp=T, robust=T, type=10)
    # # generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
    # var.ind=5:7
    # stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
    # p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
    # put together the table
    tab=cbind(est, p)
    colnames(tab)=c("HR", "P value")
    # tab=rbind(tab, "Generalized Wald Test for Markers"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
    tab
    mytex(tab, file.name=paste0("CoR_itxnBD1BD29_",a), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
    # itxn.pvals=c(itxn.pvals, last(getFixedEf(fit)[,"p.value"]))
  }
  
  # names(itxn.pvals)=config$interaction
  # itxn.pvals=itxn.pvals[!contain(config$interaction, "ICS4AnyEnv")] # remove the ones with ICS4AnyEnv
  # itx.pvals.adj.fdr=p.adjust(itxn.pvals, method="fdr")
  # itx.pvals.adj.hol=p.adjust(itxn.pvals, method="holm")
  # tab=cbind(itxn.pvals, itx.pvals.adj.hol, itx.pvals.adj.fdr)
  # colnames(tab)=c("interaction P value", "FWER", "FDR")
  # mytex(tab, file.name="CoR_itxn_multitesting", align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
}


# # joint distribution of BD1 and BD29 markers
# par(mfrow=c(1,2))
# with(dat.naive, plot(BD1bindSpike_BA.1, BD29bindSpike_BA.1, col=ifelse(EventIndPrimary, 2, 1),  main=paste0("Naive, cor ",round(cor(BD1bindSpike_BA.1, BD29bindSpike_BA.1, use="complete.obs"),2))))
# abline(0,1)
# with(dat.nnaive, plot(BD1bindSpike_BA.1, BD29bindSpike_BA.1, col=ifelse(EventIndPrimary, 2, 1), main=paste0("Naive, cor ",round(cor(BD1bindSpike_BA.1, BD29bindSpike_BA.1, use="complete.obs"),2))))
# abline(0,1)



###################################################################################################
# multivariate_assays models
# pooled over naive and nnaive

if (!is.null(config$multivariate_assays)) {
  if(verbose) print("Multiple regression")
  
  dat.ph1 = subset(dat, ph1.BD29)
  design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
  with(dat.ph1, table(Wstratum, ph2, useNA="ifany"))
  
  for (a in config$multivariate_assays) {
    for (i in 1:2) {
      # 1: per SD; 2: per 10-fold
      a.tmp=a
      aa=trim(strsplit(a, "\\+")[[1]])
      for (x in aa[!contain(aa, "\\*")]) {
        # replace every variable with scale(x) when i==1
        a.tmp=gsub(x, paste0(if(i==1) "scale","(",x,")"), a.tmp) 
      }
      f= update(form.0, as.formula(paste0("~.+", a.tmp)))
      fit=svycoxph(f, design=design.1) 
      var.ind=length(coef(fit)) - length(aa):1 + 1
      
      fits=list(fit)
      est=getFormattedSummary(fits, exp=T, robust=T, rows=var.ind, type=1)
      ci= getFormattedSummary(fits, exp=T, robust=T, rows=var.ind, type=13)
      est = paste0(est, " ", ci)
      p=  getFormattedSummary(fits, exp=T, robust=T, rows=var.ind, type=10)
      
      #generalized Wald test for whether the set of markers has any correlation (rejecting the complete null)
      stat=coef(fit)[var.ind] %*% solve(vcov(fit)[var.ind,var.ind]) %*% coef(fit)[var.ind] 
      p.gwald=pchisq(stat, length(var.ind), lower.tail = FALSE)
      
      tab=cbind(est, p)
      colnames(tab)=c(paste0("HR per ",ifelse(i==1,"sd","10 fold")," incr."), "P value")
      tab
      tab=rbind(tab, "Generalized Wald Test"=c("", formatDouble(p.gwald,3, remove.leading0 = F)))
      
      mytex(tab, file.name=paste0("CoR_multivariable_svycoxph_pretty", match(a, config$multivariate_assays), if(i==2) "_per10fold", study_name), align="c", include.colnames = T, save2input.only=T, 
            input.foldername=save.results.to)
    }
  }
  
}




print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
