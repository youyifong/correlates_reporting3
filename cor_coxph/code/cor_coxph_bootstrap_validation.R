# checking if cove.boost.collapse.strata is sensible

Sys.setenv(TRIAL = "moderna_boost"); Sys.setenv(VERBOSE = 1)

renv::activate(project = here::here(".."))     
source(here::here("..", "_common.R"))
source("code/params.R")

library(survey)
library(plotrix) # weighted.hist
library(forestplot)
library(Hmisc) # wtd.quantile, cut2

time.start=Sys.time()
TRIAL=Sys.getenv("TRIAL")
myprint(TRIAL)
myprint(verbose)
begin=Sys.time()
print(date())

# need this function b/c svycoxh may error due to singularity if, e.g. all cases have the same marker value
run.svycoxph=function(f, design) {
  fit=try(svycoxph(f, design=design), silent=T)
  if (class(fit)[1]=="try-error") NA else fit
}

# read analysis ready data
dat = read.csv(config$data_cleaned)

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

# define subsets of data    
dat.vac.naive=subset(dat,  Trt==1 & naive & ph1.BD29)
dat.pla.naive=subset(dat,  Trt==0 & naive & ph1.BD29)
dat.vac.nnaive=subset(dat, Trt==1 & !naive & ph1.BD29)
dat.pla.nnaive=subset(dat, Trt==0 & !naive & ph1.BD29)

idat=4
myprint(idat)
if (idat==1) {dat.ph1 = dat.vac.naive;  ilabel="vac_naive"}
if (idat==2) {dat.ph1 = dat.pla.naive;  ilabel="pla_naive"}
if (idat==3) {dat.ph1 = dat.vac.nnaive; ilabel="vac_nnaive"}
if (idat==4) {dat.ph1 = dat.pla.nnaive; ilabel="pla_nnaive"}


dat.ph1$ph1=dat.ph1$ph1.BD29
dat.ph1$ph2=dat.ph1$ph2.BD29
dat.ph1$wt=dat.ph1$wt.BD29
dat.ph1$EventIndPrimary =dat.ph1$EventIndOmicronBD29
dat.ph1$EventTimePrimary=dat.ph1$EventTimeOmicronBD29

dat.ph1$yy=dat.ph1$EventIndPrimary

dat.ph2 = subset(dat.ph1, ph2)


# table of ph1 and ph2 cases
tab=with(dat.ph1, table(ph2, EventIndPrimary))
names(dimnames(tab))[2]="Event Indicator"
print(tab)
mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)

design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)

f=Surv(EventTimePrimary, EventIndPrimary) ~ MinorityInd + HighRiskInd + risk_score + BD29bindSpike

fit=svycoxph(f, design=design.1); fit

dat.ph2=subset(dat.ph1, ph2)
coxph(f, dat.ph2, weights=dat.ph2$wt) # same est as fit


time.start=Sys.time()
out=mclapply(1:100, mc.cores = 1, FUN=function(seed) {  
  myprint(seed) 
  # dat.b = try(bootstrap.cove.boost(dat.ph1, seed))
  # if (inherits (dat.b, "try-error")) return (NULL)
  dat.b = bootstrap.cove.boost(dat.ph1, seed)
  
  dat.ph2.b=subset(dat.b, ph2)
  fit=try(do.call("coxph", list(formula=f, data=dat.ph2.b, weights=dat.ph2.b$wt)))
  if(!inherits(fit,"try-error")) fit$coefficients else NULL
})
boot=do.call(cbind, out)
print("run time: "%.%format(Sys.time()-time.start, digits=1)) # 50 sec for 100 replicates



# without calling cove.boost.collapse.strata. does not work well
bootstrap.cove.boost.3=function(dat.ph1, seed) {
  
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
  
  ret=dat.b
  
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



time.start=Sys.time()
out.2=mclapply(1:100, mc.cores = 1, FUN=function(seed) {  
  myprint(seed) 
  # dat.b = try(bootstrap.cove.boost(dat.ph1, seed))
  # if (inherits (dat.b, "try-error")) return (NULL)
  dat.b = bootstrap.cove.boost.2(dat.ph1, seed)
  
  dat.ph2.b=subset(dat.b, ph2)
  fit=try(do.call("coxph", list(formula=f, data=dat.ph2.b, weights=dat.ph2.b$wt)))
  if(!inherits(fit,"try-error")) fit$coefficients else NULL
})
boot.2=do.call(cbind, out.2)
print("run time: "%.%format(Sys.time()-time.start, digits=1))# 19 sec for 100 replicates


time.start=Sys.time()
out.3=mclapply(1:100, mc.cores = 1, FUN=function(seed) {  
  myprint(seed) 
  # dat.b = try(bootstrap.cove.boost(dat.ph1, seed))
  # if (inherits (dat.b, "try-error")) return (NULL)
  dat.b = bootstrap.cove.boost.3(dat.ph1, seed)
  
  dat.ph2.b=subset(dat.b, ph2)
  fit=try(do.call("coxph", list(formula=f, data=dat.ph2.b, weights=dat.ph2.b$wt)))
  if(!inherits(fit,"try-error")) fit$coefficients else NULL
})
boot.3=do.call(cbind, out.3)
print("run time: "%.%format(Sys.time()-time.start, digits=1))# 9 sec for 100 replicates


# comparing analytical stderr and bootstrap ones
# version 1 and 2 work similary
# version 3 does not work as well
fit
sd(boot["BD29bindSpike",])
summary(boot["BD29bindSpike",])
sd(boot.2["BD29bindSpike",])
summary(boot.2["BD29bindSpike",])
sd(boot.3["BD29bindSpike",])
summary(boot.3["BD29bindSpike",])


