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
dat.naive =subset(dat, naive==1)
dat.nnaive=subset(dat, naive==0)

# compute overall risks
prev.naive  = get.marginalized.risk.no.marker(form.0, dat.naive,  tfinal.tpeak)
prev.nnaive = get.marginalized.risk.no.marker(form.0, dat.nnaive, tfinal.tpeak)
prev.pooled = get.marginalized.risk.no.marker(form.0, dat, tfinal.tpeak)
overall.risks=list(prev.naive, prev.nnaive, prev.pooled)
myprint(prev.naive, prev.nnaive, prev.pooled)



###################################################################################################
cat("Univariate analyses")

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

  # loop through (1) naive, (2) nnaive, (3) pooled
  for (iNaive in 1:3) {
    # iNaive=3
    
    myprint(iNaive)
    if (iNaive==1) {dat.ph1 = dat.naive;  save.results.to = glue("{save.results.to.0}/obj{iObj}_naive/")}
    if (iNaive==2) {dat.ph1 = dat.nnaive; save.results.to = glue("{save.results.to.0}/obj{iObj}_nnaive/")}
    if (iNaive==3) {dat.ph1 = dat; save.results.to = glue("{save.results.to.0}/obj{iObj}_pooled/")}
    
    # dat.ph1 = subset(dat.ph1, Ptid!="US3302292")
    # dat.ph1 = subset(dat.ph1, !Ptid %in% c("US3302292", "US3302397", "US3492199", "US3632155"))
    # dat.ph1 = subset(dat.ph1, !Ptid %in% c("US3302292", "US3302397", "US3492199", "US3632155", "US3642188"))
    # dat.ph1 = subset(dat.ph1, Ptid!="US3642188")
    # dat.ph1[dat.ph1$Ptid=="US3302292","ph2"] =F
    

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
    design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
    with(dat.ph1, table(Wstratum, ph2, useNA="ifany"))
    
    dat.ph2 = subset(dat.ph1, ph2)
    
    # table of ph1 and ph2 cases
    tab=with(dat.ph1, table(ph2, EventIndPrimary))
    names(dimnames(tab))[2]="Event Indicator"
    print(tab)
    mytex(tab, file.name="tab1", save2input.only=T, input.foldername=save.results.to)

    # prepare for Cox models runs and margialized risks
    tpeak=29
    # the origin of followup days, may be different from tpeak, e.g., D43start48
    tpeak1 = 29
    
    # Cox regression
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
    COR=iNaive # only used in table figure labels
    
    source(here::here("code", "cor_coxph_marginalized_risk_no_marker.R"))
    
    source(here::here("code", "cor_coxph_marginalized_risk_bootstrap.R"))
    
    source(here::here("code", "cor_coxph_marginalized_risk_plotting.R"))
    
    # source(here::here("code", "cor_coxph_samplesizeratio.R"))

  } 
}



###################################################################################################
cat("Obj 3: To assess whether the CoR in 1. or 2. is modified by SARS-CoV-2 naive/non-naive status by fitting interaction models")

dat.ph1 = dat

# create data objects
design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
with(dat.ph1, table(Wstratum, ph2, useNA="ifany"))

# prepare for Cox models runs and margialized risks
tpeak=29
# the origin of followup days, may be different from tpeak, e.g., D43start48
tpeak1 = 29


# first, assess whether any demo variables are modified by naive

f= update(form.0, ~.+naive*HighRiskInd)
summary(svycoxph(f, design=design.1))

f= update(form.0, ~.+naive*MinorityInd)
summary(svycoxph(f, design=design.1))

f= update(form.0, ~.+naive*risk_score)
summary(svycoxph(f, design=design.1))

# repeat twice, with and without adjusting for naive*risk_score
for (ind in 1:2) {
  
for (iObj in 1:2) {
# iObj=1
  
  if (iObj==1) {
    all.markers = paste0("BD29", obj.assays)
  } else if (iObj==2) {
    all.markers = paste0("DeltaBD29overBD1", obj.assays)
  }
  names(all.markers)=all.markers

  save.results.to = glue("{save.results.to.0}/obj3_{iObj}/")
  if (!dir.exists(save.results.to))  dir.create(save.results.to)
  print(paste0("save results to ", save.results.to))
  
  # save info
  write(tfinal.tpeak, file=paste0(save.results.to, "timepoints_cum_risk"))
  
  # fit the interaction model and save regression results to a table
  for (a in all.markers) {
    if (ind==1) {
      f= update(form.0, as.formula(paste0("~.+naive*", a)))
    } else{
      f= update(form.0, as.formula(paste0("~.+naive*", a, "+naive:HighRiskInd")))
    }
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
    mytex(tab, file.name=paste0("CoR_itxnnaive_",a,"_",ind), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
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
  
}


###################################################################################################
cat("Obj 4: To assess whether the CoR in 1. is modified by the BD1 antibody value")

for (iNaive in 1:2) {
  # iNaive=1
  
  myprint(iNaive)
  if (iNaive==1) {dat.ph1 = dat.naive;  save.results.to = glue("{save.results.to.0}/obj4_naive/")}
  if (iNaive==2) {dat.ph1 = dat.nnaive; save.results.to = glue("{save.results.to.0}/obj4_nnaive/")}
  
  dat.ph2 = subset(dat.ph1, ph2)

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
  for (ab in config$interaction) {
    tmp=trim(strsplit(ab, " *\\* *")[[1]]); a=tmp[1]; b=tmp[2]            
    f= update(form.0, as.formula(paste0("~.+scale(", a ,") * scale(", b, ")")))
    
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
    mytex(tab, file.name=paste0("CoR_itxn_",a,"_",b), align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)
    # itxn.pvals=c(itxn.pvals, last(getFixedEf(fit)[,"p.value"]))
  }
  
  # names(itxn.pvals)=config$interaction
  # itxn.pvals=itxn.pvals[!contain(config$interaction, "ICS4AnyEnv")] # remove the ones with ICS4AnyEnv
  # itx.pvals.adj.fdr=p.adjust(itxn.pvals, method="fdr")
  # itx.pvals.adj.hol=p.adjust(itxn.pvals, method="holm")
  # tab=cbind(itxn.pvals, itx.pvals.adj.hol, itx.pvals.adj.fdr)
  # colnames(tab)=c("interaction P value", "FWER", "FDR")
  # mytex(tab, file.name="CoR_itxn_multitesting", align="c", include.colnames = T, save2input.only=T, input.foldername=save.results.to)

  
  # bootstrap risk curves
  if(!file.exists(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))) {    
    cat("make itxn.marginalized.risk\n")
    
    risks.itxn=list()      
    for (ab in config$interaction) {
      myprint(ab)
      tmp=trim(strsplit(ab, " *\\* *")[[1]]); a=tmp[1]; b=tmp[2]            
      f= update(form.0, as.formula(paste0("~.+scale(", a ,") * scale(", b, ")")))
      
      fit=svycoxph(f, design=design.1)

      # inner.id 1: treat a as the x axis variable, 2: treat b as the x axis variable
      # only do 2 because 2 is the BD29 or fold change variables
      for (inner.id in 2:2) {
        if (inner.id == 1) {vx=a; vthree=b} else {vx=b; vthree=a}        
        
        # compute risks at three values of vthree
        three.val=wtd.quantile(dat.ph1[[vthree]], dat.ph1$wt, c(.15, .5, .85))                    
        # compute risks at a sequence of vx values for each of the three vthree values
        ss=sort(c(
          wtd.quantile(dat.ph1[[vx]], dat.ph1$wt, c(0.025,0.05,0.95,0.975)), # will be included in the table
          seq(min(dat.ph1[[vx]], na.rm=TRUE), max(dat.ph1[[vx]], na.rm=TRUE), length=100) # equally spaced between min and max so that the curves look good
        ))    
        
        # estimate marginalized risks, return a matrix
        prob.ls=sapply (three.val, function(val) {
          marginalized.risk.cont.2(fit, marker.name  =vx, data=dat.ph2, weights=dat.ph2$wt, t=tfinal.tpeak, ss=ss, 
                                   marker.name.2=vthree, s.2=val)
        })
        
        #### bootstrap
        # store the current rng state
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }         
        
        seeds=1:B; names(seeds)=seeds
        out=mclapply(seeds, mc.cores = numCores, FUN=function(seed) {   
          seed=seed+560
          if (verbose>=2) myprint(seed)
          
          dat.b = bootstrap.cove.boost.2(dat.ph1, seed)
          dat.b.ph2=subset(dat.b, ph2==1)     
          with(dat.b, table(Wstratum, ph2))     
          
          # inline design object b/c it may also throw an error
          fit.b=try(svycoxph(f, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
          
          if (!inherits (fit.b, "try-error")) {
            probs=sapply (three.val, function(val) {
              marginalized.risk.cont.2(fit.b, marker.name  =vx, data=dat.b.ph2, weights=dat.b.ph2$wt, t=tfinal.tpeak, ss=ss, 
                                       marker.name.2=vthree, s.2=val)
            })
          } else {
            matrix(NA, length(ss), length(three.val))
          }
        })
        
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)    
        
        # organize bootstrap results into a list of 3, each element of which is a matrix of n.s by n.seeds
        res.ls=lapply (1:length(three.val), function(i) {
          res=sapply(out, function (x) x[,i])
          res[,!is.na(res[1,])] # remove NA's
        })
        if (verbose) str(res.ls)
        # put lb and ub into matrices
        lb.ls=sapply(res.ls, function (res) t(apply(res, 1, function(x) quantile(x, c(.025)))) )
        ub.ls=sapply(res.ls, function (res) t(apply(res, 1, function(x) quantile(x, c(.975)))) )
        
        risks.itxn[[paste0(vx,vthree)]]=list(marker=ss, prob=prob.ls, boot=res.ls, lb=lb.ls, ub=ub.ls, marker.2=three.val)
      } # end inner.id
    }
    
    save(risks.itxn, file=paste0(save.results.to, "itxn.marginalized.risk.Rdata"))
    
  } else {
    load(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))
  }  

  
  ####
  cat("make plot interaction risk curves\n")
  
  for (ab in config$interaction) {
    tmp=trim(strsplit(ab, " *\\* *")[[1]]); a=tmp[1]; b=tmp[2]            
    
    for (inner.id in 2:2) {
      if (inner.id == 1) {vx=a; vthree=b} else {vx=b; vthree=a}        
      vx.assayname = sub("BD[[0123456789]+", "", sub("DeltaBD[[0123456789]+overBD1", "", vx))
      vthree.assayname = sub("BD[[0123456789]+", "", sub("DeltaBD[[0123456789]+overBD1", "", vthree))
      
      risks=risks.itxn[[paste0(vx,vthree)]]
      
      mypdf(oma=c(0,0,0,0), onefile=F, file=paste0(save.results.to, "itxn_marginalized_risks_",ifelse(inner.id==1,a,b),"_",ifelse(inner.id==1,b,a)), mfrow=.mfrow)
      
      par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
      lwd=2
      
      shown=risks$marker>=wtd.quantile(dat.ph1[[vx]], dat.ph1$wt, 2.5/100) & 
            risks$marker<=wtd.quantile(dat.ph1[[vx]], dat.ph1$wt, 1-2.5/100)
      
      # hard code ylim if needed to make the plot look better
      ylim=range(risks$prob[shown,], 0) # [shown] so that there is not too much empty space
      xlim=range(risks$marker[shown]) 
      if(verbose) myprint(xlim, ylim)
      
      # set up an empty plot
      plot(risks$marker[shown], risks$prob[shown,1], 
           xlab=paste0("BD29/BD1", assay_labels_short[vx.assayname], " (=s)"), 
           ylab=paste0("Probability* of ",config$txt.endpoint," by Day ", tfinal.tpeak), 
           lwd=lwd, xlim=xlim, ylim=ylim, type="n", main="", xaxt="n")    
      draw.x.axis.cor(xlim, NA, NA) # no lloq for fold change
      
      # draw risk lines and confidence bands
      
      for (i in 1:length(risks$marker.2)) {
        lines(risks$marker[shown], risks$prob[shown,i], lwd=lwd, col=i, lty=1)# use dashed line for the middle so that overlaps can be seen
        col <- c(col2rgb(c("black","red","green")[i]))
        col <- rgb(col[1], col[2], col[3], alpha=255*0.2, maxColorValue=255)
        lines(risks$marker[shown], risks$lb[shown,i],   lwd=lwd, col=col, lty=1)
        lines(risks$marker[shown], risks$ub[shown,i],   lwd=lwd, col=col, lty=1)    
      }
      
      # legend for the three lines
      legend.txt=c("(15th percentile)","(median)","(85th percentile)")
      mylegend(x=3, legend=paste(signif(10**risks$marker.2, 3), legend.txt), col=1:3, lty=c(1,2,1), title=paste0("BD1 ", assay_labels_short[vthree.assayname]), lwd=lwd)
      
      # # placebo prevalance lines
      # abline(h=prev.plac[1], col="gray", lty=c(1,3,3), lwd=lwd)
      # text(x=par("usr")[2]-diff(par("usr")[1:2])/5, y=prev.plac[1]+diff(par("usr")[3:4])/30, "placebo arm "%.%formatDouble(prev.plac[1],3,remove.leading0=F))        

      # add histogram
      par(new=TRUE) 
      col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
      col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
      tmp.x=dat.ph1[[vx]][dat.ph1$ph2]
      tmp.w=dat.ph1$wt[dat.ph1$ph2]
      tmp=get.marker.histogram(tmp.x, tmp.w, attr(config,"config"))
      if (is.nan(tmp$density)[1]) tmp=hist(tmp.x, plot=F)
      plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F, xlim=xlim, ylim=c(0,max(tmp$density*1.25)))
      
      dev.off()            
    }
  }  
  
  
} # end if iNaive


# # joint distribution of BD1 and BD29 markers
# par(mfrow=c(1,2))
# with(dat.naive, plot(BD1bindSpike_BA.1, BD29bindSpike_BA.1, col=ifelse(EventIndPrimary, 2, 1),  main=paste0("Naive, cor ",round(cor(BD1bindSpike_BA.1, BD29bindSpike_BA.1, use="complete.obs"),2))))
# abline(0,1)
# with(dat.nnaive, plot(BD1bindSpike_BA.1, BD29bindSpike_BA.1, col=ifelse(EventIndPrimary, 2, 1), main=paste0("Naive, cor ",round(cor(BD1bindSpike_BA.1, BD29bindSpike_BA.1, use="complete.obs"),2))))
# abline(0,1)



###################################################################################################
# multivariate_assays models in naive + nnaive

if (!is.null(config$multivariate_assays)) {
  if(verbose) print("Multiple regression")
  
  dat.ph1 = subset(dat, ph1.BD29)
  design.1<-twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.ph1)
  with(dat.ph1, table(Wstratum, ph2, useNA="ifany"))
  
  for (ind in 1:2) { # repeat twice, with and without naive::HighRiskInd
  for (a in config$multivariate_assays) {
  for (i in 1:2) { # 1: per SD; 2: per 10-fold
      
      a.tmp=a
      if (i==1) {
        # it is key not to trim aa and it is key that a is, e.g. x1 + x2 and not x1+x2, otherwise it won't work
        aa=strsplit(a, "\\+")[[1]]
        for (x in aa[!contain(aa, "\\*")]) {
          # replace every variable with scale(x) when i==1
          a.tmp=gsub(x, paste0(if(i==1) "scale","(",x,")"), a.tmp) 
        }
      }
      
      if (ind==1) {
        f= update(form.0, as.formula(paste0("~.+naive+", a.tmp)))
      } else{
        f= update(form.0, as.formula(paste0("~.+naive*HighRiskInd+", a.tmp)))
      }
      fit=svycoxph(f, design=design.1) 
      var.ind = which(sapply(names(coef(fit)), function(x) any(sapply(aa, function (a) contain(x,kyotil::trim(a))))))
      
      # svycoxph(Surv(EventTimePrimary, EventIndPrimary) ~ MinorityInd + HighRiskInd + 
      #            risk_score + scale(BD1pseudoneutid50_BA.1) * scale(BD29pseudoneutid50_BA.1), design=design.1)
      # 
      # svycoxph(Surv(EventTimePrimary, EventIndPrimary) ~ MinorityInd + HighRiskInd + 
      #            risk_score + scale(BD1bindSpike_BA.1) * scale(BD29bindSpike_BA.1), design=design.1)
      
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
      
      mytex(tab, file.name=paste0("CoR_multivariable_svycoxph_pretty", match(a, config$multivariate_assays), if(i==2) "_per10fold", "_", ind), align="c", include.colnames = T, save2input.only=T, 
            input.foldername=save.results.to.0)
  }
  }
  }
}



print(date())
print("cor_coxph run time: "%.%format(Sys.time()-time.start, digits=1))
