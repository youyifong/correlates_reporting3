## these results are close to bootstrap results. they are not used later and only for sanity check
## compute overall risk regardless of markers in both arms by integrating over form.0. 
## the point estimate matches the results from bootstrap
## the variance is asymptotic and still needs to be figured out
#prevs=sapply (c(placebo=0, vaccine=1), function(i) {
#    dat.ph1=subset(dat.mock, Trt==i & Bserostatus==0 & ph1)
#    fit.tmp = coxph(form.0, dat.ph1, model=T) # model=T to make predict possible
#    dat.ph1[[config.cor$EventTimePrimary]]=tfinal.tpeak
#    pred.tmp=predict(fit.tmp, newdata=dat.ph1, type="expected", se.fit=T)    
#    sd.tmp=exp(mean(log(pred.tmp$se.fit)))
#    prev=c(est=NA, "2.5%"=NA, "97.5%"=NA)
#    prev[1] = mean (1 - exp(-pred.tmp$fit))    
#    #prev[2:3] = prev[1] + c(-1,1)*1.96*sd.tmp
#    prev        
#})
#prevs

if(!file.exists(paste0(save.results.to, "marginalized.risk.no.marker.Rdata"))) {    
    if (verbose) print("bootstrap marginalized.risk.no.marker Rdata")

    prob=get.marginalized.risk.no.marker(form.0, dat.ph1, tfinal.tpeak)
    
    # bootstrapping
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
    
    # prepare for bootstrapping
    if(config$sampling_scheme=="case_cohort") {
      ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (dat.ph1) 
    }
    
    # if mc.cores is >1 here, the process will be stuck in coxph for some unknown reason
    out=mclapply(1:B, mc.cores = 1, FUN=function(seed) {  
      if (verbose>=2) myprint(seed) 
      
      if (TRIAL=="moderna_boost") {
        dat.b = bootstrap.cove.boost.2(dat.ph1, seed)
        
      } else if(config$sampling_scheme=="case_cohort") {
        dat.b = get.bootstrap.data.cor(dat.ph1, ptids.by.stratum, seed)
            
      } else stop("not sure which bootstrap function to use")

      get.marginalized.risk.no.marker(form.0, dat.b, tfinal.tpeak)
    })
    boot=do.call(cbind, out)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    marginalized.risk.no.marker = list(est=c(prob, quantile(boot, c(.025,.975) )), boot=boot)      

    save(marginalized.risk.no.marker, file=paste0(save.results.to, "marginalized.risk.no.marker.Rdata"))
    
} else {
    load(paste0(save.results.to, "marginalized.risk.no.marker.Rdata"))
}
