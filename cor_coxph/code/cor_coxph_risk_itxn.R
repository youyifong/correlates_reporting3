# bootstrap and plotting for interaction models

# inner.ids is a subset of (1,2): if (inner.id == 1) {vx=a; vthree=b} else {vx=b; vthree=a}        

################################################################################

if(!file.exists(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))) {    
  cat ("bootstrap itxn model risk curves")
  
  risks.itxn=list()      
  for (ab in config$interaction) {
    myprint(ab)
    tmp=trim(strsplit(ab, " *\\* *")[[1]]); a=tmp[1]; b=tmp[2]            
    f= update(form.0, as.formula(paste0("~.+scale(", a ,") * scale(", b, ")")))
    
    fit=svycoxph(f, design=design.1)
    
    # inner.id 1: treat a as the x axis variable, 2: treat b as the x axis variable
    # only do 2 because 2 is the BD29 or fold change variables
    for (inner.id in inner.ids) {
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

  cat ("load itxn.marginalized.risk.Rdata")
  load(paste0(save.results.to, "itxn.marginalized.risk.Rdata"))
}  

write(ncol(risks.itxn[[1]]$boot[[1]]), file=paste0(save.results.to, "bootstrap_replicates"))


################################################################################
cat("make plot interaction risk curves\n")

for (ab in config$interaction) {
  tmp=trim(strsplit(ab, " *\\* *")[[1]]); a=tmp[1]; b=tmp[2]            
  
  for (inner.id in inner.ids) {
    if (inner.id == 1) {vx=a; vthree=b} else {vx=b; vthree=a}      
    # commented out b/c hard to generalize
    # vx.assayname = sub("BD[[0123456789]+", "", sub("DeltaBD[[0123456789]+overBD1", "", vx))
    # vthree.assayname = sub("BD[[0123456789]+", "", sub("DeltaBD[[0123456789]+overBD1", "", vthree))
    
    risks=risks.itxn[[paste0(vx,vthree)]]
    
    mypdf(oma=c(0,0,0,0), onefile=F, file=paste0(save.results.to, "itxn_marginalized_risks_",ifelse(inner.id==1,a,b),"_",ifelse(inner.id==1,b,a)), mfrow=.mfrow)
    
    par(las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    lwd=2
    
    shown=risks$marker>=wtd.quantile(dat.ph1[[vx]], dat.ph1$wt, 2.5/100) & 
      risks$marker<=wtd.quantile(dat.ph1[[vx]], dat.ph1$wt, 1-2.5/100)
    
    # hard code ylim if needed to make the plot look better
    ylim=range(risks$prob[shown,], 0) # [shown] so that there is not too much empty space
    xlim=range(risks$marker[shown]) 

    # set up an empty plot
    plot(risks$marker[shown], risks$prob[shown,1], 
         xlab=paste0(vx, " (=s)"), 
         # xlab=paste0("BD29/BD1", assay_labels_short[vx.assayname], " (=s)"), 
         ylab=paste0("Probability* of ",config$txt.endpoint," by Day ", tfinal.tpeak), 
         lwd=lwd, xlim=xlim, ylim=ylim, type="n", main="", xaxt="n")    
    draw.x.axis.cor(xlim, NA, NA) # no lloq for fold change
    
    # draw risk lines and confidence bands
    for (i in 1:length(risks$marker.2)) {
      lines(risks$marker[shown], risks$prob[shown,i], lwd=lwd, col=i, lty=1)# use dashed line for the middle so that overlaps can be seen
      col <- c(col2rgb(c("black","red","green")[i]))
      col <- rgb(col[1], col[2], col[3], alpha=255*0.5, maxColorValue=255)
      lines(risks$marker[shown], risks$lb[shown,i],   lwd=lwd, col=col, lty=3)
      lines(risks$marker[shown], risks$ub[shown,i],   lwd=lwd, col=col, lty=3)    
    }
    
    # legend for the three lines
    legend.txt=c("(15th percentile)","(median)","(85th percentile)")
    mylegend(x=3, legend=paste(signif(10**risks$marker.2, 3), legend.txt), col=1:3, lty=c(1,1,1), lwd=lwd, 
             title=vthree)
             # title=paste0("BD1 ", assay_labels_short[vthree.assayname]))

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


if(!file.exists(paste0(save.results.to, "itxn.coef.Rdata"))) {    
  cat ("bootstrap itxn model coef")
  
  coef.itxn=list()      
  for (ab in config$interaction) {
    myprint(ab)
    tmp=trim(strsplit(ab, " *\\* *")[[1]]); a=tmp[1]; b=tmp[2]            
    f= update(form.0, as.formula(paste0("~.+scale(", a ,") * scale(", b, ")")))
    fit=svycoxph(f, design=design.1)
    
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
        coef(fit.b)
      } else {
        rep(NA, length(coef(fit)))
      }
    })
      
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    coef.itxn[[ab]] = do.call(cbind, out)
  }
  
  save(coef.itxn, file=paste0(save.results.to, "itxn.coef.Rdata"))
  
} else {
  
  cat ("load itxn.coef.Rdata")
  load(paste0(save.results.to, "itxn.coef.Rdata"))
}  


apply(coef.itxn[[1]], 1, quantile, c(.025, .095))
