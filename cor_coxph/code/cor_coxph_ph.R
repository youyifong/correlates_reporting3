###################################################################################################
if(verbose) print("Regression for continuous markers")

# Report point and 95\% confidence interval estimates for the hazard ratio per 10-fold change in the antibody marker, 
# for the entire baseline negative vaccine cohort


fits=list()
for (a in all.markers) {
    f= update(form.0, as.formula(paste0("~.+", a)))
    fits[[a]]=svycoxph(f, design=design.1) 
}
# scaled marker
fits.scaled=list()
for (a in all.markers) {
    f= update(form.0, as.formula(paste0("~.+scale(", a, ")")))
    fits.scaled[[a]]=svycoxph(f, design=design.1) 
}

# put coxph model coef together to save
fits.cont.coef.ls = lapply(fits, function (fit) getFixedEf(fit, robust=T))

natrisk=nrow(dat.ph1)
nevents=sum(dat.ph1$yy==1)

# make pretty table
rows=length(coef(fits[[1]]))
est=getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=1)
ci= getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=13)
p=  getFormattedSummary(fits, exp=T, robust=T, rows=rows, type=10)
est.scaled=getFormattedSummary(fits.scaled, exp=T, robust=T, rows=rows, type=1)
ci.scaled= getFormattedSummary(fits.scaled, exp=T, robust=T, rows=rows, type=13)

pvals.cont = sapply(fits, function(x) {
    tmp=getFixedEf(x)
    p.val.col=which(startsWith(tolower(colnames(tmp)),"p"))
    tmp[nrow(tmp),p.val.col]
})




###################################################################################################
if(verbose) print("regression for trichotomized markers")

fits.tri=list()
for (a in all.markers) {
    if(verbose) myprint(a)
    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
    fits.tri[[a]]=run.svycoxph(f, design=design.1) 
}
fits.tri=fits.tri

fits.tri.coef.ls= lapply(fits.tri, function (fit) getFixedEf(fit, robust=T))


rows=length(coef(fits.tri[[1]]))-1:0
# get generalized Wald p values
overall.p.tri=sapply(fits.tri, function(fit) {
    if (length(fit)==1) NA else {
        stat=coef(fit)[rows] %*% solve(vcov(fit,robust=T)[rows,rows]) %*% coef(fit)[rows]
        pchisq(stat, length(rows), lower.tail = FALSE)
    }
})
#
overall.p.0=formatDouble(c(rbind(overall.p.tri, NA,NA)), digits=3, remove.leading0 = F);   overall.p.0=sub("0.000","<0.001",overall.p.0)




###################################################################################################
if(verbose) print("# multitesting adjustment for continuous and trichotomized markers together")

# If primary_assays is not defined in config, multitesting adjustment is over all assays. 
# If primary_assays defined, multitesting adjustment is only over this subset. If this set is empty, then no multitesting adjustment is done

p.unadj=c(cont=pvals.cont, tri=overall.p.tri)
# save a copy for later use
p.unadj.1 = p.unadj 
# pick out a subset based on config
if (!is.null(config$primary_assays)) {
    if (length(config$primary_assays)>0) {
        p.unadj = p.unadj[c("cont.Day"%.%tpeak%.%primary_assays, "tri.Day"%.%tpeak%.%primary_assays)]
    } else {
        p.unadj=c()
    }
}

# Holm and FDR adjustment
pvals.adj.fdr=p.adjust(p.unadj, method="fdr")
pvals.adj.hol=p.adjust(p.unadj, method="holm")

if (length(p.unadj)>1) {
        
    #### Westfall and Young permutation-based adjustment
    perm.file.name=paste0(save.results.to,"pvals.perm.Rdata")
    if(!file.exists(perm.file.name)) {
        
        dat.ph2 = design.1$phase1$sample$variables
        design.1.perm=design.1
        #design.1.perm$phase1$full$variables
    
    #    # if want to only do multitesting when liveneutmn50 is included
    #    if (!"liveneutmn50" %in% assays) numPerm=5
        
        # TODO: there is no need to permutate all.markers
        out=mclapply(1:numPerm, mc.cores = numCores, FUN=function(seed) {   
            # store the current rng state 
            save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
            if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }          
            set.seed(seed)        
            
            # permute markers in design.1.perm
            new.idx=sample(1:nrow(dat.ph2))
            tmp=dat.ph2
            for (a in all.markers) {
                tmp[[a]]=tmp[[a]][new.idx]
                tmp[[a%.%"cat"]]=tmp[[a%.%"cat"]][new.idx]
            }
            design.1.perm$phase1$sample$variables = tmp
            
            # rename all.markers so that out has the proper names. this is only effective within permutation
            names(all.markers)=all.markers
            out=c(
                cont=sapply (all.markers, function(a) {
                    f= update(form.0, as.formula(paste0("~.+", a)))
                    fit=run.svycoxph(f, design=design.1.perm) 
                    if (length(fit)==1) NA else last(c(getFixedEf(fit)))
                })        
                ,    
                tri=sapply (all.markers, function(a) {
                    f= update(form.0, as.formula(paste0("~.+", a, "cat")))
                    fit=run.svycoxph(f, design=design.1.perm) 
                    if (length(fit)==1) NA else last(c(getFixedEf(fit)))
                })
            )
            
            # restore rng state 
            assign(".Random.seed", save.seed, .GlobalEnv)    
            
            out
        })
        pvals.perm=do.call(rbind, out)
        save(pvals.perm, file=perm.file.name)
        
    } else {
        load(file=perm.file.name)
    }
    # save number of permutation replicates
    write(nrow(pvals.perm), file=paste0(save.results.to, "permutation_replicates"))
    
    
    if(any(is.na(p.unadj))) {
        pvals.adj = cbind(p.unadj=p.unadj, p.FWER=NA, p.FDR=NA)
    } else {
        #print(colnames(pvals.perm))
        #print(p.unadj)
        pvals.adj = p.adj.perm (p.unadj, pvals.perm[,names(p.unadj)], alpha=1)  
    }
    if(verbose) print(pvals.adj)
    
} else {
    print("not doing Westfall and Young")
    pvals.adj=cbind(p.unadj, p.FWER=pvals.adj.hol, p.FDR=pvals.adj.fdr)
    write(NA, file=paste0(save.results.to, "permutation_replicates"))     # so the rmd file can compile
}



# since we take ID80 out earlier, we may need to add it back for the table and we do it with the help of p.unadj.1
pvals.adj = cbind(p.unadj=p.unadj.1, pvals.adj[match(names(p.unadj.1), rownames(pvals.adj)),2:3, drop=F])


###################################################################################################
# make continuous markers table

p.1=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F); p.1=sub("0.000","<0.001",p.1)
p.2=formatDouble(pvals.adj["cont."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F); p.2=sub("0.000","<0.001",p.2)


tab.1=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p), p.2, p.1)
# rownames(tab.1)=all.markers.names.short
tab.1
mytex(tab.1, file.name="CoR_univariable_svycoxph_pretty", align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", study_name, "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per 10-fold incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    "),
    longtable=T, 
    label=paste0("tab:CoR_univariable_svycoxph_pretty"), 
    caption.placement = "top", 
    caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", config$txt.endpoint, " in the vaccine group: Hazard ratios per 10-fold increment in the marker*")
)
tab.cont=tab.1

tab.1.nop12=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est), t(ci), t(p))
# rownames(tab.1.nop12)=all.markers.names.short

# scaled markers
tab.1.scaled=cbind(paste0(nevents, "/", format(natrisk, big.mark=",")), t(est.scaled), t(ci.scaled), t(p), p.2, p.1)
# rownames(tab.1.scaled)=all.markers.names.short
tab.1.scaled
mytex(tab.1.scaled, file.name="CoR_univariable_svycoxph_pretty_scaled", align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", study_name, "} & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{2}{c}{HR per SD incr.}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{q-value}   & \\multicolumn{1}{c}{FWER} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{***} & \\multicolumn{1}{c}{} \\\\ 
         \\hline\n 
    "),
    longtable=T, 
    label=paste0("tab:CoR_univariable_svycoxph_pretty_scaled"), 
    caption.placement = "top", 
    caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", config$txt.endpoint, " in the vaccine group: Hazard ratios per SD increment in the marker*")
)
tab.cont.scaled=tab.1.scaled



###################################################################################################
# make trichotomized markers table

overall.p.1=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FWER"], 3, remove.leading0=F);   overall.p.1=sub("0.000","<0.001",overall.p.1)
overall.p.2=formatDouble(pvals.adj["tri."%.%names(pvals.cont),"p.FDR" ], 3, remove.leading0=F);   overall.p.2=sub("0.000","<0.001",overall.p.2)

# add space
overall.p.1=c(rbind(overall.p.1, NA,NA))
overall.p.2=c(rbind(overall.p.2, NA,NA))


# if "Delta"%.%tpeak%.%"overB" is included, nevents have a problem because some markers may have only two category in the cases

# n cases and n at risk
natrisk = round(c(sapply (all.markers%.%"cat", function(a) aggregate(subset(dat.ph1,ph2==1)        [["wt"]], subset(dat.ph1,ph2==1        )[a], sum, na.rm=T, drop=F)[,2] )))
nevents = round(c(sapply (all.markers%.%"cat", function(a) aggregate(subset(dat.ph1,yy==1 & ph2==1)[["wt"]], subset(dat.ph1,yy==1 & ph2==1)[a], sum, na.rm=T, drop=F)[,2] )))
natrisk[is.na(natrisk)]=0
nevents[is.na(nevents)]=0
colSums(matrix(natrisk, nrow=3))
# regression parameters
est=c(rbind(1.00,  sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=1))  ))
ci= c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=13)) ))
p=  c(rbind("N/A", sapply(fits.tri, function (fit) if(length(fit)==1) rep(NA,2) else getFormattedSummary(list(fit), exp=T, robust=T, rows=rows, type=10)) ))

tab=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0, overall.p.2, overall.p.1
)
# rownames(tab)=c(rbind(all.markers.names.short, "", ""))
tab
tab.cat=tab[1:(nrow(tab)),]
#cond.plac=dat.pla.seroneg[[config$EventTimePrimary]]<=tfinal.tpeak # not used anymore

# use longtable because this table could be long, e.g. in hvtn705second
mytex(tab[1:(nrow(tab)),], file.name="CoR_univariable_svycoxph_cat_pretty", align="c", include.colnames = F, save2input.only=T, input.foldername=save.results.to,
    col.headers=paste0("\\hline\n 
         \\multicolumn{1}{l}{", study_name, "} & \\multicolumn{1}{c}{Tertile}   & \\multicolumn{1}{c}{No. cases /}   & \\multicolumn{1}{c}{Attack}   & \\multicolumn{2}{c}{Haz. Ratio}                     & \\multicolumn{1}{c}{P-value}   & \\multicolumn{1}{c}{Overall P-}      & \\multicolumn{1}{c}{Overall q-}   & \\multicolumn{1}{c}{Overall} \\\\ 
         \\multicolumn{1}{l}{Immunologic Marker}            & \\multicolumn{1}{c}{}          & \\multicolumn{1}{c}{No. at-risk**} & \\multicolumn{1}{c}{rate}   & \\multicolumn{1}{c}{Pt. Est.} & \\multicolumn{1}{c}{95\\% CI} & \\multicolumn{1}{c}{(2-sided)} & \\multicolumn{1}{c}{value***} & \\multicolumn{1}{c}{value $\\dagger$} & \\multicolumn{1}{c}{FWER} \\\\ 
         \\hline\n 
    "),        
    # add.to.row=list(list(nrow(tab)), # insert at the beginning of table, and at the end of, say, the first table
    #     c(paste0(" \n \\multicolumn{8}{l}{} \\\\ \n", 
    #               "\n \\multicolumn{2}{l}{Placebo} & ", 
    #              paste0(sum(dat.pla.seroneg$yy), "/", format(nrow(dat.pla.seroneg), big.mark=",")), "&",  
    #              formatDouble(sum(dat.pla.seroneg$yy)/nrow(dat.pla.seroneg), digit=4, remove.leading0=F), "&",  
    #              "\\multicolumn{4}{l}{}  \\\\ \n")
    #       #"\\hline\n \\multicolumn{4}{l}{Standard Deviation 1 mcg/mL}\\\\ \n"
    #      )
    # ),
    longtable=T, 
    label=paste0("tab:CoR_univariable_svycoxph_cat_pretty"), 
    caption.placement = "top", 
    caption=paste0("Inference for Day ", tpeak, "antibody marker covariate-adjusted correlates of risk of ", config$txt.endpoint, " in the vaccine group: Hazard ratios for Middle vs. Upper tertile vs. Lower tertile*")
)




tab.nop12=cbind(
    rep(c("Lower","Middle","Upper"), length(p)/3), 
    paste0(nevents, "/", format(natrisk, big.mark=",",digit=0, scientific=F)), 
    formatDouble(nevents/natrisk, digit=4, remove.leading0=F),
    est, ci, p, overall.p.0
)
# rownames(tab.nop12)=c(rbind(all.markers.names.short, "", ""))




###################################################################################################
# multivariate_assays models

save(fits.cont.coef.ls, fits.tri.coef.ls, file=paste0(save.results.to, "coxph_fits.Rdata"))

save (tab.cont, tab.cat, tab.cont.scaled, pvals.adj, file=paste0(save.results.to, "coxph_slopes.Rdata"))
