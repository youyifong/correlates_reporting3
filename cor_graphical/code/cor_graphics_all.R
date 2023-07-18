#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
library(GGally)
library(stringr)
require(devtools)
install_version("dummies", version = "1.5.6", repos = "http://cran.us.r-project.org")
library(cowplot) # for function plot_grid
library(grid)
library(gridExtra)
install.packages("wCorr", repos = "http://cran.us.r-project.org") # weighted correlation
library(wCorr)

# There is a bug on Windows that prevents renv from working properly. The following code provides a workaround:
#if (.Platform$OS.type == "windows") .libPaths(c(paste0(Sys.getenv ("R_HOME"), "/library"), .libPaths()))

source(here::here("code", "cor_graphics_functions.R"))
source(here::here("..", "_common.R"))
#colnames(assay_metadata) = gsub("X[.]+","", colnames(assay_metadata))
# add panel for bindN in assay_metadata
assay_metadata[which(assay_metadata$assay=="bindN"), "panel"] = "bindN"


# for the order of figure panels
assay_order = assay_metadata %>% dplyr::arrange(panel, order_in_panel) %>% select(assay_label_short) %>% pull()
assay_metadata = assay_metadata %>%
    mutate(assay_label_short = factor(assay_label_short,
                                      levels = assay_order
))

dat.longer.cor.subset.plot1 <- readRDS(here::here("data_clean", "longer_cor_data_plot1.rds"))
dat.cor.subset.plot3 <- readRDS(here::here("data_clean", "cor_data_plot3.rds"))
dat.cor.subset.plot3$all_one <- 1 # as a placeholder for strata values

# path for figures and tables etc
save.results.to = here::here("output")
if (!dir.exists(save.results.to))  dir.create(save.results.to)
save.results.to = paste0(here::here("output"), "/", attr(config,"config"),"/");
if (!dir.exists(save.results.to))  dir.create(save.results.to)
#save.results.to = paste0(save.results.to, "/", COR,"/");
#if (!dir.exists(save.results.to))  dir.create(save.results.to)
print(paste0("save.results.to equals ", save.results.to))

###### Set 1 plots: BD1 and BD29 Ab distributions by case/non-case
set1_times <- c("BD1","BD29","DeltaBD29overBD1")
# ID50 614G and BA.1, at BD1, BD29, BD29-BD1
# bindSpike 614G and BA.1, at BD1, BD29, BD29-BD1
# bindRBD 614G and BA.1, at BD1, BD29, BD29-BD1
for (panel in c("pseudoneutid50", "bindSpike","bindRBD")){
    
    f_1 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1,
        assays = paste0(panel,c("","_BA.1")),
        times = set1_times,
        axis.x.text.size = 18,
        strip.x.text.size = 18)
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_2_strain_by_case_non_case_at_", set1_times[i], ".pdf")
        ggsave(plot = f_1[[i]], filename = paste0(save.results.to, file_name), width = 16, height = 16)
    }
}

# bindSpike D614, Gamma, Alpha, Beta, Delta, BA.1, at BD1, BD29, BD29-BD1
for (panel in c("bindSpike")){
    
    f_1 <- f_case_non_case_by_time_assay(
        dat = dat.longer.cor.subset.plot1,
        assays = paste0(panel,c("", "_Gamma", "_Alpha", "_Beta", "_Delta", "_BA.1")),
        times = set1_times,
        axis.x.text.size = 10,
        strip.x.text.size = 10)
    
    for (i in 1:length(set1_times)){
        
        file_name <- paste0(panel, "_6_strain_by_case_non_case_at_", set1_times[i], ".pdf")
        ggsave(plot = f_1[[i]], filename = paste0(save.results.to, file_name), width = 16, height = 16)
    }
}

###### Set 2 plots: Longitudinal plots BD1 to BD29 (and to DD1)
set2_assays = c("bindSpike","bindSpike_BA.1","pseudoneutid50","pseudoneutid50_BA.1")
# bindSpike 614G and BA.1 non-cases and cases, at BD1, BD29, DD1
# ID50 614G and BA.1 non-cases and cases, at BD1, BD29, DD1
f_2 <- f_longitude_by_assay(
    dat = dat.longer.cor.subset.plot1,
    assays = set2_assays,
    times = c("BD1","BD29","DD1")
)

for (i in 1:length(c("bindSpike","pseudoneutid50"))){
    
    panel = c("bindSpike","pseudoneutid50")[i]

    file_name <- paste0(panel, "_longitudinal_by_case_non_case.pdf")
    ggsave(plot = f_2[[i]], filename = paste0(save.results.to, file_name), width = 16, height = 11)
}


###### Set 3 plots: Correlation plots across markers at a given time point
# 9 markers, three timepoints
for (t in c("BD1","BD29","DeltaBD29overBD1")) {
    covid_corr_pairplots(
        plot_dat = dat.cor.subset.plot3,
        time = t,
        assays = assays,
        strata = "all_one",
        weight = "wt.BD29",
        plot_title = paste0(
            "Correlations of 9 ", t, " antibody markers, Corr = Weighted Spearman Rank Correlation."
        ),
        column_labels = paste(t, assay_metadata$assay_label_short),
        height = max(1.3 * length(assays) + 0.1, 5.5),
        width = max(1.3 * length(assays), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata$assay_label_short)))>40, 3.3, 3.8),
        filename = paste0(
            save.results.to, "/pairs_by_time_", t,
            "_9_markers.pdf"
        )
    )


    # 4 markers, three timepoints
    assay_metadata_sub <- subset(assay_metadata, assay %in% c("bindSpike", "bindSpike_BA.1",
                                                              "pseudoneutid50", "pseudoneutid50_BA.1"))
    covid_corr_pairplots(
        plot_dat = dat.cor.subset.plot3,
        time = t,
        assays = assay_metadata_sub$assay,
        strata = "all_one", 
        # currently strata is hard-coded to 1's and hard-coded not being used in the correlation calculation
        # strata-based model is currently commented out in the function, ggally_statistic_resample
        weight = "wt.BD29",
        plot_title = paste0(
            "Correlations of 4 ", t, " antibody markers, Corr = Weighted Spearman Rank Correlation."
        ),
        column_labels = paste(t, assay_metadata_sub$assay_label_short),
        height = max(1.3 * length(assay_metadata_sub$assay) + 0.1, 5.5),
        width = max(1.3 * length(assay_metadata_sub$assay), 5.5),
        column_label_size = ifelse(max(nchar(paste(t, assay_metadata_sub$assay_label_short)))>40, 3.3, 4.3),
        filename = paste0(
            save.results.to, "/pairs_by_time_", t,
            "_4_markers.pdf"
        )
    )
}


###### Set 4 plots: Correlation plots for a given marker across time points
# 9 markers
for (a in assays){
    panels_set <- list()
    i <- 1
    dat.cor.subset.plot3$Trt_nnaive = with(dat.cor.subset.plot3, 
                                           case_when(Trt == 1 & nnaive == 0 ~ "Vaccine Naive", 
                                                     Trt == 1 & nnaive == 1 ~ "Vaccine Non-naive", 
                                                     Trt == 0 & nnaive == 0 ~ "Placebo Naive", 
                                                     Trt == 0 & nnaive == 1 ~ "Placebo Non-naive"))
    for (tn in c("Vaccine Naive", "Vaccine Non-naive", "Placebo Naive", "Placebo Non-naive")){
        for (ce in c("Non-Cases", "Omicron Cases")){
            times_sub = c("BD1","BD29",if(ce=="Omicron Cases") "DD1")
            
            panels_set[[i]] = covid_corr_pairplots(
                plot_dat = dat.cor.subset.plot3 %>% filter(Trt_nnaive == tn & cohort_event == ce),
                time = times_sub,
                assays = a,
                strata = "all_one",
                weight = "wt.BD29",
                plot_title = "",
                column_labels = times_sub,
                height = 5.5,
                width = 5.5,
                column_label_size = 10,
                write_to_file = F
            ) 
            i = i + 1
        }
    }
    
    y.grob.1 <- textGrob("Vaccine     \nNaive", gp=gpar(fontface="bold", col="black", fontsize=9))
    y.grob.2 <- textGrob("Vaccine   \nNon-naive", gp=gpar(fontface="bold", col="black", fontsize=9))
    y.grob.3 <- textGrob("Placebo     \nNaive", gp=gpar(fontface="bold", col="black", fontsize=9))
    y.grob.4 <- textGrob("Placebo   \nNon-naive", gp=gpar(fontface="bold", col="black", fontsize=9))
    
    #x.grob <- textGrob(paste0(
    #    "Correlations of ", a, " Levels at BD1, BD29 (and DD1) by Cases/Non-cases, Corr = Weighted Spearman Rank Correlation."
    #), gp=gpar(fontface="bold", col="black", fontsize=9))
   
    #add to plot
    combined_p <- grid.arrange(
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[1]]), ggmatrix_gtable(panels_set[[2]])), left = y.grob.1), nrow=1),
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[3]]), ggmatrix_gtable(panels_set[[4]])), left = y.grob.2), nrow=1),
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[5]]), ggmatrix_gtable(panels_set[[6]])), left = y.grob.3), nrow=1),
        grid.arrange(arrangeGrob(plot_grid(
            ggmatrix_gtable(panels_set[[7]]), ggmatrix_gtable(panels_set[[8]])), left = y.grob.4), nrow=1),
        #bottom = x.grob,
        ncol = 1
    )
    
    ggsave(filename = paste0(
        save.results.to, "/pairs_by_timepoints_", a, ".pdf"), plot = combined_p, width = 8, height = 10, units="in")
}

