#Sys.setenv(TRIAL = "moderna_boost")
#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))
source(here::here("..", "_common.R"))
#colnames(assay_metadata) = gsub("X[.]+","", colnames(assay_metadata))
source(here::here("code/cor_graphics_functions.R"))
#-----------------------------------------------

times=c("BD1","BD29","DD1","DeltaBD29overBD1","DeltaDD1overBD1")
uloqs=assay_metadata$uloq; names(uloqs)=assays
pos.cutoffs=assay_metadata$pos.cutoff; names(pos.cutoffs)=assays
loqs=assay_metadata$lloq; names(loqs)=assays

#dat$EventIndPrimary=dat$EventIndOmicronBD29
#dat$EventTimePrimary=dat$EventTimeOmicronBD29
dat.cp <- dat

# assign values above the uloq to the uloq
for (a in assays){
  for (t in c("BD1","BD29","DD1")){
    if (paste0(t, a) %in% colnames(dat.cp)){
      dat.cp[, paste0(t, a)] = ifelse(dat.cp[, paste0(t, a)] > log10(uloqs[a]), log10(uloqs[a]), dat.cp[, paste0(t, a)])
    } else {print(paste0(t, a, " doesn't exist"))}
  }
}

# recalculate DeltaBD29overBD1 after the censoring above in order to calculate the response rate for cases at DD1
for (a in assays){
  if (paste0("BD29", a) %in% colnames(dat.cp) & paste0("BD1", a) %in% colnames(dat.cp)){
    dat.cp[, paste0("DeltaBD29overBD1", a)] = dat.cp[, paste0("BD29", a)] - dat.cp[, paste0("BD1", a)]
  } else {print(paste0("BD1", a, "or BD29", a, " doesn't exist"))}
}

# calculate DeltaBD57overBD1 in order to calculate the response rate for cases at DD1
for (a in assays){
  if (paste0("DD1", a) %in% colnames(dat.cp) & paste0("BD1", a) %in% colnames(dat.cp)){
    dat.cp[, paste0("DeltaDD1overBD1", a)] = dat.cp[, paste0("DD1", a)] - dat.cp[, paste0("BD1", a)]
  } else {print(paste0("BD1", a, "or DD1", a, " doesn't exist"))}
}

immuno_vars <- colnames(dat.cp)[grepl("^BD1|^BD29|^DD1|^Delta", colnames(dat.cp))]

# variable selection
dat <- dat.cp %>% select(Ptid, Trt, CalendarBD1Date, CalendarBD1Interval, nnaive,
                         BDPerprotocolIncludeSeroPos, Stage2SamplingInd, 
                      EventTimeOmicronBD1:EventIndOmicronBD29, 
                      EventTimePrimaryOmicronBD1:EventIndPrimaryOmicronBD29, 
                      #EventIndPrimary, EventTimePrimary, 
                      ph2.BD29, wt.BD29, wt.DD1,
                      all_of(immuno_vars))

# create case vs non-case indicators
# Case = COVID-19 endpoint in the interval [≥ 7 days post BD29 and ≥ Dec 1 2021, May 2022 data base lock date]. As described in the appendix the COVID-19 endpoint is documented to be Omicron BA.1 if possible whereas for some non-naive COVID-19 endpoints there was not lineage data available to document the case to be Omicron BA.1.
# Non-case = Did not acquire COVID-19 (of any strain) in the interval [BD1, data base lock date].
if (study_name=="COVEBoost") {
  dat <- dat %>%
    mutate(cohort_event = factor(
      case_when(BDPerprotocolIncludeSeroPos==1 & BD29window==1 & EventIndPrimaryOmicronBD29==1 & EventTimePrimaryOmicronBD29>=7 & Stage2SamplingInd==1 
                ~ "Omicron Cases",
                BDPerprotocolIncludeSeroPos==1 & BD29window==1 & EventIndPrimaryOmicronBD1==0 & Stage2SamplingInd==1  ~ "Non-Cases"),
      levels = c("Omicron Cases", "Non-Cases"))
      )
}
dat <- dat[!is.na(dat$cohort_event),]

write.csv(subset(dat, ph2.BD29==1), file = here::here("data_clean", "cor_data_plot3.csv"), row.names=F)
saveRDS(subset(dat, ph2.BD29==1), file = here::here("data_clean", "cor_data_plot3.rds"))

## arrange the dataset in the long form, expand by assay types
## dat.long.subject_level is the subject level covariates;
## dat.long.assay_value is the long-form time variables that differ by the assay type
dat.long.subject_level <- dat %>%
  replicate(length(assays),., simplify = FALSE) %>%
  bind_rows()

dat.long.assay_value.names <- times
dat.long.assay_value <- as.data.frame(matrix(
  nrow = nrow(dat) * length(assays),
  ncol = length(dat.long.assay_value.names)
))
colnames(dat.long.assay_value) <- dat.long.assay_value.names

for (tt in seq_along(times)) {
  dat_mock_col_names <- paste(times[tt], assays, sep = "")
  dat.long.assay_value[, dat.long.assay_value.names[tt]] <- unlist(lapply(
    # B, Day29, Delta29overB
    dat_mock_col_names,
    # BbindSpike, BbindRBD
    function(nn) {
      if (nn %in% colnames(dat)) {
        dat[, nn]
      } else {
        rep(NA, nrow(dat))
      }
    }
  ))
}

dat.long.assay_value$assay <- rep(assays, each = nrow(dat))
dat.long <- cbind(dat.long.subject_level, dat.long.assay_value)

## change the labels of the factors for plot labels
dat.long$Trt <- factor(dat.long$Trt, levels = c(0, 1), labels = c("Placebo", "Vaccine"))
dat.long$nnaive <- factor(dat.long$nnaive,
  levels = c(0, 1),
  labels = c("Naive", "Non-naive")
)
dat.long$assay <- factor(dat.long$assay, levels = assays, labels = assays)

# add label = LLoQ, uloq values to show in the plot
dat.long$LLoD = with(dat.long, log10(lods[as.character(assay)]))
dat.long$pos.cutoffs = with(dat.long, log10(pos.cutoffs[as.character(assay)]))
dat.long$LLoQ = with(dat.long, log10(loqs[as.character(assay)]))
dat.long$lb = with(dat.long, #ifelse(grepl("bind", assay), 
                   #"Pos.Cut", "LoD")
                   "LoQ")
dat.long$lbval = with(dat.long, #ifelse(grepl("bind", assay), 
                   #pos.cutoffs, LLoD)
                    LLoQ)

dat.long$ULoQ = with(dat.long, log10(uloqs[as.character(assay)]))
dat.long$lb2 = "ULoQ" #with(dat.long, ifelse(grepl("bind", assay), "ULoQ", ""))
dat.long$lbval2 =  dat.long$ULoQ #with(dat.long, ifelse(grepl("bind", assay), ULoQ, -99))

# Here, only filter based on ph2.BD29==1. Filtering by ph2.DD1 will occur downstream,
# since it should only happen for DD1-related figures.
# two timepoints study: ph2.tinterm
dat.long.cor.subset <- dat.long %>%
  dplyr::filter(ph2.BD29==1)

write.csv(dat.long.cor.subset, file = here::here("data_clean", "scatter_rug_data_plot5.csv"), row.names=F)
saveRDS(dat.long.cor.subset, file = here::here("data_clean", "scatter_rug_data_plot5.rds"))

# long to longer format by time
dat.longer.cor.subset <- dat.long.cor.subset %>%
  tidyr::pivot_longer(cols = all_of(times), names_to = "time", values_to = "value")

# phase 2 filters: 
#    include +++, ++- at BD1  for Omicron Cases
#    include +++, ++- at BD29 for Omicron Cases
#    include +++      at DD1  for Omicron Cases
#    include ++ at both BD1 and BD29 for Non-Cases
if(length(times[!grepl("Delta", times)]) > 2) {
    
  dat.longer.cor.subset <- dat.longer.cor.subset %>% 
    filter(!(time == "DD1" & cohort_event!="Omicron Cases"))  # only keep those +++ at DD1
  # table(dat.longer.cor.subset$cohort_event, dat.longer.cor.subset$time, dat.longer.cor.subset$assay)
}

# define response rates
resp <- getResponder(dat.cp, post_times = times[!times %in% c("BD1","DeltaBD29overBD1","DeltaDD1overBD1")], 
               assays=assays, pos.cutoffs = pos.cutoffs)

resp_by_time_assay <- resp[, c("Ptid", colnames(resp)[grepl("Resp", colnames(resp))])] %>%
  tidyr::pivot_longer(!Ptid, names_to = "category", values_to = "response")

# 256 unique ids, three timepoints, 15 assays
dat.longer.cor.subset <- dat.longer.cor.subset %>%
  mutate(category=paste0(time, assay, "Resp")) %>%
  select(Ptid, time, assay, category, Trt, nnaive, cohort_event, value, wt.BD29, wt.DD1,
         lbval,lbval2,
         lb,lb2) %>%
  left_join(resp_by_time_assay, by=c("Ptid", "category"))

# subsets for violin/line plots
#### figure specific data prep
# 1. define response rate:
# 2. make subsample datasets such that the violin plot only shows <= 100 non-case data points

#### for figures 1, 2, 3, 4: case vs non-case, by naive/non-naive, vaccine/placebo, (Day 1), Day 29 Day 57
groupby_vars1=c("Trt", "nnaive", "cohort_event", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1 <- get_desc_by_group(dat.longer.cor.subset, groupby_vars1)
write.csv(dat.longer.cor.subset.plot1, file = here::here("data_clean", "longer_cor_data_plot1.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1, file = here::here("data_clean", "longer_cor_data_plot1.rds"))

#### for figures 1.2, 2.2: case vs non-case, pooling all, by (Day 1), Day 29 Day 57
groupby_vars1.2=c("cohort_event", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.2 <- get_desc_by_group(dat.longer.cor.subset, groupby_vars1.2)
write.csv(dat.longer.cor.subset.plot1.2, file = here::here("data_clean", "longer_cor_data_plot1.2.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.2, file = here::here("data_clean", "longer_cor_data_plot1.2.rds"))

#### for figures 1.3, 2.3: case vs non-case, by naive vs non-naive, (Day 1), Day 29 Day 57
groupby_vars1.3=c("cohort_event", "nnaive", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.3 <- get_desc_by_group(dat.longer.cor.subset, groupby_vars1.3)
write.csv(dat.longer.cor.subset.plot1.3, file = here::here("data_clean", "longer_cor_data_plot1.3.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.3, file = here::here("data_clean", "longer_cor_data_plot1.3.rds"))

#### for figures 2.4: BD1 vs BD29, by naive vs non-naive, case/non-case, (Day 1), Day 29 Day 57
#groupby_vars1.4=c("cohort_event", "nnaive", "time", "assay")

# define response rate
dat.longer.cor.subset.plot1.4 <- dat.longer.cor.subset.plot1.3 
write.csv(dat.longer.cor.subset.plot1.4, file = here::here("data_clean", "longer_cor_data_plot1.4.csv"), row.names=F)
saveRDS(dat.longer.cor.subset.plot1.4, file = here::here("data_clean", "longer_cor_data_plot1.4.rds"))
