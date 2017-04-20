library(saas)
library(tidyverse)
library(stringr)
library(xtable)
library(cowplot)
### path to functions.r file
source('../functions.r')
### path to plasmodium folder
path = './'

## MSGF+ search
###############
### MSGF output is already available in a zip file in the msgfoutput folder.
###Please unzip if you don't want to run the msgf+ database search yourself.

## MGF file is zipped to save space. Unzip.
unzip(paste0(path,'mgf/','Dd2_ingel_soluble_Br02.mgf.zip'), exdir = paste0(path,'mgf/'))

### msgf+ search parameters
param =' -t 30ppm -ti 0,1 -tda 1 -m 1 -inst 2 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 30 -minCharge 2 -maxCharge 4 -n 1'

### msgf+ search
do_msgf_runs(msgfpluspath = '../MSGFPlus/MSGFPlus.jar',param,paste0(path,'fastas'),paste0(path,'mgf'),
             out = paste0(path,'msgfoutput'), modfile = paste0(path,'modification.txt'))

## preprocessing data
#####################
name_complete = "complete.fasta"
name_subsets = list.files(paste0(path,'fastas'),pattern = 'fasta')
name_subsets = name_subsets[name_subsets != "complete.fasta"]

dat_all = parse_msgf_mzid(paste0(path,"msgfoutput/complete.fasta.mzid"))

dat = lapply(name_subsets, function(name){
  cat(name,'\n')
  ## remove specids that are decoy and target in complete or subset search
  ## Collapse spectra that are assigned to multiple proteins to one entry in the dataframe.
  ## also indicate if PSM matched to subset proteins, non subset proteins or both
  dat_all = saas::preprocess(dat_all,is_subset = paste0(path, 'fastas/', name))

  dat_sub = parse_msgf_mzid(paste0(path,"msgfoutput/",name,".mzid")) %>%
                      saas::preprocess(is_subset = paste0(path, 'fastas/',name))

  bind_rows(dat_all,dat_sub)
})

### remove all spectra from the dataset that got removed in at least one of the subsets
id_remove_all = dat_all$spec_id[! dat_all$spec_id %in% saas::preprocess(dat_all)$spec_id]

id_remove_sub = lapply(name_subsets, function(name){
  cat(name,'\n')
  dat_sub = parse_msgf_mzid(paste0(path,"msgfoutput/",name,".mzid"))
  dat_sub$spec_id[! dat_sub$spec_id %in% saas::preprocess(dat_sub)$spec_id]
}) %>% unlist %>% unique

id_remove = unique(c(id_remove_all, id_remove_sub))

dat =  lapply(dat, function(d){filter(d,!spec_id %in% id_remove)})
save(dat,file = paste0(path,'preprocessing.rdata'))

## make overview table of database search results
#################################################
load(file = paste0(path, 'preprocessing.rdata'))

dt = make_overview_table(dat)

### Print as latex table
print(xtable(dt), include.rownames=FALSE)
sink(file = paste0(path,'overviewtable.txt'))
print(xtable(dt), include.rownames=FALSE)
sink()

### Score histograms according 3 search strategies
##################################################
load(file = paste0(path, 'preprocessing.rdata'))

hist_path = paste0(path,'histogram_PSMs/')
dir.create(hist_path)

lapply(dat,function(df){
  subset_name = filter(df,database != 'complete')$database %>% first
  cat(subset_name, '\n')
### make dataframe of PSMs for the three methods
  df = filter(df, subset, database == 'complete') %>%
    mutate(method = 'all_sub') %>%
    bind_rows(mutate(df, method = ifelse(database == 'complete', 'all_all', 'sub_sub')))

### Calculate classical FDR
  df = group_by(df,method) %>%
    arrange(score) %>%
    mutate(FDR = ifelse((subset == 1) & (decoy == 0),
                        cummax(cumsum(decoy == 1) /
                               cumsum(decoy != 1)), NA)) %>%
    ungroup

  df = mutate(df, score = -log(score))

  jpeg(paste0(hist_path,"hist_",subset_name,"3.jpg"),width = 3750,height = 2500,res = 300, quality = 100,units = 'px')
  makeplots(df)
  dev.off()
})

## diagnostic plots for all-sub method
######################################
load(file = paste0(path, 'preprocessing.rdata'))

plot_path = paste0(path, 'diagnostic_plots_all_sub/')
dir.create(plot_path)
lapply(dat,function(df){
  d = filter(df, database == 'complete', !(!decoy & !subset)) %>%
    mutate(score = -score)
  p = plot_diag(d)
  save_plot(file = paste0(plot_path, last(unique(df$database)),'.pdf'),p$all,
           ,base_width = 12, base_height = 8)
})
