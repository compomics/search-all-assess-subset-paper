library(saas)
library(tidyverse)
library(stringr)
library(xtable)
library(cowplot)
library(Biostrings)
### path to functions.r file
source('../functions.r')
### path to pyrococcus folder
path = './'

## make a fasta file for every GO subset
########################################
### make dataframe for every protein and their GO term according uniprot
dt = read.delim(paste0(path,'uniprot_goterms.txt')) %>% as_data_frame

gosid = strsplit(as.character(dt$Gene.ontology.IDs),'\\; ')
gosid_n = data_frame(go = unlist(gosid)) %>% count(go) %>% arrange(desc(n))
gosdes = strsplit(as.character(dt$Gene.ontology..GO.),'\\; ')
gosdes_n = data_frame(go_description= unlist(gosdes)) %>%
  count(go_description) %>% arrange(desc(n))

## make fasta for each GO terms (with at least  15 proteins)
dir.create(path, 'fastas/')

id= gosdes_n %>% filter(n > 14)

sapply(1:nrow(id),function(n){
  i = gosdes_n$go_description[n]

  dt_f = filter(dt,grepl(i,Gene.ontology..GO.,fixed = TRUE))

  fast = AAStringSet(dt_f$Sequence)
  names(fast) = paste0(dt_f$Entry,'|',dt_f$Protein.names)
  i = str_replace_all(i,' ','_')
  fasname = paste0(path, 'fastas/', i, '.fasta')
  writeXStringSet(fast,fasname)
})

fast = AAStringSet(dt$Sequence)
names(fast) = paste0(dt$Entry,'|',dt$Protein.names)
fasname = paste0(path, 'fastas/', 'complete.fasta')
writeXStringSet(fast,fasname)

## MSGF+ search
###############
### msgf+ search parameters
param =' -t 10ppm -ti 0,1 -tda 1 -m 3 -inst 1 -e 1 -protocol 0 -ntt 2 -minLength 6 -maxLength 30 -minCharge 2 -maxCharge 4 -n 1'

### msgf+ search
do_msgf_runs(param,paste0(path,'fastas'),paste0(path,'/mgf'),out = paste0(path,'msgfoutput'),
             modfile = paste0(path,'modification.txt'))

## preprocessing data
#####################
name_complete = "complete.fasta"
name_subsets = list.files(paste0(path,'fastas'),pattern = 'fasta')
name_subsets = name_subsets[name_subsets != "complete.fasta"]

dat_all = parse_msgf_mzid(paste0(path,"msgfoutput/complete.fasta.mzid"))

dat = lapply(name_subsets, function(name){
  cat(name,'\n')
  ## remove specids that are decoy and target in complete or subset search
  ## Collapse spectra that are assigned to multiple proteins to 1 entry in dataframe.
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

## Print as latex table
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

  jpeg(paste0(hist_path,"hist_",subset_name,".jpg"),width = 3750,height = 2500,res = 300, quality = 100,units = 'px')
  makeplots(df)
  dev.off()
})

## Diagnostic plots for all-sub method
######################################
load(file = paste0(path, 'preprocessing.rdata'))

plot_path = paste0(path, 'diagnostic_plots_all_sub/')
dir.create(plot_path)
lapply(dat,function(df){
  d = filter(df, database == 'complete', !(!decoy & !subset)) %>%
    mutate(score = -score)
  p = plot_diag(d)
  cowplot::save_plot(file = paste0(plot_path, last(unique(df$database)),'.pdf'),p$all,
                    ,base_width = 12, base_height = 8)
})

## Boxplot of fractions swapped PSMs in all-sub strategy
########################################################
load(file = paste0(path, 'preprocessing.rdata'))

plot.swapped.box = function(dat,fdr){
  df = purrr::map_df(dat,function(d){
    cat(last(d$database),'\n')
    d = group_by(d,database) %>%
      arrange(score) %>%
      mutate(FDR = ifelse((subset == 1) & (decoy == 0),
                          cummax(cumsum(decoy == 1) /
                                 cumsum(decoy != 1)), NA)) %>%
      ungroup
    filter(d,!((database != 'complete') & (decoy |( FDR > fdr)))) %>%
      group_by(spec_id) %>% filter(n() == 2) %>%
      summarise(swapped = (length(unique(sequence)) == 2)) %>%
      summarise(swapped_frac = mean(swapped))
  })

  boxplot(df$swapped_frac,horizontal = TRUE,main = paste0("Returned at ",round(fdr*100,1),'% FDR'),
          outline = FALSE,xlab = "% swapped PSM's")
  points(df$swapped_frac,jitter(rep(1,nrow(df)),amount = .1),pch = 19)
  abline(v = fdr, col = 'grey')
}

jpeg(paste0(path,"boxplot_swapped_labels_all_sub.jpg"),width = 3750,height = 2500,res = 300, quality = 100,units = 'px')
par(mfrow = c(2,1))
plot.swapped.box(dat,fdr = .01)
mtext('a',adj = 0,line = 2,cex = 1.5)
plot.swapped.box(dat,fdr = .05)
mtext('b',adj = 0,line = 2,cex = 1.5)
dev.off()
