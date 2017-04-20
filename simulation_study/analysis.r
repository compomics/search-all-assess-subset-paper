library(tidyverse)
library(saas)
library(parallel)
### path to simulation_study folder
path = './'

## simulate data under several settings
#######################################
sims.res = function(sim_n= 1, H0_mean, H1_mean, H0_sd, H1_sd, decoy_large_n, pi0, n){
  lapply(1:sim_n,function(i){
    p = saas::simulate_subset(n,pi0)
    fdrs = saas::sample_dataset(H1_n = p$H1_n, H0_n = p$H0_n, decoy_n = p$decoy_n, decoy_large_n = decoy_large_n,
                                H0_mean = H0_mean, H1_mean = H1_mean,H0_sd = H0_sd, H1_sd = H1_sd) %>%
      saas::calculate_fdr()

    fdrs = filter(fdrs,!decoy)
    if (nrow(fdrs) == 0) fdrs = data_frame(FDR = 1, FDR_stable = 1, FDR_BH = 1,H0 = TRUE)
    sumfun = function(d,alpha){
      summarise(
        fdrs,
        FDP = mean(H0[FDR <= alpha]),
        FDP_stable = mean(H0[FDR_stable <= alpha]),
        FDP_BH = mean(H0[FDR_BH <= alpha])
      ) %>% mutate(alpha)
    }
    sum = bind_rows(sumfun(fdrs,.01),sumfun(fdrs,.05))
    bind_cols(bind_rows(p,p),sum)
  }) %>% bind_rows}

set.seed(2017)
res = mclapply(c(10 , 100, 500, 1000
                 ),function(s){
                   cat('n: ',s,'\npi0: ')
                   r = lapply(c(.1 , .3, .5, .7, .9
                                ),function(pi0){
                                  cat( pi0, ' ')
                                  sims.res(sim_n = 10000,
                                           ## parameters chosen based Plasmodium dataset (see supplementary methods)
                                           H0_mean = 2.75278,
                                           H1_mean = 3.31319,
                                           H0_sd = 0.1172843,
                                           H1_sd = 0.2817541,
                                           decoy_large_n = 2000,
                                           pi0 = pi0,
                                           n = s)
                                }) %>% bind_rows
                   cat('\n')
                   r
                 }) %>% bind_rows
## when no subset decoys FDR cannot be calculated, set this FDP to zero
res = mutate(res, missing = is.na(FDP))
res[is.na(res)] = 0

save(res,file =  paste0(path,'simulation.rdata'))

## Make boxplots of the FDP in the several simulated settings
##############################################################
load( paste0(path,'simulation.rdata'))

res_FDP = res %>% filter(!missing) %>%
  gather(method,FDP, starts_with('FDP'))%>%
  mutate(method = gsub('FDP_*','',method),
         method = ifelse(method == '','TDA', method),
         method = ifelse(method == 'stable','ours', method),
         method = ifelse(method == 'BH','BH', method),
         n = as.factor(n),
         method = as.factor(method)
         ) %>%
  mutate( method = factor(method,levels = c('TDA','ours','BH'))) %>%
  transmute(pi0,n,FDR_level = alpha,method,FDP)

p = ggplot(res_FDP, aes(n, FDP, colour = method)) +
  geom_boxplot() +
  stat_summary(fun.y= mean, geom="point",
               shape=18, size=4,show.legend = FALSE, position=position_dodge(width=0.75)) +
  facet_grid(FDR_level~pi0,labeller = label_both)  +
  geom_hline(aes(yintercept=FDR_level)) +
  scale_color_discrete(name = "FDR method") +
  scale_x_discrete(name = "Number of subset PSMs")
p = p +
  theme(strip.text = element_text(size = rel(1.5))
        , axis.title = element_text(size = rel(1.5))
        , axis.text = element_text(size = rel(1.2))
        , legend.text = element_text(size = rel(1.5))
        , legend.title = element_text(size = rel(1.5))
        )
p2 = p + coord_cartesian(ylim =c(0,.1))

ggsave(p,filename = paste0(path,'/FDP.pdf'),height = 8, width = 16)
ggsave(p2,filename = paste0(path,'/FDP_zoom.pdf'),height = 8, width = 16)
