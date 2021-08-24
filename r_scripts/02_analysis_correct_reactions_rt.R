# Title     : Analysis of RT on correct trials
# Objective : Test effect of cue-probe incongrunecy
#             on reaction times after correct reactions
# Created by: Jose C. Garcia Alanis
# Created on: 24.08.21
# R version : 4.0.2 (2020-06-22), Taking Off Again

# --- set working directory ---
# use this if working from r-studio
curr_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirname(curr_path))

# source function for fast loading and installing of packages
source('./r_functions/getPacks.R')
source('./r_functions/spr2.R')

# 1) define the path to behavioral data directory -----------------------------

# get system and user information
host <- Sys.info()

# set default path or ask user for other path
if (grepl('Jo|jo|ma', host['nodename'])) {

  # default path in project structure
  path_to_data <- '../data'

  } else {

  path_to_data <- readline('Please provide path to subject data: ')
  
}

# 2) import in the data ------------------------------------------------------

# this part requires the package 'dplyr'
getPacks('dplyr')

# get RT files names
rt_files <- list.files(path = paste(path_to_data, 'derivatives/results/rt',
                                    sep = '/'),
                       full.names = T)

# read in the files
rt_list <- lapply(rt_files, read.table, sep = '\t', header = T)

# put the in a data frame
rt_df <- bind_rows(rt_list, .id = "column_label")

# recode block variable
rt_df <- rt_df %>%
  mutate(block =  ifelse(block == 'Single', 'Baseline', 'Regulation'))

# read in the ppi data
ppi_data <- read.table(paste0(path_to_data, '/subject_data/ppi_data.tsv'),
                       header = T, sep = '\t')
# remove subject 9 (no group assignment)
ppi_data <- ppi_data %>% 
  filter(!group_pattern == 0) %>%
  mutate(ix = row_number(ix))

# add descriptive group name
ppi_data <- ppi_data %>%
  mutate(pp_group = ifelse(group_pattern == 1, 'Low',
                           ifelse(group_pattern == 2, 'High', NA)))

# merge the two data frames
ppi_scales <- ppi_data %>% 
  select(ix, pp_group, total_pr, si_pr, sp_pr, f_pr)

rt_df <- rt_df %>%
  left_join(., ppi_scales, by = c('subject' = 'ix')) %>%
  select(-c(column_label, group)) %>%
  select(subject, pp_group, block, trial, cue, probe, run, rt,
         reaction_cues, reaction_probes, total_pr:f_pr) %>%
  filter(!is.na(pp_group))

# 3) exploratory analyses correct reactions -----------------------------------

# extract trials with correct responses
corrects <- rt_df %>%
  filter(reaction_cues == 'Correct' & reaction_probes == 'Correct')

# plausibility checks
# e.g., rt min should be > 0.0, max < 0.750
# e.g., only probes AX, AY, BX and BY should be present
summary(corrects); unique(corrects$probe)

# get rid of extreme values, e.g., rt values very close to 0
# using winsorsed scores
getPacks('psych')
corrects <- corrects %>%
  group_by(subject, probe) %>%
  mutate(w_rt = winsor(rt, trim = 0.05))

# set up data frames for modelling and plots
# rt to milliseconds
corr_mod <- corrects %>%
  arrange(subject) %>%
  mutate(pp_group = factor(pp_group, levels = c('High', 'Low')),
         probe = factor(probe, levels = c('AX', 'AY', 'BX', 'BY')),
         subject = factor(subject, levels = sort(unique(corrects$subject))),
         rt = w_rt * 1000)
# means
mean_rt <- corr_mod %>%
  ungroup() %>%
  group_by(subject, pp_group, probe) %>%
  summarise(mean_rt = mean(rt))

# 4) plot effect of probe -----------------------------------------------------
getPacks('ggplot2')

# create aplot to compare the distributions
png('../data/derivatives/results/figures/rt_dist.png',
    height = 3.5, width = 7, units="in", res = 600)
par(mfrow=c(1, 2))
# check distribution of rt
hist(corrects$w_rt,
     main ='5% Trimmed RT', xlab = 'RT',
     breaks = 150, xlim = c(0, 0.8), ylim = c(0, 500))
rug(corrects$w_rt)

# check distribution of rt
hist(corrects$rt,
     main ='Raw RT', xlab = 'RT',
     breaks = 150, xlim = c(0, 0.8), ylim = c(0, 500))
rug(corrects$rt)
dev.off()

# plot mean rt by cue-probe combination
pj <- position_jitter(0.10, seed = 52)
plot_mean_rt <- ggplot(data = mean_rt,
       aes(y = mean_rt,
           x = probe,
           fill = pp_group,
           shape = pp_group,
           color = pp_group,
           group = subject)) +
  geom_boxplot(aes(group = probe), position = position_nudge(x = 0.3),
               width = 0.2, outlier.color = NA, show.legend = F) +
  geom_line(size = 0.4, alpha = 0.3, position = pj) +
  geom_point(color = 'black', size = 1.5, stroke = 0.8, alpha = 1.0, position = pj) +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_manual(values = c('#B40F20FF', '#46ACC8FF')) +
  scale_color_manual(values = c('#B40F20FF', '#46ACC8FF')) +
  scale_y_continuous(limits = c(100.0, 700.0), breaks = seq(100, 700, 100)) +
  geom_segment(aes(x = -Inf, y = 100.0, xend = -Inf, yend = 700.0),
               color = 'black', size = rel(0.5), linetype = 1) +
  geom_segment(aes(x = 'AX', y = -Inf, xend = 'BY', yend = -Inf),
               color = 'black', size = rel(0.5), linetype = 1) +
  labs(title = 'Reaction time by cue-probe combination',
       x = 'Cue-Probe combination', y = 'Mean RT [ms]',
       color = 'PPI Group', fill = 'PPI Group', shape = 'PPI Group') +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = 'gray98', size = 0.5),
        strip.background = element_blank(),
        plot.title = element_text(color = 'black',
                                  size = 14,
                                  face = 'bold',
                                  family='Mukta'),
        axis.title.x = element_text(color = 'black',
                                    size = 14,
                                    face = 'bold',
                                    margin = margin(t = 15),
                                    family='Mukta'),
        axis.title.y = element_text(color = 'black',
                                    size = 14,
                                    face = 'bold',
                                    margin = margin(r = 15),
                                    family='Mukta'),
        axis.text = element_text(color = 'black',
                                 size = 12),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.key = element_blank(),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm')); plot_mean_rt
ggsave('../data/derivatives/results/figures/rt_means_plot.png',
       plot_mean_rt, width = 5.5, height = 4.5, dpi = 600)

# 5) statistical analyses correct reactions -----------------------------------
getPacks('lme4')

# fit models to assess the effect of cue probe combination
mod_rt_0 <- lmer(data = corr_mod,
                 rt ~ (1|subject))

mod_rt_1 <- lmer(data = corr_mod,
                 rt ~  probe + (1|subject),
                 contrasts = list(probe = 'contr.treatment'))

mod_rt_2 <- lmer(data = corr_mod,
                 rt ~  probe + (1|subject/probe),
                 contrasts = list(probe = 'contr.treatment'))

# 5) check model fit ----------------------------------------------------------
getPacks(c('performance', 'partR2'))
test_bf(mod_rt_1, mod_rt_2)
compare_performance(mod_rt_0, mod_rt_1, mod_rt_2, rank = T)

# check plots
plot(compare_performance(mod_rt_1, mod_rt_2, rank = T))
plot(check_distribution(mod_rt_2))

# 6) compute effects based on the best model ----------------------------------
getPacks(c('car', 'sjPlot', 'partR2', 'emmeans'))

# F-values for main effect of cue-probe combination
# (Type III Wald F tests with Kenward-Roger df)
# *** will take sometime to run ***
probes_anova <- car::Anova(mod_rt_2, test = 'F', type = 'III')
summary(mod_rt_2)

# save F-test table
tab_df(as.data.frame(probes_anova),
       title = 'Anova probe model (RT ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/anova_probe.html'),
       digits = 3)

# compute effect sizes Semi-Partial R-squared
# for main effect probe (cf. Nakagawa & Schielzeth, 2013)
partial_R2 <- partR2(mod_rt_2, partbatch = list(probe = c('probe')),
                     R2_type = "marginal", nboot = 10000, CI = 0.95)
# save r2 table
tab_df(as.data.frame(partial_R2$R2),
       title = 'Semi-Partial R-squared (model probe)',
       file = paste0(path_to_data, '/derivatives/results/tables/R2_probes.html'),
       digits = 3)

# compute estimated marginal means
# *** will take sometime to run ***
emm_options(pbkrtest.limit = 14450)
probe_means <- emmeans(mod_rt_2, ~ probe)
# save means table
tab_df(as.data.frame(probe_means),
       title = 'Emmeans probe model (RT ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/emmeans_probe.html'),
       digits = 3)

# model variance estimates and residual degrees of freedom
mod_var <- VarCorr(mod_rt_0)
totSD <- sqrt(sum(as.data.frame(mod_var)$vcov))
edf <- df.residual(mod_rt_0)

# compute effect sizes for the pairwise contrasts
es_probes <- eff_size(probe_means, sigma = totSD, edf = edf); es_probes
# save effect sizes for probe contrast
tab_df(as.data.frame(es_probes),
       title = 'Contrasts probe model (RT ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/eff_size_contrasts_probe.html'),
       digits = 3)

# compute contrasts
contr_probes <- contrast(probe_means, 'tukey', adjust = 'fdr')
# save contrasts table
tab_df(as.data.frame(contr_probes),
       title = 'Contrasts probe model (RT ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/contrasts_probe.html'),
       digits = 3)
# compute contrast CIs
ci_probes <- confint(contr_probes)
# save CIs table
tab_df(as.data.frame(ci_probes),
       title = 'Contrasts probe model (RT ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/cis_contrast_probe.html'),
       digits = 3)

# compute descriptives for the observed values
descriptives_rt <- corr_mod %>%
  group_by(probe) %>%
  summarise(mean(rt), sd(rt))
tab_df(as.data.frame(descriptives_rt),
       title = 'Descripvives RT by Cue-Probe combination.',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/descriptives_probe_rt.html'),
       digits = 3)

# save computed objects
save(list=c("probes_anova", "probe_means"),
     file=paste0(path_to_data, "/derivatives/results/stats/probe_stats.RData"))
