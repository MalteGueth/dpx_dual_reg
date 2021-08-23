# Title     : Analysis of RT
# Objective : Test effect of cue-probe incongrunecy and block effects
# Created by: Jose C. Garcia Alanis
# Created on: 05.03.21
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

# 4) plot effect of probe -----------------------------------------------------
getPacks(c('ggplot2', 'ggbeeswarm', 'see', 'viridis'))

# set up data frame for modelling and plots
corr_mod <- corrects %>%
  arrange(subject) %>%
  mutate(pp_group = factor(pp_group, levels = c('High', 'Low')),
         probe = factor(probe, levels = c('AX', 'AY', 'BX', 'BY')),
         subject = factor(subject, levels = sort(unique(corrects$subject))),
         rt = w_rt * 1000)

mean_rt <- corr_mod %>%
  ungroup() %>%
  group_by(subject, pp_group, probe) %>%
  summarise(mean_rt = mean(rt))

pj <- position_jitter(0.10, seed = 52)
plot_mean_rt <- ggplot(data = mean_rt,
       aes(y = mean_rt,
           x = probe,
           fill = pp_group,
           shape = pp_group,
           color = pp_group,
           group = subject)) +
  geom_violinhalf(aes(group = probe), show.legend = F, alpha = 0.25,
                  color = 'black', fill = 'black', size = 0.2, flip = T) +
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
       plot_mean_rt, width = 6.0, height = 4.5, dpi = 600)

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
       title = 'Anova probe model (RT ~ probe + (1+block|subject/probe))',
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
       title = 'Emmeans probe model (RT ~ probe + (1+block|subject/probe))',
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
       title = 'Contrasts probe model (RT ~ probe + (1+block|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/eff_size_contrasts_probe.html'),
       digits = 3)

# compute contrasts
contr_probes <- contrast(probe_means, 'tukey', adjust = 'fdr')
# save contrasts table
tab_df(as.data.frame(contr_probes),
       title = 'Contrasts probe model (RT ~ probe + (1+block|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/contrasts_probe.html'),
       digits = 3)
# compute contrast CIs
ci_probes <- confint(contr_probes)
# save CIs table
tab_df(as.data.frame(ci_probes),
       title = 'Contrasts probe model (RT ~ probe + (1+block|subject/probe))',
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


# 7) compute proactive control measures ---------------------------------------
# clean up a bit
rm(list = ls()[!(ls() %in% c('corr_mod', 'corrects', 'path_to_data', 'rt_df'))])

# compute total number of trials
total <- rt_df %>%
  mutate(probe_stim = ifelse(nchar(probe) > 1, substr(probe, 2, 2), probe)) %>%
  mutate(probe = paste0(cue, probe_stim)) %>%
  group_by(subject, block, probe) %>%
  summarise(n_trials = sum(!is.na(trial))) %>%
  arrange(subject, block, probe)

# compute number of errors and error rate per condition
errors <- rt_df %>%
  filter(reaction_probes == 'Incorrect') %>%
  group_by(subject, block, probe) %>%
  summarise(n_errors = sum(!is.na(trial))) %>%
  arrange(subject, block, probe) %>%
  left_join(total, ., by = c('subject', 'block', 'probe')) %>%
  mutate(n_errors = ifelse(is.na(n_errors), 0, n_errors)) %>%
  mutate(error_rate = (n_errors + 0.5) / (n_trials + 1))

# compute number of correct responses and correct rate per condition
n_corrects <- rt_df %>%
  filter(reaction_probes == 'Correct') %>%
  mutate(probe = factor(probe)) %>%
  group_by(subject, block, probe) %>%
  mutate(n_corrects = sum(!is.na(trial))) %>%
  summarise(n_corrects = mean(n_corrects)) %>%
  arrange(subject, block, probe) %>%
  left_join(total, ., by = c('subject', 'block', 'probe')) %>%
  mutate(correct_rate = (n_corrects + 0.5) / (n_trials + 1))

# merge error and corrects rates
behavioural_performance <- n_corrects %>%
  left_join(.,  errors, by = c('subject', 'block', 'probe', 'n_trials'))

# calculate the a-cue-bias index
a_bias <-  behavioural_performance %>%
  filter(probe == 'AX' | probe == 'AY') %>%
  group_by(subject) %>%
  mutate(a_bias =
           ifelse(probe == 'AX',
                  0.5 * (qnorm(correct_rate) + qnorm(lead(error_rate))),
                  NA)) %>%
  select(subject, block, a_bias) %>%
  filter(!is.na(a_bias))
# save  file
a_bias_file <- paste0(path_to_data, '/derivatives/results/stats/a_bias.tsv')
write.table(a_bias, file = a_bias_file, sep = '\t', row.names = F)

# calculate the d-prime index
d_prime <-  behavioural_performance %>%
  filter(probe == 'AX' | probe == 'BX') %>%
  group_by(subject) %>%
  mutate(d_prime =
           ifelse(probe == 'AX',
                  qnorm(correct_rate) - qnorm(lead(error_rate)),
                  NA)) %>%
  select(subject, block, d_prime) %>%
  filter(!is.na(d_prime))
# save  file
d_prime_file <- paste0(path_to_data, '/derivatives/results/stats/d_prime.tsv')
write.table(d_prime, file = d_prime_file, sep = '\t', row.names = F)

# calculate the proactive behaviour index
pbi_errors <- behavioural_performance %>%
  group_by(subject, probe) %>%
  filter(probe == 'AY' | probe == 'BX') %>%
  group_by(subject) %>%
  mutate(pbi_errors =
           ifelse(probe == 'AY',
                  (error_rate - lead(error_rate)) / (error_rate + lead(error_rate)),
                  NA)) %>%
  select(subject, block, pbi_errors) %>%
  filter(!is.na(pbi_errors))
# save  file
pbi_errors_file <- paste0(path_to_data,
                          '/derivatives/results/stats/pbi_errors.tsv')
write.table(pbi_errors, file = pbi_errors_file, sep = '\t', row.names = F)

pbi_rt <- corr_mod %>%
  group_by(subject, block, probe) %>%
  filter(probe == 'AY' | probe == 'BX') %>%
  summarise(mean_rt =  mean(rt)) %>%
  mutate(pbi_rt =
           ifelse(probe == 'AY',
                  (mean_rt - lead(mean_rt)) / (mean_rt + lead(mean_rt)),
                  NA)) %>%
  select(subject, block, pbi_rt) %>%
  filter(!is.na(pbi_rt))


 # check distribution of rt
hist(corrects$rt, breaks = 50)
rug(corrects$rt)

# get rid of extreme values, e.g., rt values very close to 0
# using winsorsed scores
getPacks('psych')
corrects <- corrects %>%
  group_by(subject, block, probe) %>%
  mutate(w_rt = winsor(rt, trim = 0.1))

# ** some descriptive statistics **
# mean rt and sd (trial type by block)
corrects %>% group_by(probe, block) %>%
  summarise(m = mean(rt), sd = sd(rt))
# by trial type
corrects %>% group_by(probe) %>%
  summarise(m = mean(rt), sd = sd(rt))
# by block
corrects %>% group_by(block) %>%
  summarise(m = mean(rt), sd = sd(rt))
# average number of trials
corrects %>%
  group_by(subject, cue) %>%
  summarise(n = sum(!is.na(rt))) %>%
  group_by(cue) %>%
  summarise(mn = mean(n), sd = sd(n))

# *** plot distribution of reaction time ***
getPacks(c('ggplot2', 'viridis', 'Hmisc'))

# data for plot
m_rt <- corrects %>%
  group_by(subject, probe, block, pp_group) %>%
  summarise(m = mean(rt))

m_rt <- corrects %>%
  group_by(subject, probe) %>%
  summarise(m = mean(rt))

# add some space between geoms
#pjd <- position_jitterdodge(
#  jitter.width = 0.5,
#  jitter.height = 0,
#  dodge.width = 0,
#  seed = 2)
pn <- position_nudge(x = 0.4)
pd <- position_jitter(0.2)


# create the plot
rt_plot <- ggplot(data = m_rt,
                  aes(x = probe, y = m,
                      fill = probe, shape = probe)) +
  facet_grid(pp_group ~ block) +
  geom_line(aes(group = subject), position = pd, alpha = 0.1, size = 0.4) +
  geom_point(position = pd) +
  scale_shape_manual(values = c(25, 24, 23, 21)) +
  geom_boxplot(position = pn, width = 0.15, alpha = 0.8) +
  scale_fill_viridis(discrete = T, begin = 0.05, end = .95) +
  labs(x = 'Cue-Probe',
       y = 'RT (ms)',
       fill = 'Cue-Probe',
       shape = 'Cue-Probe') +
  scale_y_continuous(limits = c(0.1, 0.7),
                     breaks = seq(0.1, 0.7, 0.1),
                     labels = seq(100, 700, 100)) +
  geom_segment(aes(x = -Inf, y = 0.1, xend = -Inf, yend = 0.7),
               color = 'black', size = rel(0.5), linetype = 1) +
  geom_segment(aes(x = 'AX', y = -Inf, xend = 'BY', yend = -Inf),
               color = 'black', size = rel(0.5), linetype = 1) +
  theme(axis.title.x = element_text(color = 'black', size = 12,
                                    margin = margin(t = 10)),
        axis.title.y= element_text(color = 'black', size = 12,
                                   margin = margin(r = 10)),
        axis.text = element_text(color = 'black', size = 10),
        panel.background = element_rect(fill = 'gray97'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position='bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        panel.spacing = unit(1, "lines")); rt_plot
# save to disk
ggsave(filename = '../data/derivatives/results/figures/rt_distribution.pdf',
       plot = rt_plot, width = 4, height = 5, dpi = 300)



pjn <- position_jitterdodge(
  jitter.width = 0.50,
  jitter.height = 0,
  dodge.width = 0,
  seed = NA
)

pn <- position_nudge(x = -0.2)

probe_change <- ggplot(data = m_rt,
                     aes(x = probe, y = m,
                         fill = probe, shape = probe, color = probe)) +
  facet_wrap(~ block, ncol = 2, scales = 'free_y') +
  labs(x = 'Cue-Probe',
       y = 'RT (ms)',
       fill = 'Cue-Probe',
       shape = 'Cue-Probe') +
  scale_y_continuous(limits = c(0.1, 0.7),
                     breaks = seq(0.1, 0.7, 0.1),
                     labels = seq(100, 700, 100)) +
  geom_segment(aes(x = -Inf, y = 0.1, xend = -Inf, yend = 0.7),
               color = 'black', size = rel(0.5), linetype = 1) +
  geom_segment(aes(x = 'AX', y = -Inf, xend = 'BY', yend = -Inf),
               color = 'black', size = rel(0.5), linetype = 1) +
  scale_fill_viridis(discrete = T, begin = 0.05, end = .95, option = 'C') +
  scale_color_viridis(discrete = T, begin = 0.05, end = .95, option = 'C') +
  scale_shape_manual(values = c(25, 24, 22, 23)) +
  geom_point(position = pjn, color = 'black', alpha = 1.0, size = 1.0) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = pn) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), position = pn) + 
  stat_summary(fun= mean, geom = "point", size = 1, position = pn, color = 'black') + 
  geom_boxplot(position = pn, width = 0.15, alpha = 0.8) +
  theme(plot.title = element_text(hjust = 0, face = 'bold', size = 12, 
                                  margin = margin(b = 10)),
        axis.title.x = element_text(color = 'black', size = 12, face = 'bold',
                                    margin = margin(t = 10)),
        axis.title.y = element_text(color = 'black', size = 12, face = 'bold',
                                    margin = margin(r = 10)),
        axis.text.x = element_text(color = 'black', size = 10,
                                   angle = 0, hjust = 0.5),
        axis.text.y = element_text(color = 'black', size = 10),
        panel.background = element_rect(fill = 'gray98'),
        strip.text = element_text(color = 'black', size = 10, face = 'bold'),
        strip.background = element_blank(),
        legend.title = element_text(color = 'black', size = 10, face = 'bold'),
        legend.text = element_text(color = 'black', size = 10),
        legend.position = 'bottom',
        legend.direction = 'horizontal'); probe_change
ggsave(filename = './tmi_change_jitter_line.pdf',
       plot = thi_change, width = 10, height = 3.5)
