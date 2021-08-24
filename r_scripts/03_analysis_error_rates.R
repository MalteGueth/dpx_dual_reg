# Title     : Analysis of error rates
# Objective : Test effect of cue-probe incongrunecy
#             on sujects error rates
# Created by: Jose C. Garcia Alanis
# Created on: 24.08.21
# R version : 4.1.1 (2021-08-10), Kick Things

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

# 3) compute trials numbers and errors ----------------------------------------

# compute total number of trials
total <- rt_df %>%
  mutate(probe_stim =
           ifelse(nchar(probe) > 1, substr(probe, 2, 2), probe)) %>%
  mutate(probe = paste0(cue, probe_stim)) %>%
  group_by(subject, pp_group, probe) %>%
  summarise(n_trials = sum(!is.na(trial))) %>%
  arrange(subject, probe)

# compute number of errors and error rate per condition
errors <- rt_df %>%
  filter(reaction_probes == 'Incorrect') %>%
  group_by(subject, probe) %>%
  summarise(n_errors = sum(!is.na(trial))) %>%
  arrange(subject, probe) %>%
  left_join(total, ., by = c('subject', 'probe')) %>%
  mutate(n_errors = ifelse(is.na(n_errors), 0, n_errors)) %>%
  mutate(error_rate = (n_errors + 0.5) / (n_trials + 1))

# 4) exploratory analyses errors ----------------------------------------------
getPacks(c('ggplot2', 'ggbeeswarm', 'see', 'viridis'))

# plausibility checks
# e.g., rt min should be > 0.0, max < 0.750
# e.g., only probes AX, AY, BX and BY should be present
summary(errors); unique(errors$probe)

# plot error rates
pj <- position_jitter(0.10, seed = 52)
plot_error_rates <- ggplot(data = errors,
       aes(y = error_rate,
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
  scale_y_continuous(limits = c(0.0, 0.8), breaks = seq(0.0, 0.8, 0.1)) +
  geom_segment(aes(x = -Inf, y = 0.0, xend = -Inf, yend = 0.8),
               color = 'black', size = rel(0.5), linetype = 1) +
  geom_segment(aes(x = 'AX', y = -Inf, xend = 'BY', yend = -Inf),
               color = 'black', size = rel(0.5), linetype = 1) +
  labs(title = 'Error rate by cue-probe combination',
       x = 'Cue-Probe combination', y = 'Error rate',
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
        legend.key.width = unit(1, 'cm')); plot_error_rates
ggsave('../data/derivatives/results/figures/error_rates_plot.png',
       plot_error_rates, width = 5.5, height = 4.5, dpi = 600)

# 5) statistical analyses errors ----------------------------------------------
getPacks('lme4')

# fit the model to assess the effect of cue probe combination
# (generalised version of the model fitted to reaction time data, i.e., same
# random structure)
mod_err_1 <- glmer(data = errors,
                   n_errors / n_trials ~ probe + (1|subject/probe),
                   weights = n_trials, family = binomial, nAGQ = 1,
                   contrasts = list(probe = 'contr.treatment'))

# 6) check model fit ----------------------------------------------------------
getPacks(c('car', 'DHARMa'))

# check residuals
qqp(resid(mod_err_1, 'pearson'))
plot(predict(mod_err_1), resid(mod_err_1, 'pearson'))

# test dispersion based on simulation
simulationOutput <- simulateResiduals(fittedModel = mod_err_1)
plot(simulationOutput)
testDispersion(simulationOutput)

# test dispersion based on residuals
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type="pearson")
  Pearson.chisq <- sum(rp ^ 2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = F)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf,p = pval)
}
overdisp_fun(mod_err_1)

# 6) compute effects based on the best model ----------------------------------
getPacks(c('partR2', 'emmeans'))

# Chi-square test for main effect of cue-probe combination
# (Type III Wald chisquare tests)
# *** will take sometime to run ***
probes_anova <- car::Anova(mod_err_1, type = 'III')
summary(mod_err_1)

# save F-test table
tab_df(as.data.frame(probes_anova),
       title = 'Anova probe model (error_rate ~ probe rate + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/anova_probe_errors.html'),
       digits = 3)

# compute effect sizes Semi-Partial R-squared
# for main effect probe (cf. Nakagawa & Schielzeth, 2013)
partial_R2 <- partR2(mod_err_1, nboot = 10000, CI = 0.95)
# save r2 table
tab_df(as.data.frame(partial_R2$R2),
       title = 'Semi-Partial R-squared (model probe)',
       file = paste0(path_to_data, '/derivatives/results/tables/R2_probes.html'),
       digits = 3)

# compute estimated marginal means
# *** will take sometime to run ***
probe_means <- emmeans(mod_err_1, ~ probe, type = 'response')
# save means table
tab_df(as.data.frame(probe_means),
       title = 'Emmeans probe model (error rate ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/emmeans_probe_errors.html'),
       digits = 3)

# model variance estimates and residual degrees of freedom
mod_var <- VarCorr(mod_err_1)
totSD <- sqrt(sum(as.data.frame(mod_var)$vcov))
edf <- df.residual(mod_err_1)

# compute effect sizes for the pairwise contrasts
es_probes <- eff_size(probe_means, sigma = totSD, edf = edf); es_probes
# save effect sizes for probe contrast
tab_df(as.data.frame(es_probes),
       title = 'Contrasts probe model (error rate ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/eff_size_contrasts_probe_errors.html'),
       digits = 3)

# compute contrasts
contr_probes <- contrast(probe_means, 'tukey', adjust = 'fdr')
# save contrasts table
tab_df(as.data.frame(contr_probes),
       title = 'Contrasts probe model (error rate ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/contrasts_probe_errors.html'),
       digits = 3)
# compute contrast CIs
ci_probes <- confint(contr_probes)
# save CIs table
tab_df(as.data.frame(ci_probes),
       title = 'Contrasts probe model (error rate ~ probe + (1|subject/probe))',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/cis_contrast_probe_errors.html'),
       digits = 3)

# compute descriptives for the observed values
descriptives_errors <- errors %>%
  group_by(probe) %>%
  summarise(mean(error_rate), sd(error_rate))
tab_df(as.data.frame(descriptives_errors),
       title = 'Descripvives error rate by Cue-Probe combination.',
       file = paste0(path_to_data,
                     '/derivatives/results/tables/descriptives_probe_errors.html'),
       digits = 3)

# save computed objects
save(list=c("probes_anova", "probe_means"),
     file=paste0(path_to_data, "/derivatives/results/stats/probe_stats_errors.RData"))
