# Title     : Analysis of personality data
# Objective : Check if personality variables are distributed in the correct
#             way (i.e., high and low scoring in personality questionnaires
# Created by: Jose C. Garcia Alanis
# Created on: 10.08.21
# R version : R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"

# --- set working directory ---
# use this if working from r-studio
curr_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dirname(curr_path))

# source function for fast loading and installing of packages
source('./r_functions/getPacks.R')
source('./r_functions/bootdif.R')

# 1) define the path to behavioral data directory ----------------------------

# get system and user information
host <- Sys.info()

# set default path or ask user for other path
if (grep('Joses',  host['nodename']) || grep('ma',  host['nodename'])) {
  
  # default path in project structure
  path_subject_data <- '../data/subject_data/'
  
  } else {
  
  path_subject_data <- readline('Please provide path to subject data: ')
  
}

# 2) import subject data -----------------------------------------------------
getPacks(c('dplyr', 'tidyr'))

ppi_data <- read.table(paste0(path_subject_data, 'ppi_data.tsv'), 
                           header = T, sep = '\t')

# add descriptive group name
ppi_data <- ppi_data %>%
  mutate(pp_group = ifelse(group_pattern == 1, 'low', 
                           ifelse(group_pattern == 2, 'high', NA)))

# extract psychopathic personality inventory (PPI) data
ppi_scales <- ppi_data %>%
  select(id, 
         pp_group,
         si_pr, sp_pr, f_pr)

# check number of observations per group
ppi_data %>% group_by(pp_group) %>% tally()
ppi_data %>% group_by(pp_group, sex) %>% tally()

# mean age and sex
ppi_data %>% group_by(pp_group) %>% summarise(mean(age), sd(age))

# remove that where not assigned to a group
ppi_scales <- ppi_scales %>% 
  filter(!is.na(pp_group)) %>%
  mutate(pp_group = factor(pp_group, 
                           levels = c('high', 'low'),
                           labels = c('High', 'Low')))

# extract subscales of PPI
ppi_subscales <- ppi_scales %>%
  select(pp_group, si_pr:f_pr) %>%
  gather('scale', 'score', -pp_group) %>%
  mutate(scale = factor(scale)) %>%
  mutate(scale = recode(scale, 'si_pr' = 'SI', 'sp_pr' = 'SP', 'f_pr' = 'F')) %>%
  mutate(scale = factor(scale, levels = c('SI', 'SP', 'F')))
  # mutate(scale = recode(scale,
  #                       'blame_externalization' = 'BE',
  #                       'carefree_non_planfulness' = 'CNP',
  #                       'cold_heartedness' = 'CH',
  #                       'dishonest_responsiveness' = 'DR',
  #                       'fearlessness' = 'F',
  #                       'impulsive_nonconformity' = 'IN',
  #                       'machiavellian_egocentricity' = 'ME',
  #                       'social_potency' = 'SP',
  #                       'stress_immunity' ='SI'))

# 3) create exploratory plots ------------------------------------------------
getPacks(c('ggplot2', 'ggbeeswarm', 'see', 'viridis'))

# --- create plot of achieved overall scores in the PPI ---
# allow point to spread a little bit along the x-axis

# pd <- position_nudge(c(-0.25, 0.25))
pn <- position_nudge(0.2)

# create plot
ppi_percentiles_plot <- ggplot(data = ppi_subscales,
       aes(x = scale,
           y = score,
           fill = pp_group,
           shape = pp_group)) +
  facet_wrap(~ pp_group, scales = 'free_y') +
  geom_beeswarm(size = 2.0, alpha = 0.70, 
                color = 'black', stroke = 0.75) + 
  geom_boxplot(width = 0.2, alpha = 0.75, position = pn, 
               outlier.colour = NA, size = 0.8) +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_manual(values = c('#B40F20FF', '#46ACC8FF')) +
  scale_color_manual(values = c('#B40F20FF', '#46ACC8FF')) +
  # scale_fill_viridis(option = 'A', discrete = T, begin = 0.1, end = 0.6, direction = -1) +
  # scale_color_viridis(option = 'A', discrete = T, begin = 0.1, end = 0.6, direction = -1) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 100),
               color = 'black', size = rel(0.5), linetype = 1) +
  geom_segment(aes(x = 'SI', y = -Inf, xend = 'F', yend = -Inf),
               color = 'black', size = rel(0.5), linetype = 1) +
  labs(title = 'PPI percentile scores by group',
       x = '',
       y = 'T-value',
       fill = 'Psychopathy',
       shape = 'Psychopathy') +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = 'gray98', size = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(color = 'black',
                                  size = 12,
                                  face = 'bold',
                                  family='Mukta'),
        plot.title = element_text(color = 'black',
                                  size = 12,
                                  face = 'bold',
                                  family='Mukta'),
        axis.title.x = element_text(color = 'black',
                                    size = 12,
                                    face = 'bold',
                                    margin = margin(t = 15),
                                    family='Mukta'),
        axis.title.y = element_text(color = 'black',
                                   size = 12,
                                   face = 'bold',
                                   margin = margin(r = 15),
                                   family='Mukta'),
        axis.text = element_text(color = 'black',
                                 size = 10),
        legend.position = 'none'); ppi_percentiles_plot
ggsave('../data/derivatives/results/figures/ppi_percentiles_plot.png',
       ppi_percentiles_plot, width = 7.0, height = 3.0, dpi = 600)

# --- create plot of achieved scores in the PPI subscales ---
# allow point to spread a little bit along the x-axis


ppi_means <- ppi_subscales %>%
  group_by(pp_group, scale) %>%
  summarise(mean_cl_boot(score))

getPacks(c('Hmisc'))
# create plot
ppi_means_plot <- ggplot(data = ppi_means,
                         aes(x = scale,
                             y = y,
                             fill = pp_group,
                             shape = pp_group)) +
  geom_line(aes(group = pp_group, color = pp_group), size = 0.6) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax, color = pp_group),
                width = 0.15, size = 0.6, alpha = 0.8) +
  geom_point(size = 1.5, stroke = 0.8) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_manual(values = c('#B40F20FF', '#46ACC8FF')) +
  scale_color_manual(values = c('#B40F20FF', '#46ACC8FF')) +
  geom_segment(aes(x = -Inf, y = 0, xend = -Inf, yend = 100),
               color = 'black', size = rel(0.5), linetype = 1) +
  geom_segment(aes(x = 'SI', y = -Inf, xend = 'F', yend = -Inf),
               color = 'black', size = rel(0.5), linetype = 1) +
  labs(title = 'Bootstrap Means and CIs', 
       y = 'Percentile Score', x  = 'PPI sub-scale',
       fill = 'PPI Group', shape = 'PPI Group', color = 'PPI Group')+
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = 'gray98', size = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(color = 'black',
                                  size = 12,
                                  face = 'bold',
                                  family='Mukta'),
        plot.title = element_text(color = 'black',
                                  size = 12,
                                  face = 'bold',
                                  family='Mukta'),
        axis.title.x = element_text(color = 'black',
                                    size = 12,
                                    face = 'bold',
                                    margin = margin(t = 15),
                                    family='Mukta'),
        axis.title.y = element_text(color = 'black',
                                    size = 12,
                                    face = 'bold',
                                    margin = margin(r = 15),
                                    family='Mukta'),
        axis.text = element_text(color = 'black',
                                 size = 10),
        legend.position = 'right',
        legend.direction = 'vertical',
        legend.key = element_blank(),
        legend.key.size = unit(2, 'line')); ppi_means_plot
ggsave('../data/derivatives/results/figures/ppi_means_plot.png',
       ppi_means_plot, width = 4.0, height = 3.0, dpi = 600)

# 4) test group diffrences in ppi scores -------------------------------------
# initialize place holders for results
ppi_scales <- NULL
norm_test <- NULL
w_stat <- NULL
p_vals <- NULL
boot_diffs <- NULL
boot_cis <- NULL
# loop through scales ans compute statistics
for (var in c('age', 'total_pr', 'si_pr', 'sp_pr', 'f_pr')) {
  low <- ppi_data[ppi_data$pp_group == 'low', var]
  high <- ppi_data[ppi_data$pp_group == 'high', var]

  # tests
  norm <- shapiro.test(c(low, high))
  w_res <- wilcox.test(low, high, alternative = 'less',
                       paired = FALSE, exact = FALSE)

  # concatenate
  ppi_scales <- c(ppi_scales, var)
  norm_test <- c(norm_test, norm$p.value)
  t_stat <- c(w_stat, w_res$statistic)
  p_vals <- c(p_vals, w_res$p.value)

  diff <- bootdif(y = ppi_data[, var], g = ppi_data$pp_group, b = 10000)
  boot_diffs <- c(boot_diffs, diff['Mean'])
  boot_cis <- c(boot_cis,
                paste(round(diff['Lower'], 3),
                      round(diff['Upper'], 3),
                      sep = ', '))

}

# create summary table
pp_scales <- data.frame(ppi_scales)
pp_scales$normality <- norm_test
pp_scales$wx_stat <- t_stat
pp_scales$p_value <- p_vals
pp_scales$fdr <- p.adjust(p_vals, 'fdr')
pp_scales$sig <- ifelse(pp_scales$fdr < 0.001, '***',
                        ifelse(pp_scales$fdr < 0.01, '**',
                               ifelse(pp_scales$fdr < 0.05, '*', '')))
pp_scales$boot_diff <- boot_diffs
pp_scales$boot_ci <- boot_cis

# export to html
getPacks('sjPlot')
tab_df(title = 'Group differences ins PPI scales',
       pp_scales, digits = 3,
       file = '../data/derivatives/results/tables/ppi_scales_results.html')
