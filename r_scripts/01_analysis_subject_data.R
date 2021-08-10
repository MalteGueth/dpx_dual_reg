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
  path_subject_data <- '../data/subject_data/subject_data_dpx-dualreg.tsv'
  
  } else {
  
  path_subject_data <- readline('Please provide path to subject data: ')
  
}

# 2) import subject data -----------------------------------------------------
getPacks(c('dplyr', 'tidyr'))

subject_data <- read.table(path_subject_data, header = T, sep = '\t')
subject_data <- subject_data %>%
  drop_na() %>%
  mutate(pp_group = ifelse(group == 0, 1, group)) %>%
  mutate(pp_group = factor(pp_group, labels = c('low',
                                                'high')))

# extract psychopathic personality inventory (PPI) data
ppi_data <- subject_data %>%
  select(pp_group,
         psychopathy_score,
         impulsive_nonconformity:dishonest_responsiveness)

# extract subscales of PPI
ppi_subscales <- ppi_data %>%
  select(pp_group, impulsive_nonconformity:dishonest_responsiveness) %>%
  gather('scale', 'score', -pp_group) %>%
  mutate(scale = factor(scale)) %>%
  mutate(scale = recode(scale,
                        'blame_externalization' = 'BE',
                        'carefree_non_planfulness' = 'CNP',
                        'cold_heartedness' = 'CH',
                        'dishonest_responsiveness' = 'DR',
                        'fearlessness' = 'F',
                        'impulsive_nonconformity' = 'IN',
                        'machiavellian_egocentricity' = 'ME',
                        'social_potency' = 'SP',
                        'stress_immunity' ='SI'))

# 3) create exploratory plots ------------------------------------------------
getPacks(c('ggplot2', 'viridis'))

# --- create plot of achieved overall scores in the PPI ---
# allow point to spread a little bit along the x-axis
pj <- position_jitter(0.125)
pd <- position_nudge(c(-0.25, 0.25))
# create plot
overall_ppi_plot <- ggplot(data = ppi_data,
       aes(x = pp_group,
           y = psychopathy_score,
           fill = pp_group,
           shape = pp_group)) +
  geom_boxplot(width = 0.1, alpha = 0.5, position = pd) +
  geom_point(position = pj, size = 1.5) +
  scale_shape_manual(values = c(25, 24)) +
  scale_fill_viridis(option = 'B', discrete = T, begin = 0.1, end = 0.6) +
  scale_color_viridis(option = 'B', discrete = T, begin = 0.1, end = 0.6) +
  scale_y_continuous(limits = c(15, 85)) +
  geom_segment(aes(x = -Inf, y = 20, xend = -Inf, yend = 80),
               color = 'black', size = rel(0.5), linetype = 1) +
  geom_segment(aes(x = 'low', y = -Inf, xend = 'high', yend = -Inf),
               color = 'black', size = rel(0.5), linetype = 1) +
  labs(x = '',
       y = 'Overall PPI Score',
       fill = 'Psychopathy',
       shape = 'Psychopathy') +
  theme(panel.background = element_rect(fill = 'gray98'),
        axis.title.x = element_text(color = 'black',
                                    size = 12,
                                    face = 'bold',
                                    margin = margin(t = 15)),
        axis.title.y= element_text(color = 'black',
                                   size = 12,
                                   face = 'bold',
                                   margin = margin(r = 15)),
        axis.text = element_text(color = 'black',
                                 size = 10),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.key = element_blank()); overall_ppi_plot
ggsave('../data/derivatives/results/figures/overall_ppi_plot.pdf',
       overall_ppi_plot, width = 3.0, height = 5.0)

# --- create plot of achieved scores in the PPI subscales ---
# allow point to spread a little bit along the x-axis
pj <- position_jitterdodge(jitter.width = 0.25, dodge.width = 0.25)
pd <- position_nudge(c(-0.20, 0.20))
# create plot
subscale_ppi_plot <- ggplot(data = ppi_subscales,
       aes(x = scale,
           y = score,
           fill = pp_group,
           shape = pp_group)) +
  geom_boxplot(width = 0.08, alpha = 0.5, position = pd, outlier.color = NA) +
  geom_point(position = pj, size = 1.5) +
  scale_shape_manual(values = c(25, 24)) +
  scale_fill_viridis(option = 'B', discrete = T, begin = 0.1, end = 0.6) +
  scale_color_viridis(option = 'B', discrete = T, begin = 0.1, end = 0.6) +
  scale_y_continuous(limits = c(15, 85)) +
  geom_segment(aes(x = -Inf, y = 20, xend = -Inf, yend = 80),
               color = 'black', size = rel(0.5), linetype = 1) +
  geom_segment(aes(x = 'BE', y = -Inf, xend = 'SI', yend = -Inf),
               color = 'black', size = rel(0.5), linetype = 1) +
  labs(x = 'PPI Subscale',
       y = 'Score',
       fill = 'Psychopathy',
       shape = 'Psychopathy') +
  theme(panel.background = element_rect(fill = 'gray98'),
        axis.title.x = element_text(color = 'black',
                                    size = 12,
                                    face = 'bold',
                                    margin = margin(t = 15)),
        axis.title.y= element_text(color = 'black',
                                   size = 12,
                                   face = 'bold',
                                   margin = margin(r = 15)),
        axis.text.x = element_text(color = 'black',
                                   size = 10),
        axis.text.y = element_text(color = 'black',
                                   size = 10),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.key = element_blank()); subscale_ppi_plot
ggsave('../data/derivatives/results/figures/subscale_ppi_plot.pdf',
       subscale_ppi_plot, width = 10.0, height = 5.0)

# 4) test group diffrences in ppi scores -------------------------------------
# initialize place holders for results
ppi_scales <- NULL
norm_test <- NULL
w_stat <- NULL
p_vals <- NULL
boot_diffs <- NULL
boot_cis <- NULL
# loop through scales ans compute statistics
for (var in names(ppi_data)[-1]) {
  low <- ppi_data[ppi_data$pp_group == 'low', var]
  high <- ppi_data[ppi_data$pp_group == 'high', var]

  norm <- shapiro.test(c(low, high))
  w_res <- t.test(low, high,
                  paired = FALSE, alternative = 'less')
  ppi_scales <- c(ppi_scales, var)
  norm_test <- c(norm_test, norm$p.value)
  w_stat <- c(w_stat, w_res$statistic)
  p_vals <- c(p_vals, w_res$p.value)

  diff <- bootdif(y = ppi_data[, var], g = ppi_data$pp_group)
  boot_diffs <- c(boot_diffs, diff['Mean'])
  boot_cis <- c(boot_cis,
                paste(round(diff['Lower'], 3),
                      round(diff['Upper'], 3),
                      sep = ', '))

}

# create summary table
pp_scales <- data.frame(ppi_scales)
pp_scales$normality <- norm_test
pp_scales$statistic <- w_stat
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
