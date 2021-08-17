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
if (grep('Joses', host['nodename']) || grep('ma', host['nodename'])) {

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
  mutate(block =  ifelse(block == 2, 0, block)) %>%
  mutate(block = factor(block, labels = c('Baseline', 'Regulation')))

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

# plausibility checks
# e.g., rt min should be > 0.0, max < 0.750
# e.g., only probes AX, AY, BX and BY should be present
summary(corrects); unique(corrects$probe)

# get rid of extreme values, e.g., rt values very close to 0
# using winsorsed scores
getPacks(c('psych'))
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

# 4) statistical analyses correct reactions -----------------------------------
getPacks(c('lme4', 'lmerTest', 'car'))

corr_mod <- corrects %>%
  arrange(subject) %>%
  mutate(pp_group = factor(pp_group, levels = c('High', 'Low')),
         probe = factor(probe, levels = c('AX', 'AY', 'BX', 'BY')),
         subject = factor(subject, levels = sort(unique(corrects$subject))),
         rt = w_rt * 1000) 
# %>%
#   group_by(subject, pp_group, block, probe) %>%
#   summarise(mean_rt = mean(w_rt))

mod_rt_0 <- lmer(data = corr_mod, 
                 rt ~  probe + (1+block|subject/probe), 
                 contrasts = list(probe = 'contr.sum'))
anova(mod_rt_0)
# car::Anova(mod_rt_0, test = 'F', type = 'III')
summary(mod_rt_0)

getPacks(c('emmeans'))
# emm_options(pbkrtest.limit = 14450)
probe_means <- emmeans(mod_rt_0, ~ probe)
contrast(probe_means, 'tukey', adjust = 'fdr')

getPacks(c('sjPlot', 'performance'))
plot_model(mod_rt_0, 'int')
check_model(mod_rt_0)

# 5) plot effect of probe -----------------------------------------------------
getPacks(c('ggplot2', 'ggbeeswarm', 'see', 'viridis'))

mean_rt <- corr_mod %>%
  ungroup() %>%
  group_by(subject, pp_group, probe) %>%
  summarise(mean_rt = mean(rt))

pj = position_jitter(0.10, seed = 52)
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
