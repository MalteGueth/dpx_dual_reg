# Title     : Analysis of RT
# Objective : Test effect of cue-probe incongrunecy and block effects
# Created by: Jose C. Garcia Alanis
# Created on: 05.03.21
# R version : 4.0.2 (2020-06-22), Taking Off Again

setwd('~/Documents/projects/dpx_dual_reg/code_dpx_dual-reg/')

# source function for fast loading and installing of packages
source('./r_functions/getPacks.R')
source('./r_functions/spr2.R')

# 1) define the path to behavioral data directory -----------------------------

# get system and user information
host <- Sys.info()

# set default path or promt user for other path
if (host['nodename'] == "josealanis-desktop") {
  
  # defaut path in project structure
  path_to_rt <- '../data/derivatives/results'
  
} else {
  
  path_to_rt <- readline('Please provide path to behavioral data: ')
  
}

# 2) import in the data -------------------------------------------------------
# this part requires the package 'dplyr'
getPacks('dplyr')

# files in dir
rt_files <- list.files(path = paste(path_to_rt, 'rt', sep = '/'),
                       full.names = T)

# read in the files
rt_list <- lapply(rt_files, read.table, sep = '\t', header = T)

# put the in a data frame
rt_df <- bind_rows(rt_list, .id = "column_label")

# recode block variable
rt_df <- rt_df %>% mutate(block =  ifelse(block == 2, 0, block)) %>%
  mutate(block = factor(block, labels = c('Baseline', 'Regulation')))

# # exclude subject XX if necessary
# rt_df <- rt_df %>% filter(!subject == 51)

# 3) exploratory analyses correct reactions -----------------------------------

# extract trials with correct reponses
corrects <- rt_df %>%
  filter(reaction_cues == 'Correct' & reaction_probes == 'Correct')

# plausibility checks
# e.g., rt min should be > 0.0, max < 0.750
# e.g., only probes AX, AY, BX and BY should be present
summary(corrects); unique(corrects$probe)

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
  group_by(subject, probe, block) %>%
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
  facet_wrap(~ block, ncol = 2) +
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
