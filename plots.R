#plotting and such
library(tidyverse)
library(scales)
library(ggpubr)
library(multcomp)
library(broom)
library(cowplot)

##### Aesthetic variables
fitnessColors = c("#666666",
                  "#4B5899", 
                  "#5181A6", 
                  "#60AEB3", 
                  "#5EBF9C", 
                  "#64CC7B",
                  "#999999")

fitnessLabels = c("-1",
                  "-0.014",
                  "-0.012",
                  "-0.010",
                  "-0.008",
                  "-0.006",
                  " 0")

custom_theme <- function (base_size = 16, base_family = "") {
  theme_bw(base_size, base_family, base_line_size = base_size / 30) %+replace%
    theme(text = element_text(size = base_size,
                              colour = "black",
                              family = "Geneva"),
          axis.text = element_text(family = "Courier",
                                   size = base_size / 1.5),
          legend.text = element_text(family = "Courier"),
          axis.text.x = element_text(angle = 90, hjust = 0),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,
                                      color = "black"),
          panel.grid = element_line(color = "#999999"),
          strip.background = element_rect(color = "black",
                                          fill = "#EEEEEE"),
          plot.background = element_rect(fill = NA,
                                         color = NA)
    )
}

label_newline <- function(labels){
  return(label_both(labels, sep = "\n"))
}

##### Load, mutate, filter CSV files
p = 1/6000
#fixation probabilities
f2 <- read_csv("Data/fixation_probs.csv") %>% mutate(Benefit = b, 
                                                  percentFixed = p_fix,
                                                  `Recom. Rate` = r,
                                                  fitness = factor(-s))
f2$Cost <- f2$c
f2f <- f2[,c("Benefit", "Cost", "Recom. Rate", "fitness", 
             "percentFixed", "CI_0.025", "CI_0.975")]
f2Fine <- filter(f2f, Benefit > 0.19 & Benefit < 0.21) %>%
  filter(`Recom. Rate` == 1e-8)
f2Big <- f2f %>% filter(Benefit <= 0.19 | Benefit >= 0.21 | Benefit == 0.2) %>%
  filter(`Recom. Rate` == 1e-8)
f2Sup <- f2f %>% filter(Benefit == 0.19 | Benefit == 0.21 | Benefit == 0.2)

#variation over benefit
singleStats <- read_csv("Data/singleInteractionsMil.csv") %>% 
  mutate(Benefit = benefit, Cost = cost)

singleFine <- singleStats %>% filter(benefit > 0.19 & benefit < 0.21) %>%
                              mutate(fitness = factor(fitness))
singleBig <- singleStats %>% filter(benefit <= 0.19 | benefit >= 0.21 | 
                                      benefit == 0.2) %>%
                              mutate(fitness = factor(fitness))

#variation over recom. rate & benefit
sRecomStats <- read_csv("Data/singleRecombinationFixed.csv") %>%
  mutate(Benefit = benefit,
         `Recom. Rate` = recom_rate,
         fitness = factor(fitness))

#generation vs pi by fitness/BGS
pibgs <- read_csv("Data/piBGS.csv")

##### Graphs
#benefits 0, 0.19, 0.2, 0.21
bp <- ggplot(singleBig,
#ggplot(singleBig,
       aes(x = fitness,
           y = percentFixed,
           group = factor(Benefit),
           fill = fitness)
       ) +
  scale_fill_manual(name = "Fitness",
                    values = fitnessColors,
                    labels = fitnessLabels) +
  scale_x_discrete(labels = fitnessLabels) +
  scale_y_continuous(labels = label_number(accuracy=0.00001)) +
  geom_col() +
  facet_wrap(~Benefit+Cost, 
             scales = "free_y", 
             labeller = label_bquote(cols = paste("b = ", .(Benefit),
                                                  ", c = ", .(Cost),
                                                  sep = "")),
             drop = T
  ) +
  geom_pointrange(data = f2Big,
                  mapping = aes(ymin = CI_0.025, ymax = CI_0.975),
                  shape = "diamond",
                  alpha = 0.8) +
  geom_hline(yintercept = p, 
             color = "black", 
             #size = 1.5,
             alpha = 0.7,
             linetype = "dashed") +
  xlab("Fitness of deleterious allele") +
  ylab("Fixation probability of altruistic allele") + 
  guides(fill = "none") +
  custom_theme()

#benefits 0.197, 0.198, 0.199, 0.201
fp <- ggplot(singleFine %>% filter(Benefit != 0.196 & Benefit != 0.2),
       aes(x = fitness,
           y = percentFixed,
           group = factor(Benefit),
           fill = fitness)) +
  scale_fill_manual(name = "Fitness",
                    values = fitnessColors,
                    labels = fitnessLabels) +
  scale_x_discrete(labels = fitnessLabels) +
  geom_col() +
  geom_pointrange(data = f2Fine %>% 
                    filter(Benefit != 0.196 & Benefit != 0.2),
                mapping = aes(ymin = CI_0.025, ymax = CI_0.975),
                shape = "diamond",
                alpha = 0.8) +
  scale_y_continuous(labels = label_number(accuracy=0.00001)) +
  facet_wrap(~Benefit+Cost, 
             scales = "free_y", 
             labeller = label_bquote(cols = paste("b = ", .(Benefit),
                                                  ", c = ", .(Cost),
                                                  sep = ""))
             ) +
  geom_hline(yintercept = p, 
             color = "black", 
             #size = 1.5,
             alpha = 0.7,
             linetype = "dashed") +
  xlab("Fitness of deleterious allele") +
  ylab("Fixation probability of altruistic allele") + 
  guides(fill = "none") +
  custom_theme()

sp <- plot_grid(bp, fp, labels = "AUTO", rel_widths = c(1.1, 1))

#generation vs pi by fitness/BGS
pp <- 
  ggplot(pibgs %>% mutate(fitness = factor(fitness)),
       aes(x = generation,
           y = pi_neutral,
           color = fitness)) +
  geom_smooth(se = F) +
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +
  scale_color_manual(name = "Fitness of\ndeleterious allele",
                     values = fitnessColors,
                     labels = fitnessLabels) +
  ylab(bquote(paste("Mean ", pi))) +
  xlab("Generation") + 
  custom_theme() %+replace%
  theme(panel.grid.major.x = element_line(color = "#dddddd"),
        panel.grid.minor.x = element_line(color = "#eeeeee"),
        panel.grid.major.y = element_line(color = "#dddddd"),
        panel.grid.minor.y = element_line(color = "#eeeeee"),
        axis.text.x = element_text())


#variation over recom. rate & benefit
rp <- ggplot(sRecomStats,
       aes(x = fitness, 
           y  = percentFixed,
           fill = fitness)) +
  scale_fill_manual(name = "Fitness",
                    values = fitnessColors,
                    labels = fitnessLabels) +
  scale_x_discrete(labels = fitnessLabels) +
  scale_y_continuous(limits = c(0, NA),
                     labels = label_number(accuracy=0.00001)) +
  geom_col() +
  geom_pointrange(data = f2Sup,
                mapping = aes(ymin = CI_0.025, ymax = CI_0.975),
                shape = "diamond",
                alpha = 0.8) +
  facet_grid(scales = "free",
             rows = vars(`Benefit`),
             cols = vars(`Recom. Rate`),
             labeller = label_newline) +
  geom_hline(yintercept = p, 
             color = "black", 
             alpha = 0.7,
             linetype = "dashed") +
  xlab("Fitness of deleterious allele") + 
  ylab("Fixation probability of altruistic allele") + 
  guides(fill = "none") +
  custom_theme()


##### Statistics

#fixation prob. against fitness, a factor
#linear model
singleForAnalysis <- singleStats %>% filter(benefit != 0) %>%
  mutate(fitness = factor(fitness))
model <- glm(percentFixed ~ benefit * fitness,
             data = singleForAnalysis,
             weights = numTotal,
             family = binomial())
View(tidy(glht(model)))

#fixation prob. against BGS, a boolean
singleBGS <- singleStats %>% filter(benefit != 0)
singleBGS$BGS <- (singleBGS$fitness %% 1 != 0)

bgsmodel <- glm(percentFixed ~ benefit * BGS,
                data = singleBGS,
                weights = numTotal,
                family = binomial())
summary(bgsmodel)
View(tidy(glht(bgsmodel)))
