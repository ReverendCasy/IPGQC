library(dplyr)
library(ggplot2)
library(magrittr)
library(outliers)



## Initial stats
### Assembly num
Freq10 <- read.delim(file.path('current_stats',
                               'species_stats_sorted.tsv'), 
                     sep = '\t', 
                     header = FALSE, 
                     col.names = c('Species', 'Frequency'))

### IPG_stats
ipg_singleplus <- read.delim(file.path('current_stats',
                               'with_singletons',
                               'PerAssemblyStats_singletons.tsv'),
                     header = FALSE,
                     stringsAsFactors = FALSE,
                     sep = '\t') %>% 
  `colnames<-`(c('Species', 'Assembly', 'Source', 'NoIPG'))

serr <- function(x) return(sd(x) / sqrt(length(x)))

assembly_ints_plus <- ipg_singleplus %>% aggregate(NoIPG ~ Species,
                                              data = .,
                                              function(x) c(Mean = x %>% mean() %>% round,
                                                            Std.Err = x %>% serr %>% round)) %>%
  transmute(Species = Species, Mean = NoIPG[,1], StdErr = NoIPG[,2]) %>% 
  as.data.frame() %>% 
  `colnames<-`(c('Species', 'Mean', 'StdErr'))


## Singleton-free stats
ipg_singleminus <- read.delim(file.path('current_stats',
                                       'no_singletons',
                                       'PerAssemblyStats_nosingletons.tsv'),
                             header = FALSE,
                             stringsAsFactors = FALSE,
                             sep = '\t') %>% 
  `colnames<-`(c('Species', 'Assembly', 'Source', 'NoIPG'))

assembly_ints_minus <- ipg_singleminus %>% aggregate(NoIPG ~ Species,
                                                   data = .,
                                                   function(x) c(Mean = x %>% mean() %>% round,
                                                                 Std.Err = x %>% serr %>% round)) %>%
  transmute(Species = Species, Mean = NoIPG[,1], StdErr = NoIPG[,2]) %>% 
  as.data.frame() %>% 
  `colnames<-`(c('Species', 'Mean', 'StdErr'))

assembly_stats_all <- merge(assembly_ints_plus,
                            assembly_ints_minus,
                            by = 'Species') %>% 
                      `colnames<-`(c('Species',
                                   'Mean_singleplus',
                                   'StdErr_singleplus',
                                   'Mean_singleminus',
                                   'StdErr_singleminus'))

assembly_stats_all$PercOfSingletons <- sapply(1:nrow(assembly_stats_all),
                                              function(x) (assembly_stats_all[x, 2] - assembly_stats_all[x, 4]) /
                                                assembly_stats_all[x, 2] * 100)

## Outlier removal
hfilter <- function(vec) {
  lbow <- median(vec) - 3 * mad(vec, constant = 1)
  ubow <- median(vec) + 3 * mad(vec, constant = 1)
  ind <- which(vec < ubow & vec > lbow)
  return(ind)
}



## Figures

### Total per species content
g1 <- ggplot(Freq10, 
             aes(x = reorder(Species, -Frequency),
                 y = Frequency)) +
  geom_bar(stat = 'identity', color = '#95999b', fill = '#dbdcdd') +
  labs(x = 'Species',
       y = 'No. of genomes', 
       tag = 'A') + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(color = '#383e42', size = 0.3, linetype = 'dashed'),
        panel.grid.minor.y = element_line(color = '#383e42', size = 0.15, linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.tag = element_text(size = 16, face = 'bold')) + 
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
g1

### Top ten sequenced species
g2 <- ggplot(Freq10 %>% 
               head(10) %>% 
               mutate(Species = sapply(Species, function(x) x %>% 
                                         regexPipes::gsub('([A-Z])[a-z]+ (.*)', '\\1. \\2',
                                                          perl = TRUE))), 
             aes(x = reorder(Species, -Frequency),
                 y = Frequency)) +
  geom_bar(stat = 'identity', color = '#95999b', fill = '#dbdcdd') +
  labs(x = 'Species', 
       y = 'No. of genomes',
       tag = 'B') + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 50, hjust = 1, face = 'italic'),
        axis.text.y = element_text(size = 12, color = 'black'),
        panel.grid.major.y = element_line(color = '#383e42', size = 0.3, linetype = 'dashed'),
        panel.grid.minor.y = element_line(color = '#383e42', size = 0.15, linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.tag = element_text(size = 16, face = 'bold')) + 
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

### Per species IPG distribution (singletons included)
g3 <- ggplot(assembly_ints_plus,
             aes(x = reorder(Species, -Mean),
                 y = Mean)) +
  geom_bar(stat = 'identity', color = '#95999b', fill = '#dbdcdd') +
  geom_errorbar(aes(ymin = Mean - StdErr, ymax = Mean + StdErr), color = '#383e42') +
  labs(x = 'Species', 
       y = 'Mean No. of IPGS',
       tag = 'C') + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(color = '#383e42', size = 0.3, linetype = 'dashed'),
        panel.grid.minor.y = element_line(color = '#383e42', size = 0.15, linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.tag = element_text(size = 16, face = 'bold')) + 
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

### Top ten IPG content species (singletons included)
g4 <- ggplot(assembly_ints_plus[order(assembly_ints_plus$Mean, decreasing = TRUE), ] %>%
               head(10) %>% 
               mutate(Species = sapply(Species, function(x) x %>% 
                                         regexPipes::gsub('([A-Z])[a-z]+ (.*)', '\\1. \\2',
                                                          perl = TRUE))),
             aes(x = reorder(Species, -Mean),
                 y = Mean)) +
  geom_bar(stat = 'identity', 
           color = '#95999b', 
           fill = '#dbdcdd') +
  geom_errorbar(aes(ymin = Mean - StdErr, 
                    ymax = Mean + StdErr), 
                color = '#383e42') +
  labs(x = 'Species', 
       y = 'Mean No. of IPGS',
       tag = 'D') + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 50, hjust = 1, face = 'italic'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(color = '#383e42', size = 0.3, linetype = 'dashed'),
        panel.grid.minor.y = element_line(color = '#383e42', size = 0.15, linetype = 'dashed'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.tag = element_text(size = 16, face = 'bold')) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
g4

##
