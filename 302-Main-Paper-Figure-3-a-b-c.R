
library(tidyverse)
library(gridExtra)
library(cowplot)


Geno_long <- readRDS("data/Geno_long.rds")

# now, we compute the proportion of each genotype in each interval, grouped by sex
# and we compute the SE around that estimate.
geno_summary <- Geno_long %>%
  filter(!is.na(LENGTH), !is.na(geno), !is.na(GenSex)) %>%
  group_by(GenSex, len_interval, locus) %>%
  mutate(ave_length = mean(LENGTH)) %>%
  group_by(GenSex, len_interval, ave_length, locus, geno) %>%    # include ave_length so it is retained in the output
  tally() %>%
  rename(counts = n) %>%
  mutate(total = sum(counts),
         ppn = counts / total,
         se = sqrt(ppn * (1 - ppn) / total)) %>%
  ungroup



# OK, now we can whittle that down to just the omy5locus:
omy5locus <- geno_summary %>%
  filter(locus == "SH114448.87") %>%
  mutate(`MAR genotype` = factor(plyr::revalue(geno, c("0" = "RR", "1" = "AR", "2" = "AA")), levels = c("RR", "AR", "AA")),
         Sex = GenSex)

# and we can also summarize the different length classes and look at the proportion
# of the different sexes of each
sex_ppns <- Geno_long %>%
  filter(!is.na(GenSex), !is.na(LENGTH), locus == "SexID") %>%  #   we just filter it down to the SexID locus...
  mutate(len_interval = cut_number(LENGTH, 16)) %>%   # cut this into 16 equal sized classes, so we have roughly 170 fish per size class still
  group_by(len_interval) %>%
  mutate(ave_length = mean(LENGTH)) %>%
  group_by(len_interval, ave_length, GenSex) %>%
  tally() %>%
  rename(counts = n) %>%
  mutate(total = sum(counts),
         ppn = counts / total,
         se = sqrt(ppn * (1 - ppn) / total)) %>%
  ungroup %>%
  rename(Sex = GenSex) 




#### Make the Plots  ####

# set colors to use:
AA_col <- "royalblue4"
RR_col <- "red"
AR_col <- "#FF7F00"
our_cols <- c(AA = AA_col, AR = AR_col, RR = RR_col)
#our_cols <- dichromat::dichromat(c(AA = AA_col, AR = AR_col, RR = RR_col), type = "deutan")


geno_by_sex <- ggplot(omy5locus, aes(x = ave_length, y = ppn, colour = `MAR genotype`, linetype = Sex, shape = `MAR genotype`)) +
  geom_line() +
  geom_linerange(aes(ymin = ppn + 1 * se, ymax = ppn - 1 * se), position = position_dodge(width = 1.8)) +
  geom_point(size = 3) +
  ylim(0, NA) +
  xlim(50, 200) +
  theme_bw() +
  theme(legend.key.width = grid::unit(2.0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Genotype frequency") +
  xlab("Average length (mm) of fish in size class") +
  guides(colour = guide_legend(reverse=T, title = "Omy05\nrearrangement"),
         shape = guide_legend(reverse=T, title = "Omy05\nrearrangement")) +
  scale_colour_manual(values = our_cols)

ggsave(geno_by_sex, filename = "length_vs_mar_geno_by_sex.pdf", width = 6.5, height = 5)



length_vs_sex <- ggplot(sex_ppns, aes(x = ave_length, y = ppn, linetype = Sex)) +
  geom_line() +
  geom_linerange(aes(ymin = ppn + 1 * se, ymax = ppn - 1 * se), position = position_dodge(width = 1.8)) +
  # geom_point(size = 3) +
  ylim(0, NA) +
  xlim(50, 200) +
  theme_bw() +
  theme(legend.key.width = grid::unit(2.0, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Relative frequency of sex") +
  xlab("Average length (mm) of fish in size class")

ggsave(length_vs_sex, filename = "length_vs_sex.pdf", width = 6.5, height = 5)


# if we want to put those on the same page we can do:
combo <- grid.arrange(geno_by_sex, length_vs_sex, nrow = 2)

ggsave(combo, width = 6.5, height = 8, filename = "two-stacked-figs.pdf")

#### Now, get Steve's data and use it:
lindley <- readRDS("data/lindley-results.rds") %>%
  mutate(geno = str_match(fishtype, "ale([AR][AR])")[,2]) %>%
  mutate(sex = str_match(fishtype, "([MF].*ale)")[,2]) %>%
  mutate(geno = ifelse(geno == "RA", "AR", geno)) %>%
  mutate(geno = factor(geno, levels = c("AA", "AR", "RR"))) %>%
  mutate(sex = factor(sex)) %>%
  mutate(xpos = (as.integer(sex) - 1) + (0.1 * as.integer(geno)) - 0.4 * (sex == "Male"))

# so now we can make something that looks like this plot 
# using ggplot
lindley_plot <- ggplot(lindley, aes(colour = geno, linetype = sex)) +
  geom_segment(aes(x = xpos, xend = xpos, y = plusSE, yend = minusSE)) +
  geom_point(aes(x = xpos, y = emig.frac, shape = geno), size = 3) +
  scale_colour_manual(values = our_cols) +
  ylab("Detection probability") +
  ylim(0, 1) + 
  xlim(0, 1) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0, yend = 1) + 
 # annotate("") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  theme(legend.position="none") +
  annotate("text", x = 0.23, y = 0.9, label = "Females", hjust = 0.5) +
  annotate("text", x = 0.76, y = 0.9, label = "Males", hjust = 0.5) +
  scale_shape_manual(values = c(AA = "square", AR = "triangle", RR = "circle"))




# or we can use cowplot::plot_grid and then get labels on there
# too...
left_col <- cowplot::plot_grid(lindley_plot, labels = "A")

right_col <- cowplot::plot_grid(geno_by_sex, length_vs_sex, labels = c("B", "C"), align = "v", ncol = 1)

full_panel <- cowplot::plot_grid(left_col, right_col, rel_widths = c(0.45, 1))
# OK, that is pretty slick...

ggsave(full_panel, filename = "fig3-abc-sex-geno-length-detection.pdf", width = 7, height = 6)






