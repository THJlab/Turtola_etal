##title: poly(A)tail length control
##authors: Manfred Schmid, Pawel Krawczyk, Matti Turtola
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(viridis)
library(readr)

#load data
(pA_df <- readRDS("data/pA.rds"))

distinct(pA_df, temp, rapa, plasmid, pab1DRRM4)
head(pA_df)

pA_df2 <- pA_df %>%
  separate(contig, c("name", "chr"), "::")

#read numbers in each replicate
persample_totals <- pA_df2 %>%
  group_by(temp, rapa, plasmid, pab1DRRM4, sample_nr) %>%
  summarize(cnt = n())

#read numbers in each experimental condition
percondition_totals <- pA_df2 %>%
  group_by(temp, rapa, plasmid, pab1DRRM4) %>%
  summarize(cnt = n())

#reads per transcript (contig) in each replicate
persample_per_contig_totals <- pA_df2 %>%
  group_by(name, temp, rapa, plasmid, pab1DRRM4, sample_nr) %>%
  summarize(cnt = n())

#read per transcript (contig) in each experimental condition
percondition_per_contig_totals <- pA_df2 %>%
  group_by(name, temp, rapa, plasmid, pab1DRRM4) %>%
  summarize(cnt = n())


#violin plots of pA length distributions for all transcripts
plot_violin <- function(df,
                        repeats_split = F,
                        cnt_df = percondition_totals) {
  if (repeats_split) {
    p <- df %>%
      ggplot(., aes(x = sample_nr, y = polya_length, fill = plasmid)) +
      geom_violin() +
      geom_boxplot(width = .1,
                   fill = 'lightgray',
                   outlier.shape = NA)
  } else{
    p <- df %>%
      ggplot(., aes(x = 1, y = polya_length, fill = plasmid)) +
      geom_violin() +
      geom_boxplot(width = .1,
                   fill = 'lightgray',
                   outlier.shape = NA)
  }
  p <- p +
    facet_grid(. ~ pab1DRRM4 + temp + rapa + plasmid, scales = 'free') +
    coord_cartesian(ylim = c(0, 250)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0, 'mm'),
      strip.text.x = element_text(
        angle = 90,
        hjust = 0,
        size = 8
      ),
      strip.background = element_blank(),
      panel.border = element_rect(size = .1),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  if (repeats_split)
  {
    p <- p +
      geom_text(
        data = cnt_df,
        aes(
          x = sample_nr,
          y = 250,
          label = paste0('n=', cnt)
        ),
        angle = 90,
        hjust = 'right',
        vjust = 'top'
      )
  } else{
    p <- p +
      geom_text(
        data = cnt_df,
        aes(
          x = 1,
          y = 250,
          label = paste0('n=', cnt)
        ),
        angle = 90,
        hjust = 'right',
        vjust = 'top'
      )
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  return(p)
}

#violin plot of pA length distributions for each experimental condition
##WARNING: takes very long to calculate....
plot_violin(pA_df2, repeats_split = F, cnt_df = percondition_totals)

#violin plot of pA length distributions for each replicate
##WARNING: takes very long to calculate....
plot_violin(pA_df2, repeats_split = T, cnt_df = persample_totals)


##violin plot of pA length distributions for individual genes
##replace gene name to generate different gene profiles
plot_violin_for_gene <- function(plotting_function = plot_violin,
                                 gene = 'YLL026W_HSP104',
                                 ...) {
  p <- plotting_function(filter(pA_df2, grepl(gene, name)),
                         cnt_df = filter(percondition_per_contig_totals, grepl(gene, name)),
                         ...)
  p <- p +
    ggtitle(gene)
  return(p)
}

plot_violin_for_gene(plot_violin, gene = 'YLL026W_HSP104')


##histograms of pA lengths for individual genes
plot_histogram <- function(df,
                           cnt_df = percondition_totals) {
  p <- df %>%
    ggplot(., aes(x = polya_length, fill = plasmid)) +
    geom_histogram(bins = 250) +
    coord_cartesian(xlim = c(0, 250)) +
    facet_grid(pab1DRRM4 + temp + rapa + plasmid ~ ., scales =
                 'free') +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      panel.spacing.y = unit(0, 'mm'),
      strip.text.y = element_text(
        angle = 0,
        hjust = 0,
        size = 8
      ),
      strip.background = element_blank(),
      panel.border = element_rect(size = .1),
      axis.text.y = element_blank()
    )
  p <- p +
    facet_grid(pab1DRRM4 + temp + rapa + plasmid ~ ., scales = 'free') +
    geom_text(
      data = cnt_df,
      aes(
        x = 250,
        y = Inf,
        label = paste0('n=', cnt)
      ),
      hjust = 'right',
      vjust = 'top',
      size = 3,
      nudge_y = .4
    )
  return(p)
}

plot_histogram_for_gene <- function(plot_fun = plot_histogram,
                                    gene = 'YLL026W_HSP104',
                                    ...) {
  p <- plot_fun(filter(pA_df2, grepl(gene, name)),
                cnt_df = filter(percondition_per_contig_totals,
                                grepl(gene, name)),
                ...)
  p <- p +
    ggtitle(gene)
  return(p)
}

plot_histogram_for_gene(plot_histogram, 'YLL026W_HSP104')



#define subclass for heat-induced transcripts:
#edgeR analysis of expression changes between 38C vs 25C in 'empty;no rapa' samples
# this is done with function from R package nanotail (https://github.com/smaegol/nanotail) which can be installed like this:
# install.packages("devtools")
# devtools::install_github('smaegol/nanotail')
library(nanotail)

diff_expr_temp <- nanotail::calculate_diff_exp_binom(
  pA_df %>%
    dplyr::filter(plasmid == "empty", rapa == "no", pab1DRRM4 == "") %>%
    mutate(temp = factor(temp)) %>%
    dplyr::rename(transcript = contig),
  grouping_factor = "temp",
  condition1 = "y25",
  condition2 = "y38"
)


#histogram of differential expression at 38C vs. 25C, for "MEX67-AA; pPgal-empty; no rapamycin" samples
hist(log2(diff_expr_temp$fold_change), breaks = 100)

#differential expression [38C vs. 25C, for "MEX67-AA; pPgal-empty; no rapamycin" samples] plotted with mean expression,
plot(
  log2(diff_expr_temp$fold_change) ~ diff_expr_temp$mean_expr,
  log = 'x',
  pch = ifelse(diff_expr_temp$fold_change > 2, 19, 16),
  col = case_when(
    diff_expr_temp$padj < .05 &
      diff_expr_temp$fold_change > 2 ~ 'darkred',
    diff_expr_temp$padj < .05 &
      diff_expr_temp$fold_change > 1 ~ 'red',
    diff_expr_temp$padj < .05 &
      diff_expr_temp$fold_change < 1 ~ 'blue',
    diff_expr_temp$padj >= .05 ~ 'gray'
  ),
  cex = ifelse(diff_expr_temp$padj < .05, .2, .1)
)
abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
abline(v = c(25, 50, 100), lty = 2)

#define hsRNAs; fold change > 3 [38 vs 25 in empty no rapa], mean_expr >50
heat_shock_induced <- diff_expr_temp %>%
  filter(fold_change != Inf,
         significance != "NotSig",
         mean_expr > 50,
         fold_change > 3)

heat_shock_induced$transcript

#SSA3,AAR2 and tT(UGU)G2,YGR256W are annotated together with neighboring ncRNAs but their
#pA reads originate from the mRNA 3´ends; label them later as SSA3 and GND2

hs_ids <- sub('::.*', '', heat_shock_induced$transcript)
hs_df <- pA_df2 %>% filter(name %in% hs_ids)



#hsRNA read numbers in each replicate
persample_totals_hs <- hs_df %>%
  group_by(temp, rapa, plasmid, pab1DRRM4, sample_nr) %>%
  summarize(cnt = n())

#hsRNA read numbers in each experimental condition
percondition_totals_hs <- hs_df %>%
  group_by(temp, rapa, plasmid, pab1DRRM4) %>%
  summarize(cnt = n())

#hsRNA reads per transcript (contig) in each replicate
persample_per_contig_totals_hs <- hs_df %>%
  group_by(name, temp, rapa, plasmid, pab1DRRM4, sample_nr) %>%
  summarize(cnt = n())

#hsRNA read per transcript (contig) in each experimental condition
percondition_per_contig_totals_hs <- hs_df %>%
  group_by(name, temp, rapa, plasmid, pab1DRRM4) %>%
  summarize(cnt = n())



#violin plots, hsRNAs per condition
plot_violin(hs_df, cnt_df = percondition_totals_hs)

#violin plots, hsRNAs per replicate
plot_violin(hs_df, repeats_split = T, cnt_df = persample_totals_hs)



#pA length distributions, binned, hsRNAs
binnify <-
  function(df,
           splits = c(Inf, 100, 70, 50, 30, 10, 0),
           bin_names = c('>=100', '70-99', '50-69', '30-49', '10-29', '<10'),
           totals = percondition_totals_hs,
           repeats_split = T) {
    df %<>%
      mutate(
        polya_length_bin = cut(
          polya_length,
          breaks = splits,
          labels = bin_names,
          include.lowest = T
        )
      )
    if (repeats_split) {
      df %<>%
        group_by(pab1DRRM4, temp, rapa, plasmid, sample_nr, polya_length_bin)
    } else{
      df %<>%
        group_by(pab1DRRM4, temp, rapa, plasmid, polya_length_bin)
    }
    df %<>%
      summarize(reads_in_bin = n()) %>%
      left_join(., totals) %>%
      mutate(freq_in_bin = reads_in_bin / cnt) %>%
      left_join(., totals) %>%
      ungroup %>%
      mutate(
        freq_in_bin = reads_in_bin / cnt,
        polya_length_bin = factor(polya_length_bin,
                                  levels = bin_names)
      )
  }

percondition_binned <- binnify(
  hs_df,
  splits = c(0, 10, 30, 50, 70, 100, Inf),
  bin_names = c('<10', '10-29', '30-49', '50-69', '70-99', '>=100'),
  totals = percondition_totals_hs
)

persample_binned <- binnify(
  hs_df,
  splits = c(0, 10, 30, 50, 70, 100, Inf),
  bin_names = c('<10', '10-29', '30-49', '50-69', '70-99', '>=100'),
  totals = persample_totals_hs,
  repeats_split = T
)

plot_binbars <- function(df,
                         repeats_split = T,
                         show_bins = c('>=100', '70-99', '50-69', '30-49', '10-29', '<10'),
                         bin_colors = viridis(6)[1:6],
                         cnt_df = percondition_totals_hs) {
  df <- df %>% dplyr::filter(polya_length_bin %in% show_bins) %>%
    mutate(polya_length_bin = factor(polya_length_bin, show_bins))
  if (repeats_split) {
    p <- df %>%
      ggplot +
      geom_bar(aes(x = sample_nr, y = freq_in_bin, fill = polya_length_bin),
               stat = 'identity') +
      facet_grid(. ~ pab1DRRM4 + temp + rapa + plasmid, scales = 'free')
  } else{
    p <- df %>%
      ggplot +
      geom_bar(aes(x = plasmid, y = freq_in_bin, fill = polya_length_bin),
               stat = 'identity') +
      facet_grid(. ~ pab1DRRM4 + temp + rapa, scales = 'free')
  }
  p <- p +
    scale_fill_manual(values = bin_colors) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    coord_cartesian(ylim = c(0, 1.6)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0, 'mm'),
      strip.text.x = element_text(
        angle = 90,
        hjust = 0,
        size = 8
      ),
      strip.background = element_blank(),
      panel.border = element_rect(size = .1),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  if (repeats_split) {
    p <- p +
      geom_text(
        data = cnt_df,
        aes(
          x = sample_nr,
          y = 1.01,
          label = paste0('n=', cnt)
        ),
        color = 'black',
        angle = 90,
        hjust = 'left',
        vjust = 'center',
        size = 3,
        nudge_x = .1
      )
  } else{
    p <- p +
      geom_text(
        data = cnt_df,
        aes(
          x = plasmid,
          y = 1.01,
          label = paste0('n=', cnt)
        ),
        color = 'black',
        angle = 90,
        hjust = 'left',
        vjust = 'bottom',
        size = 3,
        nudge_x = .1
      ) +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  return(p)
}

#bins per condition
plot_binbars(percondition_binned, repeats_split = F)

#bins per replicate
plot_binbars(persample_binned,
             repeats_split = T,
             cnt_df = persample_totals_hs)



#pA length density plots for hsRNAs
#plot pA-length densities for all reads within "heat_shock_induced" genes
hs_df %>%
  ggplot(., aes(
    x = polya_length,
    color = plasmid,
    fill = plasmid,
    alpha = .5
  )) +
  geom_density() +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 0.03)) +
  facet_grid(plasmid + rapa + pab1DRRM4 ~ temp, scales = 'free') +
  theme_bw() +
  theme(panel.grid = element_blank())

#overlay densities for MEX67-AA and MEX67-AA/pab1deltaRRM4 strains
hs_df %>%
  ggplot(.,
         aes(
           x = polya_length,
           color = pab1DRRM4,
           fill = plasmid,
           alpha = .5
         )) +
  geom_density() +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 0.03)) +
  facet_grid(plasmid + rapa ~ temp, scales = 'free') +
  theme_bw() +
  theme(panel.grid = element_blank())

#overlay densities for plasmids
hs_df %>%
  ggplot(., aes(
    x = polya_length,
    color = plasmid,
    fill = plasmid,
    alpha = .5
  )) +
  geom_density() +
  coord_cartesian(xlim = c(0, 250), ylim = c(0, 0.03)) +
  facet_grid(rapa + pab1DRRM4 ~ temp, scales = 'free') +
  theme_bw() +
  theme(panel.grid = element_blank())






#Analysis of 4tU RNA pA tail lengths
#load data
(tu_df <- readRDS("data/4tU_pA.rds"))

head(tu_df)
distinct(tu_df, plasmid, id, rapa, rep)

#sample id 18 was wrongly annotated as 'empty' and 'minusrapa' although it is 'pNAB2' and 'rapa'
#rename sample 18
#also fname is incorrect, remove fname as it is redundant
tu_2_df <-
  tu_df %>% mutate(plasmid = ifelse(id == 18, 'pNAB2', plasmid),
                   rapa = ifelse(id == 18, "rapa", rapa)) %>% select(-c(fname))

#add column for 4tU=17, 18, 19, 20, no4tU=21
tu_2_df$labelling <-
  recode(
    tu_2_df$id,
    "17" =  "4tU",
    "18" =  "4tU",
    "19" =  "4tU",
    "20" =  "4tU",
    "21" =  "no_labelling"
  )

#separate contig into systematic name and gene name+chr locus
tu_3_df <- tu_2_df %>%
  separate(contig, c("name", "locus"), "::")


#read numbers in each replicate
perlib_totals_tu <- tu_3_df  %>%
  group_by(rapa, plasmid, id, rep, labelling) %>%
  summarize(cnt = n())

#read numbers in each condition
percondition_totals_tu <- tu_3_df %>%
  group_by(rapa, plasmid, id, labelling) %>%
  summarize(cnt = n())

#read numbers for each transcript in each condition
percondition_pergene_totals_tu <- tu_3_df %>%
  group_by(name, rapa, plasmid, id, labelling) %>%
  summarize(cnt = n(), median_pA = median(polya_length))



#heatmaps for 4tU RNA pA length distributions
perpA_cnts <- tu_3_df %>%
  mutate(polya_length_bin = 10 * round(polya_length / 10)) %>% #round to 10nt bins
  group_by(name, strain, labelling, plasmid, rapa, id, polya_length_bin) %>%
  summarize(cnt = n()) %>%
  mutate(norm_read_cnt = cnt / sum(cnt))

##select median pA from a given condition for sorting; id=18 --> pNAB2 + rapa
median_pAlength_per_gene <- tu_3_df %>%
  filter(id == 18) %>%
  group_by(name) %>%
  summarize(median_pA = median(polya_length),
            nreads = n())

perpA_cnts %>%
  left_join(., median_pAlength_per_gene) %>%
  filter(nreads > 20) %>%
  ggplot(.,
         aes(
           x = polya_length_bin,
           y = reorder(name, median_pA),
           fill = norm_read_cnt
         )) +
  geom_tile() +
  facet_grid(labelling + rapa + plasmid + id ~ strain,
             scales = 'free',
             space = 'free') +
  scale_fill_gradient2(
    low = 'white',
    mid = 'red',
    high = 'black',
    midpoint = 0.5,
    na.value = 'white'
  ) +
  coord_cartesian(xlim = c(0, 200)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(
      size = 4,
      angle = 0,
      hjust = 0
    )
  )


#violin plots of 4tU RNAs
plot_violin_tu <- function(df,
                           cnt_df = percondition_totals_tu) {
  p <- df %>%
    ggplot(., aes(x = 1, y = polya_length, fill = plasmid)) +
    geom_violin() +
    geom_boxplot(width = .1,
                 fill = 'lightgray',
                 outlier.shape = NA)
  p <- p +
    facet_grid(. ~ rapa + plasmid + labelling, scales = 'free') +
    coord_cartesian(ylim = c(0, 250)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0, 'mm'),
      strip.text.x = element_text(
        angle = 90,
        hjust = 0,
        size = 8
      ),
      strip.background = element_blank(),
      panel.border = element_rect(size = .1),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  p <- p +
    geom_text(
      data = cnt_df,
      aes(
        x = 1,
        y = 250,
        label = paste0('n=', cnt)
      ),
      angle = 90,
      hjust = 'right',
      vjust = 'top'
    )
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  return(p)
}

#all 4tU RNAS
plot_violin_tu(tu_3_df, cnt_df = percondition_totals_tu)

#individual gene violin
plot_violin_for_gene_tu <-
  function(plotting_function = plot_violin_tu,
           gene = 'YJR048W_CYC1',
           ...) {
    p <- plotting_function(filter(tu_3_df, grepl(gene, name)),
                           cnt_df = filter(percondition_pergene_totals_tu, grepl(gene, name)),
                           ...)
    p <- p +
      ggtitle(gene)
    return(p)
  }

plot_violin_for_gene_tu(plot_violin_tu, gene = 'YJR048W_CYC1')


##pPgal plasmids also contain a CYC1 terminator;
# reads assigned to CYC1 are from endogenous CYC1 plus plasmid-derived RNAs
# to separate CYC1 reads into genomic and plasmid-derived transcripts we need the CYC1 read info
# the nanopolish pA length output contains readname but no details on the alignments
# therefore we got detailed alignment info for reads aligned to CYC1 from raw bam files
# samtools view ${bam} chrX:526272-526848 | awk '{print $1"\t"$4"\t"length($10)}' | sort -k2,2n > ${bam/.bam/_CYC1_start_end.coord}
# these were then copied to data/*_genome_mapping_CYC1_start_end.coord files

##load CYC1 read info
##reads selected to overlap CYC1 gene, ie region chrX:526272-526848 got readname, start of alignment and read length from bam file
##Manfred: can you please check this, so that the first line of code can be run:

bam_dir <-
  'data/CYC1_alignments/'
fnames <- dir(bam_dir, '_CYC1_start_end.coord')
(
  cyc1_reads <- lapply(fnames, function(fname)
    read_tsv(
      paste0(bam_dir, fname),
      col_names = c('readname', 'start', 'length')
    ) %>%
      mutate(fname = sub(
        '_genome_mapping_CYC1_start_end.coord', '', fname
      )) %>%
      tidyr::separate(fname, c('strain', 'rapa', 'plasmid', 'id', 'rep'), remove =
                        F)) %>%
    bind_rows
)

saveRDS(cyc1_reads, file = 'data/CYC1_classes.rds')

ggplot(cyc1_reads, aes(x = start)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(526272, 526848))


ggplot(cyc1_reads, aes(x = start)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(526700, 526750))


## short (pPgal-plasmids derived) reads start around 526735-526738
## this is the start of the 3'UTR and we call these CYC1term
## reads from endogenous CYC1 will extend further to 5'
## and we call those CYC1long
(
  cyc1_class <- cyc1_reads %>%
    mutate(
      cyc1_class = case_when(
        start >= 526735 & start <= 526738 ~ 'CYC1term',
        start < 526735 ~ 'CYC1long',
        start > 526738 ~ 'NA'
      )
    ) %>%
    dplyr::select(readname, strain, rapa,  plasmid, id,  rep, cyc1_class)
)


## fix labelling and plasmid, rapa info as explained for tu_2_df
cyc1_class$labelling <-
  recode(
    cyc1_class$id,
    "17" =  "4tU",
    "18" =  "4tU",
    "19" =  "4tU",
    "20" =  "4tU",
    "21" =  "no_labelling"
  )

cyc1_class %<>%
  mutate(plasmid = ifelse(id == 18, 'pNAB2', plasmid),
         rapa = ifelse(id == 18, "rapa", rapa))


(cyc1_pA <- tu_3_df %>%
    dplyr::filter(grepl('CYC1', name)))

unique(cyc1_pA$name)
distinct(cyc1_pA, strain, rapa, plasmid, id, rep, labelling)
distinct(cyc1_class, strain, rapa, plasmid, id, rep, labelling)


cyc1_pA_class <- left_join(cyc1_class, cyc1_pA) %>%
  dplyr::filter(!is.na(polya_length))

(
  cyc1_cnts <- cyc1_pA_class %>%
    group_by(strain, rapa, plasmid, id, cyc1_class, labelling) %>%
    summarize(cnt = n())
)

distinct(cyc1_pA_class, strain, rapa, plasmid, id, labelling)

cyc1_pA_class %>%
  filter(!is.na(cyc1_class), cyc1_class != 'NA') %>%
  mutate(
    xlab = paste0(labelling, '_', rapa, '_', plasmid),
    xlab = factor(
      xlab,
      levels = c(
        'no_labelling_minusrapa_empty',
        '4tU_minusrapa_empty',
        '4tU_minusrapa_pNAB2',
        '4tU_rapa_empty',
        '4tU_rapa_pNAB2'
      )
    )
  ) %>%
  ggplot(., aes(x = xlab, y = polya_length, fill = plasmid)) +
  geom_violin() +
  geom_boxplot(width = .1,
               fill = 'lightgray',
               outlier.shape = NA) +
  ggtitle('CYC1 split into classes') +
  facet_grid(. ~ cyc1_class) +
  coord_cartesian(ylim = c(0, 250)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##defining a group of hyperadenylated transcripts
##filter genes that have <20 reads and median pA higher in rapa vs control (pA>30)
reads_above_twenty <- tu_3_df %>%
  group_by(name, rapa, plasmid, id, labelling) %>%
  summarize(cnt = n(), median_pA = median(polya_length)) %>%
  filter(cnt > 20)

reads_above_twenty_pivoted <- pivot_wider(
  reads_above_twenty,
  id_cols = name,
  names_from = id,
  values_from = median_pA
)

rapa_vs_norapa <- reads_above_twenty_pivoted %>%
  select('17', '19') %>%
  na.omit %>%
  rename(pA_empty_rapa = '17') %>%
  rename(pA_empty_norapa = '19')

genes_with_pA_above <- rapa_vs_norapa %>%
  filter(pA_empty_rapa > pA_empty_norapa) %>%
  filter(pA_empty_rapa > 30)

longer_pA_byrapa_above30As <- tu_3_df %>%
  filter(name %in% genes_with_pA_above$name)

longer_pA_byrapa_above30As %>%
  mutate(
    xlab = paste0(labelling, '_', rapa, '_', plasmid),
    xlab = factor(
      xlab,
      levels = c(
        'no_labelling_minusrapa_empty',
        '4tU_minusrapa_empty',
        '4tU_minusrapa_pNAB2',
        '4tU_rapa_empty',
        '4tU_rapa_pNAB2'
      )
    )
  ) %>%
  ggplot(., aes(x = xlab, y = polya_length, fill = plasmid)) +
  geom_violin() +
  geom_boxplot(width = .1,
               fill = 'lightgray',
               outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 250)) +
    ggtitle('hyperadenylated transcripts only') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
