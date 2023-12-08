#Visualisations
library(ggplot2)

################################
#Figure 1: Basic Model Dynamics
################################

#colour schemes
theme1a <- c("#BFBBC9", "#4666FF", "#F2003C")
theme1b <- c("#100C08",  "#00CC99", "#FFDD1F")

#extract quantiles from the data
simulation_results_1_quantiles <- plyr::ddply(simulation_results_1$analysis, ~bdd+r0, summarise, 
                      t05=quantile(Duration, 0.05),
                      t5=quantile(Duration, 0.5),
                      t95=quantile(Duration, 0.95),
                      c05=quantile(Cases,0.05), 
                      c5=quantile(Cases,0.5),
                      c95=quantile(Cases,0.95),
                      g05=quantile(Generations, 0.05),
                      g5=quantile(Generations, 0.5),
                      g95=quantile(Generations, 0.95))

#Panel 1a: duration by r0 and outbreak ending
panel_1a <- ggplot(simulation_results_1$analysis) +
  geom_point(mapping=aes(x=r0, y=Duration, col=label), cex=2, alpha=0.05) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=t5), lwd=1) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=t05), lty="dashed", lwd=0.5) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=t95), lty="dashed", lwd=0.5) +
  scale_x_log10() +
  scale_y_continuous() +
  labs(x=bquote("R"[0]), y="Outbreak duration (days)", col="Outbreak size ") + theme_bw(); panel_1a

#Panel 1b: outbreak size by r0 and by outbreak ending
panel_1b <- ggplot(simulation_results_1$analysis) + 
  geom_point(mapping=aes(x=r0, y=Cases, col=label), cex=1, alpha=0.05) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=c5), lwd=0.75) + 
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=c05), lty="dashed", lwd=0.25) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=c95), lty="dashed", lwd=0.25) + 
  scale_x_log10() +
  theme_bw() + labs(y="Outbreak size", x=bquote(R[0])); panel_1b 

#Panel 1c: outbreak generation number by r0 and by outbreak ending
panel_1c <- ggplot(simulation_results_1$analysis) +
  geom_point(mapping=aes(x=r0, y=(Generations-1), col=label), cex=1.5, alpha=0.01) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=g5-1), lwd=0.75) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=g05-1), lty="dashed", lwd=0.25) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=g95-1),  lty="dashed", lwd=0.25) + 
  scale_x_log10() +
  theme_bw() + labs(y="Transmission\ngenerations", x=bquote(R[0])); panel_1c

#Panel 1d: introduction risk by journey time and by r0
introduction_risk_2$r0 <- factor(introduction_risk_2$r0, levels=unique(introduction_risk_2$r0))
panel_1d <- ggplot(introduction_risk_2) + geom_line(mapping=aes(x=time, y=p_introduction, col=r0), lwd=2) +
  labs(x="Journey time (days)", y="Introduction risk", col=bquote(R[0])) +
  theme_bw() + scale_color_manual(values=theme1b) +
  theme(legend.title=element_text(size=9), legend.position='right') +
  theme(legend.position = 'none')
panel_1d

#Panel 1e: outbreak duration by r0 and by pathogen
simulation_results_3$r0 <- factor(simulation_results_3$r0, levels= rev(r0_3))

panel_1e <- ggplot(simulation_results_3, mapping=aes(x=Duration, y=Pathogen, col=r0, group=r0)) +
  geom_point(size=2, alpha=0.05, position=position_jitterdodge(jitter.width = 0.05, dodge.width=0.6)) +
  theme_bw() +
  labs(cx=element_blank(), x="Outbreak duration (days)", y=element_blank(), col=bquote(R[0])) +
  guides(col=guide_legend(reverse=TRUE,
                          override.aes = list(alpha = 1))) +
  scale_y_discrete(labels=rev(c("Influenza<br>\u03bc<sub>E</sub> = 2<br>\u03bc<sub>I</sub> = 3",
                                "Measles<br>\u03bc<sub>E</sub> = 12<br>\u03bc<sub>I</sub> = 8",
                                "Smallpox<br>\u03bc<sub>E</sub> = 12<br>\u03bc<sub>I</sub> = 17.5"))) +
  theme(axis.text.y = ggtext::element_markdown(size=8)) +
  scale_color_manual(values=theme1b); panel_1e

#Figure 1 assembly

#make the top half of the figure
layout_1_top <- "
AB
AC
"

figure_1_top <- patchwork::wrap_plots(
                      A = panel_1a + 
                        theme(legend.title=element_text(size=9, hjust=1, vjust=1), 
                              legend.position = c(0.99, 0.995), 
                              legend.direction='vertical',
                              legend.justification=c(1,1), 
                              legend.text.align = 0,
                              legend.key.size=unit(14, 'pt'),
                              legend.text=element_text(size=7.5, margin=margin(t=3, b=3, unit='pt'))) +
                        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                        scale_color_manual(values=theme1a, labels=c("Single-generation", "Below herd\nimmunity",
                                                                    "At or above\nherd immunity")),
                      B = panel_1b + 
                        theme(axis.title.y = element_text(size=9, hjust=0.25, margin=margin(r=-5, unit='pt')), legend.position="none") +
                        scale_color_manual(values=theme1a),
                      C = panel_1c + 
                        theme(axis.title.y = element_text(size=9, hjust=0.25, margin=margin(r=-5, unit='pt'))) +
                        scale_color_manual(values=theme1a) +
                        theme(legend.position="none"),
                      design = layout_1_top,
                      widths=c(0.7, 0.3)); figure_1_top

layout_1_bottom <- "
DE
"
figure_1_bottom <- patchwork::wrap_plots(
  D = panel_1d + 
    theme(legend.position = c(0.99, 0.99),legend.justification = c(1,1)),
  E = panel_1e +
    theme(legend.position='none'),
  design = layout_1_bottom,
  widths=c(0.6, 0.4)); figure_1_bottom 

figure_1 <- figure_1_top / figure_1_bottom +
  patchwork::plot_layout(heights=c(1, 0.55)) +
  patchwork::plot_annotation(tag_levels='A')

################################################################
#Figure 2: Incorporating Population Size and Susceptibility
################################################################

#Figure 2a:
#varying S, with constant re0

#colour schemes
theme2a <- c("#BFBBC9", "#4666FF", "#F2003C")
theme2b <- c("#4C4C47", "#848FA5", "#F4B886", "#8FCB9B")

#get quantiles for plotting
simulation_results_4_quantiles <- plyr::ddply(simulation_results_4, ~re0+S, summarise, 
                                              t5=median(Duration),
                                              t05=quantile(Duration, 0.05), 
                                              t95=quantile(Duration, 0.95))

#discretise re0
simulation_results_4_quantiles$re0 <- factor(simulation_results_4_quantiles$re0, 
                                             levels=unique(simulation_results_4_quantiles$re0))
simulation_results_4$re0 <- factor(simulation_results_4$re0, 
                                   levels=unique(simulation_results_4$re0))

panel_2a <- ggplot(simulation_results_4) + 
  geom_point(mapping=aes(x=S, y=Duration, col=label), alpha=0.02) +
  facet_wrap(vars(re0), nrow=1, label=label_bquote(R[e](0)==~.(as.character(re0)))) +
  geom_line(data=simulation_results_4_quantiles, mapping=aes(x=S, y=t5), lwd=0.75) +
  geom_line(data=simulation_results_4_quantiles, mapping=aes(x=S, y=t05), lwd=0.5, lty="dashed") +
  geom_line(data=simulation_results_4_quantiles, mapping=aes(x=S, y=t95), lty="dashed") +
  scale_x_log10() + 
  ylim(0,125) +
  theme_bw() +
  scale_color_manual(values=theme2a, 
                     labels=c("Single -\ngeneration",
                              "Below herd\nimmunity",
                              "At or above\nherd immunity")) +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11, hjust=0.25),
        legend.text = element_text(margin=margin(t=0.2, unit='cm'), vjust=0.5),
        axis.text.x=element_text(angle=30, hjust=1, vjust=1),
        legend.position='none', aspect.ratio = 0.7/1,
        plot.margin = margin(b=10, unit='pt')) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(x="Initial number susceptible",
       y="Outbreak\nduration (days)",
       col="Outbreak\nsize"); panel_2a

#Figure 2b:
#median outbreak duration by N, S/N, and density-dependence

#get median outbreak durations
simulation_results_5_median <- plyr::ddply(simulation_results_5$analysis, ~N+Susceptible+bdd+bfd+q+r0, summarise,
                                           median_duration = median(Duration))

#add label for plotting
labels_figure_2b <- c(
  expression(R[0]==2),
  expression(R[0]==5),
  expression(R[0]==10),
  expression(R[0]==0.5*sqrt('N')),
  expression(R[0]==sqrt(0.5)*sqrt('N')),
  expression(R[0]==sqrt('N')),
  expression(R[0]==0.05*'N'),
  expression(R[0]==0.1*'N'),
  expression(R[0]==0.2*'N')
)

simulation_5_q_labels <- dplyr::distinct(dplyr::select(simulation_results_5_median, bdd, bfd, q))
simulation_5_q_labels$label <- labels_figure_2b
simulation_results_5_median <- dplyr::left_join(simulation_results_5_median, simulation_5_q_labels)
simulation_results_5_median$label <- factor(simulation_results_5_median$label,
                                            levels=labels_figure_2b)

#discretise N
simulation_results_5_median$N <- factor(simulation_results_5_median$N, 
                                        levels=unique(simulation_results_5_median$N))

#identify critical S/N values for each analysis
simulation_5_critical_susceptibility <- dplyr::distinct(dplyr::select(simulation_results_5_median, N, r0, label))

#re0 = S/N r0
#so when re0=1, S/N = 1/r0
simulation_5_critical_susceptibility$Susceptibility <- 1/simulation_5_critical_susceptibility$r0

panel_2b <- ggplot(simulation_results_5_median) +
  geom_line(mapping=aes(x=Susceptible, y=median_duration, col=N), lwd=1.5) +
  facet_wrap(vars(label), labeller=label_parsed) +
  scale_x_log10() +
  labs(x="Initial proportion susceptible", y="Median outbreak duration (days)", col="") +
  geom_vline(data=simulation_5_critical_susceptibility, mapping=aes(xintercept=Susceptibility, col=N), lty="dashed") +
  theme_bw() +
  theme(axis.text.x=element_text(hjust=c(0,0.5,1)),
        legend.background = element_rect(fill=NA),
        aspect.ratio=1/2) +
  scale_color_manual(values=theme2b,
                     labels=paste("N =", c(50, 100, 200, 500)))


#Figure 2 assembly

legend_panel_2b <- cowplot::get_legend(panel_2b + theme(legend.justification = c(1,0.5), 
                                           legend.box.margin = margin(b=-10, r=25, unit='pt')))

figure_2 <- (panel_2a + legend_panel_2b + patchwork::plot_layout(widths=c(0.75, 0.25))) / 
  patchwork::wrap_elements(full=panel_2b + guides(col='none')) + 
  patchwork::plot_layout(heights=c(0.6,3.2), guides=
                'collect') +
  patchwork::plot_annotation(tag_levels=list(c('A', '', 'B'))); figure_2


######################################
#Figure 3: San Francisco Port Arrivals
######################################

theme3  <- c( "#07f49e", "#42047e")

#Figure 3b: journey summaries 
data_SF <- read.csv("San_Francisco_arrivals.csv")

Longitude_order <- c(1,7,5,2,6,8,3,4)
data_SF$From_code <- factor(data_SF$From_code, levels=unique(data_SF$From_code)[Longitude_order],
                            labels=c("Liverpool", "New York City", "Panama", "Valparaíso",
                                     "Oregon", "Hawai'i", "Sydney", "Hong Kong"))

#summarise data; this becomes supplementary table 2
data_SF_summary <- plyr::ddply(data_SF,
                                ~From_code+Steam, summarise,
                                n=length(Voyage_days),
                                median=median(Voyage_days, na.rm=TRUE),
                                dmin=min(Voyage_days, na.rm=TRUE),
                                dmax=max(Voyage_days, na.rm=TRUE),
                                d1=quantile(Voyage_days, 0.1, na.rm=TRUE),
                                d9=quantile(Voyage_days, 0.9, na.rm=TRUE),
                                p_median=median(N_passengers, na.rm=TRUE),
                                pmin=min(N_passengers, na.rm=TRUE),
                                pmax=max(N_passengers, na.rm=TRUE),
                                p1=quantile(N_passengers, 0.1, na.rm=TRUE),
                                p9 = quantile(N_passengers, 0.9, na.rm=TRUE))



panel_3a <- ggplot(data_SF) +
  geom_jitter(mapping=aes(x=Voyage_days, y=From_code, col=Steam), height=0.1, alpha=0.3) +
  theme_bw() +
  labs(x="Journey Time (days)", y="Origin Port", col="Technology") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(breaks=seq(0, 250, by=50), limits = c(0,250)) +
  scale_color_manual(values=theme3, labels=c("Steam", "Sail")) +
  theme(axis.title.y=element_text(margin=margin(r=-10, unit='pt')),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_text(size=10.2),
        legend.position='none',
        aspect.ratio=1.5/1)

panel_3b <- ggplot(data_SF) + 
  geom_jitter(mapping=aes(x=N_passengers, y=From_code, col=Steam), height=0.1, alpha=0.3) +
  theme_bw() +
  labs(x="Passengers", y="", col="Propulsion\ntype") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_color_manual(values=theme3, labels=c("Sail", "Steam")) +
  theme(axis.text.y = element_blank(),
        axis.title.x=element_text(margin=margin(t=10, unit='pt')),
        legend.position="none",
        aspect.ratio = 1.5/1) + scale_x_log10()

panel_3c <- ggplot(data_SF) +
  geom_bar(mapping=aes(x=From_code, fill=Steam), show.legend = FALSE) +
  geom_point(mapping=aes(x=From_code, y=1, col=Steam), alpha=0) +
  theme_bw() +
  labs(x="Origin Port", y="Arrivals", col="") +
  scale_fill_manual(values=theme3, labels=c("Steam", "Sail")) +
  scale_color_manual(values=theme3, labels=c("Steam", "Sail"),
                     guide=guide_legend(override.aes = list(alpha=1, size=4))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(margin=margin(t=-5, unit='pt')),
        axis.title.y=element_text(margin=margin(r=-5, b=-10, unit='pt')),
        legend.position=c(0.75, 0.90),
        legend.text = element_text(margin=margin(l=-4, unit='pt')),
        legend.background=element_rect(fill=NA),
        aspect.ratio = 1.25/1)
