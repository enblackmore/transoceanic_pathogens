source('functions.R')
library(ggplot2)
set.seed(1492)

#color schemes for plots
theme1a <- c("#BFBBC9", "#4666FF", "#F2003C")
theme1b <- c("#100C08",  "#00CC99", "#FFDD1F")

####################################################################
## Panels 1a--c: Basic Dynamics                                   ##
## In a fully-susceptible population (N=100) with one index case, ##
## using a hypothetical pathogen ke=ki=3, mue=mui=5               ##
## what happens when we vary r0?                                  ##
####################################################################

### 1a-1b simulation (tracking 2 transmission generations) ###

#constant parameters
N_1abc <- 100
S_1abc <- 99
e0_1abc <- 1

mue_1abc <- mui_1abc <- 5
ke_1abc <- ki_1abc <- 3

#varying parameters
r0_1ab <- 10^seq(-1, 2, by=0.01)

#since the model function requires bdd, bfd, and q
#we set q=1 and bdd=r0 / (mue*N) 
#mathematically this is equivalent to q=0 and bfd=r0/mue

bdd_1ab <- r0_1ab/(mue_1abc*N_1abc)
q_1abc <- 1
bfd_1abc <- 0

simulation_results_1ab <-run_analysis(
  N=N_1abc,
  S=S_1abc,
  e0=e0_1abc,
  ke=ke_1abc,
  ki=ki_1abc,
  mue=mue_1abc,
  mui=mui_1abc,
  bdd=bdd_1ab,
  bfd=bfd_1abc,
  q=q_1abc,
  runs=500,
  generation_tracking = TRUE,
  generation_max = 2)

#add r0 column to analysis
simulation_results_1ab$r0 <- simulation_results_1ab$bdd * mui_1abc * N_1abc

#label for plotting
simulation_results_1ab$S <- S_1abc
simulation_results_1ab <- label_outbreaks(
  df=simulation_results_1ab,
  N=N_1abc)
simulation_results_1ab$label[which(simulation_results_1ab$label == "At or above herd immunity")] <- "Below herd immunity"

#Save
saveRDS(simulation_results_1ab, file = "simulation_results/simulation_results_1ab.RDS")

simulation_results_1ab$label[which(simulation_results_1ab$label == "At or above herd immunity" & simulation_results_1ab$r0 == 1)] <- "Below herd immunity"

#get quantiles
simulation_results_1ab_quantiles <- plyr::ddply(simulation_results_1ab, ~bdd+r0, plyr::summarise, 
                                              t05=quantile(Duration, 0.05),
                                              t5=quantile(Duration, 0.5),
                                              t95=quantile(Duration, 0.95),
                                              c05=quantile(Cases,0.05), 
                                              c5=quantile(Cases,0.5),
                                              c95=quantile(Cases,0.95))

### 1a-1b visualisation ###

#Panel 1a: duration by r0 and outbreak ending
panel_1a <- ggplot(simulation_results_1ab) +
  geom_point(mapping=aes(x=r0, y=Duration, col=label), cex=2, alpha=0.01) +
  geom_line(data=simulation_results_1ab_quantiles , mapping=aes(x=r0, y=t5), lwd=1) +
  geom_line(data=simulation_results_1ab_quantiles , mapping=aes(x=r0, y=t05), lty="dashed", lwd=0.5) +
  geom_line(data=simulation_results_1ab_quantiles , mapping=aes(x=r0, y=t95), lty="dashed", lwd=0.5) +
  scale_x_log10() +
  scale_y_continuous() +
  labs(x=bquote("R"[0]), y="Outbreak duration (days)", col="Outbreak size ") +
  theme_bw(); panel_1a

#Panel 1b: outbreak size by r0 and by outbreak ending
panel_1b <- ggplot(simulation_results_1ab) + 
  geom_point(mapping=aes(x=r0, y=Cases, col=label), cex=1, alpha=0.01) +
  geom_line(data=simulation_results_1ab_quantiles , mapping=aes(x=r0, y=c5), lwd=0.75) + 
  geom_line(data=simulation_results_1ab_quantiles , mapping=aes(x=r0, y=c05), lty="dashed", lwd=0.25) +
  geom_line(data=simulation_results_1ab_quantiles , mapping=aes(x=r0, y=c95), lty="dashed", lwd=0.25) + 
  scale_x_log10() +
  theme_bw() +
  labs(y="Outbreak size", x=bquote(R[0])); panel_1b 


### 1c simulation (lower resolution, but tracking 20 generations) ###

r0_1c <- 10^seq(-1, 2, by=0.05)
bdd_1c <- r0_1c/(mue_1abc*N_1abc)

simulation_results_1c <-run_analysis(
  N=N_1abc,
  S=S_1abc,
  e0=e0_1abc,
  ke=ke_1abc,
  ki=ki_1abc,
  mue=mue_1abc,
  mui=mui_1abc,
  bdd=bdd_1c,
  bfd=bfd_1abc,
  q=q_1abc,
  runs=250,
  generation_tracking = TRUE,
  generation_max = 20)

#add r0 column to analysis
simulation_results_1c$r0 <- simulation_results_1c$bdd * mui_1abc * N_1abc
simulation_results_1c$S <- S_1abc
  
#label for plotting
simulation_results_1c <- label_outbreaks(
  df=simulation_results_1c,
  N=N_1abc)

#save
saveRDS(simulation_results_1c, file = "simulation_results/simulation_results_1c.RDS")

#extract quantiles from the data
simulation_results_1c_quantiles <- plyr::ddply(simulation_results_1c, ~bdd+r0, plyr::summarise, 
                                                g05=quantile(Generations, 0.05),
                                                g5=quantile(Generations, 0.5),
                                                g95=quantile(Generations, 0.95))

### 1c visualisation ###

panel_1c <- ggplot(simulation_results_1c) +
  geom_point(mapping=aes(x=r0, y=(Generations-1), col=label), cex=1.5, alpha=0.05) +
  geom_line(data=simulation_results_1c_quantiles, mapping=aes(x=r0, y=g5-1), lwd=0.75) +
  geom_line(data=simulation_results_1c_quantiles, mapping=aes(x=r0, y=g05-1), lty="dashed", lwd=0.25) +
  geom_line(data=simulation_results_1c_quantiles, mapping=aes(x=r0, y=g95-1),  lty="dashed", lwd=0.25) + 
  scale_x_log10() +
  theme_bw() +
  labs(y="Transmission\ngenerations", x=bquote(R[0])); panel_1c

###########################################################
## Figure (1d): more runs, for selected r0 values        ##
## so we can bootstrap introduction risk by journey time ##
###########################################################

### 1d simulation ###

#Constant parameters
N_1d <- 100
S_1d <- 99
e0_1d <- 1
ke_1d <- ki_1d <- 3
mue_1d <- mui_1d <- 5
bfd_1d <- 0
q_1d <- 1

#selected r0 values and corresponding values of bdd
r0_1d <- c(0.5, 2, 8)
bdd_1d <- r0_1d/(mue_1d*N_1d)

simulation_results_1d <- run_analysis(
  N=N_1d,
  S=S_1d,
  e0=e0_1d,
  ke=ke_1d,
  ki=ki_1d,
  mue=mue_1d,
  mui=mui_1d,
  bdd=bdd_1d,
  bfd=bfd_1d,
  q=q_1d,
  runs=5000,
  generation_tracking = FALSE)
saveRDS(simulation_results_1d, file = "simulation_results/simulation_results_1d.RDS")

#bootstrap cumulative introduction risk:
simulation_results_1d$r0 <- simulation_results_1d$bdd * mue_1d * N_1d
simulation_results_1d$r0 <- factor(simulation_results_1d$r0, levels=r0_1d)
tmax_1d <- max(simulation_results_1d$Duration)

introduction_risk_1d <- expand.grid("r0" = unique(simulation_results_1d$r0),
                                    "time" = seq(1, tmax_1d, by=1),
                                    "p_introduction"=NA)

for(i in 1:nrow(introduction_risk_1d)){
  subset_1d <- dplyr::filter(simulation_results_1d, r0 == introduction_risk_1d$r0[i])
  introduction_risk_1d$p_introduction[i] <- length(which(subset_1d$Duration>=introduction_risk_1d$time[i]))/nrow(subset_1d)
}

### 1d visualisation ###

panel_1d <- ggplot(introduction_risk_1d) + geom_line(mapping=aes(x=time, y=p_introduction, col=r0), lwd=2) +
  labs(x="Journey time (days)", y="Introduction risk", col=bquote(R[0])) +
  theme_bw() + scale_color_manual(values=theme1b) +
  theme(legend.title=element_text(size=9), legend.position='right') +
  theme(legend.position = 'none')
panel_1d


##############################################
## Panel 1e: same r0 values as simulation 2 ##
## using mue and mui of real-life pathogens ##
##############################################

### 1e simulation ###

#Constants
N_1e <- 100
S_1e <- 99
e0_1e <- 1
ke_1e <- 3
ki_1e <- 3
q_1e <- 1
bfd_1e <- 0 #as above, since q=1 the value of bfd is irrelevant

#Variables
r0_1e <- c(0.5, 2, 8)

#bdd depends on mue, which varies by pathogen;
#calculate each separately
mue_1e_influenza <- 2
mui_1e_influenza <- 3
bdd_1e_influenza <- r0_1e/(mui_1e_influenza*N_1e)

mue_1e_measles <- 12
mui_1e_measles <- 8
bdd_1e_measles <- r0_1e/(mui_1e_measles*N_1e)

mue_1e_smallpox <- 12
mui_1e_smallpox <- 17.5
bdd_1e_smallpox <- r0_1e/(mui_1e_smallpox*N_1e)

simulation_results_1e_influenza <- run_analysis(
  N=N_1e,
  S=S_1e,
  e0=e0_1e,
  ke=ke_1e,
  ki=ki_1e,
  mue=mue_1e_influenza,
  mui=mui_1e_influenza,
  bdd=bdd_1e_influenza,
  bfd=bfd_1e,
  q=q_1e,
  runs=500,
  generation_tracking=FALSE)

simulation_results_1e_measles <- run_analysis(
  N=N_1e,
  S=S_1e,
  e0=e0_1e,
  ke=ke_1e,
  ki=ki_1e,
  mue=mue_1e_measles,
  mui=mui_1e_measles,
  bdd=bdd_1e_measles,
  bfd=bfd_1e,
  q=q_1e,
  runs=500,
  generation_tracking=FALSE)

simulation_results_1e_smallpox <- run_analysis(
  N=N_1e,
  S=S_1e,
  e0=e0_1e,
  ke=ke_1e,
  ki=ki_1e,
  mue=mue_1e_smallpox,
  mui=mui_1e_smallpox,
  bdd=bdd_1e_smallpox,
  bfd=bfd_1e,
  q=q_1e,
  runs=500,
  generation_tracking=FALSE)

#add pathogen names manually
simulation_results_1e_influenza$Pathogen <-"Influenza"
simulation_results_1e_measles$Pathogen <-"Measles"
simulation_results_1e_smallpox$Pathogen <-"Smallpox"

#combine analysis in a single data frame
simulation_results_1e <- dplyr::full_join(simulation_results_1e_smallpox,
                                         dplyr::full_join(simulation_results_1e_measles, 
                                                          simulation_results_1e_influenza))

#rewrite bdd as r0 for plotting
simulation_results_1e$mui <- NA
simulation_results_1e$mui[which(simulation_results_1e$Pathogen == "Influenza")] <- mui_1e_influenza
simulation_results_1e$mui[which(simulation_results_1e$Pathogen == "Measles")] <- mui_1e_measles
simulation_results_1e$mui[which(simulation_results_1e$Pathogen == "Smallpox")] <- mui_1e_smallpox
simulation_results_1e$r0 <- simulation_results_1e$bdd*simulation_results_1e$mui*N_1e
simulation_results_1e$r0 <- factor(simulation_results_1e$r0, levels= rev(r0_1e))
saveRDS(simulation_results_1e, file = "simulation_results/simulation_results_1e.RDS")

### 1e visualisation ###
simulation_results_1e$Pathogen <- factor(simulation_results_1e$Pathogen, levels=c("Smallpox", "Measles", "Influenza"))
panel_1e <- ggplot(simulation_results_1e, mapping=aes(x=Duration, y=Pathogen, col=r0, group=r0)) +
  geom_point(size=2, alpha=0.05, position=position_jitterdodge(jitter.width = 0.05, dodge.width=0.6)) +
  theme_bw() +
  labs(cx=element_blank(), x="Outbreak duration (days)", y=element_blank(), col=bquote(R[0])) +
  guides(col=guide_legend(reverse=TRUE,
                          override.aes = list(alpha = 1))) +
  scale_y_discrete(labels=rev(c("Influenza<br>\u03bc<sub>E</sub> = 2<br>\u03bc<sub>I</sub> = 3",
                                "Measles<br>\u03bc<sub>E</sub> = 12<br>\u03bc<sub>I</sub> = 8",
                                "Smallpox<br>\u03bc<sub>E</sub> = 12<br>\u03bc<sub>I</sub> = 17.5"))) +
  theme(axis.text.y = ggtext::element_markdown(size=8)) +
  xlim(0,400) +
  scale_color_manual(values=rev(theme1b)); panel_1e


#######################
## Figure 1 Assembly ##
#######################

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

pdf(file = "figures/figure_1.pdf", 
    width = 9, # The width of the plot in inches
    height = 7) # The height of the plot in inches
figure_1_top / figure_1_bottom +
  patchwork::plot_layout(heights=c(1, 0.55)) +
  patchwork::plot_annotation(tag_levels='A')
dev.off()
