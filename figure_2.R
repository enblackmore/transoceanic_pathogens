#Figure 2: Incorporating POpulation Size and Susceptibility
source(functions.R)
set.seed(1492)

#colour schemes
theme2a <- c("#BFBBC9", "#4666FF", "#F2003C")
theme2b <- c("#4C4C47", "#848FA5", "#F4B886", "#8FCB9B")

###################################################
#Panel 2a: same r0 values as simulation 2
#constant re0, with varying S
#across three re0 values
#corresponds with figure 2a
###################################################

#### PANEL 2a SIMULATION ###

#constants
N_2a <- 1001
S_2a <- unique(round(10^seq(0, 3, by=0.1),0))
e0_2a <- 1
ke_2a <- ki_2a <- 3
mue_2a <- mui_2a <- 5
q_2a <- 1
bfd_2a <- 0

#variables
re0_2a <- c(1, 2, 8)

#simulation parameters
runs_2a <- 500

#run simulation with generation tracking, and
#set generation_max=2, because we only need to distinguish 
#between single-generation outbreaks and everything else
generation_tracking_2a <- TRUE
generation_max_2a <- 2 

#we maintain constant re0 across different values of S
#by setting q=1
#and back-calculating bdd for each analysis:
#bdd = re0_2a/(mui_2a*S_2a) 

mui_2a <- 5
bdd_2a_re0_1 <- re0_2a[1]/(mui_2a*S_2a)
bdd_2a_re0_2 <- re0_2a[2]/(mui_2a*S_2a)
bdd_2a_re0_8 <- re0_2a[3]/(mui_2a*S_2a)

simulation_results_2a_re0_1 <- run.analysis2(
  N=N_2a,
  S=S_2a,
  e0=e0_2a,
  ke=ke_2a,
  ki=ki_2a,
  mue=mue_2a,
  mui=mui_2a,
  bdd=bdd_2a_re0_1,
  bfd=bfd_2a,
  q=q_2a,
  runs=runs_2a,
  generation_tracking = generation_tracking_2a,
  generation_max = generation_max_2a
)

simulation_results_2a_re0_2 <- run.analysis2(
  N=N_2a,
  S=S_2a,
  e0=e0_2a,
  ke=ke_2a,
  ki=ki_2a,
  mue=mue_2a,
  mui=mui_2a,
  bdd=bdd_2a_re0_2,
  bfd=bfd_2a,
  q=q_2a,
  runs=runs_2a,
  generation_tracking = generation_tracking_2a,
  generation_max = generation_max_2a
)

simulation_results_2a_re0_8 <- run.analysis2(
  N=N_2a,
  S=S_2a,
  e0=e0_2a,
  ke=ke_2a,
  ki=ki_2a,
  mue=mue_2a,
  mui=mui_2a,
  bdd=bdd_2a_re0_8,
  bfd=bfd_2a,
  q=q_2a,
  runs=runs_2a,
  generation_tracking = generation_tracking_2a,
  generation_max = generation_max_2a
)

#manually add re0 values
simulation_results_2a_re0_1$analysis$re0 <- 1
simulation_results_2a_re0_2$analysis$re0 <- 2
simulation_results_2a_re0_8$analysis$re0 <- 8

#combine analyses in a single data frame
simulation_results_2a <- dplyr::full_join(simulation_results_2a_re0_1$analysis,
                                         dplyr::full_join(simulation_results_2a_re0_2$analysis,
                                                          simulation_results_2a_re0_8$analysis))
#manually add r0 values, to identify whether outbreaks reach herd immunity
#r0 = n*re0 / s
simulation_results_2a$r0 <- N_2a*as.numeric(simulation_results_2a$re0) / simulation_results_2a$S

#label for plotting
simulation_results_2a <- label.outbreaks(df=simulation_results_2a, N=N_2a, r0=simulation_results_2a$r0)

#save output
saveRDS(simulation_results_2a, file = "simulation_results_2a.RDS")

### PANEL 2a VISUALISATION ###

#get quantiles for plotting
simulation_results_2a_quantiles <- plyr::ddply(simulation_results_2a, ~re0+S, plyr::summarise,
                                              t5=quantile(Duration, 0.5),
                                              t05=quantile(Duration, 0.05), 
                                              t95=quantile(Duration, 0.95))

#discretise re0
simulation_results_2a_quantiles$re0 <- factor(simulation_results_2a_quantiles$re0, 
                                             levels=unique(simulation_results_4_quantiles$re0))
simulation_results_2a$re0 <- factor(simulation_results_2a$re0, 
                                   levels=unique(simulation_results_2a$re0))

panel_2a <- ggplot(simulation_results_2a) + 
  geom_point(mapping=aes(x=S, y=Duration, col=label), alpha=0.02) +
  facet_wrap(vars(re0), nrow=1, label=label_bquote(R[e](0)==~.(as.character(re0)))) +
  geom_line(data=simulation_results_2a_quantiles, mapping=aes(x=S, y=t5), lwd=0.75) +
  geom_line(data=simulation_results_2a_quantiles, mapping=aes(x=S, y=t05), lwd=0.5, lty="dashed") +
  geom_line(data=simulation_results_2a_quantiles, mapping=aes(x=S, y=t95), lty="dashed") +
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

#############################################################
#Panel (2b): introducing density and frequency dependence
#varying S/N
#for three values of q
#for three values of bfd (q=0)
#for three values of bdd (q=1)
#and for three values of bdd with constant bfd (q=0.5)
#corresponds with figure 2b
#############################################################

### PANEL 2b SIMULATION ###

#fixed values
e0_2b <- 1
mue_2b <- 5
mui_2b <- 5
ke_2b <- 3
ki_2b <- 3

#variables
N_2b <- c(50, 100, 200, 500)
proportion_S_2b <- 10^seq(-2, -0.05, by=0.05)

#case 1: q=0
q_2b_q0 <- 0
bfd_2b_q0 <- c(0.4, 1, 2)
bdd_2b_q0 <- 0

input_2b_q0 <- expand.grid(N=N_2b, pp=proportion_S_2b, bfd=bfd_2b_q0)
input_2b_q0$S <- round(input_2b_q0$N * input_2b_q0$pp, 0)
any(input_2b_q0$S > input_2b_q0$N - e0_2b)

#case 2: q=0.5
q_2b_q05 <- 0.5
bfd_2b_q05 <- 1
bdd_2b_q05 <- c(0.05, 0.1, 0.2)

input_2b_q05 <- expand.grid(N=N_2b, pp=proportion_S_2b, bdd=bdd_2b_q05)
input_2b_q05$S <- round(input_2b_q05$N * input_2b_q05$pp, 0)

#case 3: q=1
q_2b_q1 <- 1
bfd_2b_q1 <- 0
bdd_2b_q1 <- c(0.05, 0.1, 0.2)

input_2b_q1 <- expand.grid(N=N_2b, pp=proportion_S_2b, bdd=bdd_2b_q1)
input_2b_q1$S <- round(input_2b_q1$N * input_2b_q1$pp, 0)

#simulation parameters
runs_2b <- 200

#run analysis
simulation_results_2b_q0 <- run.analysis2(
  N=input_2b_q0$N,
  S=input_2b_q0$S,
  e0=e0_2b,
  ke=ke_2b,
  ki=ki_2b,
  mue=mue_2b,
  mui=mui_2b,
  bdd=bdd_2b_q0,
  bfd=input_2b_q0$bfd,
  q=q_2b_q0,
  runs=runs_2b,
  generation_tracking=FALSE
)
simulation_2b_q0_analysis <- as.data.frame(simulation_results_2b_q0$analysis)
simulation_2b_q0_analysis$label <- factor(simulation_2b_q0_analysis$bfd, 
                                         levels=unique(input_2b_q0$bfd),
                                         labels=c(1:3))

simulation_results_2b_q05 <- run.analysis2(
  N=input_2b_q05$N,
  S=input_2b_q05$S,
  e0=e0_2b,
  ke=ke_2b,
  ki=ki_2b,
  mue=mue_2b,
  mui=mui_2b,
  bdd=input_2b_q05$bdd,
  bfd=bfd_2b_q05,
  q=q_2b_q05,
  runs=runs_2b,
  generation_tracking=FALSE
)
simulation_2b_q05_analysis <- as.data.frame(simulation_results_2b_q05$analysis)
simulation_2b_q05_analysis$label <- factor(simulation_2b_q05_analysis$bdd, 
                                          levels=unique(input_2b_q05$bdd),
                                          labels=c(4:6))

simulation_results_2b_q1 <- run.analysis2(
  N=input_2b_q1$N,
  S=input_2b_q1$S,
  e0=e0_2b,
  ke=ke_2b,
  ki=ki_2b,
  mue=mue_2b,
  mui=mui_2b,
  bdd=input_2b_q1$bdd,
  bfd=bfd_2b_q1,
  q=q_2b_q1,
  runs=runs_2b,
  generation_tracking=FALSE
)
simulation_2b_q1_analysis <- as.data.frame(simulation_results_2b_q1$analysis)
simulation_2b_q1_analysis$label <- factor(simulation_2b_q1_analysis$bdd, 
                                          levels=unique(input_2b_q1$bdd),
                                          labels=c(7:9))

simulation_results_2b <- dplyr::full_join(simulation_2b_q0_analysis,
                                   dplyr::full_join(simulation_2b_q05_analysis,
                                                    simulation_2b_q1_analysis))
simulation_results_2b$Susceptible <- simulation_results_2b$S / simulation_results_2b$N
simulation_results_2b$label <- as.character(simulation_results_2b$label)

#get median outbreak durations
simulation_results_2b_median <- plyr::ddply(simulation_results_2b, ~N+Susceptible+bdd+bfd+label, plyr::summarise,
                                           median_duration = median(Duration))

#add label for plotting
labels_2b <- c(
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


simulation_results_2b_median$label <- factor(simulation_results_2b_median$label,
                                             levels=c(1:9),
                                             labels=labels_2b)

#discretise N
simulation_results_2b_median$N <- factor(simulation_results_2b_median$N, 
                                        levels=unique(simulation_results_2b_median$N))

#identify critical S/N values for each analysis
simulation_2b_critical_susceptibility <- dplyr::distinct(dplyr::select(simulation_results_2b_median, N, r0, label))

#re0 = S/N r0
#so when re0=1, S/N = 1/r0
simulation_2b_critical_susceptibility$Susceptibility <- 1/simulation_2b_critical_susceptibility$r0

### PANEL 2b VISUALIZATION ###

panel_2b <- ggplot(simulation_results_2b_median) +
  geom_line(mapping=aes(x=Susceptible, y=median_duration, col=N), lwd=1.5) +
  facet_wrap(vars(label), labeller=label_parsed) +
  scale_x_log10() +
  labs(x="Initial proportion susceptible", y="Median outbreak duration (days)", col="") +
  geom_vline(data=simulation_2b_critical_susceptibility, mapping=aes(xintercept=Susceptibility, col=N), lty="dashed") +
  theme_bw() +
  theme(axis.text.x=element_text(hjust=c(0,0.5,1)),
        legend.background = element_rect(fill=NA),
        aspect.ratio=1/2) +
  scale_color_manual(values=theme2b,
                     labels=paste("N =", c(50, 100, 200, 500)))

#########################
### Figure 2 Assembly ###
#########################

legend_panel_2b <- cowplot::get_legend(panel_2b + theme(legend.justification = c(1,0.5), 
                                                        legend.box.margin = margin(b=-10, r=25, unit='pt')))

figure_2 <- (panel_2a + legend_panel_2b + patchwork::plot_layout(widths=c(0.75, 0.25))) / 
  patchwork::wrap_elements(full=panel_2b + guides(col='none')) + 
  patchwork::plot_layout(heights=c(0.6,3.2), guides=
                           'collect') +
  patchwork::plot_annotation(tag_levels=list(c('A', '', 'B'))); figure_2
