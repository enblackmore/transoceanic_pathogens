## Susceptibility Analyses
source('functions.R')
library(ggplot2)
library(RColorBrewer)
set.seed(1492)
library(tidyverse)

#figure s1a

t_s1a <- 29
N_s1a <- 163

mue_s1a_smallpox <- 12
mui_s1a_smallpox <- 17.5

r0_s1a <- 10^seq(0.5, 2, by=0.05)
bfd_s1a <- r0_s1a/mui_s1a_smallpox

S_s1a <- c(0, 1, seq(2, 162, by=2))
N_rep_s1a <- rep(N_s1a, times=length(S_s1a))
t_rep_s1a <- rep(t_s1a, times=length(S_s1a))

simulation_results_s1a <- as.data.frame(matrix(ncol=length(r0_s1a)+1, nrow=length(S_s1a)))
simulation_results_s1a[,1] <- S_s1a

for(i in 1:length(r0_s1a)){
  simulation_results_s1a[,(i+1)] <- get_ship_risk(
    t=t_rep_s1a,
    N=N_rep_s1a,
    S=S_s1a,
    bdd=0,
    bfd=bfd_s1a[i],
    mue=mue_s1a_smallpox,
    mui=mui_s1a_smallpox,
    ke=3,
    ki=3,
    q=0,
    generation_tracking=FALSE,
    runs=250
  ) 
  print(paste(i,"of",length(r0_s1a)))
}

colnames(simulation_results_s1a) <- c("S", r0_s1a)
simulation_results_s1a_plot <- pivot_longer(simulation_results_s1a, 2:ncol(simulation_results_s1a),
                                           names_to="R0", 
                                           values_to="p_introduction")
simulation_results_s1a_plot$R0 <- as.numeric(simulation_results_s1a_plot$R0)

ggplot(simulation_results_s1a_plot) +
  geom_contour_filled(mapping=aes(x=S, y=R0, z=p_introduction),
                                                         breaks=c(0, 0.05, 0.1, 0.2, 0.5, 0.75, 0.90, 1.001)) +
  scale_y_log10() 

##############################

#figure s1b

t_s1b <- 90
N_s1b <- 211

mue_s1b_measles <- 12
mui_s1b_measles <- 8

r0_s1b <- 10^seq(0.5, 2, by=0.05)
bdd_s1b <- 0
bfd_s1b <- r0_s1b / mui_s1b_measles
  
S_s1b <- c(0, 1, seq(2, 210, by=2))
N_rep_s1b <- rep(N_s1b, times=length(S_s1b))
t_rep_s1b <- rep(t_s1b, times=length(S_s1b))

simulation_results_s1b <- as.data.frame(matrix(ncol=length(r0_s1b)+1, nrow=length(S_s1b)))
simulation_results_s1b[,1] <- S_s1b

for(i in 1:length(r0_s1b)){
  simulation_results_s1b[,(i+1)] <- get_ship_risk(
    t=t_rep_s1b,
    N=N_rep_s1b,
    S=S_s1b,
    bdd=0,
    bfd=bfd_s1b[i],
    mue=mue_s1b_measles,
    mui=mui_s1b_measles,
    ke=3,
    ki=3,
    q=0,
    generation_tracking=FALSE,
    runs=250
  ) 
  print(paste(i,"of",length(r0_s1b)))
}

colnames(simulation_results_s1b) <- c("S", r0_s1b)
simulation_results_s1b_plot <- pivot_longer(simulation_results_s1b, 2:ncol(simulation_results_s1b),
                                           names_to="R0", 
                                           values_to="p_introduction")
simulation_results_s1b_plot$R0 <- as.numeric(simulation_results_s1b_plot$R0)

simulation_results_s1a_plot$ship <- 1
simulation_results_s1b_plot$ship <- 2

simulation_results_all <- full_join(simulation_results_s1a_plot,
                                    simulation_results_s1b_plot)


simulation_results_all$ship <- factor(simulation_results_all$ship,
                  levels=c(1,2),
                  labels=c("Gold Hunter smallpox introduction\n 161 passengers, 29 days",
                           "Sir Charles Napier measles extinction\n211 passengers, 90 days"))

sim_lines <- tibble(ship=unique(simulation_results_all$ship), x=NA, y0=NA, y05=NA, y1=NA,
                    xlabel=NA, y0label=NA, y05label=NA, y1label=NA)
sim_lines$x[1] <- 8
sim_lines$x[2] <- 11

sim_lines$y0[1] <- 7
sim_lines$y05[1] <- ((163*(bfd_s1a_smallpox/75))^0.5)*((bfd_s1_smallpox)^0.5)*17.5
sim_lines$y1[1] <- ((163*bfd_s1a_smallpox/75))*17.5
sim_lines$y0label[1] <- "q=0: R0=7"
sim_lines$y05label[1] <- "q=0.5: R0=10.3"
sim_lines$y1label[1] <- "q=1: R0=15.2"
sim_lines$xlabel[1] <- "S(0)=8"

sim_lines$y0[2] <- 15
sim_lines$y05[2] <- ((211*(bfd_s1b_measles/75))^0.5)*((15/8)^0.5)*8
sim_lines$y1[2] <- 211*(bfd_s1b_measles/75)*8
sim_lines$y0label[2] <- "q=0: R0=15"
sim_lines$y05label[2] <- "q=0.5: R0=25.1"
sim_lines$y1label[2] <- "q=1: R0=42.2"
sim_lines$xlabel[2] <- "S(0)=11"

ggplot(simulation_results_all) + 
  geom_contour_filled(mapping=aes(x=S, y=R0, z=p_introduction),
                      breaks=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.95, 1.0001)) +
  labs(y=bquote("R"[0]), x="Initial number of susceptible people, S(0)",
       fill="Introduction\nRisk") +
  geom_vline(data=sim_lines, mapping=aes(xintercept=x), lty="dashed", lwd=0.5) +
  geom_hline(data=sim_lines, mapping=aes(yintercept=y0), lty="dashed", lwd=0.5) +
  geom_hline(data=sim_lines, mapping=aes(yintercept=y05), lty="dashed", lwd=0.5) +
  geom_hline(data=sim_lines, mapping=aes(yintercept=y1), lty="dashed", lwd=0.5) +
  geom_label(data=sim_lines, mapping=aes(label=xlabel, x=x+2.5, y=115),
            cex=3,
            hjust=0) +
  geom_label(data=sim_lines, mapping=aes(label=y0label, x=100, y=y0*0.9),
            cex=3,
            hjust=0) +
  geom_label(data=sim_lines, mapping=aes(label=y05label, x=100, y=y05*1.1),
             cex=3,
             hjust=0) +
  geom_label(data=sim_lines, mapping=aes(label=y1label, x=100, y=y1*1.1),
            cex=3,
            hjust=0) +
  guides(fill=guide_colorsteps(show.limits=TRUE, frame.colour="black", frame.linewidth=0.5, frame.linetype=1)) +
  scale_y_log10() +
  scale_fill_brewer(type="seq", palette="Blues", direction = 1) +
  facet_wrap(vars(ship)) + theme_bw() 

write_csv(simulation_results_all, 'sim_results_s1ab.csv')



