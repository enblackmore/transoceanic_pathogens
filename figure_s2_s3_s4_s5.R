#figure 4 sensitivity analysis

### Figure S2: Influenza

S_proportion_s2 <- 0.05
e0_s2 <- 1

N_s2 <- unique(round(10^seq(1,3.2,by=0.01),0))
S_s2 <- round(S_proportion_s2*N_s2,0)
any(N_s2 < S_s2 + e0_s2) #sanity check

ke_s2 <- 3
ki_s2 <- 3

mue_s2_influenza <- 2
mui_s2_influenza <- 3
bfd_s2_influenza <- 1.5/mui_s2_influenza

q_s2 <- c(0, 0.25, 0.5, 0.75, 1)
c_s2 <- c(50, 100, 150)

runs <- 500

parm_combinations <- expand.grid(q=q_s2, c=c_s2)
results <- array(NA, dim=c((length(N_s2)*runs), 7, length(q_s2)*length(c_s2)))

for(i in 1:nrow(parm_combinations)){
  q_i <-parm_combinations[i,1]
  c_i <- parm_combinations[i,2]
  results <- run_analysis2(
    N=N_s2,
    S=S_s2,
    e0=e0_s2,
    ke=ke_s2,
    ki=ki_s2,
    mue=mue_s2_influenza,
    mui=mui_s2_influenza,
    bdd=(bfd_s2_influenza/c_i),
    bfd=bfd_s2_influenza,
    q=q_i,
    runs=runs,
    generation_tracking=FALSE)
  results$q <- q_i
  results$c <- c_i
  assign(paste0("results_", i), results)
}

df <- full_join(results_1, results_2)
for(i in 3:nrow(parm_combinations)){
  df <- full_join(df, get(paste0("results_", i)))
}

simulation_results_s2 <- df

#bootstrap introduction risk as a function of pathogen, N and journey time
#identify longest outbreak duration in dataset
time_max <- 250

#list of times to assess introduction risk
times <- seq(0, time_max, by=1)

introduction_risk_s2 <- expand.grid(N=N_s2, time=times, q=q_s2, c=c_s2, p_introduction=NA)

for(i in 1:nrow(parm_combinations)){
  q_i <- parm_combinations$q[i]
  c_i <- parm_combinations$c[i]
  for(j in 1:length(N_s2)){
    subset_2 <- intersect(intersect(which(simulation_results_s2$q==q_i), which(simulation_results_s2$c==c_i)),
                          which(simulation_results_s2$N==N_s2[j]))
    subset_intro <- intersect(intersect(which(introduction_risk_s2$N == N_s2[j]),
                                        which(introduction_risk_s2$q == q_i)), which(introduction_risk_s2$c==c_i))
    vec_j <- numeric(length(times))
    for(k in 1:length(times)){
      vec_j[k] <- length(which(simulation_results_s2$Duration[subset_2] > times[k]))/length(subset_2)
    }
    introduction_risk_s2$p_introduction[subset_intro] <- vec_j
  }
  print(paste(i, "of", nrow(parm_combinations)))
}

saveRDS(introduction_risk_s2, file="simulation_results/introduction_risk_s2.RDS")

### Figure S2 Visualisation ###
theme4_shape <- c(22,23,25, 8, 4,21,10,24)
introduction_risk_s2$p_introduction[which(introduction_risk_s2$p_introduction == 1)] <- 0.99

introduction_risk_s2$q <- factor(introduction_risk_s2$q, levels=q_s2, 
                                 labels=c("q = 0", "q = 0.25", "q = 0.5", "q = 0.75", "q = 1"))

introduction_risk_s2$c <- factor(introduction_risk_s2$c, levels=c_s2, 
                                 labels=c("c = 50", "c = 100", "c = 150"))

figure_s2_contours <- ggplot(introduction_risk_s2) +
  geom_contour_filled(bins=6, mapping=aes(y=N, x=time, z=p_introduction), breaks=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1)) + 
  labs(x="Journey time (days)", y="Total population (N)", fill="Introduction Risk") + 
  facet_grid(q~c) +
  theme_bw() +
  scale_fill_manual(values=alpha(grey.colors(6, start=1, end=0.2), c(0,1,1,1,1,1))) +
  xlim(0,200)

#Panel S2: overplotting selected historical voyages
data_panel_s2 <- read.csv("data/Selected_voyages.csv")

#factor column 'Ship' so ships appear in chronological order
data_panel_s2$Ship <- factor(data_panel_s2$Ship,
                              levels=data_panel_s2$Ship)

figure_s2_influenza <- figure_s2_contours +
  geom_point(data=data_panel_s2, mapping=aes(x=t, y=N, pch=Ship), col="black", fill="white",  cex=3) + 
  theme(legend.text = element_text(margin = margin(b=2, t=2, unit='pt')), legend.justification='top') +
  guides(fill=guide_colorsteps(show.limits=TRUE, frame.colour="black", frame.linewidth=0.5, frame.linetype=1)) +
  labs(col="Voyage") +
  scale_shape_manual(values=theme4_shape) +
  scale_y_log10(limits=c(10,1500)); figure_s2_influenza


pdf(file = "figures/figure_s2.pdf", 
    width = 9.32, # The width of the plot in inches
    height = 6.21) # The height of the plot in inches
figure_s2_influenza
dev.off()




  
### S3: measles

S_proportion_s3 <- 0.05
e0_s3 <- 1

N_s3 <- unique(round(10^seq(1,3.2,by=0.01),0))
S_s3 <- round(S_proportion_s3*N_s3,0)
any(N_s3 < S_s3 + e0_s3) #sanity check

ke_s3 <- 3
ki_s3 <- 3

mue_s3_measles <- 12
mui_s3_measles <- 8
bfd_s3_measles <- 15/mui_s3_measles

q_s3 <- c(0, 0.25, 0.5, 0.75, 1)
c_s3 <- c(50, 100, 150)

runs <- 1000

parm_combinations <- expand.grid(q=q_s3, c=c_s3)
results <- array(NA, dim=c((length(N_s3)*runs), 7, length(q_s3)*length(c_s3)))

for(i in 1:nrow(parm_combinations)){
  q_i <-parm_combinations[i,1]
  c_i <- parm_combinations[i,2]
  results <- run_analysis2(
    N=N_s3,
    S=S_s3,
    e0=e0_s3,
    ke=ke_s3,
    ki=ki_s3,
    mue=mue_s3_measles,
    mui=mui_s3_measles,
    bdd=(bfd_s3_measles/c_i),
    bfd=bfd_s3_measles,
    q=q_i,
    runs=runs,
    generation_tracking=FALSE)
  results$q <- q_i
  results$c <- c_i
  assign(paste0("results_", i), results)
}

df <- full_join(results_1, results_2)
for(i in 3:nrow(parm_combinations)){
  df <- full_join(df, get(paste0("results_", i)))
}

simulation_results_s3 <- df

#bootstrap introduction risk as a function of pathogen, N and journey time
#identify longest outbreak duration in dataset
time_max <- 250

#list of times to assess introduction risk
times <- seq(0, time_max, by=1)

introduction_risk_s3 <- expand.grid(N=N_s3, time=times, q=q_s3, c=c_s3, p_introduction=NA)

for(i in 1:nrow(parm_combinations)){
  q_i <- parm_combinations$q[i]
  c_i <- parm_combinations$c[i]
  for(j in 1:length(N_s3)){
    subset_2 <- intersect(intersect(which(simulation_results_s3$q==q_i), which(simulation_results_s3$c==c_i)),
                          which(simulation_results_s3$N==N_s3[j]))
    subset_intro <- intersect(intersect(which(introduction_risk_s3$N == N_s3[j]),
                                        which(introduction_risk_s3$q == q_i)), which(introduction_risk_s3$c==c_i))
    vec_j <- numeric(length(times))
    for(k in 1:length(times)){
      vec_j[k] <- length(which(simulation_results_s3$Duration[subset_2] > times[k]))/length(subset_2)
    }
    introduction_risk_s3$p_introduction[subset_intro] <- vec_j
  }
  print(paste(i, "of", nrow(parm_combinations)))
}

saveRDS(introduction_risk_s3, file="simulation_results/introduction_risk_s3.RDS")

### Figure s3 Visualisation ###
theme4_shape <- c(22,23,25, 8, 4,21,10,24)
introduction_risk_s3$p_introduction[which(introduction_risk_s3$p_introduction == 1)] <- 0.99

introduction_risk_s3$q <- factor(introduction_risk_s3$q, levels=q_s3, 
                                 labels=c("q = 0", "q = 0.25", "q = 0.5", "q = 0.75", "q = 1"))

introduction_risk_s3$c <- factor(introduction_risk_s3$c, levels=c_s3, 
                                 labels=c("c = 50", "c = 100", "c = 150"))



figure_s3_contours <- ggplot(introduction_risk_s3) +
  geom_contour_filled(bins=6, mapping=aes(y=N, x=time, z=p_introduction), breaks=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1)) + 
  labs(x="Journey time (days)", y="Total population (N)", fill="Introduction Risk") + 
  facet_grid(q~c) +
  theme_bw() +
  scale_fill_manual(values=alpha(grey.colors(6, start=1, end=0.2), c(0,1,1,1,1,1))) +
  xlim(0,200)

#Panel s3_b: overplotting selected historical voyages
data_panel_s3 <- read.csv("data/Selected_voyages.csv")

#factor column 'Ship' so ships appear in chronological order
data_panel_s3$Ship <- factor(data_panel_s3$Ship,
                              levels=data_panel_s3$Ship)


figure_s3_measles <- figure_s3_contours +
  geom_point(data=data_panel_s3, mapping=aes(x=t, y=N, pch=Ship), col="black", fill="white",  cex=3) + 
  theme(legend.text = element_text(margin = margin(b=2, t=2, unit='pt')), legend.justification='top') +
  guides(fill=guide_colorsteps(show.limits=TRUE, frame.colour="black", frame.linewidth=0.5, frame.linetype=1)) +
  labs(col="Voyage") +
  scale_shape_manual(values=theme4_shape) +
  scale_y_log10(limits=c(10,1500)); figure_s3_measles

pdf(file = "figures/figure_s3.pdf", 
    width = 9.32, # The width of the plot in inches
    height = 6.21) # The height of the plot in inches
figure_s3_measles
dev.off()

#figure 4 sensitivity analysis

### 1: Smallpox

S_proportion_s4 <- 0.05
e0_s4 <- 1

N_s4 <- unique(round(10^seq(1,3.2,by=0.01),0))
S_s4 <- round(S_proportion_s4*N_s4,0)
any(N_s4 < S_s4 + e0_s4) #sanity check

ke_s4 <- 3
ki_s4 <- 3

mue_s4_smallpox <- 8
mui_s4_smallpox <- 17.5
bfd_s4_smallpox <- 7/mui_s4_smallpox

q_s4 <- c(0, 0.25, 0.5, 0.75, 1)
c_s4 <- c(50, 100, 150)

runs <- 500

parm_combinations <- expand.grid(q=q_s4, c=c_s4)
results <- array(NA, dim=c((length(N_s4)*runs), 7, length(q_s4)*length(c_s4)))

for(i in 1:nrow(parm_combinations)){
  q_i <-parm_combinations[i,1]
  c_i <- parm_combinations[i,2]
  results <- run_analysis2(
    N=N_s4,
    S=S_s4,
    e0=e0_s4,
    ke=ke_s4,
    ki=ki_s4,
    mue=mue_s4_smallpox,
    mui=mui_s4_smallpox,
    bdd=(bfd_s4_smallpox/c_i),
    bfd=bfd_s4_smallpox,
    q=q_i,
    runs=runs,
    generation_tracking=FALSE)
  results$q <- q_i
  results$c <- c_i
  assign(paste0("results_", i), results)
}

df <- full_join(results_1, results_2)
for(i in 3:nrow(parm_combinations)){
  df <- full_join(df, get(paste0("results_", i)))
}

simulation_results_s4 <- df

#bootstrap introduction risk as a function of pathogen, N and journey time
#identify longest outbreak duration in dataset
time_max <- 250

#list of times to assess introduction risk
times <- seq(0, time_max, by=1)

introduction_risk_s4 <- expand.grid(N=N_s4, time=times, q=q_s4, c=c_s4, p_introduction=NA)

for(i in 1:nrow(parm_combinations)){
  q_i <- parm_combinations$q[i]
  c_i <- parm_combinations$c[i]
  for(j in 1:length(N_s4)){
    subset_2 <- intersect(intersect(which(simulation_results_s4$q==q_i), which(simulation_results_s4$c==c_i)),
                          which(simulation_results_s4$N==N_s4[j]))
    subset_intro <- intersect(intersect(which(introduction_risk_s4$N == N_s4[j]),
                                        which(introduction_risk_s4$q == q_i)), which(introduction_risk_s4$c==c_i))
    vec_j <- numeric(length(times))
    for(k in 1:length(times)){
      vec_j[k] <- length(which(simulation_results_s4$Duration[subset_2] > times[k]))/length(subset_2)
    }
    introduction_risk_s4$p_introduction[subset_intro] <- vec_j
  }
  print(paste(i, "of", nrow(parm_combinations)))
}

saveRDS(introduction_risk_s4, file="simulation_results/introduction_risk_s4.RDS")

### Figure s4 Visualisation ###
theme4_shape <- c(22,23,25, 8, 4,21,10,24)
introduction_risk_s4$p_introduction[which(introduction_risk_s4$p_introduction == 1)] <- 0.99

introduction_risk_s4$q <- factor(introduction_risk_s4$q, levels=q_s4, 
                                 labels=c("q = 0", "q = 0.25", "q = 0.5", "q = 0.75", "q = 1"))

introduction_risk_s4$c <- factor(introduction_risk_s4$c, levels=c_s4, 
                                 labels=c("c = 50", "c = 75", "c = 100"))



figure_s4_contours <- ggplot(introduction_risk_s4) +
  geom_contour_filled(bins=6, mapping=aes(y=N, x=time, z=p_introduction), breaks=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1)) + 
  labs(x="Journey time (days)", y="Total population (N)", fill="Introduction Risk") + 
  facet_grid(q~c) +
  theme_bw() +
  scale_fill_manual(values=alpha(grey.colors(6, start=1, end=0.2), c(0,1,1,1,1,1))) +
  xlim(0,200)

#Panel S2_b: overplotting selected historical voyages
data_panel_s4 <- read.csv("data/Selected_voyages.csv")

#factor column 'Ship' so ships appear in chronological order
data_panel_s4$Ship <- factor(data_panel_s4$Ship,
                              levels=data_panel_s4$Ship)

figure_s4_smallpox <- figure_s4_contours +
  geom_point(data=data_panel_s4, mapping=aes(x=t, y=N, pch=Ship), col="black", fill="white",  cex=3) + 
  theme(legend.text = element_text(margin = margin(b=2, t=2, unit='pt')), legend.justification='top') +
  guides(fill=guide_colorsteps(show.limits=TRUE, frame.colour="black", frame.linewidth=0.5, frame.linetype=1)) +
  labs(col="Voyage", title="Sensitivity of smallpox introduction risk to q and c") +
  scale_shape_manual(values=theme4_shape) +
  scale_y_log10(limits=c(10,1500)); figure_s4_smallpox

pdf(file = "figures/figure_s4.pdf", 
    width = 9.32, # The width of the plot in inches
    height = 6.21) # The height of the plot in inches
figure_s4_smallpox
dev.off()

## sensitivity analysis on p

#figure 4 sensitivity analysis
S_proportion_s5 <- c(0.01, 0.02, 0.05, 0.1, 0.25, 0.5)
pathogen <- c("influenza", "measles", "smallpox")
e0_s5 <- 1

N_s5 <- unique(round(10^seq(1,3.2,by=0.01),0))


ke_s5 <- 3
ki_s5 <- 3

mue_s5 <- c(2,12,12)
mui_s5 <- c(3, 8, 17.5)
bfd_s5<- c(1.5, 15, 7)/mui_s5
bdd_s5 <- bfd_s5/100
q <- 0.5

runs <- 500

parm_combinations <- expand.grid(c(1:3), S=S_proportion_s5)
results <- array(NA, dim=c((length(N_s5)*runs), 7, length(pathogen)*length(S_proportion_s5)))

for(i in 1:nrow(parm_combinations)){
  pathogen_index <- parm_combinations[i,1]
  S_i <-  parm_combinations[i,2]
  results <- run_analysis2(
    N=N_s5,
    S=round(N_s5*S_i, 0),
    e0=e0_s5,
    ke=ke_s5,
    ki=ki_s5,
    mue=mue_s5[pathogen_index],
    mui=mui_s5[pathogen_index],
    bdd=bdd_s5[pathogen_index],
    bfd=bfd_s5[pathogen_index],
    q=q,
    runs=runs,
    generation_tracking=FALSE)
  results$pathogen <- pathogen[pathogen_index]
  results$prop_S <- S_i
  assign(paste0("results_", i), results)
  print(paste(i, "of", nrow(parm_combinations)))
}

df <- full_join(results_1, results_2)
for(i in 1:nrow(parm_combinations)){
  df <- full_join(df, get(paste0("results_", i)))
}

simulation_results_s5 <- df

#bootstrap introduction risk as a function of pathogen, N and journey time
#identify longest outbreak duration in dataset
time_max <- 250

#list of times to assess introduction risk
times <- seq(0, time_max, by=1)

introduction_risk_s5 <- expand.grid(N=N_s5, S = S_proportion_s5, time=times, 
                                    pathogen= pathogen, p_introduction=NA)

for(i in 1:nrow(parm_combinations)){
  pathogen_index <- parm_combinations[i,1]
  S_i <-  parm_combinations[i,2]
  for(j in 1:length(N_s5)){
    subset_2 <- intersect(intersect(which(simulation_results_s5$prop_S==S_i), which(simulation_results_s5$pathogen==pathogen[pathogen_index])),
                          which(simulation_results_s5$N==N_s5[j]))
    subset_intro <- intersect(intersect(which(introduction_risk_s5$N == N_s5[j]),
                                        which(introduction_risk_s5$S == S_i)), which(introduction_risk_s5$pathogen==pathogen[pathogen_index]))
    vec_j <- numeric(length(times))
    for(k in 1:length(times)){
      vec_j[k] <- length(which(simulation_results_s5$Duration[subset_2] > times[k]))/length(subset_2)
    }
    introduction_risk_s5$p_introduction[subset_intro] <- vec_j
  }
  print(paste(i, "of", nrow(parm_combinations)))
}

saveRDS(introduction_risk_s5, file="simulation_results/introduction_risk_s5.RDS")

### Figure s5 Visualisation ###
theme4_shape <- c(22,23,25, 8, 4,21,10,24)
introduction_risk_s5$p_introduction[which(introduction_risk_s5$p_introduction == 1)] <- 0.99
introduction_risk_s5$pathogen <- factor(introduction_risk_s5$pathogen, levels=pathogen,
                                        labels=c("Influenza", "Measles", "Smallpox"))

figure_s5_contours <- ggplot(introduction_risk_s5) +
  geom_contour_filled(bins=6, mapping=aes(y=N, x=time, z=p_introduction), breaks=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1)) + 
  labs(x="Journey time (days)", y="Total population (N)", fill="Introduction Risk") + 
  facet_grid(pathogen~S) +
  theme_bw() +
  scale_fill_manual(values=alpha(grey.colors(6, start=1, end=0.2), c(0,1,1,1,1,1))) +
  xlim(0,200)

#Panel s5_b: overplotting selected historical voyages
data_panel_s5 <- read.csv("data/Selected_voyages.csv")

#factor column 'Ship' so ships appear in chronological order
data_panel_s5$Ship <- factor(data_panel_s5$Ship,
                              levels=data_panel_s5$Ship)


panel_s5 <- figure_s5_contours +
  geom_point(data=data_panel_s5, mapping=aes(x=t, y=N, pch=Ship), col="black", fill="white",  cex=3) + 
  theme(legend.text = element_text(margin = margin(b=2, t=2, unit='pt')), legend.justification='top') +
  guides(fill=guide_colorsteps(show.limits=TRUE, frame.colour="black", frame.linewidth=0.5, frame.linetype=1)) +
  labs(col="Voyage") +
  scale_shape_manual(values=theme4_shape) +
  scale_y_log10(limits=c(10,1500)); panel_s5

pdf(file = "figures/figure_s5.pdf", 
    width = 9.32, # The width of the plot in inches
    height = 6.21) # The height of the plot in inches
   panel_s5
dev.off()



