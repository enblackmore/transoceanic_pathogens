#Simulations and Analysis

set.seed(1492)

################################
#Simulation (1): Basic Dynamics
#corresponds with figures 1a-1c
################################

#In a fully-susceptible population (N=100) with one index case,
#using a hypothetical pathogen ke=ki=3, mue=mui=5
#what happens when we vary r0?

#constant parameters
N_1 <- 100
S_1 <- 99
e0_1 <- 1

mue_1 <- mui_1 <- 5
ke_1 <- ki_1 <- 3

#varying parameters
r0_1 <- 10^seq(-1, 2, by=0.01)

#since the modeul function requires bdd, bfd, and q
#we set q=1 and bdd=r0 / (mue*N) 
#mathematically this is equivalent to q=0 and bfd=r0/mue
bdd_1 <- r0_1/(mue_1*N_1)
q_1 <- 1
bfd_1 <- 0

simulation_results_1 <- run.analysis(
  N=N_1,
  S=S_1,
  e0=e0_1,
  ke=ke_1,
  ki=ki_1,
  mue=mue_1,
  mui=mui_1,
  bdd=bdd_1,
  bfd=bfd_1,
  q=q_1,
  runs=5,
  generation_tracking = TRUE,
  generation_max=20)

#check for apporpriate generation_max
check.generation.max(simulation_results_1)

#add r0 column to analysis
simulation_results_1$analysis$r0 <- simulation_results_1$analysis$bdd*mue_1*N_1

#label simulations with:
#(1) single-generation transmission
#(2) transmission below herd immunity
#(3) transmission at or above herd immunity

simulation_results_1$analysis$label <- NA
simulation_results_1$analysis$label[which(simulation_results_1$analysis$Generations==1)] <- 1
simulation_results_1$analysis$label[which(simulation_results_1$analysis$Generations > 1 
                                          & simulation_results_1$analysis$Cases/N_1 < max(1-(1/simulation_results_1$analysis$r0),0))] <- 2
simulation_results_1$analysis$label[which(simulation_results_1$analysis$Cases/N_1 >= 1-(1/simulation_results_1$analysis$r0))] <- 3
simulation_results_1$analysis$label <- factor(simulation_results_1$analysis$label, levels=c(1, 2, 3), 
                                              labels=c("Single-generation", "Below herd immunity", "At or above herd immunity"))


##########################################################
#Simulation (2): more runs, for selected r0 values
#so we can bootstrap cdf for introduction risk across time
#corresponds with figure 1d
##########################################################

#Constant parameters
N_2 <- 100
S_2 <- 99
e0_2 <- 1
ke_2 <- ki_2 <- 3
mue_2 <- mui_2 <- 5
bfd_2 <- 0
q_2 <- 1

#selected r0 values and corresponding values of bdd
r0_2 <- c(0.5, 2, 8)
bdd_2 <- r0_2/(mue_2*N_2)

simulation_results_2 <- run.analysis(
  N=N_2,
  S=S_2,
  e0=e0_2,
  ke=ke_2,
  ki=ki_2,
  mue=mue_2,
  mui=mui_2,
  bdd=bdd_2,
  bfd=bfd_2,
  q=q_2,
  runs=15,
  generation_tracking = FALSE)

#bootstrap cumulative introduction risk:
introduction_risk_2 <- introduction.risk.time(simulation_results_2$analysis)

#rewrite bdd as r0 for plotting
introduction_risk_2$r0 <- introduction_risk_2$bdd*mue_2*N_2

###################################################
#Simulation (3): same r0 values as simulation 2
#using mue and mui of real-life pathogens
#corresponds with figure 1e
###################################################
#Constants
N_3 <- 100
S_3 <- 99
e0_3 <- 1
ke_3 <- 3
ki_3 <- 3
q_3 <- 1
bfd_3 <- 0 #as above, since q=1 the value of bfd is irrelevant

#Variables
r0_3 <- c(0.5, 2, 8)

#bdd depends on mue, which varies by pathogen;
#calculate each separately
mue_3_influenza <- 2
mui_3_influenza <- 3
bdd_3_influenza <- r0_3/(mui_3_influenza*N_3)


mue_3_measles <- 12
mui_3_measles <- 8
bdd_3_measles <- r0_3/(mui_3_measles*N_3)

mue_3_smallpox <- 12
mui_3_smallpox <- 17.5
bdd_3_smallpox <- r0_3/(mui_3_smallpox*N_3)

simulation_results_3_influenza <- run.analysis(
  N=N_3,
  S=S_3,
  e0=e0_3,
  ke=ke_3,
  ki=ki_3,
  mue=mue_3_influenza,
  mui=mui_3_influenza,
  bdd=bdd_3_influenza,
  bfd=bfd_3,
  q=q_3,
  runs=500,
  generation_tracking=FALSE)

simulation_results_3_measles <- run.analysis(
  N=N_3,
  S=S_3,
  e0=e0_3,
  ke=ke_3,
  ki=ki_3,
  mue=mue_3_measles,
  mui=mui_3_measles,
  bdd=bdd_3_measles,
  bfd=bfd_3,
  q=q_3,
  runs=500,
  generation_tracking=FALSE)

simulation_results_3_smallpox <- run.analysis(
  N=N_3,
  S=S_3,
  e0=e0_3,
  ke=ke_3,
  ki=ki_3,
  mue=mue_3_smallpox,
  mui=mui_3_smallpox,
  bdd=bdd_3_smallpox,
  bfd=bfd_3,
  q=q_3,
  runs=500,
  generation_tracking=FALSE)

#add pathogen names manually
simulation_results_3_influenza$analysis$Pathogen <-"Influenza"
simulation_results_3_measles$analysis$Pathogen <-"Measles"
simulation_results_3_smallpox$analysis$Pathogen <-"Smallpox"

#combine analysis in a single data frame
simulation_results_3 <- dplyr::full_join(simulation_results_3_smallpox$analysis,
                                         dplyr::full_join(simulation_results_3_measles$analysis, 
                                                          simulation_results_3_influenza$analysis))


###################################################
#Simulation (4): same r0 values as simulation 2
#constant re0, with varying S
#across three re0 values
#corresponds with figure 2a
###################################################

#constants
N_4 <- 1001
S_4 <- unique(round(10^seq(0, 3, by=0.1),0))
e0_4 <- 1
ke_4 <- ki_4 <- 3
mue_4 <- mui_4 <- 5
q_4 <- 1
bfd_4 <- 0

#variables
re0_4 <- c(1, 2, 8)

#we maintain constant re0 across different values of S
#by setting q=1
#and back-calculating bdd for each analysis:
#bdd = re0_4/(mui_4*S_4) 

mui_4 <- 5
bdd_4_re0_1 <- re0_4[1]/(mui_4*S_4)
bdd_4_re0_2 <- re0_4[2]/(mui_4*S_4)
bdd_4_re0_8 <- re0_4[3]/(mui_4*S_4)

simulation_results_4_re0_1 <- run.analysis2(
  N=N_4,
  S=S_4,
  e0=e0_4,
  ke=ke_4,
  ki=ki_4,
  mue=mue_4,
  mui=mui_4,
  bdd=bdd_4_re0_1,
  bfd=bfd_4,
  q=q_4,
  runs=500,
  generation_tracking = FALSE
)

simulation_results_4_re0_2 <- run.analysis2(
  N=N_4,
  S=S_4,
  e0=e0_4,
  ke=ke_4,
  ki=ki_4,
  mue=mue_4,
  mui=mui_4,
  bdd=bdd_4_re0_2,
  bfd=bfd_4,
  q=q_4,
  runs=500,
  generation_tracking = FALSE
)

simulation_results_4_re0_8 <- run.analysis2(
  N=N_4,
  S=S_4,
  e0=e0_4,
  ke=ke_4,
  ki=ki_4,
  mue=mue_4,
  mui=mui_4,
  bdd=bdd_4_re0_8,
  bfd=bfd_4,
  q=q_4,
  runs=500,
  generation_tracking = FALSE
)

#manually add re0 values
simulation_results_4_re0_1$analysis$re0 <- 1
simulation_results_4_re0_2$analysis$re0 <- 2
simulation_results_4_re0_8$analysis$re0 <- 8

#combine analyses in a single data frame
simulation_results_4 <- dplyr::full_join(simulation_results_4_re0_1$analysis,
                                         dplyr::full_join(simulation_results_4_re0_2$analysis,
                                                          simulation_results_4_re0_8$analysis))

#############################################################
#Simulation (5): introducing density and frequency dependence
#varying S/N
#for three values of q
#for three values of bfd (q=0)
#for three values of bdd (q=1)
#and for three values of bdd with constant bfd (q=0.5)
#corresponds with figure 2b
#############################################################

#fixed values
e0_5 <- 1
mue_5 <- 5
mui_5 <- 5
ke_5 <- 3
ki_5 <- 3

#variables
N_5 <- c(50, 100, 200, 500)
proportion_S_5 <- 10^seq(-2, -0.05, by=0.05)

q_5 <- c(0,0.5,1)
bfd_varying_5 <- c(2, 5, 10)/mui_5
bfd_fixed_5 <- 5/mui_5
bdd_5 <- c(0.05, 0.1, 0.20)/mue_5

#reformat as input for run.analysis2()
all_variables_5 <- expand.grid(N=N_5, proportion=proportion_S_5, q=q_5, b_index=c(1,2,3))
all_variables_5$S <-  round(all_variables_5$proportion*all_variables_5$N,0)
any(all_variables_5$S > all_variables_5$N+e0_5) #check to make sure rounding doesn't introduce impossible S values

#add bfd and bdd values, which are only used for some analyses
#is a bdd or bfd value is mathematically irrelevant (e.g. q=0 or q=1), use placeholder value 0
all_variables_5$bfd <- all_variables_5$bdd <- 0

#for q=0, bdd is irrelevant and bfd varies
all_variables_5$bfd[which(all_variables_5$q==0)] <- bfd_varying_5[all_variables_5$b_index[which(all_variables_5$q==0)]]

#for q=1, bfd is irrelevant and bdd varies
all_variables_5$bdd[which(all_variables_5$q==1)] <- bdd_5[all_variables_5$b_index[which(all_variables_5$q==1)]]

#for q=0.5, bdd varies and bfd is constant
all_variables_5$bdd[which(all_variables_5$q==0.5)] <- bdd_5[all_variables_5$b_index[which(all_variables_5$q==0.5)]]
all_variables_5$bfd[which(all_variables_5$q==0.5)] <- bfd_fixed_5

#run analysis

simulation_results_5 <- run.analysis2(
  N=all_variables_5$N,
  S=all_variables_5$S,
  e0=e0_5,
  ke=ke_5,
  ki=ki_5,
  mue=mue_5,
  mui=mui_5,
  bdd=all_variables_5$bdd,
  bfd=all_variables_5$bfd,
  q=all_variables_5$q,
  runs=500,
  generation_tracking=FALSE
)

#############################################################
#Simulation (6): cumulative introduction risk
#varying pathogens and ship size (N)
#############################################################

S_proportion_6 <- 0.05
e0_6 <- 1

N_6 <- round(10^seq(1,3.2,by=0.01),0)
S_6 <- round(S_proportion_6*N_6,0)
any(N_6 < S_6 + e0_6) #sanity check

ke_6 <- 3
ki_6 <- 3

q_6 <- 0.5
bdd_6 <- 0.05


#bdd depends on mue, which varies by pathogen;
#calculate each separately
mue_6_influenza <- 2
mui_6_influenza <- 3
bfd_6_influenza <- 1.5/mui_6_influenza


mue_6_measles <- 12
mui_6_measles <- 8
bfd_6_measles <- 15/mui_6_measles

mue_6_smallpox <- 12
mui_6_smallpox <- 17.5
bfd_6_smallpox <- 7/mui_6_smallpox


simulation_results_6_influenza <- run.analysis2(
  N=N_6,
  S=S_6,
  e0=e0_6,
  ke=ke_6,
  ki=ki_6,
  mue=mue_6_influenza,
  mui=mui_6_influenza,
  bdd=bdd_6,
  bfd=bfd_6_influenza,
  q=q_6,
  runs=5,
  generation_tracking =FALSE)

simulation_results_6_measles <- run.analysis2(
  N=N_6,
  S=S_6,
  e0=e0_6,
  ke=ke_6,
  ki=ki_6,
  mue=mue_6_measles,
  mui=mui_6_measles,
  bdd=bdd_6,
  bfd=bfd_6_measles,
  q=q_6,
  runs=5,
  generation_tracking=FALSE)


simulation_results_6_smallpox <- run.analysis2(
  N=N_6,
  S=S_6,
  e0=e0_6,
  ke=ke_6,
  ki=ki_6,
  mue=mue_6_smallpox,
  mui=mui_6_smallpox,
  bdd=bdd_6,
  bfd=bfd_6_smallpox,
  q=q_6,
  runs=5,
  generation_tracking=FALSE)

#add pathogen indicators manually
simulation_results_6_influenza$analysis$Pathogen <- "Influenza"
simulation_results_6_measles$analysis$Pathogen <- "Measles"
simulation_results_6_smallpox$analysis$Pathogen <- "Smallpox"

#Combine analysis in a data frame
simulation_results_6 <- dplyr::full_join(simulation_results_6_influenza$analysis,
                                         dplyr::full_join(simulation_results_6_measles$analysis,
                                                          simulation_results_6_smallpox$analysis))

#bootstrap introduction risk as a function of pathogen, N and journey time
introduction_risk_6 <- introduction.risk.time(df)

