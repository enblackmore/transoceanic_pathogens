source('functions.R')
set.seed(1492)

############################################
## Table 1: introduction risk estimates   ##
## for selected voyages to San Francisco  ##  
############################################

# Parameters for tables 1 and 2
S_proportion_tables <- 0.05
e0_tables <- 1

ke_tables <- 3
ki_tables <- 3

q_tables <- 0.5
bdd_tables <- 0.02

#bdd depends on mue, which varies by pathogen;
#calculate each separately
mue_tables_influenza <- 2
mui_tables_influenza <- 3
bfd_tables_influenza <- 1.5/mui_tables_influenza

mue_tables_measles <- 12
mui_tables_measles <- 8
bfd_tables_measles <- 15/mui_tables_measles

mue_tables_smallpox <- 12
mui_tables_smallpox <- 17.5
bfd_tables_smallpox <- 7/mui_tables_smallpox

#Get ship data
SF_data <- read.csv('data/San_Francisco_arrivals.csv')
SF_selected_ships_simulation_tables <- c(448, 397, 25, 64, 168, 218, 341, 373, 545, 522, 132, 136, 494, 492)
results_table1 <- SF_data[SF_selected_ships_simulation_tables,c(4,5,9,10,15)]
results_table1$S <- round(S_proportion_tables*results_table1$N_passengers, 0)

#Bootstrap introduction risk
results_table1$risk_influenza <- get_ship_risk(
  N=results_table1$N_passengers,
  S=results_table1$S,
  t=results_table1$Voyage_days,
  e0=1,
  bdd=bdd_tables,
  bfd=bfd_tables_influenza,
  mue=mue_tables_influenza,
  mui=mui_tables_influenza,
  ke=ke_tables,
  ki=ki_tables,
  q=q_tables,
  generation_tracking = FALSE,
  runs=5000
)

results_table1$risk_measles <- get_ship_risk(
  N=results_table1$N_passengers,
  S=results_table1$S,
  t=results_table1$Voyage_days,
  e0=1,
  bdd=bdd_tables,
  bfd=bfd_tables_measles,
  mue=mue_tables_measles,
  mui=mui_tables_measles,
  ke=ke_tables,
  ki=ki_tables,
  q=q_tables,
  generation_tracking = FALSE,
  runs=5000
)

results_table1$risk_smallpox <- get_ship_risk(
  N=results_table1$N_passengers,
  S=results_table1$S,
  t=results_table1$Voyage_days,
  e0=1,
  bdd=bdd_tables,
  bfd=bfd_tables_smallpox,
  mue=mue_tables_smallpox,
  mui=mui_tables_smallpox,
  ke=ke_tables,
  ki=ki_tables,
  q=q_tables,
  generation_tracking = FALSE,
  runs=5000
)

saveRDS(results_table1, file = "simulation_results/results_table1.RDS")

##########################################
## Table 2: introduction risk estimates ##
## for selected voyages 1492-1918       ##
##########################################

#Get data
results_table2 <- read.csv('data/selected_voyages.csv')
results_table2$S <- round(results_table2$N*S_proportion_tables, 0)

#Bootstrap introduction risk
results_table2$risk_influenza <- get_ship_risk(
  N=results_table2$N,
  S=results_table2$S,
  t=results_table2$t,
  e0=1,
  bdd=bdd_tables,
  bfd=bfd_tables_influenza,
  mue=mue_tables_influenza,
  mui=mui_tables_influenza,
  ke=ke_tables,
  ki=ki_tables,
  q=q_tables,
  generation_tracking = FALSE,
  runs=5000
)

results_table2$risk_measles <- get_ship_risk(
  N=results_table2$N,
  S=results_table2$S,
  t=results_table2$t,
  e0=1,
  bdd=bdd_tables,
  bfd=bfd_tables_measles,
  mue=mue_tables_measles,
  mui=mui_tables_measles,
  ke=ke_tables,
  ki=ki_tables,
  q=q_tables,
  generation_tracking = FALSE,
  runs=5000
)

results_table2$risk_smallpox <- get_ship_risk(
  N=results_table2$N,
  S=results_table2$S,
  t=results_table2$t,
  e0=1,
  bdd=bdd_tables,
  bfd=bfd_tables_smallpox,
  mue=mue_tables_smallpox,
  mui=mui_tables_smallpox,
  ke=ke_tables,
  ki=ki_tables,
  q=q_tables,
  generation_tracking = FALSE,
  runs=5000
)

saveRDS(results_table2, file = "simulation_results/results_table2.RDS")


##############################################
## Figure 4: cumulative introduction risk   ##
## varying pathogens and ship size (N)      ##
##############################################

### Figure 4 simulation ###

S_proportion_4 <- 0.05
e0_4 <- 1

N_4 <- round(10^seq(1,3.2,by=0.01),0)
S_4 <- round(S_proportion_4*N_4,0)
any(N_4 < S_4 + e0_4) #sanity check

ke_4 <- 3
ki_4 <- 3

q_4 <- 0.5
bdd_4 <- 0.02


#bdd depends on mue, which varies by pathogen;
#calculate each separately
mue_4_influenza <- 2
mui_4_influenza <- 3
bfd_4_influenza <- 1.5/mui_4_influenza

mue_4_measles <- 12
mui_4_measles <- 8
bfd_4_measles <- 15/mui_4_measles

mue_4_smallpox <- 12
mui_4_smallpox <- 17.5
bfd_4_smallpox <- 7/mui_4_smallpox

pathogens_4 <- c("Influenza", "Measles", "Smallpox")

simulation_results_4_influenza <- run_analysis2(
  N=N_4,
  S=S_4,
  e0=e0_4,
  ke=ke_4,
  ki=ki_4,
  mue=mue_4_influenza,
  mui=mui_4_influenza,
  bdd=bdd_4,
  bfd=bfd_4_influenza,
  q=q_4,
  runs=500,
  generation_tracking=FALSE)

simulation_results_4_measles <- run_analysis2(
  N=N_4,
  S=S_4,
  e0=e0_4,
  ke=ke_4,
  ki=ki_4,
  mue=mue_4_measles,
  mui=mui_4_measles,
  bdd=bdd_4,
  bfd=bfd_4_measles,
  q=q_4,
  runs=500,
  generation_tracking=FALSE)

simulation_results_4_smallpox <- run_analysis2(
  N=N_4,
  S=S_4,
  e0=e0_4,
  ke=ke_4,
  ki=ki_4,
  mue=mue_4_smallpox,
  mui=mui_4_smallpox,
  bdd=bdd_4,
  bfd=bfd_4_smallpox,
  q=q_4,
  runs=500,
  generation_tracking=FALSE)

#add pathogen indicators manually
simulation_results_4_influenza$Pathogen <- "Influenza"
simulation_results_4_measles$Pathogen <- "Measles"
simulation_results_4_smallpox$Pathogen <- "Smallpox"

#Combine analysis in a data frame
simulation_results_4 <- dplyr::bind_rows(simulation_results_4_influenza,
                                         dplyr::bind_rows(simulation_results_4_measles,
                                                          simulation_results_4_smallpox))
saveRDS(simulation_results_4, file = "simulation_results/simulation_results_4.RDS")

#bootstrap introduction risk as a function of pathogen, N and journey time
#identify longest outbreak duration in dataset
time_max <- round(max(simulation_results_4$Duration),0)

#list of times to assess introduction risk
times <- seq(0, time_max, by=1)

#create output data frame to store results
introduction_risk_4 <- expand.grid(pathogens_4, N_4, times, NA)
colnames(introduction_risk_4) <- c("Pathogen", "N", 'time', 'p_introduction')

for(i in 1:nrow(introduction_risk_4)){
  subset <- dplyr::filter(simulation_results_4, 
                          N==introduction_risk_4$N[i] & Pathogen==introduction_risk_4$Pathogen[i])
  introduction_risk_4$p_introduction[i] <- length(which(subset$Duration > introduction_risk_4$time[i]))/nrow(subset)
  print(i/nrow(introduction_risk_4))
}

saveRDS(introduction_risk_4, file="simulation_results/introduction_risk_4.RDS")

### Figure 4 Visualisation ###

#aesthetics
theme4_color <- c("#CE1FFF","#710F0F", "#09f04a", "#FFCA09",
                  "#0038E0", "#FF0A3F", "#0cbcff")
theme4_shape <- c(22,23,25, 8, 4,21,10,24)

#convert p_introduction values of 1 to 0.99
#as no contour bin for p >= 1
introduction_risk_4$p_introduction[which(introduction_risk_4$p_introduction == 1)] <- 0.99

figure_4_contours <- ggplot(introduction_risk_4) +
  geom_contour_filled(bins=6, mapping=aes(y=N, x=time, z=p_introduction), breaks=c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1)) + 
  labs(x="Journey time (days)", y="Total population (N)", fill="Introduction Risk") + 
  facet_wrap(vars(Pathogen)) +
  theme_bw() +
  scale_fill_manual(values=alpha(grey.colors(6, start=1, end=0.2), c(0,1,1,1,1,1))) +
  xlim(0,200)

#Panel 4a: overplotting San Francisco journey data
#Adding journey time data
data_SF <- read.csv('data/San_Francisco_arrivals.csv')

#Reduce number of origins for easier plotting
data_panel_4a <- dplyr::filter(data_SF, From_code != "LVP")
data_panel_4a$From_code <- factor(data_panel_4a$From_code, levels=c("NYC", "VAL", "PAN", "ORE", "HAI", "SYD", "HKG"),
                                  labels=c("New York City", "Valparaíso", "Panama", "Oregon", "Hawai'i", "Sydney", "Hong Kong"))

#split voyage origins across facets for visual clarity
data_panel_4a$Pathogen <- "Influenza"
data_panel_4a$Pathogen[which(data_panel_4a$From_code %in% c("Hong Kong", "Sydney", "Hawai'i"))] <- "Measles"
data_panel_4a$Pathogen[which(data_panel_4a$From_code %in% c("Panama", "New York City"))] <- "Smallpox"

#create a data frame of ships for which numerical estimates are given in table 1
#to overplot crosses
results_table1$From_code <- factor(results_table1$From_code, levels=c("NYC", "VAL", "PAN", "ORE", "HAI", "SYD", "HKG"),
                                  labels=c("New York City", "Valparaíso", "Panama", "Oregon", "Hawai'i", "Sydney", "Hong Kong"))
results_table1$Pathogen <- "Influenza"
results_table1$Pathogen[which(results_table1$From_code %in% c("Hong Kong", "Sydney", "Hawai'i"))] <- "Measles"
results_table1$Pathogen[which(results_table1$From_code %in% c("Panama", "New York City"))] <- "Smallpox"

panel_4a <- figure_4_contours +
  geom_point(data=data_panel_4a, mapping=aes(x=Voyage_days, y=N_passengers, pch=Steam, col=From_code), cex=2) +
  geom_point(data=results_table1, mapping=aes(x=Voyage_days, y=N_passengers), col='black', pch=4, cex=3) +
  guides(fill=guide_colorsteps(frame.colour="black", frame.linewidth=0.5, frame.linetype=1)) +
  labs(col="Origin port", pch="Technology") + theme(legend.justification = 'top') +
  scale_color_manual(values=theme4_color) +
  scale_shape_manual(values=c(16,17), labels=c("Sail", "Steam")) +
  scale_y_log10(limits=c(10,1500)); panel_4a

#Panel 4b: overplotting selected historical voyages
data_panel_4b <- read.csv("data/Selected_voyages.csv")

#factor column 'Ship' so ships appear in chronological order
data_panel_4b$Ship <- factor(data_panel_4b$Ship,
                             levels=data_panel_4b$Ship)

panel_4b <- figure_4_contours +
  geom_point(data=data_panel_4b, mapping=aes(x=t, y=N, pch=Ship), col="black", fill="white",  cex=3) + 
  theme(legend.text = element_text(margin = margin(b=2, t=2, unit='pt')), legend.justification='top') +
  guides(fill=guide_colorsteps(show.limits=TRUE, frame.colour="black", frame.linewidth=0.5, frame.linetype=1)) +
  labs(col="Voyage") +
  scale_shape_manual(values=theme4_shape) +
  scale_y_log10(limits=c(10,1500)); panel_4b

### Figure 4 assembly ###

#Extract legends for figure rearrangement
fig4a_col <- cowplot::get_legend(panel_4a + 
                                   guides(fill='none', 
                                          shape='none'
                                   ))

fig4a_shape <-  cowplot::get_legend(panel_4a +
                                      guides(col='none',
                                             fill='none'))

fig4a_fill <- cowplot::get_legend(panel_4a + 
                                    guides(col='none',
                                           shape='none', 
                                           fill=guide_colorsteps(show.limits=TRUE, 
                                                                 frame.colour="black",
                                                                 frame.linewidth=0.5,
                                                                 frame.linetype=1,
                                                                 title.position='top')) &
                                    theme(legend.position = 'bottom',
                                          legend.box = 'horizontal',
                                          legend.box.margin=margin(t=-12, unit='pt'),
                                          legend.key.width=unit(22, unit='pt'),
                                          legend.text=element_text(size=unit(8, 'pt'))))

fig4b_shape <-  cowplot::get_legend(panel_4b + guides(fill='none') +
                                      theme(legend.box.margin = margin(l=35, unit='pt'),
                                            legend.text = element_text(margin = margin(b=2, t=2, unit='pt'))))


layout_4 <- "
AAAEC
BBBF#
D####
"

pdf(file = "figures/figure_4.pdf", 
    width = 11, # The width of the plot in inches
    height = 7) # The height of the plot in inches
patchwork::wrap_plots(
  A = panel_4a + guides(shape='none', fill='none', col='none'),
  B = panel_4b + guides(fill='none', shape='none') + scale_y_log10(limits=c(10,1500)),
  C = fig4a_shape,
  D = fig4a_fill,
  E = fig4a_col,
  'F' = fig4b_shape,
  design = layout_4,
  widths = c(1,1,1, 0.5, 0.25),
  heights = c(1,1,0.25)
) + patchwork::plot_annotation(tag_levels = list(c('A', 'B','', '', '', '')))
dev.off()



