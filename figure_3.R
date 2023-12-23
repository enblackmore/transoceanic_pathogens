
####################################################
## Figure 3, Table S2: San Francisco Port Arrivals #
####################################################

theme3  <- c( "#07f49e", "#42047e")

#######################################################
### Table S2: San Francisco port arrival statistics ###
#######################################################

#read and summarise arrivals data
data_SF <- read.csv("data/San_Francisco_arrivals.csv")

Longitude_order <- c(1,7,5,2,6,8,3,4)
data_SF$From_code <- factor(data_SF$From_code, levels=unique(data_SF$From_code)[Longitude_order],
                            labels=c("Liverpool", "New York City", "Panama", "Valparaíso",
                                     "Oregon", "Hawai'i", "Sydney", "Hong Kong"))

#summarise data; this becomes supplementary table 2
data_SF_summary <- plyr::ddply(data_SF,
                               ~From_code+Steam, plyr::summarise,
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

########################################
### Map 1: journeys to San Francisco ###
########################################

#Map 1: median journey time to San Francisco
figure_3_map_data <- read_csv('data/figure_3_map_data.csv')
pacific_centered_map <- ggplot2::map_data("world", wrap=c(0,360))


map1 <- ggplot(pacific_centered_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(col='#B1B1AA', fill='#B1B1AA') + 
  theme_bw() + 
  coord_cartesian(xlim=c(120,330)) + 
  labs(x="", y="") +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill='#F0F8FF'),
        axis.text=element_blank(), 
        axis.ticks = element_blank(),
        aspect.ratio = 4/5.25)

#geom_curve does not allow curvature and angle as a mapped aesthetics
#so pass through a 'for' loop
for(i in (1:nrow(figure_3_map_data))[-5]){
  map1 <- map1 + geom_curve(data=figure_3_map_data[i,], mapping=aes(x=x, y=y, xend=xend, yend=yend),
                            curvature=figure_3_map_data$curvature[i], 
                            angle=figure_3_map_data$angle[i], ncp=10)
}

#add points and labels
map1 <- map1 + geom_point(data=figure_3_map_data, mapping=aes(x=x, y=y), col='#28464B') +
  geom_label(data=figure_3_map_data[-5,], mapping=aes(x=labelx, y=labely, label=labelname), col='#28464B',
             cex=2.5, hjust=0, fill='#FEFEFA'); map1

########################################
### Panel 3a: journey time by origin ###
########################################
data_SF$Steam <- factor(data_SF$Steam, levels=c(TRUE, FALSE), labels=c("Steam", "Sail"))

panel_3a <- ggplot(data_SF) +
  geom_jitter(mapping=aes(x=Voyage_days, y=From_code, col=Steam), height=0.1, alpha=0.3) +
  theme_bw() +
  labs(x="Journey Time (days)", y="Origin Port", col="Technology") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(breaks=seq(0, 250, by=50), limits = c(0,250)) +
  scale_color_manual(values=rev(theme3), labels=c("Sail", "Steam")) +
  theme(axis.title.y=element_text(margin=margin(r=-10, unit='pt')),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        axis.title.x=element_text(size=10.2),
        legend.position='none',
        aspect.ratio=1.5/1)

##############################################
### Panel 3b: journey passengers by origin ###
##############################################

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

##########################################
### Panel 3c: journey number by origin ###
##########################################

panel_3c <- ggplot(data_SF) +
  geom_bar(mapping=aes(x=From_code, fill=Steam), show.legend = FALSE) +
  geom_point(mapping=aes(x=From_code, y=1, col=Steam), alpha=0) +
  theme_bw() +
  labs(x="Origin Port", y="Arrivals", col="") +
  scale_fill_manual(values=theme3) +
  scale_color_manual(values=theme3,
                     guide=guide_legend(override.aes = list(alpha=1, size=4))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(margin=margin(t=-5, unit='pt')),
        axis.title.y=element_text(margin=margin(r=0, b=-10, unit='pt')),
        legend.position=c(0.75, 0.90),
        legend.text = element_text(margin=margin(l=-4, unit='pt')),
        legend.background=element_rect(fill=NA),
        aspect.ratio = 1.25/1)

#########################
### Figure 3 assembly ###
#########################


pdf(file = "figures/figure_3.pdf", 
    width = 13.75, # The width of the plot in inches
    height = 4.75) # The height of the plot in inches
patchwork::wrap_elements(full=map1) + panel_3a + panel_3b + panel_3c +
  patchwork::plot_layout(widths=c(3,1,1,1)) +
  patchwork::plot_annotation(tag_levels='A')
dev.off()