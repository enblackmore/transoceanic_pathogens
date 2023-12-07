#Visualisations
library(ggplot2)

#Figure 1

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


#Panel 1a: duration by r0 and outbreak size
panel_1a <- ggplot(simulation_results_1$analysis) +
  geom_point(mapping=aes(x=r0, y=Duration, col=label), cex=2, alpha=0.05) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=t5), lwd=1) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=t05), lty="dashed", lwd=0.5) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=t95), lty="dashed", lwd=0.5) +
  scale_x_log10() +
  scale_y_continuous() +
  labs(x=bquote("R"[0]), y="Outbreak duration (days)", col="Outbreak size ") + theme_bw(); panel_1a

panel_1b <- ggplot(simulation_results_1$analysis) + 
  geom_point(mapping=aes(x=r0, y=Cases, col=label), cex=1, alpha=0.05) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=c5), lwd=0.75) + 
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=c05), lty="dashed", lwd=0.25) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=c95), lty="dashed", lwd=0.25) + 
  scale_x_log10() +
  theme_bw() + labs(y="Outbreak size", x=bquote(R[0])); panel_1b 

panel_1c <- ggplot(simulation_results_1$analysis) +
  geom_point(mapping=aes(x=r0, y=(Generations-1), col=label), cex=1.5, alpha=0.01) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=g5-1), lwd=0.75) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=g05-1), lty="dashed", lwd=0.25) +
  geom_line(data=simulation_results_1_quantiles, mapping=aes(x=r0, y=g95-1),  lty="dashed", lwd=0.25) + 
  scale_x_log10() +
  theme_bw() + labs(y="Transmission\ngenerations", x=bquote(R[0])); panel_1c

introduction_risk_2$r0 <- factor(introduction_risk_2$r0, levels=unique(introduction_risk_2$r0))
panel_1d <- ggplot(introduction_risk_2) + geom_line(mapping=aes(x=time, y=p_introduction, col=r0), lwd=2) +
  labs(x="Journey time (days)", y="Introduction risk", col=bquote(R[0])) +
  theme_bw() + scale_color_manual(values=theme1b) +
  theme(legend.title=element_text(size=9), legend.position='right') +
  theme(legend.position = 'none')
panel_1d

