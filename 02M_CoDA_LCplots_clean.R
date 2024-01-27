
#plotting the forecast

#reset the plot (if needed)
dev.off()
graphics.off()
#plot the forecast
col.pro <- rainbow(16)
col.cause <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#28E2E5")

CoD.label <- c("Other non-cardio", 
               "Rheumatic","Hypertension", 
               "Hypertensive disease","Myocardial infarction","Other IHD",
               "Pulmonary", "Non-rheumatic", "Cardiac arrest", "Heart failure", "Other cardio")

plot_obs <- pro.dea.obs

colnames(plot_obs) = years
rownames(plot_obs) = CoD.label

datLong_obs <- melt(data          = plot_obs,
                    id.vars       = c(CoD.label),
                    measure.vars  = c(years),
                    variable.name = "year",
                    value.name    = "proportion")
datLong_obs
colnames(datLong_obs) = c("CoD", "year", "proportion")

g <- ggplot(data = datLong_obs, aes(x = CoD, y = proportion, group = year, colour = year))

theme_set(theme_bw())

observed_m <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  #scale_x_continuous(label = CoD.label)
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Observed") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none")

plot_data_clr <- pro.dea.clr.LC[,1:16]

colnames(plot_data_clr) = years
rownames(plot_data_clr) = CoD.label

datLong_clr <- melt(data          = plot_data_clr,
                    id.vars       = c(CoD.label),
                    measure.vars  = c(years),
                    variable.name = "year",
                    value.name    = "proportion")
datLong_clr
colnames(datLong_clr) = c("CoD", "year", "proportion")

g <- ggplot(data = datLong_clr, aes(x = CoD, y = proportion, group = year, colour = year))

theme_set(theme_bw())

clr_m <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("CLR") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none")


plot_data_ilr <- pro.dea.ilr.LC[,1:16]

colnames(plot_data_ilr) = years
rownames(plot_data_ilr) = CoD.label

datLong_ilr <- melt(data          = plot_data_ilr,
                    id.vars       = c(CoD.label),
                    measure.vars  = c(years),
                    variable.name = "year",
                    value.name    = "proportion")
datLong_ilr
colnames(datLong_ilr) = c("CoD", "year", "proportion")

g <- ggplot(data = datLong_ilr, aes(x = CoD, y = proportion, group = year, colour = year))

theme_set(theme_bw())

ilr_m <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("ILR") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none")


plot_data_alpha <- pro.dea.alpha.LC[,1:16]

colnames(plot_data_alpha) = years
rownames(plot_data_alpha) = CoD.label

datLong_alpha <- melt(data          = plot_data_alpha,
                      id.vars       = c(CoD.label),
                      measure.vars  = c(years),
                      variable.name = "year",
                      value.name    = "proportion")
datLong_alpha
colnames(datLong_alpha) = c("CoD", "year", "proportion")

g <- ggplot(data = datLong_alpha, aes(x = CoD, y = proportion, group = year, colour = year))

theme_set(theme_bw())

alpha_m <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Alpha") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14))

savepdf("comparison_m", width = 32, height = 12, toplines = 0.8)
par(mar=c(5, 4, 4, 8), xpd=TRUE)
observed_m + clr_m + ilr_m + alpha_m + 
  plot_layout(ncol = 4)
dev.off()

# plot 2: time series by cause comparing alpha and observed

# CLR Data
plot2_data_clr <- pro.dea.clr.LC
colnames(plot2_data_clr) = total_years
rownames(plot2_data_clr) = CoD.label
datLong2_clr <- melt(data          = plot2_data_clr,
                     id.vars       = c(CoD.label),
                     measure.vars  = c(total_years),
                     variable.name = "year",
                     value.name    = "proportion")
colnames(datLong2_clr) = c("CoD", "year", "proportion")
# ILR Data
plot2_data_ilr <- pro.dea.ilr.LC
colnames(plot2_data_ilr) = total_years
rownames(plot2_data_ilr) = CoD.label
datLong2_ilr <- melt(data          = plot2_data_ilr,
                     id.vars       = c(CoD.label),
                     measure.vars  = c(total_years),
                     variable.name = "year",
                     value.name    = "proportion")
colnames(datLong2_ilr) = c("CoD", "year", "proportion")
# Alpha Data
plot2_data_alpha <- pro.dea.alpha.LC
colnames(plot2_data_alpha) = total_years
rownames(plot2_data_alpha) = CoD.label
datLong2_alpha <- melt(data          = plot2_data_alpha,
                       id.vars       = c(CoD.label),
                       measure.vars  = c(total_years),
                       variable.name = "year",
                       value.name    = "proportion")
colnames(datLong2_alpha) = c("CoD", "year", "proportion")

# Select causes to include in plot
include.causes <- c("Rheumatic","Hypertension", "Hypertensive disease","Myocardial infarction","Other IHD","Pulmonary", "Non-rheumatic", "Cardiac arrest", "Heart failure", "Other cardio")

# Create plot
theme_set(theme_bw())

clr_proj_m <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong2_clr, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_obs, CoD %in% include.causes)) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  scale_color_brewer(palette = "Paired") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(2001, 2026, 2)) +
  ggtitle("CLR") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none")

ilr_proj_m <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong2_ilr, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_obs, CoD %in% include.causes)) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  scale_color_brewer(palette = "Paired") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(2001, 2026, 2)) +
  ggtitle("ILR") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none")

alpha_proj_m <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong2_alpha, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_obs, CoD %in% include.causes)) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  scale_color_brewer(palette = "Paired") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(2001, 2026, 2)) +
  ggtitle("Alpha") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14))

# Group plots together
savepdf("all_proj_m", width = 32, height = 12, toplines = 0.8)
par(mar=c(2, 2, 2, 2), xpd=TRUE)
clr_proj_m + ilr_proj_m + alpha_proj_m + 
  plot_layout(ncol = 3)
dev.off()