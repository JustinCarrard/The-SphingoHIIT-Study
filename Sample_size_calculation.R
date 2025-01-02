#==========================================================================================
#This R-code was designed to calculate the sample size of the SphingoHIIT study.
#More information is available on ClinicalTrials.gov (ID NCT05390866) and 
#in the study protocol (https://doi.org/10.12688/f1000research.128978.1).

#==========================================================================================
#------------------------------------------------------------------------------------------
# Set paths
#------------------------------------------------------------------------------------------

setwd(" ") # Main path

dat_path <- "./data"
graphics_path <- "./output/graphics"
text_path <- "./output/text"
scripts_path <- "./scripts"

#------------------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------------------

library(emmeans)
library(lme4)
library(nlme)
library(ggplot2)
library(tidyverse)
library(visreg)
library(readxl)

#------------------------------------------------------------------------------------------
# Sample size
#------------------------------------------------------------------------------------------

dat <- read_xlsx(paste0(dat_path, "/Ceramide Data.xlsx"), sheet = 2)

# Set drop outs
dropout <- 0.0

#log2_transformation_lipid_subclasses_and_clinical_lipids
(variables_to_transform <- names(dat)[2:5]) # Name of variables to transform with log2

for (i in variables_to_transform) {
  dat[, paste0(i, "_log2")] <- log2(dat[, i])
}

names(dat)

head(dat)

# Calculate standard deviations
names(dat)

sigmas <- dat %>%
  summarise(across(.cols = names(dat)[6:9], .fns = sd)) %>%
  unlist()

assumed_sd <- mean(sigmas) # Mean standard deviation

diff_normalscale <- log2(1.194)

rho <- seq(0.05, 0.95, 0.05)
pow <- c(0.8)

required_n_ancova <- sapply(rho, FUN = function(x) {
  sapply(pow, function(y) {
    ceiling(power.t.test(
      delta = diff_normalscale
      , sd = sqrt(assumed_sd^2*(1 - x^2))
      , sig.level = 0.05
      , power = y
      , type = "two.sample"
    )$n*(1/(1 - dropout)))
  })
})

plot_dat <- data.frame(
  rho = rep(rho, times = length(pow))
  , reqN = as.vector(t(required_n_ancova))
  , power = factor(paste0(100*rep(pow, each = length(rho)), "%"))
)

theme_set(theme_bw())
p <- ggplot(data = plot_dat, aes(x = rho, y = reqN, group = power)) +
  geom_point(aes(colour = power), size = 4) +
  geom_line(aes(colour = power)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(
    y = "Required N per group"
    , x = expression("Correlation coefficient "*rho*" between pre and post measurements")
    , title = paste0("Assumed standard deviation: ", signif(assumed_sd, 3), "; assumed geometric mean ratio: ", round(2^(diff_normalscale), 2))
  ) +
  scale_colour_manual(name = "Power", values = c("#008FD0", "#F07E00")) +
  theme(
    axis.title.y=element_text(colour = "black", size = 17, hjust = 0.5, margin=margin(0,12,0,0)),
    axis.title.x=element_text(colour = "black", size = 17, margin=margin(12,0,0,0)),
    axis.text.x=element_text(colour = "black", size=15),
    axis.text.y=element_text(colour = "black", size=15),
    legend.position="right",
    legend.text = element_text(size=15),
    legend.title = element_text(size = 15),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(face = "bold"),
    strip.background=element_rect(fill="white"),
    strip.text.x=element_text(size=12)
  )

p

ggsave(paste0(graphics_path, "/samplesize_", round(assumed_sd, 3), ".tiff"), p, width = 19*0.6, height = 12*0.6, units = "in", dpi = 300)












