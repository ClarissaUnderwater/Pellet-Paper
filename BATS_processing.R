
# Aug 21 2025 

# contains data plotting and statistics for respiration, sinking, images
# input: ST1_BATS_pellets.xlsx

##### libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggpmisc)
library(forcats)
library(readxl)
library(patchwork)

#######################################################################################################
# load and prepare data
#######################################################################################################

##### choose working directory
setwd("your/location")

# read the data sheet
Pellets <- read_excel("ST1_BATS_pellets.xlsx", sheet = "all_pellets"
                      , na = c("NA", "na", "", " "))
Pellets <- Pellets[-1, ]  # Remove the second row containing units


# convert numbers to numeric values
Pellets[] <- lapply(Pellets, function(x) {
  if (is.character(x) || is.factor(x)) {
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    if (all(!is.na(x_num) | is.na(x))) {
      return(x_num)
    } else {
      return(x)  # keep original if not safely numeric
    }
  } else {
    return(x)
  }
})

# Convert to numeric (non-numeric entries become NA)
Pellets$sinking_numeric <- as.numeric(Pellets$`sinking velocity`)
Pellets$sinking_numeric[is.infinite(Pellets$sinking_numeric)] <- NA
#######################################################################################################
### graphics settings
#######################################################################################################

taxon_colors <- c(
  "Salp big" = "#5C351B",
  "Salp" = "#B79600",
  "Copepod" = "#00BFF2",
  "Pteropod" = "#69C380",
  "Snail" = "#F48100",
  "Euphausid" = "purple3",
  "Shrimp" = "#C10000"
)


taxon_shapes <- c(
  "Salp big" = 8,     # asterisk
  "Salp" = 8,         # asterisk
  "Copepod" = 16,     # filled circles
  "Pteropod" = 15,    # filled square
  "Snail" = 22,     # square with cross
  "Euphausid" = 17,   # triangle
  "Shrimp" = 3        # cross
)


# Assign colors based on Taxon
Pellets$txcolor <- taxon_colors[Pellets$Taxon]
Pellets$txshape <- taxon_shapes[Pellets$Taxon]

# add taxon as factor and order pellets by taxon
Pellets$Taxon <- factor(Pellets$Taxon, levels = c("Salp big", "Salp", "Snail", "Pteropod","Copepod", "Euphausid", "Shrimp"))
Pellets$animal_ID <- with(Pellets,fct_reorder(animal_ID, as.numeric(Taxon)))

# Define coarse taxa
Taxon2 <- c("Salp big" = "Salp", "Salp" = "Salp", "Snail" = "Gastropod", "Pteropod" = "Gastropod",
            "Copepod" = "Copepod","Euphausid" = "Eucarida", "Shrimp" = "Eucarida")

# Assign coarse taxa
Pellets$Taxon2 <- Taxon2[Pellets$Taxon]

# order degradation levels
Pellets$degradation <- factor(Pellets$degradation, levels = c("intact", "degraded", "casing"))




#######################################################################################################
# calculating porosities and modeling correlations with size
#######################################################################################################

###############
############### Constants
###############

g <- 9.81 # gravity (m/s^2)
rho_f <- 1025 # seawater density (kg/m³) at 21°C and North Atlantic conditions
mu <- 0.001 # dynamic viscosity of seawater (Pa·s) at 21°C
rho_s <- 1650 # assumed solid particle density (kg/m³), based on Plough et al.

###############
############### pellet radius and length
###############

###### Salp Durkin
Pellets$width_SP_Durkin <- 0.63 * Pellets$ESD 
Pellets$length_SP_Durkin <- (pi * (Pellets$ESD / 2)^2) / Pellets$width_SP_Durkin
Pellets$volume_SP <- Pellets$length_SP_Durkin * (Pellets$width_SP_Durkin^2) / 4
Pellets$volume_SP <- ifelse( Pellets$short != "SP", Pellets$volume, Pellets$volume_SP)

####### cylinder Laurenceau-Cornec
# Step 1: Solve for radius `r` using the quadratic formula
# 4r^2 - P*r + A = 0  -> solving for r
# Roots are: (P ± sqrt(P^2 - 16 * A)) / 8
Pellets$radius_LC <- (Pellets$perimeter - sqrt(Pellets$perimeter^2 - 16 * Pellets$`2D area`)) / 8  
# Step 2: Calculate length `L` using `L = A / (2 * r)`
Pellets$length_LC <-  Pellets$`2D area` / (2 * Pellets$radius_LC)
# Step 3: Calculate volume as `V = L * pi * r^2`
Pellets$volume_cylinder_LC <- Pellets$length_LC * pi * (Pellets$radius_LC)^2

Pellets$volume_SP <- ifelse( Pellets$shape == "cylinder", Pellets$volume_cylinder_LC, Pellets$volume_SP)

####### radius in m
Pellets$radius <- ifelse( Pellets$shape == "cylinder", Pellets$radius_LC/1000, 
                          ifelse( Pellets$shape == "salp", Pellets$width_SP_Durkin/2000, Pellets$ESD/2000) )

Pellets$length <- ifelse( Pellets$shape == "cylinder", Pellets$length_LC/1000, 
                          ifelse( Pellets$shape == "salp", Pellets$length_SP_Durkin/1000, Pellets$ESD/1000) )

Pellets$aspectratio <- Pellets$length_LC/Pellets$radius_LC/2
Pellets$aspectratio_SP <- Pellets$length_SP_Durkin/Pellets$width_SP_Durkin
Pellets$aspectratio <- ifelse( Pellets$shape == "salp", Pellets$aspectratio_SP, Pellets$aspectratio)

Pellets$circularity <- 4 * pi * Pellets$`2D area` / (Pellets$perimeter^2)
Pellets$elongation <- Pellets$perimeter^2 / Pellets$`2D area`

###############
####### traditional Stoke's law
###############

# Calculate effective particle density (rho_p)
# converting sinking speed from m/day to m/s
Pellets$rho_p <- (9 * mu * Pellets$sinking_numeric / (24 * 60 * 60)) / (2 * g * Pellets$radius^2) + rho_f
median(Pellets$rho_p[!is.na(Pellets$rho_p) & Pellets$degradation != "casing"])

# Calculate porosity (phi)
Pellets$porosity_plough <- (rho_s - Pellets$rho_p) / (rho_s - rho_f)
median(Pellets$porosity_plough[!is.na(Pellets$porosity_plough) & Pellets$degradation != "casing"])

###############
####### adapted Stoke's law based on Komar
###############
# v = 0.0790 * (1 / mu) * (rho_s - rho_f) * g * L^2 *(L/D)^(-1.664)

# Calculate effective particle density (rho_p)
Pellets$rho_p_Komar <- Pellets$sinking_numeric / (24 * 60 * 60) / 0.0790 * mu / g / (Pellets$length)^2 / ((Pellets$length) /(Pellets$radius*2))^(-1.664) + rho_f
mean(Pellets$rho_p_Komar[!is.na(Pellets$rho_p_Komar) & Pellets$degradation != "casing" ])

Pellets$porosity_komar <- (rho_s - Pellets$rho_p_Komar) / (rho_s - rho_f)
mean(Pellets$porosity_komar[!is.na(Pellets$porosity_komar)])

###############
# Model Correlations: Porosity vs volume
###############

# Filter data for modeling
df_fit <- subset(Pellets,degradation != "casing" & !is.na(porosity_komar) & !is.na(volume_SP) &
                   is.finite(porosity_komar) & is.finite(volume_SP) &
                   porosity_komar > 0 & porosity_komar < 1 &
                   volume_SP > 0)

##### Exponential model
exp_porosity <- nls(porosity_komar ~ phi_max * (1 - exp(-k * volume_SP)), data = df_fit, start = list(phi_max = 0.99, k = 1))
summary(exp_porosity)  # Equation: porosity_komar = phi_max * (1 - exp(-k * volume_SP))

df_fit$predicted_exp <- predict(exp_porosity) # Predicted values from exp model

RSS <- sum((df_fit$porosity_komar - df_fit$predicted)^2, na.rm = TRUE) # Residual Sum of Squares
TSS <- sum((df_fit$porosity_komar - mean(df_fit$porosity_komar, na.rm = TRUE))^2, na.rm = TRUE) # Total Sum of Squares
R2_nls <- 1 - RSS / TSS # Pseudo-R²
R2_nls

##### Log-linear model
log_model <- lm((porosity_komar) ~ (log(volume_SP)), data = df_fit)
summary(log_model)
df_fit$predicted_log <- predict(log_model) # Predict values from log model

##### Hill model with lower/upper bounds
hill_model <- nls(
  porosity_komar ~ phi_max * (volume_SP^n / (K^n + volume_SP^n)),
  data = df_fit,
  start = list(phi_max = 0.9, K = 5, n = 2),  # safer start values
  algorithm = "port",  # allows bounds
  lower = c(phi_max = 0.5, K = 0.01, n = 0.1),
  upper = c(phi_max = 1.0, K = 100, n = 10)
)
summary(hill_model)

df_fit$predicted_hill <- predict(hill_model) # Predict values from Hill model

###### Compare model AICs
AIC(exp_porosity, log_model, hill_model, fit_power_log, fit_power_nls)


ggplot() +
  geom_point(data = df_fit, aes(x = predicted_exp, y = porosity_komar),alpha = 0.3) +
  geom_point(data = df_fit, aes(x = predicted_log, y = porosity_komar),alpha = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Predicted porosity",
    y = "Observed porosity",
    title = "Model Fit: Observed vs. Predicted log(Sinking Velocity)"
  ) +
  theme_minimal()


ggplot() +
  geom_line(data = df_fit, aes(x = volume_SP, y = predicted_exp), color = "blue", size = 0.5) +
  geom_line(data = df_fit, aes(x = volume_SP, y = predicted_log), color = "black", size = 0.5) +
  geom_point(data = subset(Pellets, degradation != "casing"), 
             aes(x = (volume_SP), y = porosity_komar, shape = Taxon)) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  scale_x_log10() +
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 5), axis.text.y = element_text(size = 5), 
        legend.position = "none")
ggsave("porosity_1r.jpg", width=4.5, height=2, dpi=300)

# average by animal or Taxon
ggplot(data = subset(Pellets, degradation != "casing"), 
               aes(x = animal_ID, y = porosity_komar, color = Taxon, shape = Taxon)) +
  geom_boxplot(aes(color = Taxon), alpha = 0.5, size = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 1, dodge.width = 0.75), 
             size = 1) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +
  scale_fill_manual(values = taxon_colors) +  # matches box fill to line color
  ylab(expression("Respiration ("*O[2]*"/pellet day"^-1*")"))+
  xlab(NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

###############
# Model Correlations: predicting sinking velocity
###############
# porosity_plough # from particle sinking based on Stokes law and Plougs solid density
# porosity_komar # from particle sinking based on Komars law and Plougs solid density
# porosity_tx # taxon2 averages
# predicted_exp # exponential function of porosity with size

# averages by taxon
Pellets$porosity_tx <- ave(Pellets$porosity_komar, Pellets$Taxon2, FUN = function(x) mean(x, na.rm = TRUE))
Pellets$rho_p_tx <- ave(Pellets$rho_p, Pellets$Taxon2, FUN = function(x) mean(x, na.rm = TRUE))
df_fit$porosity_tx <- ave(df_fit$porosity_komar, df_fit$Taxon2, FUN = function(x) mean(x, na.rm = TRUE))
df_fit$rho_p_tx <- ave(df_fit$rho_p, df_fit$Taxon2, FUN = function(x) mean(x, na.rm = TRUE))

# predicted from porosity-size correlation exponential version
df_fit$predicted_exp
df_fit$rho_p_exp <- -((df_fit$predicted_exp * (rho_s - rho_f))-rho_s) #porosity_komar <- (rho_s - rho_p_Komar) / (rho_s - rho_f)

# calculate sinking velocities from these porosities and each pellet's geometry
df_fit$komar_velocity_exp <- 0.0790 * (1 / mu) *
  (df_fit$rho_p_exp - rho_f) * g * (df_fit$length/1000)^2 *
  ((df_fit$length/1000) / (df_fit$radius/500))^(-1.664)

df_fit$komar_velocity_tx <- 0.0790 * (1 / mu) *
  (df_fit$rho_p_tx - rho_f) * g * (df_fit$length/1000)^2 *
  ((df_fit$length/1000) / (df_fit$radius/500))^(-1.664)


## models

# 1. baseline model: size only
fit_size <- lm(log(sinking_numeric) ~ log(`2D area`), data = df_fit)
summary(fit_size)
df_fit$sinking_size <- predict(fit_size)

# 2. taxon-specific porosity
fit_size_tx <- lm(log(sinking_numeric) ~ log(`2D area`) + porosity_tx, data = df_fit)
df_fit$sinking_size_tx <- predict(fit_size_tx)

# 3. exponential fit for porosity
fit_size_exp <- lm(log(sinking_numeric) ~ log(`2D area`) + predicted_exp, data = df_fit)
df_fit$sinking_size_exp <- predict(fit_size_exp)

AIC(fit_size, fit_size_tx, fit_size_exp)


# 4. add shape factors
fit_shape_all <- lm(log(sinking_numeric) ~ log(`2D area`) + circularity + aspectratio + predicted_exp, data = df_fit)
df_fit$sinking_size_shape <- NA  # Create a new column with all NAs
df_fit$sinking_size_shape[!is.na(df_fit$sinking_numeric) &
                            !is.na(df_fit$`2D area`) &
                            !is.na(df_fit$circularity) &
                            !is.na(df_fit$aspectratio) &
                            !is.na(df_fit$predicted_exp)] <- predict(fit_shape_all)
summary(fit_shape_all)
AIC(fit_shape_all)


ggplot() +
  geom_point(data = df_fit, aes(x = log(`2D area`), y = log(sinking_numeric),  shape = Taxon)) +
  geom_line(data = df_fit, aes(x = log(`2D area`), y = sinking_size,  shape = Taxon), color = "red") +
  geom_point(data = df_fit, aes(x = log(`2D area`), y = sinking_size_exp,  shape = Taxon), color = "purple") +
  geom_point(data = df_fit, aes(x = log(`2D area`), y = sinking_size_tx,  shape = Taxon), color = "blue") +
  geom_point(data = df_fit, aes(x = log(`2D area`), y = sinking_size_shape,  shape = Taxon), color = "green") +
  geom_point(shape = 21, size = 1) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  ylab(expression("log(Sinking velocity (m day"^-1*"))"))+
  xlab(expression("log(2D Area (mm"^2*"))"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5f.jpg", width=3, height=3, dpi=300)


ggplot() +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = log(sinking_numeric),  shape = Taxon)) +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = sinking_size,  shape = Taxon), color = "red") +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = sinking_size_exp,  shape = Taxon), color = "purple") +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = sinking_size_tx,  shape = Taxon), color = "blue") +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = sinking_size_shape,  shape = Taxon), color = "green") +
  geom_point(shape = 21, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ylab(expression("log(Sinking velocity (m day"^-1*"))"))+
  xlab(expression("log(Measured sinking velocity (m day"^-1*"))"))+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5f.jpg", width=3, height=3, dpi=300)

#######################################################################################################
# Respiration means
#######################################################################################################
mean(subset(Pellets$Respiration, Pellets$degradation != "casing" & !is.na(Pellets$Respiration)))
sd(subset(Pellets$Respiration, Pellets$degradation != "casing" & !is.na(Pellets$Respiration)))

mean(subset(Pellets$`Carbon-specific respiration`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon-specific respiration`)))
sd(subset(Pellets$`Carbon-specific respiration`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon-specific respiration`)))

mean(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`)))/12
sd(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`)))/12

### carbon content means
mean(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon == "Salp big"))
sd(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon == "Salp big"))

mean(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon == "Salp"))
sd(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon == "Salp"))

mean(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Eucarida"))
sd(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Eucarida"))

mean(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Gastropod"))
sd(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Gastropod"))

mean(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Copepod"))
sd(subset(Pellets$`Carbon content`, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Copepod"))

### carbon concentration means   
mean(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon == "Salp"))
sd(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon == "Salp"))

mean(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon == "Salp big"))
sd(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon == "Salp big"))

mean(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Eucarida"))
sd(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Eucarida"))

mean(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Gastropod"))
sd(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Gastropod"))

mean(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Copepod"))
sd(subset(Pellets$`Carbon content`/Pellets$volume_SP, Pellets$degradation != "casing" & !is.na(Pellets$`Carbon content`) & Pellets$Taxon2 == "Copepod"))


#######################################################################################################
# Figure 4 - by taxon 
#######################################################################################################
FT1 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
              aes(x = Taxon, y = `Carbon-specific respiration`)) +
  geom_boxplot(aes(fill = Taxon), alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("Respiration ("*O[2]*"/C day"^-1*")"))+
  scale_color_manual(values = taxon_colors) +
    scale_fill_manual(values = taxon_colors) +
    scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FT2 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
               aes(x = Taxon, y = Respiration/volume_SP)) +
  geom_boxplot(aes(fill = Taxon),alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  #ylim(0,8)+
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Respiration ("*O[2]*"/mm"^3*" day"^-1*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FT5 <- ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
       aes(x = Taxon, y = `Carbon content` / volume_SP)) +
  geom_boxplot(aes(fill = Taxon), alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(color = Taxon, shape = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.5), size = 1.25) +
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Carbon concentration ("*mu*"g mm"^-3*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5c.jpg", width=3, height=3, dpi=300)

FT7 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
       aes(x = Taxon, y = sinking_numeric)) +
  geom_boxplot(aes(fill = Taxon), alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(color = Taxon, shape = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.5), size = 1.25) +
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Sinking velocity (m day"^-1*")"))+
  scale_fill_manual(values = taxon_colors) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.x = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5e.jpg", width=3, height=3, dpi=300)

FT6 <-
ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
       aes(x = Taxon, y = sinking_numeric * `Carbon content`)) +
  geom_boxplot(aes(fill = Taxon), alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Carbon flux ("*mu*"g C day"^-1*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FT3 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Taxon, y = porosity_komar)) +
  geom_boxplot(aes(fill = Taxon), alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("Porosity (%)"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)


FT4 <- 
ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
       aes(x = Taxon, y = aspectratio)) +
  geom_boxplot(aes(fill = Taxon), alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  #ylim(0,0.17)+
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Elongation (aspect ratio)"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)


FT8 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Taxon , y = `Carbon content`)) +
  geom_boxplot(aes(fill = Taxon), alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  #ylim(0,0.17)+
  xlab(NULL)+
  scale_y_log10()+
    ylab(expression("Carbon content ("*mu*"g C)"))+
   scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FTall <- (FT8 | FT7 | FT6 | FT5 | FT1 | FT2 | FT4 | FT3) + plot_annotation(tag_levels = 'a') + plot_layout(ncol = 4, nrow = 2)
print(FTall) # Show the combined figure
ggsave("FTall.jpg", FTall, width = 7.1, height = 6, dpi = 900) # Save the combined figure
#######################################################################################################
# Figure 7 and S8 - by degradatation
#######################################################################################################

# making this not plot a boxplot if there are <3 data points
Pellets_filtered <- Pellets %>% filter(!is.na(`Carbon content`)) %>% group_by(degradation) %>% filter(n() >= 3)  # keep only groups with 5 or more data points

FD1 <- 
  ggplot(data = subset(Pellets), 
         aes(x = degradation, y = `Carbon-specific respiration`)) +
  geom_boxplot(data = Pellets_filtered, alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("Respiration ("*O[2]*"/C day"^-1*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FD2 <- 
  ggplot(data = subset(Pellets), 
         aes(x = degradation, y = Respiration/volume_SP)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  #ylim(0,8)+
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Respiration ("*O[2]*"/mm"^3*" day"^-1*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)


FD3 <- 
  ggplot(data = subset(Pellets), 
              aes(x = degradation, y = `Carbon content` / volume_SP)) +
  geom_boxplot(data = Pellets_filtered, alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(color = Taxon, shape = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.5), size = 1.25) +
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Carbon concentration ("*mu*"g mm"^-3*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.x = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5c.jpg", width=3, height=3, dpi=300)


FD4 <- 
  ggplot(data = subset(Pellets), 
         aes(x = degradation, y = porosity_komar)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("Porosity (%)"))+
  #facet_wrap(~Taxon2)+
  #scale_y_log10()+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "right")
#ggsave("Supplement_degr_taxa.jpg", width=5.5, height=5.5, dpi=300)


#FD5 <- 
  ggplot(data = subset(Pellets), 
         aes(x = degradation, y = sinking_numeric)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("Porosity (%)"))+
  facet_wrap(~Taxon2)+
  #scale_y_log10()+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "right")
#ggsave("Supplement_degr_taxa.jpg", width=5.5, height=5.5, dpi=300)


FTD <- (FD1 | FD3 | FD4 ) + plot_annotation(tag_levels = 'a') + plot_layout(ncol = 3)
print(FTD) # Show the combined figure
ggsave("FTD.jpg", FTD, width = 7.1, height = 3, dpi = 900) # Save the combined figure



FDS1 <- 
  ggplot(data = subset(Pellets), 
         aes(x = degradation, y = porosity_komar)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("Porosity (%)"))+
  facet_wrap(~Taxon2, ncol=4)+
  #scale_y_log10()+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Supplement_degr_taxa.jpg", width=5.5, height=5.5, dpi=300)


FSD2 <- 
ggplot(data = subset(Pellets), 
       aes(x = degradation, y = sinking_numeric)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("sinking velocity (m day"^-1*")"))+
  facet_wrap(~Taxon2, ncol=4)+
  #scale_y_log10()+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Supplement_degr_taxa.jpg", width=5.5, height=5.5, dpi=300)

FSD <- (FDS1 | FSD2 ) + plot_annotation(tag_levels = 'a') + plot_layout(nrow = 2)
print(FSD) # Show the combined figure
ggsave("FSD.jpg", FSD, width = 7.1, height = 6, dpi = 900) # Save the combined figure

#######################################################################################################
# Figure S7 - by station
#######################################################################################################
FS1 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Station, y = `Carbon-specific respiration`)) +
  geom_boxplot( alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("Respiration ("*O[2]*"/C day"^-1*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FS2 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Station, y = Respiration/volume_SP)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  #ylim(0,8)+
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Respiration ("*O[2]*"/mm"^3*" day"^-1*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FS5 <- ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
              aes(x = Station, y = `Carbon content` / volume_SP)) +
  geom_boxplot( alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(color = Taxon, shape = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.5), size = 1.25) +
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Carbon concentration ("*mu*"g mm"^-3*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.x = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5c.jpg", width=3, height=3, dpi=300)

FS7 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Station, y = sinking_numeric)) +
  geom_boxplot(, alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(color = Taxon, shape = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 0.5), size = 1.25) +
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Sinking velocity (m day"^-1*")"))+
  scale_fill_manual(values = taxon_colors) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.x = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5e.jpg", width=3, height=3, dpi=300)

FS6 <-
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Station, y = sinking_numeric * `Carbon content`)) +
  geom_boxplot( alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Carbon flux ("*mu*"g C day"^-1*")"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FS3 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Station, y = porosity_komar)) +
  geom_boxplot( alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  xlab(NULL)+
  ylab(expression("Porosity (%)"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)


FS4 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Station, y = aspectratio)) +
  geom_boxplot( alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  #ylim(0,0.17)+
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Elongation (aspect ratio)"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)


FS8 <- 
  ggplot(data = subset(Pellets, Pellets$degradation!="casing"), 
         aes(x = Station , y = `Carbon content`)) +
  geom_boxplot( alpha = 0.3, outlier.shape = NA, size = 0.3) +
  geom_point(aes(shape = Taxon, color = Taxon), position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 1.25) +
  #ylim(0,0.17)+
  xlab(NULL)+
  scale_y_log10()+
  ylab(expression("Carbon content ("*mu*"g C)"))+
  scale_color_manual(values = taxon_colors) +
  scale_fill_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        #axis.text.y = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("C_spec_resp.pdf", width=5.5, height=3, dpi=300)

FSall <- (FS8 | FS7 | FS6 | FS5 | FS1 | FS2 | FS4 | FS3) + plot_annotation(tag_levels = 'a') + plot_layout(ncol = 4, nrow = 2)
print(FSall) # Show the combined figure
ggsave("FSall.jpg", FSall, width = 7.1, height = 6, dpi = 900) # Save the combined figure
#######################################################################################################
# Statistics: carbon-specific respiration by taxon, station, degradation
#######################################################################################################

### Respiration by taxa

X1 <- subset(Pellets$`Carbon-specific respiration`, Pellets$Taxon == "Snail" & !is.na(Pellets$`Carbon-specific respiration`))
X2 <- subset(Pellets$`Carbon-specific respiration`, Pellets$Taxon != "Snail" & !is.na(Pellets$`Carbon-specific respiration`))
t.test(X1,X2)

X1 <- subset(Pellets$`Carbon-specific respiration`, Pellets$Taxon == "Shrimp" & !is.na(Pellets$`Carbon-specific respiration`))
X2 <- subset(Pellets$`Carbon-specific respiration`, Pellets$Taxon != "Shrimp" & !is.na(Pellets$`Carbon-specific respiration`))
t.test(X1,X2)

X1 <- subset(Pellets$`Carbon-specific respiration`, Pellets$Taxon2 == "Salp" & !is.na(Pellets$`Carbon-specific respiration`))
X2 <- subset(Pellets$`Carbon-specific respiration`, Pellets$Taxon2 != "Salp" & !is.na(Pellets$`Carbon-specific respiration`))
t.test(X1,X2)

X1 <- subset(Pellets$`Carbon-specific respiration`, Pellets$Taxon == "Copepod" & !is.na(Pellets$`Carbon-specific respiration`))
X2 <- subset(Pellets$`Carbon-specific respiration`, Pellets$Taxon != "Copepod" & !is.na(Pellets$`Carbon-specific respiration`))
t.test(X1,X2)

# Run pairwise t-tests between all Taxon levels
ttest_results <- pairwise.t.test(Pellets$`Carbon-specific respiration`, 
                                 Pellets$Taxon,
                                 p.adjust.method = "BH", # or "none", "bonferroni", etc.
                                 pool.sd = FALSE, # Welch t-test (unequal variance)
                                 paired = FALSE)

# Print the p-value matrix
ttest_results$p.value


### intact vs degraded

# absolute respiration
X1 <- subset(Pellets$Respiration, Pellets$degradation == "intact" & !is.na(Pellets$Respiration))
X2 <- subset(Pellets$Respiration, Pellets$degradation == "degraded" & !is.na(Pellets$Respiration))
t.test(X1,X2)

# C specific all
X1 <- subset(Pellets$`Carbon-specific respiration`, Pellets$degradation == "intact" & !is.na(Pellets$`Carbon-specific respiration`))
X2 <- subset(Pellets$`Carbon-specific respiration`, Pellets$degradation == "degraded" & !is.na(Pellets$`Carbon-specific respiration`))
t.test(X1,X2)

# C-specific only salp, pteropod, euphausiid
X1 <- subset(Pellets$`Carbon-specific respiration`, Pellets$degradation == "intact" & !is.na(Pellets$`Carbon-specific respiration`) &  (Pellets$Taxon == "Euphausid" | Pellets$Taxon2 == "Salp" | Pellets$Taxon == "Pteropod"))
X2 <- subset(Pellets$`Carbon-specific respiration`, Pellets$degradation == "degraded" & !is.na(Pellets$`Carbon-specific respiration`) & (Pellets$Taxon == "Euphausid" | Pellets$Taxon2 == "Salp" | Pellets$Taxon == "Pteropod"))
t.test(X1,X2)

# volumetric respiration
X1 <- subset(Pellets$Respiration / Pellets$volume_SP, Pellets$degradation == "intact" & !is.na(Pellets$Respiration))
X2 <- subset(Pellets$Respiration / Pellets$volume_SP, Pellets$degradation == "degraded" & !is.na(Pellets$Respiration))
t.test(X1,X2)

#######################################################################################################
# Figure S6 absolute and volumetric respiration by animal
#######################################################################################################

FS6a <- ggplot(data = subset(Pellets, degradation != "casing"), 
               aes(x = animal_ID, y = Respiration, color = Taxon, shape = Taxon)) +
  geom_boxplot(aes(color = Taxon), alpha = 0.5, size = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 1, dodge.width = 0.75), 
             size = 1) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +
  scale_fill_manual(values = taxon_colors) +  # matches box fill to line color
  ylab(expression("Respiration ("*O[2]*"/pellet day"^-1*")"))+
  xlab(NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.x = element_blank(),
        plot.tag = element_text(size = 9),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

FS6b <- ggplot(data = subset(Pellets, degradation != "casing"), 
               aes(x = animal_ID, y = volume_SP, color = Taxon, shape = Taxon)) +
  geom_boxplot(aes(color = Taxon), alpha = 0.5, size = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 1, dodge.width = 0.75), 
             size = 1) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +
  scale_fill_manual(values = taxon_colors) +  # matches box fill to line color
  ylab(expression("Volume (mm"^3*")"))+
  xlab(NULL)+
  scale_y_log10() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        plot.tag = element_text(size = 9),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

FS6c <- ggplot(data = subset(Pellets, degradation != "casing"), 
               aes(x = animal_ID, y = Respiration/volume_SP, color = Taxon, shape = Taxon)) +
  geom_boxplot(aes(color = Taxon), alpha = 0.5, size = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 1, dodge.width = 0.75), 
             size = 1) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +
  scale_fill_manual(values = taxon_colors) +  # matches box fill to line color
  ylab(expression("Respiration ("*O[2]*"/mm"^3*" day"^-1*")"))+
  #ylim(0,11)+
  xlab(NULL)+
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 9),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

FS6d <- ggplot(data = subset(Pellets, degradation != "casing"), 
               aes(x = animal_ID, y = sinking_numeric, color = Taxon, shape = Taxon)) +
  geom_boxplot(aes(color = Taxon), alpha = 0.5, size = 0.3) +
  geom_point(position = position_jitterdodge(jitter.width = 1, dodge.width = 0.75), 
             size = 1) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +
  scale_fill_manual(values = taxon_colors) +  # matches box fill to line color
  ylab(expression("Sinking Velocity (m day"^-1*")"))+
  #xlab("Animal ID") +
  #ylim(0,11)+
  xlab(NULL)+
  scale_y_log10() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 9),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

FS6 <- (FS6a | FS6b | FS6c | FS6d) + plot_annotation(tag_levels = 'a') + plot_layout(nrow = 4)

print(FS6) # Show the combined figure
#ggsave("FS6.pdf", FS6, width = 7.1, height = 3, dpi = 300) # Save the combined figure
ggsave("FS6.jpg", FS6, width = 7.5, height = 9, dpi = 1200) # Save the combined figure


#######################################################################################################
# Figure 6 respiration against carbon, volume, perimeter
#######################################################################################################

F4a <- 
  ggplot(data = subset(Pellets, degradation!="casing"), 
         aes(x = `Carbon content`, y = Respiration)) +
  geom_errorbar(data = subset(Pellets, degradation!="casing"),
                aes(x = `Carbon content`, y = Respiration, 
                    ymin = Respiration - `Standard Error`, ymax = Respiration + `Standard Error`), size = 0.1) +
  geom_point(data = subset(Pellets, degradation != "casing"), 
             aes(x =`Carbon content`, y = Respiration, color = Taxon,  shape = Taxon)) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  ylab(expression("Respiration ("*mu*"mol "*O[2]* " day"^-1*")"))+
  xlab(expression("Carbon content ("*mu*"g C)"))+
  scale_x_log10() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.tag = element_text(size = 9),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = c(0.02, 0.98),  # near top-left inside
        legend.justification = c(0, 1),    # anchor top-left corner of the legend
        legend.background = element_rect(fill = "white", color = NA),  # optional: make background solid
        legend.box.background = element_rect(color = "black", size = 0.3),  # optional: border
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank()
        #legend.title = element_text(size = 7)
  )


F4b <- 
  ggplot(data = subset(Pellets, degradation!="casing"), 
         aes(x = volume, y = Respiration)) +
  geom_errorbar(data = subset(Pellets, degradation!="casing"),
                aes(x = volume, y = Respiration, 
                    ymin = Respiration - `Standard Error`, ymax = Respiration + `Standard Error`), size = 0.1) +
  geom_point(data = subset(Pellets, degradation != "casing"), 
             aes(x =volume, y = Respiration, color = Taxon,  shape = Taxon)) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  ylab(NULL)+
  xlab(expression("Volume (mm"^3*")"))+
  scale_x_log10(breaks = c(0.01, 0.1, 1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.tag = element_text(size = 9),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")


df <- subset(Pellets, degradation != "casing")
reg_model <- lm(Respiration ~ perimeter, data = df)
summary(reg_model)
coefs <- coef(reg_model)
r2 <- summary(reg_model)$r.squared
eq_label <- bquote(y == .(round(coefs[1], 3)) + .(round(coefs[2], 3)) * x ~~ ~~ R^2 == .(round(r2, 3)))

# Build equation text manually
eq_label <- paste0(
  "y = ", round(coefs[1], 3),
  " + ", round(coefs[2], 3), "x",
  '\n',
  "R² = ", round(r2, 2)
)
 

F4c <- 
  ggplot(data = subset(Pellets, degradation!="casing"), 
         aes(x =perimeter, y = Respiration)) +
  geom_errorbar(data = subset(Pellets, degradation!="casing"),
                aes(x = perimeter, y = Respiration, 
                    ymin = Respiration - `Standard Error`, ymax = Respiration + `Standard Error`), size = 0.1) +
  geom_smooth(data = subset(Pellets, degradation!="casing"), 
              mapping = aes(x =perimeter, y = Respiration), method = "lm", size = 0.2) +
  geom_point(data = subset(Pellets, degradation != "casing"), 
             aes(x =perimeter, y = Respiration, color = Taxon,  shape = Taxon)) +
  annotate("text", x = 0.7, y = 0.15, label = deparse(eq_label), parse = TRUE, size = 2.5, hjust = 0) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  ylab(NULL)+
  xlab("Perimeter (mm)")+
  scale_x_log10() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.tag = element_text(size = 9),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

F4 <- F4a + F4b + F4c + plot_layout(ncol = 3) # Combine them into one layout (1 row, 3 columns)

F4 <- (F4a | F4b | F4c) + plot_annotation(tag_levels = 'a')

print(F4) # Show the combined figure
#ggsave("F3vol.pdf", F3, width = 7.1, height = 3, dpi = 300) # Save the combined figure
ggsave("F4.jpg", F4, width = 7.5, height = 3, dpi = 900) # Save the combined figure


ggplot(data = subset(Pellets, degradation!="casing" ), # `Carbon-specific respiration`
       aes(x =porosity_komar, y = `Carbon content`/volume_SP)) +
  geom_point(aes( color = Taxon,  shape = Taxon)) +
  #annotate("text", x = 0.7, y = 0.15, label = deparse(eq_label), parse = TRUE, size = 2.5, hjust = 0) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  #ylim(0,11)+
  #facet_wrap(~shape)+
  #scale_y_log10() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

ggplot(data = subset(Pellets, degradation!="casing" ), # `Carbon-specific respiration`
       aes(x =porosity_komar, y = Respiration/volume_SP)) +
  geom_point(aes( color = Taxon,  shape = Taxon)) +
  #annotate("text", x = 0.7, y = 0.15, label = deparse(eq_label), parse = TRUE, size = 2.5, hjust = 0) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  ylim(0,6)+
  #facet_wrap(~shape)+
  #scale_y_log10() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

ggplot(data = subset(Pellets, degradation!="casing" ), # `Carbon-specific respiration`
       aes(x =porosity_komar, y = `Carbon-specific respiration`)) +
  geom_point(aes( color = Taxon,  shape = Taxon)) +
  #annotate("text", x = 0.7, y = 0.15, label = deparse(eq_label), parse = TRUE, size = 2.5, hjust = 0) +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  #ylim(0,6)+
  #facet_wrap(~shape)+
  #scale_y_log10() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#######################################################################################################
# Figure 6 Table S2 statistics
#######################################################################################################

run_lm_stats <- function(X, Y) {
  model <- lm(Y ~ X)
  sum_model <- summary(model)
  
  slope <- coef(model)[2]
  intercept <- coef(model)[1]
  p_value <- coef(sum_model)[2,4]
  r2 <- sum_model$r.squared
  df <- sum_model$df[2]
  
  results <- c(slope = slope, intercept = intercept, p_value = p_value, r2 = r2, df = df)
  return(results)  # Return as usual
}

X <- subset(Pellets$`Carbon content`, Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$degradation!="casing")
run_lm_stats(X,Y) 

X <- subset(Pellets$`Carbon content`, Pellets$short != "SP" & Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$short != "SP" & Pellets$degradation!="casing")
run_lm_stats(X,Y) 

X <- subset(Pellets$volume_SP, Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$degradation!="casing")
run_lm_stats(X,Y) 
X <- subset(Pellets$volume_SP, Pellets$short != "SP" & Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$short != "SP" & Pellets$degradation!="casing")
run_lm_stats(X,Y) 

X <- subset(Pellets$`2D area`, Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$degradation!="casing")
run_lm_stats(X,Y) 

X <- subset(Pellets$`2D area`, Pellets$short != "SP" & Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$short != "SP" & Pellets$degradation!="casing")
run_lm_stats(X,Y) 

X <- subset(Pellets$perimeter, Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$degradation!="casing")
run_lm_stats(X,Y)

X <- subset(Pellets$perimeter, Pellets$short != "SP" & Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$short != "SP" & Pellets$degradation!="casing")
run_lm_stats(X,Y)





X <- subset(Pellets$aspectratio, Pellets$short != "SP" & Pellets$degradation!="casing")
Y <- subset(Pellets$Respiration, Pellets$short != "SP" & Pellets$degradation!="casing")
run_lm_stats(X,Y)

X <- subset(Pellets$aspectratio, Pellets$degradation!="casing" & Pellets$shape == "cylinder")
Y <- subset(Pellets$`Carbon-specific respiration`, Pellets$degradation!="casing" & Pellets$shape == "cylinder")
run_lm_stats(X,Y) 

X <- subset(Pellets$porosity_komar, Pellets$degradation!="casing" )
Y <- subset(Pellets$Respiration/Pellets$volume_SP, Pellets$degradation!="casing" )
run_lm_stats(X,Y)

#######################################################################################################
# Figure 5: contribution to CARBON FLUX - carbon against volume, sinking and porosity against size
#######################################################################################################

X <- subset(Pellets$volume_SP, Pellets$degradation!="casing")
Y <- subset(Pellets$`Carbon content`, Pellets$degradation!="casing")

X <- subset(Pellets$volume_SP, Pellets$degradation!="casing" & Pellets$short == "SP" )
Y <- subset(Pellets$`Carbon content`, Pellets$degradation!="casing"  & Pellets$short == "SP")

X <- log(X)
Y <- log(Y)

run_lm_stats(X,Y) 

# Alldredge exponential model fit to pellets
fit <- lm(log(Y) ~ log(X), data = df)
summary(fit)

# 4. Extract coefficients
a <- exp(coef(fit)[1])  # Intercept transformed back
b <- coef(fit)[2]       # Slope

# Print the result
cat("Estimated equation: C =", round(a, 3), "* V^", round(b, 3), "\n")


# Create a sequence of X values for the fitted curve
fit_data <- data.frame(X = seq(0.001, 12, length.out = 100))

# equations 
fit_data$Y_fit <- 17.7 * fit_data$X^0.166
fit_data$Y_SP <- 19.9 * fit_data$X^0.35
fit_data$Y_Alldredge <- 1.05 * fit_data$X^0.51
fit_data$Y_Peru<- 8.9 * fit_data$X^0.51


F5a <- 
  ggplot() +
  geom_point(data = subset(Pellets), 
             aes(x = volume_SP, y = `Carbon content`, color = Taxon,  shape = Taxon)) +
  geom_line(data = fit_data, aes(x = X, y = Y_fit), color = "black", size = 0.5) +  # Fitted curve
  geom_line(data = subset(fit_data, X > 0.08), aes(x = X, y = Y_Alldredge), color = "black", size = 0.5, linetype = "dashed") +  # Fitted curve
  #geom_line(data = subset(fit_data, X > 0.08), aes(x = X, y = Y_SP), color = "black", size = 0.5) +  # Fitted curve
  ylab("Carbon (µg C)")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  xlab(expression("Volume (mm"^3*")"))+
    annotate("text", x = 0.001, y = 45, label = "C == 17.7 * V^0.166", parse = TRUE, size = 3, hjust = 0)+
  annotate("text", x = 0.001, y = 33, label = "R^2 == 0.28", parse = TRUE, size = 3, hjust = 0)+
  annotate("text", x = 0.03, y = 1, label = "Alldredge", parse = TRUE, size = 3, hjust = 0)+
  scale_y_log10() +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1)) +
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        plot.tag = element_text(size = 9),
        legend.position = "none")
#ggsave("Fig5d.jpg", width=3, height=3, dpi=300)


F5b <- 
ggplot() +
  geom_point(data = df_fit, aes(x = log(`2D area`), y = log(sinking_numeric),  color = Taxon, shape = Taxon)) +
  geom_smooth(data = df_fit, aes(x = log(`2D area`), y = sinking_size), size = 0.5, color = "black") +
  #geom_point(data = df_fit, aes(x = log(`2D area`), y = sinking_size_exp,  shape = Taxon), color = "purple") +
  #geom_point(data = df_fit, aes(x = log(`2D area`), y = sinking_size_tx,  shape = Taxon), color = "blue") +
  #geom_point(data = df_fit, aes(x = log(`2D area`), y = sinking_size_shape,  shape = Taxon), color = "darkred") +
  geom_point(shape = 21, size = 1) +
  annotate("text",x = -0.75, y = 4,  label = "v == 536.4 %*% `A`^0.443", parse = TRUE, size = 3, hjust = 0)+
  annotate("text",x = -0.75, y = 3.75,  label = "R^2 == 0.43", parse = TRUE, size = 3, hjust = 0)+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  ylab(expression("log(Sinking velocity (m day"^-1*"))"))+
  xlab(expression("log(2D Area (mm"^2*"))"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        plot.tag = element_text(size = 9),
        axis.text.x = element_text(size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5f.jpg", width=3, height=3, dpi=300)


F5d <- 
  ggplot() +
  # geom_point(data = subset(Pellets, type == "pellet" & degradation != "4_old" & degradation != "6_rounded"), 
  #            aes(x = (ESD/2), y = porosity_plough, color = Taxon,  shape = tx)) +
  # geom_smooth(data = subset(Pellets, type == "pellet" & degradation != "4_old" & degradation != "6_rounded"), 
  #              mapping = aes(x = (ESD/2), y = porosity_plough), method = "lm", size = 0.2) +
  #  geom_line(data = df_fit, aes(x = volume_SP, y = predicted_hill), color = "red", size = 0.5) +
  geom_line(data = df_fit, aes(x = volume_SP, y = predicted_exp), color = "black", size = 0.5) +
  #geom_line(data = df_fit, aes(x = volume_SP, y = predicted_log), color = "black", size = 0.5) +
  geom_point(data = subset(Pellets, degradation != "casing"), 
             aes(x = (volume_SP), y = porosity_komar, color = Taxon, shape = Taxon)) +
  xlab(expression("Volume (mm"^-3*"))"))+
  ylab(expression("Porosity (%)"))+
  ylim(0.2,1)+
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1)) +
  #scale_y_log10() +
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.tag = element_text(size = 9),
        panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 5), 
        legend.position = c(0.5, 0.5),  # near top-left inside
        legend.justification = c(0, 1),    # anchor top-left corner of the legend
        legend.background = element_rect(fill = "white", color = NA),  # optional: make background solid
        legend.box.background = element_rect(color = "black", size = 0.3),  # optional: border
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank())
        #legend.position = "none")
#ggsave("porosity_1r.jpg", width=4.5, height=2, dpi=300)


F5 <- F5a + F5d + F5b + plot_layout(ncol = 3, nrow = 1) + plot_annotation(tag_levels = 'a') # Combine them into one layout (1 row, 3 columns)

print(F5) # Show the combined figure
#ggsave("F3vol.pdf", F3, width = 7.1, height = 3, dpi = 300) # Save the combined figure
ggsave("F5.jpg", F5, width = 7.5, height = 3, dpi = 900) # Save the combined figure


X <- subset(log(Pellets$volume_SP), Pellets$degradation!="casing" )
Y <- subset(Pellets$porosity_komar, Pellets$degradation!="casing" )
run_lm_stats(X,Y) 

X <- subset((Pellets$`Carbon content`/ Pellets$volume_SP), Pellets$degradation!="casing" )
Y <- subset(Pellets$porosity_komar , Pellets$degradation!="casing" )
run_lm_stats(X,Y) 

#######################################################################################################
# Figure S5: model sinking
#######################################################################################################

FSa <- ggplot() +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = log(sinking_numeric),  shape = Taxon)) +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = sinking_size,  shape = Taxon), color = "darkred") +
  geom_point(shape = 21, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ylab(expression("log(Predicted sinking velocity (m day"^-1*"))"))+
  xlab(expression("log(Measured sinking velocity (m day"^-1*"))"))+
  ggtitle("only size")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5f.jpg", width=3, height=3, dpi=300)

FSb <- ggplot() +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = log(sinking_numeric),  shape = Taxon)) +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = sinking_size_exp,  shape = Taxon), color = "purple") +
  geom_point(shape = 21, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ggtitle("size + exponential porosity eq.")+
  ylab(expression("log(Predicted sinking velocity (m day"^-1*"))"))+
  xlab(expression("log(Measured sinking velocity (m day"^-1*"))"))+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5f.jpg", width=3, height=3, dpi=300)

FSc <- ggplot() +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = log(sinking_numeric),  shape = Taxon)) +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = sinking_size_tx,  shape = Taxon), color = "blue") +
  geom_point(shape = 21, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ggtitle("size + porosity taxon average")+
  ylab(expression("log(Predicted sinking velocity (m day"^-1*"))"))+
  xlab(expression("log(Measured sinking velocity (m day"^-1*"))"))+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5f.jpg", width=3, height=3, dpi=300)

FSd <- ggplot() +
  geom_point(data = df_fit, aes(x = log(sinking_numeric), y = log(sinking_numeric),  shape = Taxon)) +
   geom_point(data = df_fit, aes(x = log(sinking_numeric), y = sinking_size_shape,  shape = Taxon), color = "turquoise") +
  geom_point(shape = 21, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ggtitle("size + exponential porosity + shape")+
  ylab(expression("log(Predicted sinking velocity (m day"^-1*"))"))+
  xlab(expression("log(Measured sinking velocity (m day"^-1*"))"))+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Fig5f.jpg", width=3, height=3, dpi=300)

FS <- FSa + FSb + FSc + FSd + plot_layout(ncol = 2, nrow = 2) + plot_annotation(tag_levels = 'a') # Combine them into one layout (1 row, 3 columns)

print(FS) # Show the combined figure
#ggsave("F3vol.pdf", F3, width = 7.1, height = 3, dpi = 300) # Save the combined figure
ggsave("FS.jpg", FS, width = 7.5, height = 7.5, dpi = 900) # Save the combined figure



#######################################################################################################
# Figure 3: diversity correlations
#######################################################################################################

DivA <-
  ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = Chao16S, y = Chao18S, color = Taxon,  shape = Taxon)) +
  geom_smooth(data = subset(Pellets, !is.na(Shannon_16S)), 
              aes(x = Chao16S, y = Chao18S), method = "lm", size = 0.2) +
  ylab("Eukarya")+
  xlab("Prokarya")+
  ggtitle("Chao 1 Diversity")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
          title = element_text(size = 10, color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
          axis.text.y = element_text(size = 8, color = "black"), 
          legend.position = c(0.02, 0.98),  # near top-left inside
          legend.justification = c(0, 1),    # anchor top-left corner of the legend
          legend.background = element_rect(fill = "white", color = NA),  # optional: make background solid
          legend.box.background = element_rect(color = "black", size = 0.3),  # optional: border
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size = 7),
          legend.title = element_blank()
    )
#ggsave("Chao1_correlation.pdf", width=5.5, height=3, dpi=300)

model <- lm(Pellets$Chao18S ~ Pellets$Chao16S)
summary(model)

DivB <-
ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = Shannon_16S, y = Shannon_18S, color = Taxon,  shape = Taxon)) +
  ylab("Eukarya")+
  xlab("Prokarya")+
  ggtitle("Shannon Diversity")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Chao1_correlation.pdf", width=5.5, height=3, dpi=300)

model <- lm(alpha_diversity_18S$Shannon_18S ~ alpha_diversity_16S$Shannon_16S)
summary(model)

Div <- DivA + DivB + plot_layout(ncol = 2, nrow = 1) + plot_annotation(tag_levels = 'a') # Combine them into one layout

print(Div) # Show the combined figure
ggsave("Div.jpg", Div, width = 7.5, height = 3, dpi = 900) # Save the combined figure


#######################################################################################################
# Figure S3: diversity all
#######################################################################################################

Div18Sh <-
  ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = animal_ID, y = Shannon_18S, color = Taxon,  shape = Taxon)) +
  ylab("Shannon Diversity 18S")+
  xlab(NULL)+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")
#ggsave("Chao1_correlation.pdf", width=5.5, height=3, dpi=300)

Div18Ch <-
  ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = animal_ID, y = Chao18S, color = Taxon,  shape = Taxon)) +
  ylab("Chao1 Diversity 16S")+
  xlab(NULL)+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

Div16Sh <-
  ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = animal_ID, y = Shannon_16S, color = Taxon,  shape = Taxon)) +
  ylab("Shannon Diversity 16S")+
  xlab(NULL)+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

Div16Ch <-
  ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = animal_ID, y = Chao16S, color = Taxon,  shape = Taxon)) +
  ylab("Chao1 Diversity 16S")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        legend.position = "none")

DivS <- Div18Sh + Div18Ch + Div16Sh + Div16Ch + plot_layout(ncol = 1, nrow = 4) + plot_annotation(tag_levels = 'a') # Combine them into one layout

print(DivS) # Show the combined figure
ggsave("DivS.jpg", DivS, width = 7.5, height = 10, dpi = 900) # Save the combined figure

#######################################################################################################
# Figure S4: diversity against size
#######################################################################################################
DivP18C <-
  ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = perimeter, y = Chao18S, color = Taxon,  shape = Taxon)) +
  ylab(NULL)+
  xlab(NULL)+
  ggtitle("Chao 1 Diversity")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.position = "none")


DivP18S <-
ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = perimeter, y = Shannon_18S, color = Taxon,  shape = Taxon)) +
  ylab("Eukarya (18S)")+
  xlab(NULL)+
  ggtitle("Shannon Diversity")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        legend.position = "none")

DivP16C <-
ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = perimeter, y = Chao16S, color = Taxon,  shape = Taxon)) +
  ylab(NULL)+
  xlab("Perimeter (mm)")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"),
        legend.position = c(0.7, 0.95),  # near top-left inside
        legend.justification = c(0, 1),    # anchor top-left corner of the legend
        legend.background = element_rect(fill = "white", color = NA),  # optional: make background solid
        legend.box.background = element_rect(color = "black", size = 0.3),  # optional: border
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank()
  )


DivP16S <-
ggplot() +
  geom_point(data = subset(Pellets, !is.na(Shannon_16S)), 
             aes(x = perimeter, y = Shannon_16S, color = Taxon,  shape = Taxon)) +
  ylab("Prokarya (16S)")+
  xlab("Perimeter (mm)")+
  scale_color_manual(values = taxon_colors) +
  scale_shape_manual(values = taxon_shapes) +  
  geom_point(shape = 21, size = 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5),
        title = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"),
        legend.position = "none")

DivP <- DivP18S + DivP18C + DivP16S + DivP16C + plot_layout(ncol = 2, nrow = 2) + plot_annotation(tag_levels = 'a') # Combine them into one layout

print(DivP) # Show the combined figure
ggsave("DivP.jpg", DivP, width = 7.5, height = 7.5, dpi = 900) # Save the combined figure


#######################################################################################################
# export data
#######################################################################################################
#write.csv(Pellets, "BATS_pellets.csv", row.names = FALSE)
#######################################################################################################

          







