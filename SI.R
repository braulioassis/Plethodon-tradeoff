### ----------------------------------------------------------
### Linear models and ANOVA with interaction effects in R
### Examining relationships between metabolic rate (VO2),
### body mass, and environmental variables (humidity, temperature)
### ----------------------------------------------------------

library(car)

# Import the dataset containing response variables (vo2, delta_vo2)
# and predictors such as mass, humidity, temperature, and respirometry index (ri)
df <- read.csv("delta_ri_and_vo2_with_pc.csv")

### ----------------------------------------------------------
### 1. VO2 models
### ----------------------------------------------------------

# Full model: includes all main effects and three-way interaction
m1 <- lm(vo2_final ~ mass_initial + ri_final * humidity * temperature, data = df)

# Reduced model A: excludes humidity (keeps interaction between ri_final and temperature)
m1a <- lm(vo2_final ~ mass_initial + ri_final * temperature, data = df)

# Reduced model B: excludes temperature (keeps interaction between ri_final and humidity)
m1b <- lm(vo2_final ~ mass_initial + ri_final * humidity, data = df)


### ----------------------------------------------------------
### 2. Variance Inflation Factor (VIF) and ANOVA
### ----------------------------------------------------------

# The VIF quantifies multicollinearity among predictors — values >5 or >10 can indicate concern
# ANOVA (Type II) tests for significance of each term in the model, accounting for others

# --- Model m1 ---
vif <- vif(m1)                   # Calculate VIFs for predictors
vif <- c(vif, NA)                # Add NA to align with ANOVA output row count
t1 <- cbind(Anova(m1, type = 2), vif)  # Combine Type II ANOVA table with VIF values

# --- Model m1a ---
vif <- vif(m1a)
vif <- c(vif, NA)
t1a <- cbind(Anova(m1a, type = 2), vif)

# --- Model m1b ---
vif <- vif(m1b)
vif <- c(vif, NA)
t1b <- cbind(Anova(m1b, type = 2), vif)


### ----------------------------------------------------------
### 3. ΔVO2 models (change in metabolic rate)
### ----------------------------------------------------------

# Full model: delta_vo2 (change in metabolic rate) as the response
m2 <- lm(delta_vo2 ~ mass_initial + delta_ri * humidity * temperature, data = df)

# Reduced model A: humidity interaction only
m2a <- lm(delta_vo2 ~ mass_initial + delta_ri * humidity, data = df)

# Reduced model B: temperature interaction only
m2b <- lm(delta_vo2 ~ mass_initial + delta_ri * temperature, data = df)


### ----------------------------------------------------------
### 4. VIF and ANOVA for ΔVO2 models
### ----------------------------------------------------------

# --- Model m2 ---
vif <- vif(m2)
vif <- c(vif, NA)
t2 <- cbind(Anova(m2, type = 2), vif)

# --- Model m2a ---
vif <- vif(m2a)
vif <- c(vif, NA)
t2a <- cbind(Anova(m2a, type = 2), vif)

# --- Model m2b ---
vif <- vif(m2b)
vif <- c(vif, NA)
t2b <- cbind(Anova(m2b, type = 2), vif)


### ----------------------------------------------------------
### 5. Combine and export all results
### ----------------------------------------------------------

# Combine ANOVA + VIF tables for all models into one large summary table
S1 <- rbind(t1, t1a, t1b, t2, t2a, t2b)

# Save combined results to a CSV file for downstream interpretation
write.csv(S1, "S1.csv")
