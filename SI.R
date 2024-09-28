library(car)
df <- read.csv("delta_ri_and_vo2_with_pc.csv")

m1 <- lm(vo2_final ~ mass_initial + ri_final * humidity * temperature, data = df)
m1a <- lm(vo2_final ~ mass_initial + ri_final * temperature, data = df)
m1b <- lm(vo2_final ~ mass_initial + ri_final * humidity, data = df)

vif <- vif(m1)
vif <- c(vif, NA)
t1 <- cbind(Anova(m1, type = 2), vif)

vif <- vif(m1a)
vif <- c(vif, NA)
t1a <- cbind(Anova(m1a, type = 2), vif)

vif <- vif(m1b)
vif <- c(vif, NA)
t1b <- cbind(Anova(m1b, type = 2), vif)

m2 <- lm(delta_vo2 ~ mass_initial + delta_ri * humidity * temperature, data = df)
m2a <- lm(delta_vo2 ~ mass_initial + delta_ri * humidity, data = df)
m2b <- lm(delta_vo2 ~ mass_initial + delta_ri * temperature, data = df)

vif <- vif(m2)
vif <- c(vif, NA)
t2 <- cbind(Anova(m2, type = 2), vif)

vif <- vif(m2a)
vif <- c(vif, NA)
t2a <- cbind(Anova(m2a, type = 2), vif)

vif <- vif(m2b)
vif <- c(vif, NA)
t2b <- cbind(Anova(m2b, type = 2), vif)

S1 <- rbind(t1,t1a,t1b,t2,t2a,t2b)
write.csv(S1, "S1.csv")
