# Read the CSV data
library(shiny)
library(rsconnect)

df <- read.csv("dataset.csv", header = TRUE)

# Convert Variant_ID to factor
df$Variant_ID <- as.factor(df$Variant_ID)

#install necessary libraries if you haven't
#install.packages(c("psych","combinat","MASS","effsize"))
library(psych) # for describeBy
library(combinat) # for permn
library(MASS) # for fitdistr
library(effsize) # for cohen.d
library(pwr) # pwr.t.test

########################################### Box Plots
boxplot(Time ~ Variant_ID, data=df) # Plotting boxplot of Time for both Variants
boxplot(Shampoos_Added ~ Variant_ID, data=df) # Plotting boxplot of Time for both Variants

describeTime <- describeBy(df$Time, df$Variant_ID) # Descriptive statistics
describeShampoosAdded <- describeBy(df$Shampoos_Added, df$Variant_ID) # Descriptive statistics

########################################### Central limit theorem

# Subset the data for Variant A and Variant B
variantA <- subset(df, Variant_ID == 1)
variantB <- subset(df, Variant_ID == 2)

# Extract the Time and Shampoos_Added columns for each variant
timeA <- variantA$Time
timeB <- variantB$Time
shampoosA <- variantA$Shampoos_Added
shampoosB <- variantB$Shampoos_Added

# Central limit theorem for Time column
sampleA_time <- replicate(100, mean(sample(timeA, 100, replace = TRUE)))
sampleB_time <- replicate(100, mean(sample(timeB, 100, replace = TRUE)))

# Central limit theorem for Shampoos_Added column
sampleA_shampoos <- replicate(100, mean(sample(shampoosA, 100, replace = TRUE)))
sampleB_shampoos <- replicate(100, mean(sample(shampoosB, 100, replace = TRUE)))

# Histograms for Time column
par(mfrow = c(2, 2))
hist(sampleA_time, main = "Variant A - Time", xlab = "Time")
hist(sampleB_time, main = "Variant B - Time", xlab = "Time")

# Density plots for Time column
plot(density(sampleA_time), main = "Variant A - Time", xlab = "Time")
plot(density(sampleB_time), main = "Variant B - Time", xlab = "Time")

# Shapiro-Wilk test for Time column
shapiro_test_time_A <- shapiro.test(sampleA_time)
shapiro_test_time_B <- shapiro.test(sampleB_time)

print(paste("Shapiro-Wilk Test - Time (Variant A): p-value =", shapiro_test_time_A$p.value))
print(paste("Shapiro-Wilk Test - Time (Variant B): p-value =", shapiro_test_time_B$p.value))

# QQ plots for Time column
par(mfrow = c(1, 2))
qqnorm(sampleA_time); qqline(sampleA_time)
qqnorm(sampleB_time); qqline(sampleB_time)

# Histograms for Shampoos_Added column
par(mfrow = c(2, 2))
hist(sampleA_shampoos, main = "Variant A - Shampoos Added", xlab = "Shampoos Added")
hist(sampleB_shampoos, main = "Variant B - Shampoos Added", xlab = "Shampoos Added")

# Density plots for Shampoos_Added column
plot(density(sampleA_shampoos), main = "Variant A - Shampoos Added", xlab = "Shampoos Added")
plot(density(sampleB_shampoos), main = "Variant B - Shampoos Added", xlab = "Shampoos Added")

# Shapiro-Wilk test for Shampoos_Added column
shapiro_test_shampoos_A <- shapiro.test(sampleA_shampoos)
shapiro_test_shampoos_B <- shapiro.test(sampleB_shampoos)

print(paste("Shapiro-Wilk Test - Shampoos Added (Variant A): p-value =", shapiro_test_shampoos_A$p.value))
print(paste("Shapiro-Wilk Test - Shampoos Added (Variant B): p-value =", shapiro_test_shampoos_B$p.value))

# QQ plots for Shampoos_Added column
par(mfrow = c(1, 2))
qqnorm(sampleA_shampoos); qqline(sampleA_shampoos)
qqnorm(sampleB_shampoos); qqline(sampleB_shampoos)

# Kolmogorov-Smirnov test for Time column
ks_test_time <- ks.test(sampleA_time, "pnorm", mean = mean(sampleB_time), sd = sd(sampleB_time))
print(paste("Kolmogorov-Smirnov Test - Time:", ks_test_time$p.value))

# Kolmogorov-Smirnov test for Shampoos_Added column
ks_test_shampoos <- ks.test(sampleA_shampoos, "pnorm", mean = mean(sampleB_shampoos), sd = sd(sampleB_shampoos))
print(paste("Kolmogorov-Smirnov Test - Shampoos Added:", ks_test_shampoos$p.value))


########################################### Z-test Time
X_var1 <- mean(df[df$Variant_ID == 1,]$Time)
X_var2 <- mean(df[df$Variant_ID == 2,]$Time)
S2_var1 <- var(df[df$Variant_ID == 1,]$Time)
S2_var2 <- var(df[df$Variant_ID == 2,]$Time)
n_var1 <- length(df[df$Variant_ID == 1,]$Time)
n_var2 <- length(df[df$Variant_ID == 2,]$Time)

X_diff = X_var2 - X_var1
(Z_diff_Time = X_diff / sqrt(S2_var1/n_var1 + S2_var2/n_var2))

########################################### Tails Time
(p_value = 1-pnorm(Z_diff_Time)) # 1-tailed, Variant 2 > Variant 1
(p_value = pnorm(Z_diff_Time)) # 1-tailed, Variant 2 < Variant 1
(p_value_Tails_Time = 1-pnorm(Z_diff_Time) + pnorm(-1 * Z_diff_Time)) # 2-tailed, Variant 2 <> Variant 1

########################################### Confidence Intervals
z_alpha_half <- qnorm(1-0.05/2)
SE_diff = sqrt(S2_var1/n_var1 + S2_var2/n_var2)
(cinterval_Time <- c(X_diff-SE_diff*z_alpha_half, X_diff+SE_diff*z_alpha_half))

########################################### Z-test Time
X_var1 <- mean(df[df$Variant_ID == 1,]$Shampoos_Added)
X_var2 <- mean(df[df$Variant_ID == 2,]$Shampoos_Added)
S2_var1 <- var(df[df$Variant_ID == 1,]$Shampoos_Added)
S2_var2 <- var(df[df$Variant_ID == 2,]$Shampoos_Added)
n_var1 <- length(df[df$Variant_ID == 1,]$Shampoos_Added)
n_var2 <- length(df[df$Variant_ID == 2,]$Shampoos_Added)

X_diff = X_var2 - X_var1
(Z_diff_Shampoos_Added = X_diff / sqrt(S2_var1/n_var1 + S2_var2/n_var2))

########################################### Tails Add to cart
(p_value = 1-pnorm(Z_diff_Shampoos_Added)) # 1-tailed, Variant 2 > Variant 1
(p_value = pnorm(Z_diff_Shampoos_Added)) # 1-tailed, Variant 2 < Variant 1
(p_value_Tails_Shampoos_Added = 1-pnorm(Z_diff_Shampoos_Added) + pnorm(-1 * Z_diff_Shampoos_Added)) # 2-tailed, Variant 2 <> Variant 1

########################################### Confidence Intervals
z_alpha_half <- qnorm(1-0.05/2)
SE_diff = sqrt(S2_var1/n_var1 + S2_var2/n_var2)
(cinterval_Shampoos_Added <- c(X_diff-SE_diff*z_alpha_half, X_diff+SE_diff*z_alpha_half))

###########################################  t-test returns the confidence interval automatically
t_test_Time <- t.test(Time ~ Variant_ID, data=df, alternative="two.sided", conf.level = 0.95, var.equal=TRUE)
t_test_Shampoos_Added <- t.test(Shampoos_Added ~ Variant_ID, data=df, alternative="two.sided", conf.level = 0.95, var.equal=TRUE)

###########################################  Effect size
effect_size_test_Time <- (cohen.d(Time ~ Variant_ID, hedges.correction=TRUE, data = df))
effect_size_test_Shampoos_Added <- (cohen.d(Shampoos_Added ~ Variant_ID, hedges.correction=TRUE, data = df))

########################################## Power analysis
nsample <- nrow(df)
nreplications <- 10000

# The power is calculated with the SUBTRACTION of the populations, not with the populations themselves
dist_null_time <- replicate(nreplications, mean(rnorm(nsample)) - mean(rnorm(nsample)))
dist_alternate_time <- replicate(nreplications, mean(rnorm(nsample, mean = 0.2)) - mean(rnorm(nsample, mean = 0.2)))

dist_null_shampoos <- replicate(nreplications, mean(rnorm(nsample)) - mean(rnorm(nsample)))
dist_alternate_shampoos <- replicate(nreplications, mean(rnorm(nsample, mean = 0.2)) - mean(rnorm(nsample, mean = 0.2)))

plot(density(dist_null_time), xlim = c(-1, 1), ylim = c(0, 5), col = "blue", main = "", xlab = "")
par(new = TRUE)
plot(density(dist_alternate_time), xlim = c(-1, 1), ylim = c(0, 5), col = "red", main = "", xlab = "")

plot(density(dist_null_shampoos), xlim = c(-1, 1), ylim = c(0, 5), col = "blue", main = "", xlab = "")
par(new = TRUE)
plot(density(dist_alternate_shampoos), xlim = c(-1, 1), ylim = c(0, 5), col = "red", main = "", xlab = "")

# Positive 1-tail for Time
mean_null_time <- fitdistr(dist_null_time, "normal")$estimate[1]
sd_null_time <- fitdistr(dist_null_time, "normal")$estimate[2]
critical_value_time <- qnorm(0.95, mean = mean_null_time, sd = sd_null_time)
abline(v = critical_value_time, col = "green")

mean_alternate_time <- fitdistr(dist_alternate_time, "normal")$estimate[1]
sd_alternate_time <- fitdistr(dist_alternate_time, "normal")$estimate[2]
cat("Power for Time:", pnorm(critical_value_time, mean = mean_alternate_time, sd = sd_alternate_time), "\n\n")

# Positive 1-tail for Shampoos_Added
mean_null_shampoos <- fitdistr(dist_null_shampoos, "normal")$estimate[1]
sd_null_shampoos <- fitdistr(dist_null_shampoos, "normal")$estimate[2]
critical_value_shampoos <- qnorm(0.95, mean = mean_null_shampoos, sd = sd_null_shampoos)
abline(v = critical_value_shampoos, col = "green")

mean_alternate_shampoos <- fitdistr(dist_alternate_shampoos, "normal")$estimate[1]
sd_alternate_shampoos <- fitdistr(dist_alternate_shampoos, "normal")$estimate[2]
cat("Power for Shampoos_Added:", pnorm(critical_value_shampoos, mean = mean_alternate_shampoos, sd = sd_alternate_shampoos), "\n\n")

########################################## Sample Size
delta_time <- 0.2  # effect size for Time variable
delta_shampoos <- 0.3  # effect size for Shampoos_Added variable
sigma_time <- 1  # standard deviation for Time variable
sigma_shampoos <- 1  # standard deviation for Shampoos_Added variable
alpha <- 0.05
beta <- 0.2

# Calculate sample size for Time variable
n_time <- ceiling(((qnorm(1-alpha) + qnorm(1-beta))^2 * (sigma_time^2 + sigma_time^2))/(delta_time^2))

# Calculate sample size for Shampoos_Added variable
n_shampoos <- ceiling(((qnorm(1-alpha) + qnorm(1-beta))^2 * (sigma_shampoos^2 + sigma_shampoos^2))/(delta_shampoos^2))

cat("Sample size for Time variable (per group): ", n_time, "\n")
cat("Sample size for Shampoos_Added variable (per group): ", n_shampoos, "\n")

########################################## Power
# Variable 1
delta_time <- 1  # effect size for Time variable
sigma_time <- 5  # standard deviation for Time variable
es_time <- delta_time / sigma_time

# Variable 2
delta_shampoos <- 0.5  # effect size for Shampoos_Added variable
sigma_shampoos <- 2  # standard deviation for Shampoos_Added variable
es_shampoos <- delta_shampoos / sigma_shampoos

# Calculate power for Variable 1
power_time <- pwr.t.test(d = es_time, sig.level = 0.05, power = 0.8,
                         type = "two.sample", alternative = "greater")

# Calculate power for Variable 2
power_shampoos <- pwr.t.test(d = es_shampoos, sig.level = 0.05, power = 0.8,
                             type = "two.sample", alternative = "greater")

cat("Power for Time variable:", power_time$power, "\n")
cat("Power for Shampoos_Added variable:", power_shampoos$power, "\n")

########################################### Power Calculation
# Variables and effect sizes
grand_average_time <- mean(df$Time)
standard_deviation_time <- sd(df$Time)
grand_average_shampoos <- mean(df$Shampoos_Added)
standard_deviation_shampoos <- sd(df$Shampoos_Added)


effect_time <- function(level) {
  switch(level, 'a' = 0, 'b' = 1)
}

effect_shampoos <- function(level) {
  switch(level, 'a' = 0, 'b' = 0.5)
}

# The model is: time ~ effect(factor) + N(grand_average, standard_deviation)
distr_group_a_time <- replicate(nreplications, 
                                (mean(effect_time('a') + rnorm(nsample, mean = grand_average_time, sd = standard_deviation_time))) -
                                  (mean(effect_time('a') + rnorm(nsample, mean = grand_average_time, sd = standard_deviation_time))))
distr_group_b_time <- replicate(nreplications, 
                                (mean(effect_time('b') + rnorm(nsample, mean = grand_average_time, sd = standard_deviation_time))) -
                                  (mean(effect_time('a') + rnorm(nsample, mean = grand_average_time, sd = standard_deviation_time))))

distr_group_a_shampoos <- replicate(nreplications, 
                                    (mean(effect_shampoos('a') + rnorm(nsample, mean = grand_average_shampoos, sd = standard_deviation_shampoos))) -
                                      (mean(effect_shampoos('a') + rnorm(nsample, mean = grand_average_shampoos, sd = standard_deviation_shampoos))))
distr_group_b_shampoos <- replicate(nreplications, 
                                    (mean(effect_shampoos('b') + rnorm(nsample, mean = grand_average_shampoos, sd = standard_deviation_shampoos))) -
                                      (mean(effect_shampoos('a') + rnorm(nsample, mean = grand_average_shampoos, sd = standard_deviation_shampoos))))

# Plot the distributions
plot(density(distr_group_a_time), xlim = c(-1, 2), ylim = c(0, 1), col = "blue", main = "", xlab = "")
par(new=TRUE)
plot(density(distr_group_b_time), xlim = c(-1, 2), ylim = c(0, 1), col = "red", main = "", xlab = "")
par(new=TRUE)

plot(density(distr_group_a_shampoos), xlim = c(-1, 2), ylim = c(0, 1), col = "blue", main = "", xlab = "")
par(new=TRUE)
plot(density(distr_group_b_shampoos), xlim = c(-1, 2), ylim = c(0, 1), col = "red", main = "", xlab = "")
par(new=TRUE)

# Power analysis for Time variable
critical_value_time <- quantile(distr_group_a_time, 0.95)
abline(v=critical_value_time, col = "green")

beta_time <- length(distr_group_b_time[distr_group_b_time <= critical_value_time]) / length(distr_group_b_time)
cat("Beta for Time variable:", beta_time, "\n")
cat(paste("Time - 95% percentile (alpha): ", critical_value <- quantile(distr_group_a_time, 0.95), "\n", sep=''))

# Power analysis for Shampoos_Added variable
abline(v=critical_value, col = "green")

beta_shampoos <- length(distr_group_b_shampoos[distr_group_b_shampoos <= critical_value_shampoos]) / length(distr_group_b_shampoos)
cat("Beta for Shampoos_Added variable:", beta_shampoos, "\n")
cat(paste("Shampoos Added - 95% percentile (alpha): ", critical_value <- quantile(distr_group_a_shampoos, 0.95), "\n", sep=''))


########################################### New Power



########################################### Run App

ui <- fluidPage(
  titlePanel("Evavuation Dashboard"),
  sidebarLayout(
    sidebarPanel(
      h4("Evaluation of Interactive Systems"),
    ),
    mainPanel(
      h2("Experiment Description"),
      p("The aim of the experiment was to investigate whether the inclusion of these labels would improve user experience and task performance on the website.\nThe research question guiding this study was:"),
      strong("Does the provision of additional labels about shampoo products improve the user experience and task performance compared to a version without these labels?"),
      div(
        br(),
        strong("Power Analysis"),
        br(),
        p(paste("Shapiro-Wilk Test - Time (Variant A): p-value =", shapiro_test_time_A$p.value)),
        p(paste("Shapiro-Wilk Test - Time (Variant B): p-value =", shapiro_test_time_B$p.value)),
        br(),
        p(paste("Shapiro-Wilk Test - Shampoos Added (Variant A): p-value =", shapiro_test_shampoos_A$p.value)),
        p(paste("Shapiro-Wilk Test - Shampoos Added (Variant B): p-value =", shapiro_test_shampoos_B$p.value)),
        br(),
        p(paste("Kolmogorov-Smirnov Test - Time:", ks_test_time$p.value)),
        p(paste("Kolmogorov-Smirnov Test - Shampoos Added:", ks_test_shampoos$p.value)),
        br(),
        p(paste("Kolmogorov-Smirnov Test - Time:", ks_test_time$p.value)),
        p(paste("Kolmogorov-Smirnov Test - Shampoos Added:", ks_test_shampoos$p.value)),
        br(),
        p(paste("Z Test - Time:", Z_diff_Time)),
        p(paste("Tail Test - Time:", Z_diff_Time)),
        p(paste("Confidence Interval Test - Time: min ", min(cinterval_Time), " - max ", max(cinterval_Time))),
        br(),
        p(paste("Z Test - Shampoos Added:", Z_diff_Shampoos_Added)),
        p(paste("Tail Test - Shampoos Added:", p_value_Tails_Shampoos_Added)),
        p(paste("Confidence Interval Test - Shampoos Added: min ", min(cinterval_Shampoos_Added), " - max ", max(cinterval_Shampoos_Added))),
        br(),
        p(paste("T Test - Time: p-value = ", t_test_Time$p.value, ", method = ", t_test_Time$method)),
        p(paste("Effect Size Test - Time: J = ", effect_size_test_Time$J, ", sd = ", effect_size_test_Time$sd, ", method = ", effect_size_test_Time$method)),
        br(),
        p(paste("T Test - Shampoos Added: p-value = ", t_test_Shampoos_Added$p.value, ", method = ", t_test_Shampoos_Added$method)),
        p(paste("Effect Size Test - Shampoos Added: J = ", effect_size_test_Shampoos_Added$J, ", sd = ", effect_size_test_Shampoos_Added$sd, ", method = ", effect_size_test_Shampoos_Added$method)),
        style = "background-color: white;border-radius:10px;color: black;padding: 10px;margin-top: 20px;margin-bottom: 20px;box-shadow: 0px 0px 10px #ccc;",
      ),
    )
  )
)
