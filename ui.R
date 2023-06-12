# Read the CSV data
library(shiny)
library(rsconnect)

df <- read.csv("dataset.csv", header = TRUE, stringsAsFactors = TRUE)

# Convert Variant_ID to factor
df$Variant_ID <- as.factor(df$Variant_ID)

#install necessary libraries if you haven't
#install.packages(c("psych","combinat","MASS","effsize"))
library(psych) # for describeBy
library(combinat) # for permn
library(MASS) # for fitdistr
library(effsize) # for cohen.d
library(pwr) # pwr.t.test
library(broom)
library(car)
library(sjstats) # eta_sq

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

###########################################
# Independent samples Wilcoxon
#wilcox.test(Time ~ Variant_ID, data=df, alternative="greater")
#wilcox.test(Shampoos_Added ~ Variant_ID, data=df, alternative="greater")

# Paired
#wilcox.test(Time ~ Variant_ID, data=df, alternative="greater", paired = TRUE)
#wilcox.test(Shampoos_Added ~ Variant_ID, data=df, alternative="greater", paired = TRUE)

# Kruskal -Wallis
(fit <- kruskal.test(Time ~ Variant_ID, data=df))
print(paste("Kruskal -Wallis = ", fit))
(fit <- kruskal.test(Shampoos_Added ~ Variant_ID, data=df))
print(paste("Kruskal -Wallis = ", fit))

# Friedman's test
# (fit <- friedman.test(Time ~ Variant_ID|Test_Id, data=df))
# (fit <- friedman.test(Shampoos_Added ~ Variant_ID|Test_Id, data=df))

# Dunn's test (post-hoc)
library(dunn.test) # dunn,test

dunn.test(df$Time, df$Variant_ID)
dunn.test(df$Shampoos_Added, df$Variant_ID)

# Cliff's delta
cliff.delta(Time ~ Variant_ID, data=df)
cliff.delta(Shampoos_Added ~ Variant_ID, data=df)

# AOV Time
fit <- aov(Time ~ Variant_ID, data = df)
summary(fit)
print(paste("Coefficients = ", fit$coefficients))

# R squared
(tidy_fit <- tidy(fit))
(R_squared <- tidy_fit$sumsq[1] / (tidy_fit$sumsq[1] + tidy_fit$sumsq[2]))
print(paste("R Squared = ", R_squared))
# Tukey's test
# Maybe to Remove
(comp <- TukeyHSD(fit))
plot(comp)

# ANOVA assumptions 
# Show diagnostic plots
plot(fit, which=c(1,2))

# Check normality visually
mean_resid <- fitdistr(resid(fit), "normal")$estimate[1]
sd_resid <- fitdistr(resid(fit), "normal")$estimate[2]
xvalue <- seq(-6, 6, length = 200)
hist(resid(fit), freq = FALSE)
lines(xvalue, dnorm(xvalue, mean = mean_resid, sd = sd_resid))

# Diagnostic tests - Normality
shapiro.test(resid(fit))
ks.test(resid(fit), "pnorm")
qqnorm(resid(fit))
qqline(resid(fit))

# Diagnostic tests - Homogeneity of variances
leveneTest(Time ~ Variant_ID, data=df)

# Effect Size
# eta_squared(fit) Deprecated

# AOV Shampoos Added
fit <- aov(Shampoos_Added ~ Variant_ID, data = df)
summary(fit)
print(paste("Coefficients = ", fit$coefficients))

# R squared
(tidy_fit <- tidy(fit))
(R_squared <- tidy_fit$sumsq[1] / (tidy_fit$sumsq[1] + tidy_fit$sumsq[2]))
print(paste("R Squared = ", R_squared))
# Tukey's test
# Maybe to Remove
(comp <- TukeyHSD(fit))
print(comp)
plot(comp)

# ANOVA assumptions
# Show diagnostic plots
plot(fit, which=c(1,2))

# Check normality visually
mean_resid <- fitdistr(resid(fit), "normal")$estimate[1]
sd_resid <- fitdistr(resid(fit), "normal")$estimate[2]
xvalue <- seq(-6, 6, length = 200)
hist(resid(fit), freq = FALSE)
lines(xvalue, dnorm(xvalue, mean = mean_resid, sd = sd_resid))

# Diagnostic tests - Normality
shapiro.test(resid(fit))
ks.test(resid(fit), "pnorm")
qqnorm(resid(fit))
qqline(resid(fit))

# Diagnostic tests - Homogeneity of variances
leveneTest(Shampoos_Added ~ Variant_ID, data=df)

# Effect Size
# eta_squared(fit) Deprecated

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
