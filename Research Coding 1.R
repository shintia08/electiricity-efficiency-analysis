# Import Data 
library("readxl")
data <- read_excel("C:/SKRIPSI/DOKUMEN dan DATA/SFA + CLUSTERING/DATA/Research Data.xlsx",sheet = 'SFA', skip = 1)

# 1.	Descriptive Analysis 
num_data <- data[, 3:6] 
summary(num_data) # statistic summaaries
sapply(num_data, sd, na.rm = TRUE) # Standard deviation (SD) as a summary statistic measuring data dispersion for each numerical variable.

# 2. Check Outliers by MD
X_MD_SFA <- data[, c(3, 4, 5)]
X_MD_SFA<- as.data.frame(lapply(X_MD_SFA, as.numeric)) # Convert selected variables to numeric data frame for Mahalanobis distance analysis.

# Compute Mahalanobis distance to detect multivariate outliers for SFA datasets.
md_SFA <- mahalanobis(
  X_MD_SFA,
  colMeans(X_MD_SFA),
  cov(X_MD_SFA)
) # Calculate Mahalanobis distance for SFA variables using their mean vector and covariance matrix.

install.packages("MVN")
library(MVN)

# Perform multivariate normality test (Henze-Zirkler), set Chi-square cutoff, and identify outliers based on Mahalanobis distance threshold.
mvn(X_MD_SFA, mvn_test = "hz")
p_SFA   <- ncol(X_MD_SFA)
cutoff_SFA <- qchisq(0.95, df = p_SFA)
outlier_flag_SFA <- md2_SFA > cutoff_SFA

# Plot Mahalanobis distance for SFA data, highlight outliers in red, and add Chi-square cutoff reference line and legend.
plot(md2_SFA,
     pch = 16,
     col = ifelse(outlier_flag_SFA, "red", "black"),
     xlab = "Observation",
     ylab = expression("Mahalanobis Distance"^2),
     main = "Outlier Detection Using Mahalanobis Distance SFA")

abline(h = cutoff_SFA, col = "blue", lwd = 2, lty = 2)

legend("topright",
       legend = c("Non-outlier", "Outlier", "Chi-square cutoff (5%)"),
       col = c("black", "red", "blue"),
       pch = c(16, 16, NA),
       lty = c(NA, NA, 2),
       lwd = c(NA, NA, 2),
       bty = "n")

# Jumlah outlier
sum(outlier_flag_SFA)

# 3. Build SFA Model

# Separate input and output variables and convert them to numeric format for analysis.

Output <- as.numeric(data[[6]])
Input1 <- as.numeric(data[[3]])
Input2 <- as.numeric(data[[4]])
Input3 <- as.numeric(data[[5]])
library(car)

df_log <- data.frame(
  Y_log = log(Output),
  X1_log = log(Input1),
  X2_log = log(Input2),
  X3_log = log(Input3))

# Check Asumption multicoliniear

vif(lm(Input1 ~  Input2 + Input3))
vif(lm(Input2 ~  Input1 + Input3))
vif(lm(Input3 ~  Input1 + Input2))

# Check Asumption Model 
library(frontier)

# 1. Check Frontier Form : CD vs TL using Likelihood Ratio Tes
#  === MODEL SFA (Cobb-Douglas)  ===
sfa_cd <- sfa(
  Y_log ~  X1_log + X2_log  + X3_log, data=df_log)
summary(sfa_cd)
#  === MODEL SFA (Translog)  ===
sfa_tl <- sfa(Y_log ~ X1_log + X2_log + X3_log + I(0.5 * X1_log^2) +
    I(0.5 * X2_log^2) +I(0.5 * X3_log^2) +I(X1_log * X2_log) + 
    I(X1_log * X3_log) + I(X2_log * X3_log), data = df_log)
summary(sfa_tl)

# Likelihood Ratio Test: CD vs TL

LL_tl <- as.numeric(logLik(sfa_tl))# TL
LL_cd <- as.numeric(logLik(sfa_cd)) # CD
LR_test <- -2 * (LL_cd - LL_tl)
df <- attr(logLik(sfa_tl), "df") - attr(logLik(sfa_cd), "df")
chi_critical <- qchisq(0.95, df)

# Decision
cat("\n===== LR TEST RESULT =====\n")
cat("LogLik Cobb-Douglas :", LL_cd, "\n")
cat("LogLik Translog     :", LL_tl, "\n")
cat("LR Statistic        :", LR_test, "\n")
cat("Chi-Square (5%)     :", chi_critical, "\n")
cat("Degrees of freedom  :", df, "\n")

if (LR_test > chi_critical) {
  cat("\n>> Kesimpulan: Tolak H0 — Model Translog secara signifikan lebih baik.\n")
} else {
  cat("\n>> Kesimpulan: Gagal tolak H0 — Model Cobb-Douglas cukup, Translog tidak memberikan perbaikan signifikan.\n")
}

# 2. Check Asumption Existences Inefficiency from Frontier that select vs Model Ordinary

# Ordinary Least Squares (restricted model: tidak ada inefficiency) 
ols <- lm(Y_log ~ X1_log + X2_log + X3_log , data = df_log)
summary(ols)

# Likelihood
LL_cd <- logLik(sfa_cd)
LL_ols <- logLik(ols)

# Likelihood Ratio Test OLS VS SFA
LR_ineff <- -2 * (LL_ols - LL_cd)
df_ineff <- attr(logLik(sfa_cd), "df") - attr(logLik(ols), "df")
df_ineff
chi_critical_ineff <- qchisq(0.95, df_ineff)

# Decision
cat("\n===== LR TEST INEFFICIENCY =====\n")
cat("LogLik OLS:", LL_ols, "\n")
cat("LogLik SFA:", LL_cd, "\n")
cat("LR Statistic:", LR_ineff, "\n")
cat("Critical Value (5%):", chi_critical_ineff, "\n")

if (LR_ineff > chi_critical_ineff) {
  cat("\n>> Kesimpulan: Tolak H0 — inefficiency signifikan → model SFA valid.\n")
} else {
  cat("\n>> Kesimpulan: Gagal tolak H0 — tidak ada inefficiency → cukup pakai OLS.\n")
}

# 3. Check Asumption Distribusi Ineff using Likelihood Ratio Test

# Distribusi Half Normal
LL_HN <- logLik(sfa_cd)

# Distribusi Truncated Normal 
sfa_TN <- sfa(Y_log ~ X1_log + X2_log + X3_log , data=df_log, ineffDecrease = TRUE,
              truncNorm = TRUE)
summary(sfa_TN)
LL_TN <- logLik(sfa_TN)

# === Likelihood Ratio Test ===
LR_dist <- -2 * (LL_HN - LL_TN)
df_dist <- attr(logLik(sfa_TN), "df") - attr(logLik(sfa_cd), "df")
chi_critical_dist <- qchisq(0.95, df_dist)

# Decision
cat("\n===== LR TEST: Half-Normal vs Truncated Normal =====\n")
cat("LogLik Half Normal:", LL_HN, "\n")
cat("LogLik Truncated Normal:", LL_TN, "\n")
cat("LR Statistic:", LR_dist, "\n")
cat("Critical Value (5%):", chi_critical_dist , "\n")

if (LR_dist > chi_critical_dist) {
  cat("\n>> Kesimpulan: Tolak H0 → Truncated Normal lebih baik.\n")
} else {
  cat("\n>> Kesimpulan: Gagal tolak H0 → Half Normal cukup → pakai model Half Normal.\n")
}


# 4. Check Endogenitas 
library(car)

# Uji Endogenitas Error U_i inefficiency
lm_u <- lm(
  u_i ~ X1_log + X2_log + X3_log,
  data = df_log
)

# Wald test parsial per variabel
wald_X1_ui <- linearHypothesis(lm_u, "X1_log = 0", test = "Chisq")
wald_X1_ui
wald_X2_ui <- linearHypothesis(lm_u, "X2_log = 0", test = "Chisq")
wald_X2_ui
wald_X3_ui <- linearHypothesis(lm_u, "X3_log = 0", test = "Chisq")
wald_X2_ui


wald_u <- linearHypothesis(
  lm_u,
  c("X1_log = 0", "X2_log = 0", "X3_log = 0"),
  test = "Chisq"
)
wald_u

# Uji Endogenitas Error V_i Noise Acak
lm_v <- lm(
  v_i ~ X1_log + X2_log + X3_log,
  data = df_log
)

# Wald test parsial per variabel
wald_X1_vi <- linearHypothesis(lm_v, "X1_log = 0", test = "Chisq")
wald_X1_vi
wald_X2_vi <- linearHypothesis(lm_u, "X2_log = 0", test = "Chisq")
wald_X2_vi
wald_X3_vi <- linearHypothesis(lm_u, "X3_log = 0", test = "Chisq")
wald_X2_vi


wald_v <- linearHypothesis(
  lm_v,
  c("X1_log = 0", "X2_log = 0", "X3_log = 0"),
  test = "Chisq"
)
wald_v



# 5. Check Heterokedastisitas

# inefficiency (u)
library(lmtest)
TE <- efficiencies(sfa_cd,log=FALSE)
TE
u_i <- -log(TE)
bptest(u_i ~ X1_log + X2_log + X3_log, data = df_log)

# Noise (v)
eps <- residuals(sfa_cd)  # ε_i = v_i + u_i
v_i <- eps + u_i
bptest(v_i ~ X1_log + X2_log + X3_log, data = df_log)

# 5. Check Normalitas Residual 
# Ekstrak residual komponen stochastic error (v_i)

v_hat <- residuals(sfa_cd, type = "noise")
shapiro.test(v_hat)

# Separate input-output variables and convert them to numeric format; Bayesian estimation is used since MLE is not appropriate due to non-normality of the inefficiency term (vi).

install.packages("rstan")
install.packages("ggplot2")
install.packages("tidyr")
library(rstan)
library(ggplot2)
library(tidyr)

# Data Structuring for Stan Model
data_stan <- list(
  N = nrow(df_log),
  K =3 ,
  y = df_log$Y_log,
  X = as.matrix(df_log[, c( "X1_log" , "X2_log", "X3_log")])
)

# prior distribution selection for the parameter gamma
# The prior γ ~ Beta(2,2) is chosen due to its stable estimation performance.
stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N, K] X;
}

parameters {
  real alpha;
  vector[K] beta;

  real<lower=0> tau;            // total SD
  real<lower=0, upper=1> gamma;   // SAME as R (variance proportion)

  vector<lower=0>[N] U_raw;       // half-normal inefficiency
  real<lower=2> nu;               // df Student-t
}

transformed parameters {
  real<lower=0> sigma_u;
  real<lower=0> sigma_v;

  vector[N] u;
  vector[N] mu;

  sigma_u = tau * sqrt(gamma);
  sigma_v = tau * sqrt(1 - gamma);

  u  = sigma_u * U_raw;
  mu = alpha + X * beta - u;
}

model {
  // =====================
  // PRIORS (R-COMPATIBLE)
  // =====================
  alpha ~ student_t(3, 0, 2.5);
  beta  ~ student_t(3, 0, 2.5);

  tau ~ normal(0, 0.5);          // total scale
  gamma ~ beta(15, 2);              // weak prior on (0,1)

  U_raw ~ normal(0, 1);            // half-normal
  nu ~ gamma(2, 0.1);

  // =====================
  // LIKELIHOOD
  // =====================
  y ~ student_t(nu, mu, sigma_v);
}

generated quantities {
  vector[N] TE;
  vector[N] log_lik;

  for (i in 1:N){
    TE[i] = exp(-u[i]);
    log_lik[i] = student_t_lpdf(y[i] | nu, mu[i], sigma_v);
    }
}
"

# Bayesian SFA Model Estimation via MCMC
fit <- stan(
  model_code = stan_code,
  data = data_stan,
  chains = 4,
  iter = 4000,      # total iterasi
  warmup = 2000,     # warmup < iter
  seed = 123,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
)

print(
  fit,
  pars = c("alpha","beta","tau", "gamma"),
  probs = c(0.025, 0.975)
)

# Calculation of the Deviance Information Criterion (DIC)

log_lik <- rstan::extract(fit)$log_lik
deviance <- -2 * rowSums(log_lik)
D_bar <- mean(deviance)
log_lik_mean <- colMeans(log_lik)
D_theta_bar <- -2 * sum(log_lik_mean)
p_D <- D_bar - D_theta_bar
DIC <- D_bar + p_D
DIC

# Posterior Sample Extraction
post <- rstan::extract(fit)

# Posterior Parameter Extraction
beta0_post <- post$alpha      # intercept (alpha)
beta_post  <- post$beta       # slope: iter x K

tau_post <- post$tau
gamma_post <- post$gamma

# Slope Parameter Labeling
colnames(beta_post) <- paste0("beta", 1:ncol(beta_post))
# Coefficient Formatting for Visualization
coef_post <- cbind(beta0 = beta0_post, beta_post)

# Posterior analysis of regression coefficients (β parameters)

df_beta <- as.data.frame(coef_post)
colnames(df_beta) <- c("beta0","beta1","beta2","beta3")
df_beta$iter <- 1:nrow(df_beta)

# ubah ke long format pakai tidyr
df_beta_long <- pivot_longer(df_beta, cols=c("beta0","beta1","beta2","beta3"),
                             names_to="variable", values_to="value")
# Density plot posterior beta
ggplot(df_beta_long, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(
    title = "Posterior Distribution of Parameters",
    x = "Value",
    y = "Density"
  ) +
  theme(legend.position = "none")
# Traceplot posterior beta
ggplot(df_beta_long, aes(x=iter, y=value, color=variable)) +
  geom_line() +
  labs(title="Traceplot of Beta", x="Iteration", y="Value") +
  theme_minimal()

# Posterior analysis of τ (tau) and γ (gamma) parameters
df_parameterisasi <- data.frame(
  tau = tau_post,
  gamma = gamma_post,
  iter = 1:length(tau_post)
)

# Ubah ke long format
df_parameterisasi_long <- pivot_longer(df_parameterisasi,
                              cols = c("tau","gamma"),
                              names_to = "variable",
                              values_to = "value")
# Density plot
ggplot(df_parameterisasi_long, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  labs(title = "Posterior Distribution of tau and gamma",
       x = "Value", y = "Density") +
  theme_minimal()

# Traceplot
ggplot(df_parameterisasi_long, aes(x = iter, y = value, color = variable)) +
  geom_line() +
  labs(title = "Traceplot of tau and gamma",
       x = "Iteration", y = "Value") +
  theme_minimal()

U_post       <- post$U        # iter x N
sigma_u_post <- post$sigma_u  # iter
u_post <- sweep(U_post, 1, sigma_u_post, "*")
TE_post <- exp(-u_post)
TE_mean <- apply(TE_post, 2, mean)
TE_median <- apply(TE_post, 2, median)
TE_mean

# === Sensitivity analysis using alternative Beta priors for gamma: Beta(5,2) to Beta(20,2) ===


# ==== CLUSTERING =====

# Create dataset for clustering analysis using mean TE  and Input3 variable.
clust_data <- data.frame(
  TE = TE_mean,Input3)

# Statistic Summaries
summary(clust_data)
sapply(clust_data, sd, na.rm = TRUE)

# Compute Spearman rank correlation between TE and Input3 variables.

# Outlier detection using squared Mahalanobis distance with Chi-square threshold

md_Clust <- mahalanobis(
  clust_data,
  colMeans(clust_data),
  cov(clust_data)
)

md2_Clust <- md_Clust^2
p_Clust   <- ncol(clust_data)

cutoff_Clust <- qchisq(0.95, df = p_Clust)
outlier_flag_Clust <- md2_Clust > cutoff_Clust

par(mar = c(5, 4, 4, 10))

plot(md2_Clust,
     pch = 16,
     col = ifelse(outlier_flag_Clust, "red", "black"),
     xlab = "Observation",
     ylab = expression("Mahalanobis Distance"^2))

# garis cutoff hanya dalam plot
segments(x0 = 1,
         x1 = length(md2_Clust),
         y0 = cutoff_Clust,
         y1 = cutoff_Clust,
         col = "blue",
         lwd = 2,
         lty = 2)

# legend di luar plot
par(xpd = TRUE)
legend(x = length(md2_Clust) + 3,
       y = max(md2_Clust),
       legend = c("Non-outlier", "Outlier", "Chi-square cutoff (5%)"),
       col = c("black", "red", "blue"),
       pch = c(16, 16, NA),
       lty = c(NA, NA, 2),
       lwd = c(NA, NA, 2),
       bty = "n")

# Normalisasi / scaling using Z-Score
clust_data_scaled <- scale(clust_data)

# Clustering using K-Medoids
install.packages("clusterCrit")
library(cluster)

# Average silhouette width method for determining the optimal number of clusters
sil_width <- sapply(2:5, function(k){
  pam_fit <- pam((clust_data_scaled), k = k,metric = "manhattan")
  pam_fit$silinfo$avg.width
})
plot(2:5, sil_width, type="b", xlab="Number of clusters", ylab="Average Silhouette Width")


# Cluster analysis using PAM (K-medoids) with Manhattan distance (k = 3).
install.packages("GGally")
library(ggplot2)
library(GGally)
set.seed(123)

pam_fit <- pam(
  clust_data_scaled,
  k = 3,
  metric = "manhattan"
)

clust_data$cluster <- factor(pam_fit$clustering)
cluster <- pam_fit$clustering
length(cluster)     # harus 33
table(cluster)
pam_fit$medoids
medoids <- clust_data[pam_fit$id.med, ]
medoids

# Cluster validation
# Silhouette score
sil <- silhouette(pam_fit$clustering, dist(clust_data_scaled))
mean_sil <- mean(sil[, 3])
mean_sil

# Total Sum of Squares
SST <- sum(scale(clust_data_scaled, scale = FALSE)^2)
# Within-cluster Sum of Squares
SSW <- 0
for (k in unique(pam_fit$clustering)) {
  idx <- which(pam_fit$clustering == k)
  center <- colMeans(clust_data_scaled[idx, ])
  SSW <- SSW + sum(scale(clust_data_scaled[idx, ], scale = FALSE)^2)
}
# R-squared clustering
R2 <- 1 - SSW / SST
R2

# Cluster Analysis Results Table
clust_data$"Satuan PLN/Provinsi" <- data$"Satuan PLN/Provinsi"
tabel_cluster <- clust_data[, c("Satuan PLN/Provinsi", "TE", "Input3", "cluster")]
tabel_cluster

