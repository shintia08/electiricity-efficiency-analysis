# Import Data 
library("readxl")

data <- read_excel("C:/SKRIPSI/DOKUMEN dan DATA/SFA + CLUSTERING/DATA/Research Data.xlsx",sheet = 'SFA', skip = 1)

# 1.	Descriptive Analysis 

num_data <- data[, 3:6]
summary(num_data)
sapply(num_data, sd, na.rm = TRUE)

# 2. Check Outliers by MD
X_MD_SFA <- data[, c(6, 3, 4, 5)]
X_MD_SFA<- as.data.frame(lapply(X_MD_SFA, as.numeric))

mvn

# Hitung Mahalanobis Distance
md_SFA <- mahalanobis(
  X_MD_SFA,
  colMeans(X_MD_SFA),
  cov(X_MD_SFA)
)

md_Clust <- mahalanobis(
  X_MD_Clust,
  colMeans(X_MD_Clust),
  cov(X_MD_Clust)
)

install.packages("MVN")
library(MVN)
mvn(X_MD_SFA, mvn_test = "hz")
mvn(X_MD_Clust, mvnTest = "hz")
p_SFA   <- ncol(X_MD_SFA)
p_Clust   <- ncol(X_MD_Clust)

cutoff_SFA <- qchisq(0.95, df = p_SFA)
cutoff_Clust <- qchisq(0.95, df = p_Clust)
outlier_flag_SFA <- md2_SFA > cutoff_SFA
outlier_flag_Clust <- md2_Clust > cutoff_Clust

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

plot(md2_Clust,
     pch = 16,
     col = ifelse(outlier_flag_Clust, "red", "black"),
     xlab = "Observation",
     ylab = expression("Mahalanobis Distance"^2),
     main = "Outlier Detection Using Mahalanobis Distance Cluster")

abline(h = cutoff_Clust, col = "blue", lwd = 2, lty = 2)

legend("topright",
       legend = c("Non-outlier", "Outlier", "Chi-square cutoff (5%)"),
       col = c("black", "red", "blue"),
       pch = c(16, 16, NA),
       lty = c(NA, NA, 2),
       lwd = c(NA, NA, 2),
       bty = "n")

# Jumlah outlier
sum(outlier_flag_SFA)
sum(outlier_flag_Clust)


# 3. Build SFA Model by MLE and Check Asumption 
library(frontier)

Output <- as.numeric(data[[6]])
Input1 <- as.numeric(data[[3]])
Input2 <- as.numeric(data[[4]])
Input3 <- as.numeric(data[[5]])
library(car)
# Check Asumption Before Build Model
# 1. Multicoliniearity
vif(lm(Input1 ~  Input2 + Input3))
vif(lm(Input2 ~  Input1 + Input3))
vif(lm(Input3 ~  Input1 + Input2))

# Check Asumption After Build Model 
# 2. Check Frontier Form : CD vs TL
df_log <- data.frame(
  Y_log = log(Output),
  X1_log = log(Input1),
  X2_log = log(Input2),
  X3_log = log(Input3))

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

# 3. Check Asumption Existences Inefficiency from Frontier that select vs Model Ordinary

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

# 4. Check Asumption Distribusi Ineff

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

# 5. Check Normalitas Residual 
# Ekstrak residual komponen stochastic error (v_i)

v_hat <- residuals(sfa_cd, type = "noise")
shapiro.test(v_hat)

# 6. Check Heterokedastisitas

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

# 7. Check Endogenitas 
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

# Ambil data yang ingin digunakan
clust_data <- data.frame(
  TE, X1, X2, X3, X4, )

# Opsional: normalisasi / scaling
clust_data_scaled <- scale(clust_data)
install.packages("GGally")

library(GGally)
library(cluster)

set.seed(123)

# K-medoids
pam_fit <- pam(scale(clust_data_scaled), k = 5, metric = "euclidean")
clust_data$cluster <- factor(pam_fit$clustering)

# Silhouette score
sil <- silhouette(pam_fit$clustering, dist(clust_data_scaled))
mean_sil <- mean(sil[, 3])
mean_sil

install.packages("clusterSim")
library(clusterSim)
index <- index.G1(clust_data_scaled, pam_fit$clustering, centrotypes="medoids")
print(index)

library(lmtest)

bptest
