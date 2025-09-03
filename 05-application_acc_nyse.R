
# nyse data ---------------------------------------------

library(tidyverse)
library(ggrepel)
library(quantmod)
library(patchwork)

Sys.setlocale("LC_TIME", "C") # para las fechas en inglés

source("scripts/sym_GH_estimation_functions.R")

# data --------------------------------------------------------------------
fecha_inicial <- as.Date("2000-01-01")
fecha_final <- as.Date("2024-03-31")

# Genera la secuencia de fechas mensuales desde la fecha inicial
secuencia_fechas <- seq(from = fecha_inicial, to = fecha_final, by = "month")
# Da formato a la secuencia de fechas
secuencia_fechas <- format(secuencia_fechas, "%b/%y")


# *S&P500 ----
getSymbols("^GSPC", src = "yahoo", from = fecha_inicial, to = fecha_final, periodicity =  "monthly")
sp500 <- monthlyReturn(GSPC) 

# *tasa libre de riesgo RF----
# fuente de datos de la tasa de interés libre de riesgo
# https://fred.stlouisfed.org/categories/22

getSymbols("TB3MS", src = "FRED", from = fecha_inicial, to = fecha_final, periodicity =  "monthly")
#tb3ms <- to.monthly(TB3MS)
tb3ms <- (1+TB3MS[,1]/100)^(1/12)-1 
tb3ms <- data.frame(tb3ms) %>% 
  tibble::rownames_to_column("fecha") %>% 
  mutate(fecha = as_datetime(fecha))


# tb_3ms <- readxl::read_excel("data/TB3MS.xls", skip = 10) %>% 
#   filter(observation_date > as.Date("1999-12-31 00:00:00")) %>% 
#   filter(observation_date < as.Date("2024-04-01 00:00:00")) %>% 
#   transmute(fecha = observation_date, 
#             tb3ms_a = (TB3MS/100), # tasa anual
#             tb3ms_m = (1+TB3MS/100)^(1/12)-1) #tasa mensual

# *stocks ----
# para obtener información historica de los precios de las acciones:

# bank of america
getSymbols("BAC", src = "yahoo", from = fecha_inicial, to = fecha_final, periodicity =  "monthly")
bac <- monthlyReturn(BAC) %>% as.data.frame()

# boing
getSymbols("BA", src = "yahoo", from = fecha_inicial, to = fecha_final, periodicity =  "monthly")
ba <- monthlyReturn(BA) %>% as.data.frame()

# ford motor company
getSymbols("F", src = "yahoo", from = fecha_inicial, to = fecha_final, periodicity =  "monthly")
ford <- monthlyReturn(F) %>% as.data.frame()

# general electric
getSymbols("GE", src = "yahoo", from = fecha_inicial, to = fecha_final, periodicity =  "monthly")
ge <- monthlyReturn(GE) %>% as.data.frame()

# microsoft
getSymbols("MSFT", src = "yahoo", from = fecha_inicial, to = fecha_final, periodicity =  "monthly")
msft <- monthlyReturn(MSFT) %>% as.data.frame()

acc_nyse <- data.frame(sp500 = sp500,
                       bac = bac,
                       ba = ba,
                       ford = ford,
                       ge = ge,
                       msft = msft) %>% 
  tibble::rownames_to_column("fecha") %>% 
  mutate(fecha = as_datetime(fecha))
colnames(acc_nyse) <- c("fecha","sp500", "bac", "ba", "ford", "ge", "msft")

acc_nyse_data <- left_join(tb3ms,acc_nyse)

# factores
SMB_HML <- read.delim("data/F-F_Research_Data_Factors.csv", sep = ",", skip = 3) %>% 
  mutate(fecha = paste0(str_sub(X,1,4),"-",str_sub(X,5,6),"-01"),
         fecha = as_datetime(fecha),
         smb = SMB/100,
         hml = HML/100,
         mkt.rf = Mkt.RF/100,
         rf = RF/100) %>% 
  #select(-c(X, SMB, HML)) %>% 
  #select(fecha, everything()) %>% 
  filter(fecha >= fecha_inicial) %>% 
  filter(fecha <= fecha_final)
head(SMB_HML)

acc_nyse_data <- left_join(acc_nyse_data, SMB_HML)


rm(sp500, GSPC, BAC, bac, ba, BA, F, ford, GE, ge, msft, MSFT, tb3ms, TB3MS, acc_nyse, SMB_HML )

#log rendimientos
X. <- acc_nyse_data %>% transmute(SP500 = log(sp500+1)*100) 
Y. <- acc_nyse_data %>% transmute(BankOfAmerica = log(bac+1)*100,
                                  Boeing = log(ba+1)*100, 
                                  Ford = log(ford+1)*100,
                                  GElectric = log(ge+1)*100,
                                  Microsoft = log(msft+1)*100)
Z. <- cbind(X.,Y.)

# retornos 
X_ <- acc_nyse_data %>% transmute(SP500 = sp500) 
Y_ <- acc_nyse_data %>% transmute(BankOfAmerica = bac,
                                  Boeing = ba, 
                                  Ford = ford,
                                  GElectric = ge,
                                  Microsoft = msft)
Z_ <- cbind(X_,Y_)

# the excess returns
X <- acc_nyse_data %>% transmute(SP500 = (sp500 - TB3MS)) 
Y <- acc_nyse_data %>% transmute(BankOfAmerica = (bac - TB3MS),
                                 Boeing = (ba - TB3MS), 
                                 Ford   = (ford - TB3MS),
                                 GElectric= (ge - TB3MS),
                                 Microsoft = (msft - TB3MS))
Z <- cbind(X,Y)
rownames(Z) <- secuencia_fechas

# save(X,Y,Z, file = "data/NYSE.RData")
# load("data/NYSE.RData")

# Descriptives ------------------------------------------------------------

descriptives_N <- function(Z) {
  Z <- as.matrix(Z)
  media <- colMeans(Z)
  sd <- apply(Z, 2, sd)
  covmat <- cov(Z)
  kurt <- apply(Z, 2, moments::kurtosis)
  skew <- apply(Z, 2, moments::skewness)
  
  if (dim(Z)[2] != 1) {
    mardia <-
      data.frame(
        MVN::mvn(Z, mvn_test = "mardia")$multivariate_normality[1, 1:3],
        MVN::mvn(Z, mvn_test = "mardia")$multivariate_normality[2, 1:3]
      ) %>%
      mutate(Statistic = as.numeric(levels(Statistic))[Statistic],
             p.value = as.numeric(levels(p.value))[p.value],
             Statistic.1 = as.numeric(levels(Statistic.1))[Statistic.1],
             p.value.1 = as.numeric(levels(p.value.1))[p.value.1]) %>% 
      transmute(Mardia_skew_s = Statistic,
                Mardia_skew_pv = ifelse(p.value < 0.009,
                                        format(p.value,
                                               scientific = TRUE,
                                               digits = 1),
                                        format(p.value,
                                               digits = 3)),
                Mardia_kurt_s = Statistic.1,
                Mardia_kurt_pv = ifelse(p.value.1 < 0.009,
                                        format(p.value.1,
                                               scientific = TRUE,
                                               digits = 1),
                                        format(p.value.1,
                                               digits = 3))
      )
    univariate <- data.frame(
      "mean" = round(media*100,3),
      "sd" = sd*100,
      #round(covmat,3),
      "skewness" = skew,
      "kurtosis" = kurt
    )
    
    tab <- data.frame(univariate, mardia) %>% rownames_to_column(var = "Variables")
    
    list(univariate = univariate,
         multivariate = mardia,
         tab = tab)
    
  } else {
    univariate <- data.frame(
      "mean" = media,
      "sd" = sd,
      #"var" = covmat,
      "skewness" = skew,
      "kurtosis" = kurt
    )
    list(univariate = univariate)
  }
}

# descriptivo de los excesos de retornos
desc_acc_nyse <- descriptives_N(Z)


library(GGally)

# grafico de los excesos de retorno.
plot_data_acc_nyse <- GGally::ggpairs(Z, 
                                     axisLabels = "none",
                                     lower = list(continuous = wrap("points", 
                                                                    color = "darkblue",
                                                                    alpha = 0.5)),
                                     upper = list(continuous = wrap("cor", 
                                                                    size = 3,
                                                                    method = "pearson", 
                                                                    color = "darkgray", 
                                                                    stars = FALSE)),
                                     diag = list(continuous = wrap("densityDiag", 
                                                                   color = "cadetblue4",
                                                                   fill = "cadetblue3",
                                                                   alpha = 0.2)  # Densidades en la diagonal
                                     )) + # Correlaciones con asteriscos) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),  # Quitar líneas principales
        panel.grid.minor = element_blank())  # Quitar líneas secundarias

# boxplot de los excesos de retornos

stocks <- c("BankOfAmerica", "Boeing", "Ford", "GElectric", "Microsoft", "SP500")

plot_bp_acc_nyse <- data.frame(Z) %>% 
  #select(stocks) %>% 
  pivot_longer(cols = SP500:Microsoft ,names_to = "serie") %>% 
  mutate(serie = factor(serie, 
                        labels = stocks, 
                        levels = stocks)) %>% 
  ggplot() +
  geom_boxplot(aes(y = value, x = serie), fill = "lightblue", color = "darkblue") +
  theme_bw() + 
  labs(x =" ", y = " ")

# Test de autocorrelacion y heteroscedasticidad ----

pvalues <- matrix(nrow = 3, ncol = 5)
X <- X$SP500  
for (i in 1:5) {
fit <- lm(Y[,i] ~ X)
par(mfrow = c(2,2))
plot(fit)
wh <- lmtest::bptest(fit, ~ X + I(X^2)) # test de white H0: No hay heteroscedasticidad
bp <- lmtest::bptest(fit) # test de breush pagan H0: No hay heteroscedasticidad
dw <- lmtest::dwtest(fit) # test de durbin watson H0: No hay autocorrelacion
pvalues[,i] <- c(wh$p.value, bp$p.value, dw$statistic)
}

colnames(pvalues) <- colnames(Y)
round(pvalues, 4)
dev.off()
# segun el test de white hay heteroscedasticidad en todos, excepto en Microsoft
# segun el test de breusch pagan solo hay heteroscedasticidad en GElectric
# segun el test de dw hay autocorrelación de primer grado en GElectric.

# el test de w indica heteroscedasticidad y el de pg no porque la heteroscedasticidad
# detectada no es lineal

# Model -------------------------------------------------------------------

# modelo normal
fit_N_acc_nyse <- MRMME::mrmme( Y= Y, X = X, method = "EM")

# nig
fit_SGH_acc_nyse_nig <- fit_SGH_MRMME_profile(X = X, Y = Y, 
                                             etas = seq(0.01,5,length = 20) , 
                                             type ="S-NIG")
#hyp
fit_SGH_acc_nyse_hyp <- fit_SGH_MRMME_profile(X = X, Y = Y, 
                                             etas = seq(0.01,5,length = 20) , 
                                             type ="S-HYP")

plot_perfilada_nig <- fit_SGH_acc_nyse_nig$plot_profile
plot_perfilada_hyp <- fit_SGH_acc_nyse_hyp$plot_profile
plot_perfilada <- plot_perfilada_nig | plot_perfilada_hyp

plot_perfilada <- plot_perfilada + 
  plot_annotation(
    tag_levels = "a",
    tag_prefix = '(',
    tag_sep = '',
    tag_suffix = ')'
  ) &
  theme(plot.tag.position = "top",
        plot.tag = element_text(
          face = "bold",
          size = 8,
          vjust = 5
        ))

# model assessment --------------------------------------------------------

# modelo NIG
plot_assess_acc_nyse_nig <- model_assess(X,Y,
                                        fit_SGH_acc_nyse_nig$logLz_i,
                                        fit_SGH_acc_nyse_nig$theta,
                                        n_sample = 500, type = "S-NIG")

plot_assess_acc_nyse_hyp <- model_assess(X,Y,
                                         fit_SGH_acc_nyse_hyp$logLz_i,
                                         fit_SGH_acc_nyse_hyp$theta,
                                         n_sample = 500, type = "S-HYP")

plot_assess_acc_nyse_N <- model_assess_N(X,Y,theta = fit_N_acc_nyse$theta, n_sample = 500)

# estimacion y se ---------------------------------------------------------



nam.tex <- c("$a_1$", "$a_2$", "$a_3$", "$a_4$", "$a_5$",
             "$\\beta_{1}$", "$\\beta_{2}$","$\\beta_{3}$","$\\beta_{4}$","$\\beta_{5}$",
             "$\\phi$", "$\\mu$", "$\\sigma$")

se_N <- fit_N_acc_nyse$se.fit$se %>% round(3)
se_nig <- se_SGH(fit_SGH_acc_nyse_nig$theta, X, Y)$se %>% round(3)
se_hyp <- se_SGH(fit_SGH_acc_nyse_hyp$theta, X, Y)$se %>% round(3)

estim_N <- paste0(round(fit_N_acc_nyse$theta,3),"(",se_N,")")
estim_nig <- paste0(round(fit_SGH_acc_nyse_nig$theta,3)[-c(14:15)],"(",se_nig,")")
estim_hyp <- paste0(round(fit_SGH_acc_nyse_hyp$theta,3)[-c(14:15)],"(",se_hyp,")")

estim_nyse <- cbind(nam.tex, c(stocks[1:5], stocks[1:5], " ", " ", " "), estim_N, estim_nig, estim_hyp)


# Grafico de distancias de mahalanobis ------------------------------------

# *modelo NIG ----
X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2

theta <- fit_SGH_acc_nyse_nig$theta
a <- matrix(theta[1:q], ncol = 1, nrow = q)
B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
phi <- theta[q + p * q + 1]
mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
matrixcalc::is.positive.definite(Sigma_x)
lambda <- theta[(d + 1)]
eta <- theta[(d + 2)]

alpha <- matrix(c(rep(0, p), a), ncol = 1)
Lambda <- rbind(diag(1, p), B)
Sigma_e <- diag(phi, r, r)
mu_z <- c(alpha + Lambda %*% mu_x)
Sigma_z <- Lambda %*% tcrossprod(Sigma_x, Lambda) + Sigma_e
Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
aux_mah <- t(t(Z) - mu_z)
delta <- diag(aux_mah %*% solve(Sigma_z, tol = 1e-30) %*% t(aux_mah))

data_mah <- data.frame(delta = delta, 
                       fecha =  as.Date(paste("01", secuencia_fechas), format = "%d %b/%y")) 

#lim_mah <- quantile(data_mah$delta, probs = c(0.98)) %>% round(0)
lim_mah <- quantile(data_mah$delta, probs = c(0.99))

data_mah_lab <- data_mah %>%
  #mutate(delta2 = round(delta,0)) %>% 
  filter(delta >= lim_mah)

plot_dist_mah_nyse_nig <- ggplot(data_mah, aes(x = fecha, y = delta)) + 
  geom_point(color="gray48") + 
  theme_bw() + 
  labs(x = "time", y = expression(hat(delta))) +
  #geom_hline(yintercept = 40, linetype = "dashed", color = "red") + 
  geom_point(data = data_mah_lab, aes(x = fecha, y = delta), col = "red") + 
  geom_text_repel(data = data_mah_lab, 
            aes(label = format(fecha, "%b/%y")), 
            # vjust = 2, 
            # hjust = -0.1, 
            color = "darkblue", size = 3) # Añadir etiquetas



# *modelo HYP ----
X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2

theta <- fit_SGH_acc_nyse_hyp$theta
a <- matrix(theta[1:q], ncol = 1, nrow = q)
B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
phi <- theta[q + p * q + 1]
mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
matrixcalc::is.positive.definite(Sigma_x)
lambda <- theta[(d + 1)]
eta <- theta[(d + 2)]

alpha <- matrix(c(rep(0, p), a), ncol = 1)
Lambda <- rbind(diag(1, p), B)
Sigma_e <- diag(phi, r, r)
mu_z <- c(alpha + Lambda %*% mu_x)
Sigma_z <- Lambda %*% tcrossprod(Sigma_x, Lambda) + Sigma_e
Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
aux_mah <- t(t(Z) - mu_z)
delta <- diag(aux_mah %*% solve(Sigma_z, tol = 1e-30) %*% t(aux_mah))

data_mah <- data.frame(delta = delta, 
                       fecha =  as.Date(paste("01", secuencia_fechas), format = "%d %b/%y")) 

lim_mah <- quantile(data_mah$delta, probs = c(0.99))
#lim_mah <- quantile(data_mah$delta, probs = c(0.98)) %>% round(0)

data_mah_lab <- data_mah %>%
  #mutate(delta2 = round(delta,0)) %>% 
  filter(delta >= lim_mah)

plot_dist_mah_nyse_hyp <-ggplot(data_mah, aes(x = fecha, y = delta)) + 
  geom_point(color="gray48") + 
  theme_bw() + 
  labs(x = "time", y = expression(hat(delta))) +
  geom_point(data = data_mah_lab, aes(x = fecha, y = delta), col = "red") + 
  geom_text_repel(data = data_mah_lab, 
                  aes(label = format(fecha, "%b/%y")), 
                  # vjust = 2, 
                  # hjust = -0.1, 
                  color = "darkblue", size = 3) # Añadir etiquetas


# *modelo Normal  ----
X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2

# theta --> (a,B,phi,mu,Sigma)
theta <- fit_N_acc_nyse$theta
a <- matrix(theta[1:q], ncol = 1, nrow = q)
B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
phi <- theta[q + p * q + 1]
mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))

# eta and psi
Lambda <- rbind(diag(1, p, p), B)
alpha <- c(rep(0,p),a)
eta <- c(alpha + Lambda %*% mu_x)
psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
psi <- 0.5 * (psi + t(psi))
delta <- stats::mahalanobis(Z,eta,psi)

data_mah <- data.frame(delta = delta, 
                       fecha =  as.Date(paste("01", secuencia_fechas), format = "%d %b/%y")) 

lim_mah <- qchisq(0.99,6)

data_mah_lab <- data_mah %>%
  filter(delta >= lim_mah)

plot_dist_mah_nyse_nor <-ggplot(data_mah, aes(x = fecha, y = delta)) + 
  geom_point(color="gray48") + 
  theme_bw() + 
  labs(x = "time", y = expression(hat(delta))) +
  geom_point(data = data_mah_lab, aes(x = fecha, y = delta), col = "red") + 
  geom_text_repel(data = data_mah_lab, 
                  aes(label = format(fecha, "%b/%y")), 
                  # vjust = 2, 
                  # hjust = -0.1, 
                  color = "darkblue", size = 3) # Añadir etiquetas



# Grafico de distancias de mahalanobis ------------------------------------

# distancias estandarizadas
X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2


# nig----

theta <- fit_SGH_acc_nyse_nig$theta
a <- matrix(theta[1:q], ncol = 1, nrow = q)
B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
phi <- theta[q + p * q + 1]
mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
matrixcalc::is.positive.definite(Sigma_x)
lambda <- theta[(d + 1)]
eta <- theta[(d + 2)]

alpha <- matrix(c(rep(0, p), a), ncol = 1)
Lambda <- rbind(diag(1, p), B)
Sigma_e <- diag(phi, r, r)
mu_z <- c(alpha + Lambda %*% mu_x)
Sigma_z <- Lambda %*% tcrossprod(Sigma_x, Lambda) + Sigma_e
Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
aux_mah <- t(t(Z) - mu_z)
delta_nig <- diag(aux_mah %*% solve(Sigma_z, tol = 1e-30) %*% t(aux_mah))

# hyp----
theta <- fit_SGH_acc_nyse_hyp$theta
a <- matrix(theta[1:q], ncol = 1, nrow = q)
B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
phi <- theta[q + p * q + 1]
mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))
matrixcalc::is.positive.definite(Sigma_x)
lambda <- theta[(d + 1)]
eta <- theta[(d + 2)]

alpha <- matrix(c(rep(0, p), a), ncol = 1)
Lambda <- rbind(diag(1, p), B)
Sigma_e <- diag(phi, r, r)
mu_z <- c(alpha + Lambda %*% mu_x)
Sigma_z <- Lambda %*% tcrossprod(Sigma_x, Lambda) + Sigma_e
Sigma_z <- 0.5 * (Sigma_z + t(Sigma_z))
aux_mah <- t(t(Z) - mu_z)
delta_hyp <- diag(aux_mah %*% solve(Sigma_z, tol = 1e-30) %*% t(aux_mah))

# nor----

theta <- fit_N_acc_nyse$theta
a <- matrix(theta[1:q], ncol = 1, nrow = q)
B <- matrix(theta[(q + 1):(q + p * q)], nrow = q, ncol = p)
phi <- theta[q + p * q + 1]
mu_x <- theta[(q + p * q + 2):(q + p * q + 1 + p)]
Sigma_x <- ks::invvech(theta[(q + p * q + 2 + p):d])
Sigma_x <- 0.5 * (Sigma_x + t(Sigma_x))

# eta and psi
Lambda <- rbind(diag(1, p, p), B)
alpha <- c(rep(0,p),a)
eta <- c(alpha + Lambda %*% mu_x)
psi <- Lambda %*% tcrossprod(Sigma_x, Lambda) + diag(phi, r, r)
psi <- 0.5 * (psi + t(psi))
delta_nor <- stats::mahalanobis(Z,eta,psi)


data_mah <- data.frame(delta_nig =  scale(delta_nig),
                       delta_hyp =  scale(delta_hyp),
                       delta_nor =  scale(delta_nor),
                       fecha =  as.Date(paste("01", secuencia_fechas), format = "%d %b/%y")) 

lim_nig <- quantile(data_mah$delta_nig, probs = c(0.99))
lim_hyp <- quantile(data_mah$delta_hyp, probs = c(0.99))
lim_nor <- quantile(data_mah$delta_nor, probs = c(0.99))

data_nig_lbl <- data_mah %>% filter(delta_nig >= lim_nig)
data_hyp_lbl <- data_mah %>% filter(delta_hyp >= lim_hyp)
data_nor_lbl <- data_mah %>% filter(delta_nor >= lim_nor)

plot_dist_nig <- ggplot(data_mah, aes(x = fecha, y = delta_nig)) + 
  geom_point(color="gray48") + 
  theme_bw() + 
  labs(x = "", y = expression(hat(delta))) +
  geom_point(data = data_nig_lbl, aes(x = fecha, y = delta_nig), col = "red") + 
  geom_text_repel(data = data_nig_lbl, 
                  aes(label = format(fecha, "%b/%y")),
                  color = "darkblue", size = 3)  + # Añadir etiquetas
  ylim(c(-1,15))

plot_dist_hyp <- ggplot(data_mah, aes(x = fecha, y = delta_hyp)) + 
  geom_point(color="gray48") + 
  theme_bw() + 
  labs(x = "", y = expression(hat(delta))) +
  geom_point(data = data_hyp_lbl, aes(x = fecha, y = delta_hyp), col = "red") + 
  geom_text_repel(data = data_hyp_lbl, 
                  aes(label = format(fecha, "%b/%y")),
                  color = "darkblue", size = 3)  + # Añadir etiquetas
  ylim(c(-1,15))

plot_dist_nor <- ggplot(data_mah, aes(x = fecha, y = delta_nor)) + 
  geom_point(color="gray48") + 
  theme_bw() + 
  labs(x = "", y = expression(hat(delta))) +
  geom_point(data = data_nor_lbl, aes(x = fecha, y = delta_nor), col = "red") + 
  geom_text_repel(data = data_nor_lbl, 
                  aes(label = format(fecha, "%b/%y")),
                  color = "darkblue", size = 3)  + # Añadir etiquetas
  ylim(c(-1,15))


mahalanobis_dist_nyse <- plot_dist_nor | plot_dist_nig | plot_dist_hyp


# Grafico de las rectas para cada activo ----

intercepts_nig <- fit_SGH_acc_nyse_nig$a
intercepts_hyp <- fit_SGH_acc_nyse_hyp$a
intercepts_nor <- fit_N_acc_nyse$a

slopes_nig <- fit_SGH_acc_nyse_nig$B
slopes_hyp <- fit_SGH_acc_nyse_hyp$B
slopes_nor <- fit_N_acc_nyse$B

assets <- colnames(Y)

plots_rec <- list()
for( i in 1:5){
  
  data <- tibble(asset = Y[,i], "SP500" = X[,1])
  
  line_data <- tibble(model = factor(c("N-MRMME", "S-NIG MRMME", "S-HYP MRMME"),
                                     levels = c("N-MRMME", "S-NIG MRMME", "S-HYP MRMME")),
                      intercept = c(intercepts_nor[i], intercepts_nig[i], intercepts_hyp[i]),
                      slope = c(slopes_nor[i], slopes_nig[i], slopes_hyp[i]))
  
  plots_rec[[i]] <- ggplot(data) +
    geom_point(aes(x = SP500, y=asset), color = "gray48") +
    geom_abline(aes(intercept = intercept,
                  slope = slope, 
                  color = model), data = line_data) +
    theme_bw() + 
    labs(x = "X", y = "Y", title = assets[i], color = "Model") +
    scale_color_manual(values = c("S-NIG MRMME" = "red", 
                                  "S-HYP MRMME" = "blue",
                                  "N-MRMME" = "green")) +
    geom_hline(yintercept = 0, color = "gray48", linetype = "dashed") + 
    geom_vline(xintercept = 0, color = "gray48", linetype = "dashed") 
    
    
}

# grafico conjunto
plot_lines <- wrap_plots(plots_rec, ncol = 3, nrow = 2) + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# Test de hipotesis -------------------------------------------------------

# H0: a=0

# *Modelo NIG -----
# modelo no restringido
mod_nr <- fit_SGH_acc_nyse_nig

# modelo restringido

mod_r_nig <- fit_SGH_MRMME_profile_R(X = X, Y = Y, 
                                     etas = seq(0.01,5,length = 20), 
                                     type ="S-NIG")
mod_r <- mod_r_nig

X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2
 
theta_hat <- mod_nr$theta
theta_tilde <- mod_r$theta
 
S_tilde <- score_matrix_SGH(theta_tilde, X = X, Y = Y, total = TRUE) %>% matrix(ncol=1)
 
# A matrix and g vector
A <- diag(c(rep(1,q), rep(0,d-q)))
A <- A[colSums(A) != 0, ]
g <- matrix(rep(0,q), ncol =1)
 
 # Fisher information
 
 # score_mat <- score_matrix_SGH(theta_tilde, X, Y)
 # fim <- crossprod(score_mat) # empirical fisher information matrix 
 # fim <- 0.5 * (fim + t(fim))
 # colnames(fim) <- nam(p, q)
 # rownames(fim) <- nam(p, q)
 # D <-  Rfast::spdinv(fim) # inversa de la fim
 # colnames(D) <- nam(p, q)
 # rownames(D) <- nam(p, q)
 
FIM_inv_tilde <-  se_SGH(theta_tilde, X = X, Y = Y)$covmat 
FIM_inv_hat <-  se_SGH(theta_hat, X = X, Y = Y)$covmat 
 
# LR
LR <- 2 * (lLz(X, Y, theta_hat, type = "S-NIG")$logL - 
            lLz(X, Y, theta_tilde, type = "S-NIG")$logL)
            
# WD
g_aux <- (A %*% theta_hat[1:d]) - g
WD <-  t(g_aux) %*% solve(A %*% FIM_inv_hat %*% t(A)) %*% g_aux
 
# SC
SC <- t(S_tilde) %*% FIM_inv_tilde %*% S_tilde
 
# GR
GR <- t(S_tilde) %*% (theta_hat[1:d] - theta_tilde[1:d])
GR
 
statistics_nig <- round(c(
   "LR" = LR,
   "WD" = WD,
   "SC" = SC,
   "GR" = GR
 ), 5)
 
 
d_nr <- d
d_r <- d - q
df <- d_nr - d_r
 
p.value_nig <- stats::pchisq(statistics_nig, df, lower.tail = FALSE)
qchisq(0.95,5)

AIC_nr_nig <- mod_nr$aic
AIC_r_nig <- mod_r$aic
 
BIC_nr_nig <- d_nr * log(n) - 2 * mod_nr$logLz
BIC_r_nig <- d_r * log(n) - 2 * mod_r$logLz
 
logL_nr_nig <- mod_nr$logLz
logL_r_nig <- mod_r$logLz

# *Modelo HYP -----
# modelo no restringido
mod_nr <- fit_SGH_acc_nyse_hyp

# modelo restringido

mod_r_hyp <-  fit_SGH_MRMME_profile_R(X = X, Y = Y, 
                                      etas = seq(0.01,5,length = 20) , 
                                      type ="S-HYP")
mod_r <- mod_r_hyp

X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2

theta_hat <- mod_nr$theta
theta_tilde <- mod_r$theta

S_tilde <- score_matrix_SGH(theta_tilde, X = X, Y = Y, total = TRUE) %>% matrix(ncol=1)

# A matrix and g vector
A <- diag(c(rep(1,q), rep(0,d-q)))
A <- A[colSums(A) != 0, ]
g <- matrix(rep(0,q), ncol =1)

# Fisher information

# score_mat <- score_matrix_SGH(theta_tilde, X, Y)
# fim <- crossprod(score_mat) # empirical fisher information matrix 
# fim <- 0.5 * (fim + t(fim))
# colnames(fim) <- nam(p, q)
# rownames(fim) <- nam(p, q)
# D <-  Rfast::spdinv(fim) # inversa de la fim
# colnames(D) <- nam(p, q)
# rownames(D) <- nam(p, q)

FIM_inv_tilde <-  se_SGH(theta_tilde, X = X, Y = Y)$covmat 
FIM_inv_hat <-  se_SGH(theta_hat, X = X, Y = Y)$covmat 

# LR
LR <- 2 * (lLz(X, Y, theta_hat, type = "S-HYP")$logL - 
             lLz(X, Y, theta_tilde, type = "S-HYP")$logL)

# WD
g_aux <- (A %*% theta_hat[1:d]) - g
WD <-  t(g_aux) %*% solve(A %*% FIM_inv_hat %*% t(A)) %*% g_aux

# SC
SC <- t(S_tilde) %*% FIM_inv_tilde %*% S_tilde

# GR
GR <- t(S_tilde) %*% (theta_hat[1:d] - theta_tilde[1:d])
GR

statistics_hyp <- round(c(
  "LR" = LR,
  "WD" = WD,
  "SC" = SC,
  "GR" = GR
), 5)


d_nr <- d
d_r <- d - q
df <- d_nr - d_r

p.value_hyp <- stats::pchisq(statistics_hyp, df, lower.tail = FALSE)

AIC_nr_hyp <- mod_nr$aic
AIC_r_hyp <- mod_r$aic

BIC_nr_hyp <- d_nr * log(n) - 2 * mod_nr$logLz
BIC_r_hyp <- d_r * log(n) - 2 * mod_r$logLz

logL_nr_hyp <- mod_nr$logLz
logL_r_hyp <- mod_r$logLz


# *Modelo Normal -----

# modelo no restringido
mod_nr <- fit_N_acc_nyse

# modelo restringido
mod_r_nor <- fit_EM_N_R(X,Y)
mod_r <- mod_r_nor

X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2

theta_hat <- mod_nr$theta
theta_tilde <- mod_r$theta

S_tilde <- score_matrix_N(theta_tilde, X = X, Y = Y, total = TRUE) %>% matrix(ncol=1)

# A matrix and g vector
A <- diag(c(rep(1,q), rep(0,d-q)))
A <- A[colSums(A) != 0, ]
g <- matrix(rep(0,q), ncol =1)

# Fisher information

FIM_inv_tilde <-  se_N(theta_tilde, X = X, Y = Y)$covmat 
FIM_inv_hat <-  se_N(theta_hat, X = X, Y = Y)$covmat 

# LR
LR <- 2 * (logL(theta_hat, X, Y) - 
           logL(theta_tilde, X, Y))
# WD
g_aux <- (A %*% theta_hat[1:d]) - g
WD <-  t(g_aux) %*% solve(A %*% FIM_inv_hat %*% t(A)) %*% g_aux

# SC
SC <- t(S_tilde) %*% FIM_inv_tilde %*% S_tilde

# GR
GR <- t(S_tilde) %*% (theta_hat[1:d] - theta_tilde[1:d])

statistics_N <- round(c(
  "LR" = LR,
  "WD" = WD,
  "SC" = SC,
  "GR" = GR
), 5)


d_nr <- d
d_r <- d - q
df <- d_nr - d_r

p.value_N <- stats::pchisq(statistics_N, df, lower.tail = FALSE)

AIC_nr_N <- mod_nr$AIC
AIC_r_N <- mod_r$AIC

BIC_nr_N <- d_nr * log(n) - 2 * mod_nr$logL
BIC_r_N <- d_r * log(n) - 2 * mod_r$logL

logL_nr_N <- mod_nr$logL
logL_r_N <- mod_r$logL

# *cuadro de resultados ----

nyse.test <- data.frame(
           "Statistic" = statistics_N, 
           "p-value" = p.value_N, 
           "Statistic" = statistics_nig, 
           "p-value" = p.value_nig,
           "Statistic" = statistics_hyp, 
           "p-value" = p.value_hyp)
rownames(nyse.test) <- paste0(c("Likelihood \n ratio","Wald","Score","Gradient"),
                              " (", names(statistics_N), ")")
nyse.test <- rownames_to_column(nyse.test,"Test") 

# *comparacion modelos nr y r ----
# para los 3 modelos, normal, s-nig y s-hyp

comp.tab.nyse <- data.frame("Under $H_0$" = c(logL_r_N, AIC_r_N, BIC_r_N),
                           "Complete" = c(logL_nr_N, AIC_nr_N, BIC_nr_N),
                           "Under $H_0$" = c(logL_r_nig, AIC_r_nig, BIC_r_nig),
                           "Complete" = c(logL_nr_nig, AIC_nr_nig, BIC_nr_nig),
                           "Under $H_0$" = c(logL_r_hyp, AIC_r_hyp, BIC_r_hyp),
                           "Complete" = c(logL_nr_hyp, AIC_nr_hyp, BIC_nr_hyp))
rownames(comp.tab.nyse) <- c("$\\mathcal{L}(\\btheta)$","AIC", "BIC")





# Diagnostico de influencias ----------------------------------------------

# Modelo NIG ----
influence.nig <- influence_SGH(fit_SGH_acc_nyse_nig$theta, X, Y)

# *plot scheme 1 ----
p1 <-
  ggplot(data.frame(influence.nig$G1, ind = names(influence.nig$G1)), aes(label = ind)) +
  #scale_y_continuous(limits = c(0,1)) +
  geom_hline(aes(yintercept = influence.nig$cut_point_G1), col = "red") +
  geom_point(aes(x = ind, y = influence.nig$G1)) +
  labs(x = "", y = "F values in Scheme I") +
  geom_text(aes(
    x = ind,
    y = influence.nig$G1,
    label = ifelse(influence.nig$G1 >= influence.nig$cut_point_G1, ind, "")
  ),
  hjust = 0,
  vjust = -1) +
  theme_bw() +
  ylim(c(0,0.45))

# *plot scheme 2----

p2 <- ggplot(data.frame(influence.nig$G2, 
                        ind = as.Date(paste("01", secuencia_fechas), format = "%d %b/%y")),
         aes(label = ind)) +
  geom_hline(aes(yintercept = influence.nig$cut_point_G2), col = "red") +
  geom_point(aes(x = ind, y = influence.nig$G2), color="gray48") +
  labs(x = "", y = "F values in Scheme II") +
  geom_text(
    aes(
      x = ind,
      y = influence.nig$G2,
      label = ifelse(round(influence.nig$G2,2) > influence.nig$cut_point_G2, format(ind,"%b/%y"), "")
    ),
    hjust = -1,
    vjust = 1,
    size = 2
  ) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  ylim(c(0,0.45))

library(patchwork)
p2t <- p1 | p2
p2t_nig <- p2t  +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = '(',
    tag_sep = '',
    tag_suffix = ')'
  ) &
  theme(plot.tag.position = "top",
        plot.tag = element_text(
          face = "bold",
          size = 8,
          vjust = 5
        ))

# Modelo HYP ----
influence.hyp <- influence_SGH(fit_SGH_acc_nyse_hyp$theta, X, Y)

# *plot scheme 1 ----
p1 <-
  ggplot(data.frame(influence.hyp$G1, ind = names(influence.hyp$G1)), aes(label = ind)) +
  #scale_y_continuous(limits = c(0,1)) +
  geom_hline(aes(yintercept = influence.hyp$cut_point_G1), col = "red") +
  geom_point(aes(x = ind, y = influence.hyp$G1)) +
  labs(x = "", y = "F values in Scheme I") +
  geom_text(aes(
    x = ind,
    y = influence.hyp$G1,
    label = ifelse(influence.hyp$G1 >= influence.hyp$cut_point_G1, ind, "")
  ),
  hjust = 0,
  vjust = -1) +
  theme_bw() +
  ylim(c(0,0.45))


# *plot scheme 2----

p2 <- ggplot(data.frame(influence.hyp$G2, 
                        ind = as.Date(paste("01", secuencia_fechas), format = "%d %b/%y")),
             aes(label = ind)) +
  geom_hline(aes(yintercept = influence.hyp$cut_point_G2), col = "red") +
  geom_point(aes(x = ind, y = influence.hyp$G2),color="gray48") +
  labs(x = "", y = "F values in Scheme II") +
  geom_text(
    aes(
      x = ind,
      y = influence.hyp$G2,
      label = ifelse(round(influence.hyp$G2,2) > influence.hyp$cut_point_G2, format(ind,"%b/%y"), "")
    ),
    hjust = -1,
    vjust = 1,
    size = 2
  ) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  ylim(c(0,0.45))

p2t <- p1 | p2
p2t_hyp <- p2t  +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = '(',
    tag_sep = '',
    tag_suffix = ')'
  ) &
  theme(plot.tag.position = "top",
        plot.tag = element_text(
          face = "bold",
          size = 8,
          vjust = 5
        ))

# Modelo Normal ----

influence.nor <- influence_N(fit_N_acc_nyse$theta, X, Y)

# *plot scheme 1 ----
p1 <-
  ggplot(data.frame(influence.nor$G1, ind = names(influence.nor$G1)), aes(label = ind)) +
  #scale_y_continuous(limits = c(0,1)) +
  geom_hline(aes(yintercept = influence.nor$cut_point_G1), col = "red") +
  geom_point(aes(x = ind, y = influence.nor$G1)) +
  labs(x = "", y = "F values in Scheme I") +
  geom_text(aes(
    x = ind,
    y = influence.nor$G1,
    label = ifelse(influence.nor$G1 >= influence.nor$cut_point_G1, ind, "")
  ),
  hjust = 0,
  vjust = -1) +
  theme_bw() +
  ylim(c(0,0.45))

# *plot scheme 2----

p2 <- ggplot(data.frame(influence.nor$G2, 
                        ind = as.Date(paste("01", secuencia_fechas), format = "%d %b/%y")),
             aes(label = ind)) +
  geom_hline(aes(yintercept = influence.nor$cut_point_G2), col = "red") +
  geom_point(aes(x = ind, y = influence.nor$G2), color="gray48") +
  labs(x = "", y = "F values in Scheme II") +
  geom_text(
    aes(
      x = ind,
      y = influence.nor$G2,
      label = ifelse(round(influence.nor$G2,2) > influence.nor$cut_point_G2, format(ind,"%b/%y"), "")
    ),
    hjust = -1,
    vjust = 1,
    size = 2
  ) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  ylim(c(0,0.45))

p2t <- p1 | p2
p2t_nor <- p2t  +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = '(',
    tag_sep = '',
    tag_suffix = ')'
  ) &
  theme(plot.tag.position = "top",
        plot.tag = element_text(
          face = "bold",
          size = 8,
          vjust = 5
        ))

# grafico de influencia de los 3 modelos

plot_influencia <- p2t_nig/ p2t_hyp / p2t_nor +
plot_annotation(
    tag_levels = "a",
    tag_prefix = '(',
    tag_sep = '',
    tag_suffix = ')'
  ) &
  theme(plot.tag.position = "top",
        plot.tag = element_text(
          face = "bold",
          size = 8,
          vjust = 5
        ))


# Guardar ----

save(acc_nyse_data,
     X, Y, Z, desc_acc_nyse, plot_bp_acc_nyse, plot_data_acc_nyse,
     fit_N_acc_nyse, fit_SGH_acc_nyse_nig, fit_SGH_acc_nyse_hyp,
     plot_perfilada,
     plot_lines,
     estim_nyse,
     plot_assess_acc_nyse_nig, plot_assess_acc_nyse_hyp, plot_assess_acc_nyse_N,
     plot_dist_mah_nyse_nig, plot_dist_mah_nyse_hyp, plot_dist_mah_nyse_nor,
     mahalanobis_dist_nyse,
     mod_r_nig, mod_r_hyp, mod_r_nor,
     nyse.test, comp.tab.nyse,
     influence.nig, influence.hyp, influence.nor,
     p2t_nig, p2t_hyp, p2t_nor, plot_influencia,
     file = "scripts/05-application_acc_nyse.RData")

load("scripts/05-application_acc_nyse.RData")

# AJUSTE SIN PUNTOS INFLUYENTES ----

X <- acc_nyse_data %>% transmute(SP500 = (sp500 - TB3MS)) 
Y <- acc_nyse_data %>% transmute(BankOfAmerica = (bac - TB3MS),
                                 Boeing = (ba - TB3MS), 
                                 Ford   = (ford - TB3MS),
                                 GElectric= (ge - TB3MS),
                                 Microsoft = (msft - TB3MS))
Z <- cbind(X,Y)
rownames(Z) <- secuencia_fechas

# para quitar abril/2009 (fila 112)
X <- X %>% slice(-112)
Y <- Y %>% slice(-112)
Z <- cbind(X,Y)
rownames(Z) <- secuencia_fechas[-112]


# Model -------------------------------------------------------------------

# modelo normal
fit_N_acc_nyse_influencia <- MRMME::mrmme( Y= Y, X = X, method = "EM")

# nig
fit_SGH_acc_nyse_nig_influencia <- fit_SGH_MRMME_profile(X = X, Y = Y, 
                                              etas = seq(0.01,5,length = 20) , 
                                              type ="S-NIG")
#hyp
fit_SGH_acc_nyse_hyp_influencia <- fit_SGH_MRMME_profile(X = X, Y = Y, 
                                              etas = seq(0.01,5,length = 20) , 
                                              type ="S-HYP")

# model assessment --------------------------------------------------------

# modelo NIG
plot_assess_acc_nyse_nig_influencia <- model_assess(X,Y,
                                                    fit_SGH_acc_nyse_nig_influencia$logLz_i,
                                                    fit_SGH_acc_nyse_nig_influencia$theta,
                                                    n_sample = 500, type = "S-NIG")

plot_assess_acc_nyse_hyp_influencia <- model_assess(X,Y,
                                                    fit_SGH_acc_nyse_hyp_influencia$logLz_i,
                                                    fit_SGH_acc_nyse_hyp_influencia$theta,
                                                    n_sample = 500, type = "S-HYP")

plot_assess_acc_nyse_N_influencia <- model_assess_N(X,Y,theta = fit_N_acc_nyse_influencia$theta, n_sample = 500)


# estimacion y se ---------------------------------------------------------

nam.tex <- c("$a_1$", "$a_2$", "$a_3$", "$a_4$", "$a_5$",
             "$\\beta_{1}$", "$\\beta_{2}$","$\\beta_{3}$","$\\beta_{4}$","$\\beta_{5}$",
             "$\\phi$", "$\\mu$", "$\\sigma$")

se_N <- fit_N_acc_nyse_influencia$se.fit$se %>% round(3)
se_nig <- se_SGH(fit_SGH_acc_nyse_nig_influencia$theta, X, Y)$se %>% round(3)
se_hyp <- se_SGH(fit_SGH_acc_nyse_hyp_influencia$theta, X, Y)$se %>% round(3)

estim_N <- paste0(round(fit_N_acc_nyse_influencia$theta,3),"(",se_N,")")
estim_nig <- paste0(round(fit_SGH_acc_nyse_nig_influencia$theta,3)[-c(14:15)],"(",se_nig,")")
estim_hyp <- paste0(round(fit_SGH_acc_nyse_hyp_influencia$theta,3)[-c(14:15)],"(",se_hyp,")")

estim_nyse_influencia <- cbind(nam.tex, c(stocks[1:5], stocks[1:5], " ", " ", " "), estim_N, estim_nig, estim_hyp)

# Test de hipotesis -------------------------------------------------------

# H0: a=0

# *Modelo NIG -----
# modelo no restringido
mod_nr <- fit_SGH_acc_nyse_nig_influencia

# modelo restringido

mod_r_nig_influencia <- fit_SGH_MRMME_profile_R(X = X, Y = Y, 
                                     etas = seq(0.01,5,length = 20) , 
                                     type ="S-NIG")
mod_r <- mod_r_nig

X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2

theta_hat <- mod_nr$theta
theta_tilde <- mod_r$theta

S_tilde <- score_matrix_SGH(theta_tilde, X = X, Y = Y, total = TRUE) %>% matrix(ncol=1)

# A matrix and g vector
A <- diag(c(rep(1,q), rep(0,d-q)))
A <- A[colSums(A) != 0, ]
g <- matrix(rep(0,q), ncol =1)

# Fisher information

# score_mat <- score_matrix_SGH(theta_tilde, X, Y)
# fim <- crossprod(score_mat) # empirical fisher information matrix 
# fim <- 0.5 * (fim + t(fim))
# colnames(fim) <- nam(p, q)
# rownames(fim) <- nam(p, q)
# D <-  Rfast::spdinv(fim) # inversa de la fim
# colnames(D) <- nam(p, q)
# rownames(D) <- nam(p, q)

FIM_inv_tilde <-  se_SGH(theta_tilde, X = X, Y = Y)$covmat 
FIM_inv_hat <-  se_SGH(theta_hat, X = X, Y = Y)$covmat 

# LR
LR <- 2 * (lLz(X, Y, theta_hat, type = "S-NIG")$logL - 
             lLz(X, Y, theta_tilde, type = "S-NIG")$logL)

# WD
g_aux <- (A %*% theta_hat[1:d]) - g
WD <-  t(g_aux) %*% solve(A %*% FIM_inv_hat %*% t(A)) %*% g_aux

# SC
SC <- t(S_tilde) %*% FIM_inv_tilde %*% S_tilde

# GR
GR <- t(S_tilde) %*% (theta_hat[1:d] - theta_tilde[1:d])
GR

statistics_nig <- round(c(
  "LR" = LR,
  "WD" = WD,
  "SC" = SC,
  "GR" = GR
), 5)


d_nr <- d
d_r <- d - q
df <- d_nr - d_r

p.value_nig <- stats::pchisq(statistics_nig, df, lower.tail = FALSE)

AIC_nr_nig <- mod_nr$aic
AIC_r_nig <- mod_r$aic

BIC_nr_nig <- d_nr * log(n) - 2 * mod_nr$logLz
BIC_r_nig <- d_r * log(n) - 2 * mod_r$logLz

logL_nr_nig <- mod_nr$logLz
logL_r_nig <- mod_r$logLz

# *Modelo HYP -----
# modelo no restringido
mod_nr <- fit_SGH_acc_nyse_hyp

# modelo restringido

mod_r_hyp <-  fit_SGH_MRMME_profile_R(X = X, Y = Y, 
                                      etas = seq(0.01,5,length = 20) , 
                                      type ="S-HYP")
mod_r <- mod_r_hyp

X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2

theta_hat <- mod_nr$theta
theta_tilde <- mod_r$theta

S_tilde <- score_matrix_SGH(theta_tilde, X = X, Y = Y, total = TRUE) %>% matrix(ncol=1)

# A matrix and g vector
A <- diag(c(rep(1,q), rep(0,d-q)))
A <- A[colSums(A) != 0, ]
g <- matrix(rep(0,q), ncol =1)

# Fisher information

# score_mat <- score_matrix_SGH(theta_tilde, X, Y)
# fim <- crossprod(score_mat) # empirical fisher information matrix 
# fim <- 0.5 * (fim + t(fim))
# colnames(fim) <- nam(p, q)
# rownames(fim) <- nam(p, q)
# D <-  Rfast::spdinv(fim) # inversa de la fim
# colnames(D) <- nam(p, q)
# rownames(D) <- nam(p, q)

FIM_inv_tilde <-  se_SGH(theta_tilde, X = X, Y = Y)$covmat 
FIM_inv_hat <-  se_SGH(theta_hat, X = X, Y = Y)$covmat 

# LR
LR <- 2 * (lLz(X, Y, theta_hat, type = "S-HYP")$logL - 
             lLz(X, Y, theta_tilde, type = "S-HYP")$logL)

fit_SGH_acc_nyse_hyp$logLz

# WD
g_aux <- (A %*% theta_hat[1:d]) - g
WD <-  t(g_aux) %*% solve(A %*% FIM_inv_hat %*% t(A)) %*% g_aux

# SC
SC <- t(S_tilde) %*% FIM_inv_tilde %*% S_tilde

# GR
GR <- t(S_tilde) %*% (theta_hat[1:d] - theta_tilde[1:d])
GR

statistics_hyp <- round(c(
  "LR" = LR,
  "WD" = WD,
  "SC" = SC,
  "GR" = GR
), 5)

d_nr <- d
d_r <- d - q
df <- d_nr - d_r

p.value_hyp <- stats::pchisq(statistics_hyp, df, lower.tail = FALSE)

AIC_nr_hyp <- mod_nr$aic
AIC_r_hyp <- mod_r$aic

BIC_nr_hyp <- d_nr * log(n) - 2 * mod_nr$logLz
BIC_r_hyp <- d_r * log(n) - 2 * mod_r$logLz

logL_nr_hyp <- mod_nr$logLz
logL_r_hyp <- mod_r$logLz

# *Modelo Normal -----

# modelo no restringido
mod_nr <- fit_N_acc_nyse

# modelo restringido
mod_r_nor <- fit_EM_N_R(X,Y)
mod_r <- mod_r_nor

X <- as.matrix(X)
Y <- as.matrix(Y)
Z <- cbind(X, Y)
n <- dim(Y)[1]
q <- dim(Y)[2]
p <- dim(X)[2]
r <- p + q
d <- (p + 1) * (2 + p + 2 * q) / 2

theta_hat <- mod_nr$theta
theta_tilde <- mod_r$theta

S_tilde <- score_matrix_N(theta_tilde, X = X, Y = Y, total = TRUE) %>% matrix(ncol=1)

# A matrix and g vector
A <- diag(c(rep(1,q), rep(0,d-q)))
A <- A[colSums(A) != 0, ]
g <- matrix(rep(0,q), ncol =1)

# Fisher information

FIM_inv_tilde <-  se_N(theta_tilde, X = X, Y = Y)$covmat 
FIM_inv_hat <-  se_N(theta_hat, X = X, Y = Y)$covmat 

# LR
LR <- 2 * (logL(theta_hat, X, Y) - 
             logL(theta_tilde, X, Y))
# WD
g_aux <- (A %*% theta_hat[1:d]) - g
WD <-  t(g_aux) %*% solve(A %*% FIM_inv_hat %*% t(A)) %*% g_aux

# SC
SC <- t(S_tilde) %*% FIM_inv_tilde %*% S_tilde

# GR
GR <- t(S_tilde) %*% (theta_hat[1:d] - theta_tilde[1:d])

statistics_N <- round(c(
  "LR" = LR,
  "WD" = WD,
  "SC" = SC,
  "GR" = GR
), 5)

d_nr <- d
d_r <- d - q
df <- d_nr - d_r

p.value_N <- stats::pchisq(statistics_N, df, lower.tail = FALSE)

AIC_nr_N <- mod_nr$AIC
AIC_r_N <- mod_r$AIC

BIC_nr_N <- d_nr * log(n) - 2 * mod_nr$logL
BIC_r_N <- d_r * log(n) - 2 * mod_r$logL

logL_nr_N <- mod_nr$logL
logL_r_N <- mod_r$logL

# *cuadro de resultados ----

nyse.test_influencia <- data.frame(
  "Statistic" = statistics_N, 
  "p-value" = p.value_N, 
  "Statistic" = statistics_nig, 
  "p-value" = p.value_nig,
  "Statistic" = statistics_hyp, 
  "p-value" = p.value_hyp)
rownames(nyse.test_influencia) <- paste0(c("Likelihood \n ratio","Wald","Score","Gradient"),
                              " (", names(statistics_N), ")")
nyse.test_influencia <- rownames_to_column(nyse.test_influencia,"Test") 

# *comparacion modelos nr y r ----
# para los 3 modelos, normal, s-nig y s-hyp

comp.tab.nyse_influencia <- data.frame("Under $H_0$" = c(logL_r_N, AIC_r_N, BIC_r_N),
                            "Complete" = c(logL_nr_N, AIC_nr_N, BIC_nr_N),
                            "Under $H_0$" = c(logL_r_nig, AIC_r_nig, BIC_r_nig),
                            "Complete" = c(logL_nr_nig, AIC_nr_nig, BIC_nr_nig),
                            "Under $H_0$" = c(logL_r_hyp, AIC_r_hyp, BIC_r_hyp),
                            "Complete" = c(logL_nr_hyp, AIC_nr_hyp, BIC_nr_hyp))

rownames(comp.tab.nyse_influencia) <- c("$\\mathcal{L}(\\btheta)$","AIC", "BIC")

fit_N_acc_nyse_influencia = fit_N_acc_nyse
fit_SGH_acc_nyse_nig_influencia = fit_SGH_acc_nyse_nig
fit_SGH_acc_nyse_hyp_influencia = fit_SGH_acc_nyse_hyp

mod_r_nig_influencia = mod_r_nig
mod_r_hyp_influencia = mod_r_hyp
mod_r_nor_influencia = mod_r_nor

# Guardar ----

save(estim_nyse_influencia,
     fit_N_acc_nyse_influencia,
     fit_SGH_acc_nyse_nig_influencia,
     fit_SGH_acc_nyse_hyp_influencia,
     plot_assess_acc_nyse_nig_influencia, plot_assess_acc_nyse_hyp_influencia,
     plot_assess_acc_nyse_N_influencia,
     nyse.test_influencia,
     comp.tab.nyse_influencia,
     mod_r_nig_influencia,
     mod_r_hyp_influencia,
     mod_r_nor_influencia,
     file = "scripts/05-application_acc_nyse_inf.RData")

#load("scripts/05-application_acc_nyse_inf.RData")

