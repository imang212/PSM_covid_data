#Načtení potřebných knihoven
library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(stats)
library(corrplot)
library(gridExtra)
library(TTR)
library(seasonal)
library(vars)

PATH = "covid_timeseries.csv"

Sys.setlocale("LC_TIME", "C")

#Načtení dat
data <- read.csv(paste(PATH))

#Kontrola prázdných sloupců
prazdne_sloupce <- sapply(data, function(x) all(is.na(x)))
prazdne_sloupce
names(data)[prazdne_sloupce]

set.seed(123)
dates <- seq(as.Date("2020-01-01"), as.Date("2023-12-31"), by = "day")
n <- length(dates)

# Simulace COVID dat s trendem a sezónností
trend <- seq(0, 500, length.out = n) + rnorm(n, 0, 10)
seasonal <- 50 * sin(2 * pi * as.numeric(dates) / 365.25) + 30 * sin(2 * pi * as.numeric(dates) / 7)
new_cases <- pmax(0, trend + seasonal + rnorm(n, 0, 20))

# Další proměnné
new_deaths <- new_cases * 0.02 + rnorm(n, 0, 2); new_deaths <- pmax(0, new_deaths)

stringency_index <- 30 + 40 * sin(2 * pi * as.numeric(dates) / 365.25) + rnorm(n, 0, 5)
stringency_index <- pmin(100, pmax(0, stringency_index))

hosp_patients <- new_cases * 0.15 + rnorm(n, 0, 5); hosp_patients <- pmax(0, hosp_patients)

data <- data.frame(date = dates, location = "Czech Republic", new_cases = new_cases, new_deaths = new_deaths, stringency_index = stringency_index, hosp_patients = hosp_patients, total_cases = cumsum(new_cases), reproduction_rate = 1 + 0.3 * sin(2 * pi * as.numeric(dates) / 365.25) + rnorm(n, 0, 0.1))

# Převod na časovou řadu
ts_data <- ts(data$new_cases, start = c(2020, 1), frequency = 365.25)
head(data)
#####################################################################################
##Grafické zobrazení řady
# Grafické zobrazení hlavní řady
p1 <- ggplot(data, aes(x = date, y = new_cases)) +
  geom_line(color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Vývoj nových případů COVID-19",
       x = "Datum", y = "Počet nových případů",
       subtitle = "Denní data s trendem (červená linie)") +
  theme_minimal()

# Boxplot podle měsíců pro identifikaci sezónnosti
data$month <- factor(month(data$date), labels = month.abb)
p2 <- ggplot(data, aes(x = month, y = new_cases)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Sezónní rozdělení nových případů",
       x = "Měsíc", y = "Počet nových případů") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p1, p2, ncol = 1)

##Dekompozice časové řady
# STL dekompozice
stl_decomp <- stl(ts_data, s.window = "periodic", t.window = 365)
plot(stl_decomp, main = "STL dekompozice časové řady nových případů")
######################################################################################
##Identifikace trendu pomocí vyhlazení
# Klouzavé průměry různých řádů
ma_7 <- SMA(data$new_cases, n = 7)   
ma_30 <- SMA(data$new_cases, n = 30) 
ma_365 <- SMA(data$new_cases, n = 365) 

# Exponenciální vyrovnání - jednodušší přístup
exp_smooth_model <- HoltWinters(ts_data, gamma = FALSE)

# Vytvoříme jednoduché exponenciální vyhlazení pomocí ETS
ets_model <- ets(ts_data, model = "AAN", damped = FALSE)
exp_smooth_fitted <- fitted(ets_model)

# Zajistíme stejnou délku (ETS vrací stejnou délku jako originální data)
n_obs <- length(data$new_cases)

# Pokud je exp_smooth_fitted kratší, doplníme NA na začátek
if(length(exp_smooth_fitted) < n_obs) {
  exp_smooth_full <- c(rep(NA, n_obs - length(exp_smooth_fitted)), 
                       as.numeric(exp_smooth_fitted))
} else {
  exp_smooth_full <- as.numeric(exp_smooth_fitted[1:n_obs])
}

# Grafické srovnání - pouze metody, které máme
trend_data <- data.frame(
  date = data$date,
  original = data$new_cases,
  ma_7 = ma_7,
  ma_30 = ma_30
)

# Přidáme ma_365 pouze pokud má dostatek dat
if(sum(!is.na(ma_365)) > 100) {
  trend_data$ma_365 <- ma_365
}

# Přidáme exponenciální vyhlazení
trend_data$exp_smooth <- exp_smooth_full

trend_plot <- trend_data %>%
  pivot_longer(cols = -date, names_to = "method", values_to = "value") %>%
  filter(!is.na(value)) %>%  # Odfiltrujeme NA hodnoty
  ggplot(aes(x = date, y = value, color = method)) +
  geom_line(alpha = 0.7, size = 0.8) +
  scale_color_manual(values = c("original" = "gray50", "ma_7" = "blue", 
                                "ma_30" = "green", "ma_365" = "red",
                                "exp_smooth" = "purple"),
                     name = "Metoda",
                     labels = c("original" = "Originální data", 
                                "ma_7" = "MA(7)", "ma_30" = "MA(30)", 
                                "ma_365" = "MA(365)", "exp_smooth" = "Exp. vyhlazení")) +
  labs(title = "Srovnání metod vyhlazení trendu",
       x = "Datum", y = "Hodnota") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(trend_plot)

# Analýza trendové složky
trend_component <- stl_decomp$time.series[,"trend"]
cat("Souhrn trendové složky:\n")
print(summary(trend_component))

# Analýza sezónní složky
seasonal_component <- stl_decomp$time.series[,"seasonal"]
cat("\nSouhrn sezónní složky:\n")
print(summary(seasonal_component))
##########################################################################################
##Hledání optimálního funkčního modelu
# Příprava dat pro regresní model
data$t <- 1:nrow(data)
data$sin_365 <- sin(2 * pi * data$t / 365.25)
data$cos_365 <- cos(2 * pi * data$t / 365.25)
data$sin_7 <- sin(2 * pi * data$t / 7)
data$cos_7 <- cos(2 * pi * data$t / 7)

# Model 1: Lineární trend + roční sezónnost
model1 <- lm(new_cases ~ t + sin_365 + cos_365, data = data)
# Model 2: Polynomiální trend + roční sezónnost
model2 <- lm(new_cases ~ poly(t, 2) + sin_365 + cos_365, data = data)
# Model 3: Lineární trend + roční + týdenní sezónnost
model3 <- lm(new_cases ~ t + sin_365 + cos_365 + sin_7 + cos_7, data = data)
# Model 4: Polynomiální trend + roční + týdenní sezónnost
model4 <- lm(new_cases ~ poly(t, 2) + sin_365 + cos_365 + sin_7 + cos_7, data = data)
# Srovnání modelů
models <- list(model1, model2, model3, model4)
model_names <- c("Lin+Roční", "Poly+Roční", "Lin+Roční+Týdenní", "Poly+Roční+Týdenní")

# AIC a BIC srovnání
comparison <- data.frame(
  Model = model_names,
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  R_squared = sapply(models, function(x) summary(x)$r.squared)
)
print(comparison)

# Nejlepší model podle AIC
best_func_model <- models[[which.min(comparison$AIC)]]
print(summary(best_func_model))

# Grafické zobrazení fitování
data$fitted_func <- fitted(best_func_model)
ggplot(data, aes(x = date)) +
  geom_line(aes(y = new_cases), color = "black", alpha = 0.6, size = 0.5) +
  geom_line(aes(y = fitted_func), color = "red", size = 1) +
  labs(title = "Nejlepší funkční model",
       x = "Datum", y = "Počet případů",
       subtitle = paste("Model:", model_names[which.min(comparison$AIC)])) +
  theme_minimal()
###########################################################################################
##Hledání optimálního SARIMA modelu
# Kontrola stacionarity
adf_test <- adf.test(ts_data)
print(paste("ADF test p-hodnota:", round(adf_test$p.value, 4)))
# KPSS test stacionarity
kpss_test <- kpss.test(ts_data)
print(paste("KPSS test p-hodnota:", round(kpss_test$p.value, 4)))

# ACF a PACF grafy
par(mfrow = c(2,1))
acf(ts_data, main = "ACF původní řady", lag.max = 100)
pacf(ts_data, main = "PACF původní řady", lag.max = 100)

# Diferenciace pokud je potřeba
if(adf_test$p.value > 0.05) {
  ts_diff <- diff(ts_data)
  adf_diff <- adf.test(ts_diff)
  print(paste("ADF test diferencované řady:", round(adf_diff$p.value, 4)))
  
  par(mfrow = c(2,1))
  acf(ts_diff, main = "ACF diferencované řady", lag.max = 100)
  pacf(ts_diff, main = "PACF diferencované řady", lag.max = 100)
} else {
  ts_diff <- ts_data
}

cat("Testování vybraných SARIMA kombinací...\n")

# Vybrané rozumné kombinace na základě ACF/PACF a zkušeností
(fit1 <- Arima(ts_data, order = c(1,0,0), seasonal = list(order = c(0,1,0))))
(fit2 <- Arima(ts_data, order = c(0,0,1), seasonal = list(order = c(0,1,0))))
(fit3 <- Arima(ts_data, order = c(1,0,1), seasonal = list(order = c(0,1,0))))
(fit4 <- Arima(ts_data, order = c(0,0,0), seasonal = list(order = c(1,1,0))))
(fit5 <- Arima(ts_data, order = c(0,0,0), seasonal = list(order = c(0,1,1))))
(fit6 <- Arima(ts_data, order = c(0,0,0), seasonal = list(order = c(1,1,1))))
(fit7 <- Arima(ts_data, order = c(1,0,0), seasonal = list(order = c(1,1,1))))
BIC(fit1); BIC(fit2); BIC(fit3); BIC(fit4); BIC(fit5); BIC(fit6); BIC(fit7)

(auto_sarima_a <- auto.arima(ts_data)) #sezonní diference Zt = 0.004 + 0156*Zt-1 + 0.09*Et-1 + 0.11*Et-2 - 0.53Zt-1
BIC(auto_sarima_a) # automaticky model rozhodne optimalni neni - je prilis komplikovany
(auto_sarima_b <- auto.arima(ts_data, ic = "bic"))

# Použijeme auto.arima výsledek jako finální
final_sarima <- auto_sarima_b

# Diagnostika residuí
checkresiduals(final_sarima)
#######################################################################################
## Analýza závislostí na jiných řadách
# Příprava dalších časových řad
other_series <- c("new_deaths", "stringency_index", "hosp_patients", "reproduction_rate")

# Křížová korelace s různými zpožděními
ccf_results <- list()
max_lag <- 30

par(mfrow = c(2,2))
for(i in 1:length(other_series)) {
  y_var <- data[[other_series[i]]]
  y_var <- y_var[!is.na(y_var)]
  x_var <- data$new_cases[1:length(y_var)]
  
  ccf_result <- ccf(x_var, y_var, lag.max = max_lag, 
                    main = paste("CCF:", other_series[i]))
  ccf_results[[other_series[i]]] <- ccf_result
}

# Identifikace významných korelací
significant_lags <- list()
for(series in names(ccf_results)) {
  ccf_vals <- ccf_results[[series]]$acf
  lags <- ccf_results[[series]]$lag
  
  # 95% interval spolehlivosti (přibližně)
  n <- length(data$new_cases)
  ci <- 1.96/sqrt(n)
  
  significant_idx <- which(abs(ccf_vals) > ci)
  if(length(significant_idx) > 0) {
    significant_lags[[series]] <- data.frame(
      lag = lags[significant_idx],
      correlation = ccf_vals[significant_idx]
    )
    print(paste("Významné korelace pro", series, ":"))
    print(significant_lags[[series]])
  }
}

# Příprava zpožděných proměnných na základě CCF analýzy
data$deaths_lag <- c(rep(NA, 1), data$new_deaths[1:(nrow(data)-1)])
data$stringency_lag <- c(rep(NA, 5), data$stringency_index[1:(nrow(data)-5)])
data$hosp_lag <- c(rep(NA, 2), data$hosp_patients[1:(nrow(data)-2)])

# Model s externími regresory
external_vars <- cbind(
  deaths_lag = data$deaths_lag,
  stringency_lag = data$stringency_lag,
  hosp_lag = data$hosp_lag
)

# Odstranění NA hodnot
complete_idx <- complete.cases(cbind(data$new_cases, external_vars))
ts_complete <- ts(data$new_cases[complete_idx], start = c(2020, 1), frequency = 365.25)
external_complete <- external_vars[complete_idx, ]

# ARIMAX model
arimax_model <- auto.arima(ts_complete, xreg = external_complete)
print(summary(arimax_model))

# VAR model pro srovnání
var_data <- data[complete_idx, c("new_cases", "new_deaths", "stringency_index", "hosp_patients")]
var_data <- var_data[complete.cases(var_data), ]

# Optimální lag pro VAR
var_select <- VARselect(var_data, lag.max = 10, type = "both")
print(var_select$selection)

optimal_lag <- var_select$selection["AIC(n)"]
var_model <- VAR(var_data, p = optimal_lag, type = "both")
print(summary(var_model))

# Srovnání modelů
models_comparison <- data.frame(
  Model = c("Funkční", "SARIMA", "ARIMAX", "VAR"),
  AIC = c(AIC(best_func_model), AIC(final_sarima), AIC(arimax_model), AIC(var_model)),
  BIC = c(BIC(best_func_model), BIC(final_sarima), BIC(arimax_model), BIC(var_model))
)
print(models_comparison)

# Diagnostika residuí pro všechny modely
models_to_check <- list(
  "Funkční" = best_func_model,
  "SARIMA" = final_sarima,
  "ARIMAX" = arimax_model
)

par(mfrow = c(2,2))
for(model_name in names(models_to_check)) {
  model <- models_to_check[[model_name]]
  
  if(inherits(model, "lm")) {
    residuals <- residuals(model)
  } else {
    residuals <- residuals(model)
  }
  
  # ACF residuí
  acf(residuals, main = paste("ACF residuí -", model_name), lag.max = 50)
  # Ljung-Box test
  lb_test <- Box.test(residuals, lag = 20, type = "Ljung-Box")
  print(paste(model_name, "- Ljung-Box test p-hodnota:", round(lb_test$p.value, 4)))
  # Normalita residuí
  shapiro_test <- shapiro.test(sample(residuals, min(5000, length(residuals))))
  print(paste(model_name, "- Shapiro test p-hodnota:", round(shapiro_test$p.value, 4)))
}

# QQ plots pro normalitu
par(mfrow = c(2,2))
for(model_name in names(models_to_check)) {
  model <- models_to_check[[model_name]]
  residuals <- residuals(model)
  qqnorm(residuals, main = paste("Q-Q plot -", model_name))
  qqline(residuals)
}
##########################################################################################
##Predikce budoucích hodnot
# Predikce na 10 období dopředu
horizon <- 10
# Funkční model predikce
future_t <- (nrow(data)+1):(nrow(data)+horizon)
future_data <- data.frame(
  t = future_t,
  sin_365 = sin(2 * pi * future_t / 365.25),
  cos_365 = cos(2 * pi * future_t / 365.25),
  sin_7 = sin(2 * pi * future_t / 7),
  cos_7 = cos(2 * pi * future_t / 7)
)
# Jednoduchá predikce bez složitého rozlišování modelů
pred_func <- tryCatch({
  predict(best_func_model, newdata = future_data, interval = "prediction")
}, error = function(e) {
  cat("Chyba při predikci funkčního modelu:", e$message, "\n")
  # Náhradní predikce - jen fit bez intervalů
  fit_values <- predict(best_func_model, newdata = future_data)
  # Vytvoříme jednoduché intervaly na základě residuální std. chyby
  residual_se <- summary(best_func_model)$sigma
  cbind(fit = fit_values, 
        lwr = fit_values - 1.96 * residual_se, 
        upr = fit_values + 1.96 * residual_se)
})

cat("Funkční model - predikce dokončena\n")

# SARIMA predikce
pred_sarima <- forecast(final_sarima, h = horizon)

# ARIMAX predikce (potřebujeme budoucí hodnoty externích proměnných)
# Pro zjednodušení použijeme poslední dostupné hodnoty
last_external <- external_complete[nrow(external_complete), , drop = FALSE]
future_external <- matrix(rep(as.numeric(last_external), horizon), 
                          nrow = horizon, byrow = TRUE)

# Zajistíme správné názvy sloupců
colnames(future_external) <- colnames(external_complete)

cat("Názvy sloupců v trénovacích datech:", colnames(external_complete), "\n")
cat("Názvy sloupců v predikčních datech:", colnames(future_external), "\n")

pred_arimax <- tryCatch({
  forecast(arimax_model, xreg = future_external, h = horizon)
}, error = function(e) {
  cat("Chyba při ARIMAX predikci:", e$message, "\n")
  # Náhradní predikce bez externích regresorů
  forecast(final_sarima, h = horizon)
})

cat("ARIMAX predikce dokončena\n")

# Grafické zobrazení predikcí
future_dates <- seq(max(data$date) + 1, by = "day", length.out = horizon)

# Příprava dat pro graf
plot_data <- data.frame(
  date = c(tail(data$date, 50), future_dates),
  actual = c(tail(data$new_cases, 50), rep(NA, horizon)),
  func_pred = c(rep(NA, 50), pred_func[,1]),
  func_lower = c(rep(NA, 50), pred_func[,2]),
  func_upper = c(rep(NA, 50), pred_func[,3]),
  sarima_pred = c(rep(NA, 50), pred_sarima$mean),
  sarima_lower = c(rep(NA, 50), pred_sarima$lower[,2]),
  sarima_upper = c(rep(NA, 50), pred_sarima$upper[,2]),
  arimax_pred = c(rep(NA, 50), pred_arimax$mean),
  arimax_lower = c(rep(NA, 50), pred_arimax$lower[,2]),
  arimax_upper = c(rep(NA, 50), pred_arimax$upper[,2])
)

# Graf predikcí
ggplot(plot_data, aes(x = date)) +
  geom_line(aes(y = actual), color = "black", size = 1, alpha = 0.8) +
  geom_line(aes(y = func_pred), color = "red", size = 1) +
  geom_ribbon(aes(ymin = func_lower, ymax = func_upper), fill = "red", alpha = 0.2) +
  geom_line(aes(y = sarima_pred), color = "blue", size = 1) +
  geom_ribbon(aes(ymin = sarima_lower, ymax = sarima_upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = arimax_pred), color = "green", size = 1) +
  geom_ribbon(aes(ymin = arimax_lower, ymax = arimax_upper), fill = "green", alpha = 0.2) +
  geom_vline(xintercept = max(data$date), linetype = "dashed", alpha = 0.5) +
  labs(title = "Predikce nových případů - srovnání modelů",
       x = "Datum", y = "Počet případů",
       subtitle = "Červená: Funkční model, Modrá: SARIMA, Zelená: ARIMAX") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Výpis číselných predikcí
predictions_summary <- data.frame(
  Datum = future_dates,
  Funkční = round(pred_func[,1], 0),
  SARIMA = round(pred_sarima$mean, 0),
  ARIMAX = round(pred_arimax$mean, 0)
)
print("Predikce na příštích 10 dní:")
print(predictions_summary)
###############################################################################################
##Závěrečné srovnání a doporučení
# Finální srovnání všech modelů
final_comparison <- data.frame(
  Model = c("Funkční (trend+sezónost)", "SARIMA", "ARIMAX"),
  AIC = c(AIC(best_func_model), AIC(final_sarima), AIC(arimax_model)),
  BIC = c(BIC(best_func_model), BIC(final_sarima), BIC(arimax_model)),
  Parametrů = c(length(coef(best_func_model)), 
                length(coef(final_sarima)),
                length(coef(arimax_model))),
  RMSE = c(
    sqrt(mean(residuals(best_func_model)^2)),
    sqrt(mean(residuals(final_sarima)^2)),
    sqrt(mean(residuals(arimax_model)^2))
  )
)

print("Finální srovnání modelů:")
print(final_comparison)

# Doporučení nejlepšího modelu
best_model_idx <- which.min(final_comparison$AIC)
best_model_name <- final_comparison$Model[best_model_idx]

cat("\n=== DOPORUČENÍ ===\n")
cat("Nejlepší model podle AIC:", best_model_name, "\n")
cat("AIC:", round(final_comparison$AIC[best_model_idx], 2), "\n")
cat("RMSE:", round(final_comparison$RMSE[best_model_idx], 2), "\n")



######################################################################################




# Načtení potřebných knihoven (bez seasonal a vars)
library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(stats)
library(corrplot)
library(gridExtra)
library(TTR)

Sys.setlocale("LC_TIME", "C")

# Načtení dat (nebo simulace)
# ... (same as in your original script up to creation of `data` and ts_data) ...
# For brevity in this snippet, assume the data creation block is identical to yours
# (dates, simulated new_cases, new_deaths, stringency_index, hosp_patients, total_cases, reproduction_rate)
# and ts_data <- ts(data$new_cases, start = c(2020, 1), frequency = 365.25)
#Načtení dat
data <- read.csv("covid_timeseries.csv")

#Kontrola prázdných sloupců
prazdne_sloupce <- sapply(data, function(x) all(is.na(x)))
prazdne_sloupce
names(data)[prazdne_sloupce]

set.seed(123)
dates <- seq(as.Date("2020-01-01"), as.Date("2023-12-31"), by = "day")
n <- length(dates)

# Simulace COVID dat s trendem a sezónností
trend <- seq(0, 500, length.out = n) + rnorm(n, 0, 10)
seasonal <- 50 * sin(2 * pi * as.numeric(dates) / 365.25) + 30 * sin(2 * pi * as.numeric(dates) / 7)
new_cases <- pmax(0, trend + seasonal + rnorm(n, 0, 20))

# Další proměnné
new_deaths <- new_cases * 0.02 + rnorm(n, 0, 2); new_deaths <- pmax(0, new_deaths)

stringency_index <- 30 + 40 * sin(2 * pi * as.numeric(dates) / 365.25) + rnorm(n, 0, 5)
stringency_index <- pmin(100, pmax(0, stringency_index))

hosp_patients <- new_cases * 0.15 + rnorm(n, 0, 5); hosp_patients <- pmax(0, hosp_patients)

data <- data.frame(date = dates, location = "Czech Republic", new_cases = new_cases, new_deaths = new_deaths, stringency_index = stringency_index, hosp_patients = hosp_patients, total_cases = cumsum(new_cases), reproduction_rate = 1 + 0.3 * sin(2 * pi * as.numeric(dates) / 365.25) + rnorm(n, 0, 0.1))

# Převod na časovou řadu
ts_data <- ts(data$new_cases, start = c(2020, 1), frequency = 365.25)
head(data)
######################################################################
##Grafické zobrazení řady
# Grafické zobrazení hlavní řady
p1 <- ggplot(data, aes(x = date, y = new_cases)) +
  geom_line(color = "steelblue", alpha = 0.7) +
  geom_smooth(method = "loess", color = "red", se = FALSE) +
  labs(title = "Vývoj nových případů COVID-19",
       x = "Datum", y = "Počet nových případů",
       subtitle = "Denní data s trendem (červená linie)") +
  theme_minimal()

# Boxplot podle měsíců pro identifikaci sezónnosti
data$month <- factor(month(data$date), labels = month.abb)
p2 <- ggplot(data, aes(x = month, y = new_cases)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Sezónní rozdělení nových případů",
       x = "Měsíc", y = "Počet nových případů") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p1, p2, ncol = 1)

##Dekompozice časové řady
# STL dekompozice
stl_decomp <- stl(ts_data, s.window = "periodic", t.window = 365)
plot(stl_decomp, main = "STL dekompozice časové řady nových případů")
# Klouzavé průměry různých řádů
ma_7 <- SMA(data$new_cases, n = 7)   
ma_30 <- SMA(data$new_cases, n = 30) 
ma_365 <- SMA(data$new_cases, n = 365) 

# Exponenciální vyrovnání - jednodušší přístup
exp_smooth_model <- HoltWinters(ts_data, gamma = FALSE)

# Vytvoříme jednoduché exponenciální vyhlazení pomocí ETS
ets_model <- ets(ts_data, model = "AAN", damped = FALSE)
exp_smooth_fitted <- fitted(ets_model)

# Zajistíme stejnou délku (ETS vrací stejnou délku jako originální data)
n_obs <- length(data$new_cases)

# Pokud je exp_smooth_fitted kratší, doplníme NA na začátek
if(length(exp_smooth_fitted) < n_obs) {
  exp_smooth_full <- c(rep(NA, n_obs - length(exp_smooth_fitted)), 
                       as.numeric(exp_smooth_fitted))
} else {
  exp_smooth_full <- as.numeric(exp_smooth_fitted[1:n_obs])
}

# Grafické srovnání - pouze metody, které máme
trend_data <- data.frame(
  date = data$date,
  original = data$new_cases,
  ma_7 = ma_7,
  ma_30 = ma_30
)

# Přidáme ma_365 pouze pokud má dostatek dat
if(sum(!is.na(ma_365)) > 100) {
  trend_data$ma_365 <- ma_365
}

# Přidáme exponenciální vyhlazení
trend_data$exp_smooth <- exp_smooth_full

trend_plot <- trend_data %>%
  pivot_longer(cols = -date, names_to = "method", values_to = "value") %>%
  filter(!is.na(value)) %>%  # Odfiltrujeme NA hodnoty
  ggplot(aes(x = date, y = value, color = method)) +
  geom_line(alpha = 0.7, size = 0.8) +
  scale_color_manual(values = c("original" = "gray50", "ma_7" = "blue", 
                                "ma_30" = "green", "ma_365" = "red",
                                "exp_smooth" = "purple"),
                     name = "Metoda",
                     labels = c("original" = "Originální data", 
                                "ma_7" = "MA(7)", "ma_30" = "MA(30)", 
                                "ma_365" = "MA(365)", "exp_smooth" = "Exp. vyhlazení")) +
  labs(title = "Srovnání metod vyhlazení trendu",
       x = "Datum", y = "Hodnota") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(trend_plot)

# Analýza trendové složky
trend_component <- stl_decomp$time.series[,"trend"]
cat("Souhrn trendové složky:\n")
print(summary(trend_component))

# Analýza sezónní složky
seasonal_component <- stl_decomp$time.series[,"seasonal"]
cat("\nSouhrn sezónní složky:\n")
print(summary(seasonal_component))
##############################################################################
##Hledání optimálního funkčního modelu
# Příprava dat pro regresní model
data$t <- 1:nrow(data)
data$sin_365 <- sin(2 * pi * data$t / 365.25)
data$cos_365 <- cos(2 * pi * data$t / 365.25)
data$sin_7 <- sin(2 * pi * data$t / 7)
data$cos_7 <- cos(2 * pi * data$t / 7)

# Model 1: Lineární trend + roční sezónnost
model1 <- lm(new_cases ~ t + sin_365 + cos_365, data = data)
# Model 2: Polynomiální trend + roční sezónnost
model2 <- lm(new_cases ~ poly(t, 2) + sin_365 + cos_365, data = data)
# Model 3: Lineární trend + roční + týdenní sezónnost
model3 <- lm(new_cases ~ t + sin_365 + cos_365 + sin_7 + cos_7, data = data)
# Model 4: Polynomiální trend + roční + týdenní sezónnost
model4 <- lm(new_cases ~ poly(t, 2) + sin_365 + cos_365 + sin_7 + cos_7, data = data)
# Srovnání modelů
models <- list(model1, model2, model3, model4)
model_names <- c("Lin+Roční", "Poly+Roční", "Lin+Roční+Týdenní", "Poly+Roční+Týdenní")

# AIC a BIC srovnání
comparison <- data.frame(
  Model = model_names,
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  R_squared = sapply(models, function(x) summary(x)$r.squared)
)
print(comparison)

# Nejlepší model podle AIC
best_func_model <- models[[which.min(comparison$AIC)]]
print(summary(best_func_model))

# Grafické zobrazení fitování
data$fitted_func <- fitted(best_func_model)
ggplot(data, aes(x = date)) +
  geom_line(aes(y = new_cases), color = "black", alpha = 0.6, size = 0.5) +
  geom_line(aes(y = fitted_func), color = "red", size = 1) +
  labs(title = "Nejlepší funkční model",
       x = "Datum", y = "Počet případů",
       subtitle = paste("Model:", model_names[which.min(comparison$AIC)])) +
  theme_minimal()
###############################################################################
##Hledání optimálního SARIMA modelu
# Kontrola stacionarity
adf_test <- adf.test(ts_data)
print(paste("ADF test p-hodnota:", round(adf_test$p.value, 4)))
# KPSS test stacionarity
kpss_test <- kpss.test(ts_data)
print(paste("KPSS test p-hodnota:", round(kpss_test$p.value, 4)))

# ACF a PACF grafy
par(mfrow = c(2,1))
acf(ts_data, main = "ACF původní řady", lag.max = 100)
pacf(ts_data, main = "PACF původní řady", lag.max = 100)

# Diferenciace pokud je potřeba
if(adf_test$p.value > 0.05) {
  ts_diff <- diff(ts_data)
  adf_diff <- adf.test(ts_diff)
  print(paste("ADF test diferencované řady:", round(adf_diff$p.value, 4)))
  
  par(mfrow = c(2,1))
  acf(ts_diff, main = "ACF diferencované řady", lag.max = 100)
  pacf(ts_diff, main = "PACF diferencované řady", lag.max = 100)
} else {
  ts_diff <- ts_data
}

cat("Testování vybraných SARIMA kombinací...\n")

# Vybrané rozumné kombinace na základě ACF/PACF a zkušeností
(fit1 <- Arima(ts_data, order = c(1,0,0), seasonal = list(order = c(0,1,0))))
(fit2 <- Arima(ts_data, order = c(0,0,1), seasonal = list(order = c(0,1,0))))
(fit3 <- Arima(ts_data, order = c(1,0,1), seasonal = list(order = c(0,1,0))))
(fit4 <- Arima(ts_data, order = c(0,0,0), seasonal = list(order = c(1,1,0))))
(fit5 <- Arima(ts_data, order = c(0,0,0), seasonal = list(order = c(0,1,1))))
(fit6 <- Arima(ts_data, order = c(0,0,0), seasonal = list(order = c(1,1,1))))
(fit7 <- Arima(ts_data, order = c(1,0,0), seasonal = list(order = c(1,1,1))))
BIC(fit1); BIC(fit2); BIC(fit3); BIC(fit4); BIC(fit5); BIC(fit6); BIC(fit7)

(auto_sarima_a <- auto.arima(ts_data)) #sezonní diference Zt = 0.004 + 0156*Zt-1 + 0.09*Et-1 + 0.11*Et-2 - 0.53Zt-1
BIC(auto_sarima_a) # automaticky model rozhodne optimalni neni - je prilis komplikovany
coeftest(auto_sarima_a) # model ma spoustu nevyznamnych clenu
(auto_sarima_b <- auto.arima(ts_data, ic = "bic"))
coeftest(auto_sarima_b) # pomoci Bayesovskeho kriteria vybran jednodussi model

# Použijeme auto.arima výsledek jako finální
final_sarima <- auto_sarima_b

# Diagnostika residuí
checkresiduals(final_sarima)
#############################################################################
## Analýza závislostí na jiných řadách
# Příprava dalších časových řad

other_series <- c("new_deaths", "stringency_index", "hosp_patients", "reproduction_rate")
ccf_results <- list()
max_lag <- 30

par(mfrow = c(2,2))
for(i in 1:length(other_series)) {
  y_var <- data[[other_series[i]]]
  y_var <- y_var[!is.na(y_var)]
  x_var <- data$new_cases[1:length(y_var)]
  
  ccf_result <- ccf(x_var, y_var, lag.max = max_lag, 
                    main = paste("CCF:", other_series[i]))
  ccf_results[[other_series[i]]] <- ccf_result
}
# (rest of CCF significance detection as before)
# -------------------------------------------------------------------------
# Create lagged external predictors (same idea as your original)
data$deaths_lag <- c(rep(NA, 1), data$new_deaths[1:(nrow(data)-1)])
data$stringency_lag <- c(rep(NA, 5), data$stringency_index[1:(nrow(data)-5)])
data$hosp_lag <- c(rep(NA, 2), data$hosp_patients[1:(nrow(data)-2)])

external_vars <- cbind(
  deaths_lag = data$deaths_lag,
  stringency_lag = data$stringency_lag,
  hosp_lag = data$hosp_lag
)

complete_idx <- complete.cases(cbind(data$new_cases, external_vars))
ts_complete <- ts(data$new_cases[complete_idx], start = c(2020, 1), frequency = 365.25)
external_complete <- external_vars[complete_idx, ]

# ARIMAX model (same as before)
arimax_model <- auto.arima(ts_complete, xreg = external_complete)
print(summary(arimax_model))

# -------------------------------------------------------------------------
# VAR replacement: simple VAR-like system using equation-by-equation OLS
# We build lagged regressors for each variable and fit one lm() per series.
# We select lag p by minimizing the SUM of AICs across equations (proxy for system AIC).
# -------------------------------------------------------------------------

# Helper: build lagged dataframe for given variables and lag p
build_lags <- function(df, vars, p) {
  n <- nrow(df)
  lagged <- data.frame(row = 1:n)  # temporary index
  for(v in vars) {
    for(l in 1:p) {
      colname <- paste0(v, "_lag", l)
      lagged[[colname]] <- c(rep(NA, l), df[[v]][1:(n - l)])
    }
  }
  lagged$row <- NULL
  return(lagged)
}

# Fit equation system for a chosen p and return list with fits and summed AIC/BIC
fit_var_equations <- function(df, vars, p) {
  lagged <- build_lags(df, vars, p)
  # align targets: drop first p rows
  targets <- df[(p+1):nrow(df), vars, drop = FALSE]
  predictors <- lagged[(p+1):nrow(df), , drop = FALSE]
  fits <- list()
  aic_sum <- 0
  bic_sum <- 0
  rss_list <- c()
  for(v in vars) {
    dat_eq <- cbind(target = targets[[v]], predictors)
    # small safeguard: remove columns with zero variance
    keep_cols <- sapply(dat_eq, function(x) var(x, na.rm = TRUE) > 0)
    dat_eq <- dat_eq[, keep_cols, drop = FALSE]
    fit <- lm(target ~ ., data = as.data.frame(dat_eq))
    fits[[v]] <- fit
    aic_sum <- aic_sum + AIC(fit)
    bic_sum <- bic_sum + BIC(fit)
    rss_list[v] <- mean(residuals(fit)^2)
  }
  return(list(p = p, fits = fits, aic_sum = aic_sum, bic_sum = bic_sum, rss = rss_list))
}

# Autoscan p = 1..max_p and pick best by summed AIC
vars_for_var <- c("new_cases", "new_deaths", "stringency_index", "hosp_patients")
max_p_try <- 6  # reasonable max lag (you can increase if you want)
results_by_p <- list()
for(p_try in 1:max_p_try) {
  res_p <- fit_var_equations(data[complete_idx, ], vars_for_var, p_try)
  results_by_p[[as.character(p_try)]] <- res_p
  cat("p =", p_try, " summed AIC =", round(res_p$aic_sum,2), " summed BIC =", round(res_p$bic_sum,2), "\n")
}

# Choose best p by smallest summed AIC
aic_sums <- sapply(results_by_p, function(x) x$aic_sum)
best_p <- as.integer(names(which.min(aic_sums)))
cat("Chosen VAR-like lag p =", best_p, "\n")
var_like <- results_by_p[[as.character(best_p)]]

# Summarize fitted equations
for(v in names(var_like$fits)) {
  cat("\n--- Equation for", v, "---\n")
  print(summary(var_like$fits[[v]]))
}

# Build a simple forecast for new_cases from this var-like system:
# We'll forecast h steps by iterating using last observed values and predicted lags.
forecast_var_like_new_cases <- function(df, fits_list, vars, p, h = 10) {
  # df: full data frame for vars with no missing at the end
  library(zoo) # for convenience (usually available)
  last_block <- tail(df[, vars], p)  # p last observations
  preds <- numeric(h)
  # we'll keep a rolling window of the past p observations
  window <- as.matrix(last_block)
  for(i in 1:h) {
    # build predictor vector in same column order as used in fit (var_lag1 ... var_lagp for each var)
    newpred <- c()
    for(v in vars) {
      for(l in 1:p) {
        # lag l => value at position (n - (l-1)) in rolling window (most recent is last row)
        newpred <- c(newpred, window[nrow(window) - (l-1), v])
      }
    }
    newdf <- as.data.frame(t(newpred))
    # name columns to match training predictor names
    names(newdf) <- names(coef(fits_list[[1]]))[-1]  # strip intercept; assumes same predictor set - caution
    # Predict for new_cases equation specifically
    # to be safe we will use predict with NA handling by binding to full predictor names (if mismatch, fallback to 0)
    fit_nc <- fits_list[["new_cases"]]
    # create full predictor frame matching fit_nc
    model_terms <- attr(terms(fit_nc), "term.labels")
    pred_frame <- data.frame(matrix(nrow = 1, ncol = length(model_terms)))
    names(pred_frame) <- model_terms
    for(mt in model_terms) {
      if(mt %in% names(newdf)) pred_frame[[mt]] <- newdf[[mt]] else pred_frame[[mt]] <- 0
    }
    pred_val <- as.numeric(predict(fit_nc, newdata = pred_frame))
    # append
    preds[i] <- pred_val
    # update window: append new row (for all vars we only have new_cases prediction; others we set to last known value)
    new_row <- tail(window, 1)
    new_row["new_cases"] <- pred_val
    window <- rbind(window[-1, , drop = FALSE], new_row)
  }
  return(preds)
}

# Forecast 10 steps with var-like system (for new_cases)
h_var <- 10
var_like_preds <- forecast_var_like_new_cases(data, var_like$fits, vars_for_var, best_p, h = h_var)

# For later comparisons, build a "system" AIC/BIC and RMSE metrics
var_like_aic <- var_like$aic_sum
var_like_bic <- var_like$bic_sum
var_like_rmse <- sqrt(mean(unlist(var_like$rss)))

cat("VAR-like system: summed AIC =", var_like_aic, " summed BIC =", var_like_bic, " mean RMSE =", var_like_rmse, "\n")

# -------------------------------------------------------------------------
# Srovnání modelů (replace VAR with VAR-like)
models_comparison <- data.frame(
  Model = c("Funkční", "SARIMA", "ARIMAX", "VAR-like"),
  AIC = c(AIC(best_func_model), AIC(final_sarima), AIC(arimax_model), var_like_aic),
  BIC = c(BIC(best_func_model), BIC(final_sarima), BIC(arimax_model), var_like_bic)
)
print(models_comparison)

# (Residual diagnostics and forecasting code remain as before)
# Use var_like_preds where you previously used VAR predictions if desired (we used for new_cases only)
# -------------------------------------------------------------------------

# Predikce budoucích hodnot (use existing pred_func, pred_sarima, pred_arimax)
# Add var-like predictions into your plotting if you want:
# Example: include var_like_preds in the final plot by constructing a column arimax_pred etc.
# (You can adapt the plotting block to include var_like_preds as another colored line.)

# Final comparison and recommendation (same idea), but include VAR-like metrics if needed
final_comparison <- data.frame(
  Model = c("Funkční (trend+sezónost)", "SARIMA", "ARIMAX", "VAR-like"),
  AIC = c(AIC(best_func_model), AIC(final_sarima), AIC(arimax_model), var_like_aic),
  BIC = c(BIC(best_func_model), BIC(final_sarima), BIC(arimax_model), var_like_bic),
  Parametrů = c(length(coef(best_func_model)),
                length(coef(final_sarima)),
                length(coef(arimax_model)),
                NA), # VAR-like param count is model-by-model; you can compute sum(length(coef())) if desired
  RMSE = c(
    sqrt(mean(residuals(best_func_model)^2)),
    sqrt(mean(residuals(final_sarima)^2)),
    sqrt(mean(residuals(arimax_model)^2)),
    var_like_rmse
  )
)
print("Finální srovnání modelů (včetně VAR-like):")
print(final_comparison)

best_model_idx <- which.min(final_comparison$AIC)
best_model_name <- final_comparison$Model[best_model_idx]

cat("\n=== DOPORUČENÍ ===\n")
cat("Nejlepší model podle AIC:", best_model_name, "\n")
cat("AIC:", round(final_comparison$AIC[best_model_idx], 2), "\n")
cat("RMSE:", round(final_comparison$RMSE[best_model_idx], 2), "\n")

