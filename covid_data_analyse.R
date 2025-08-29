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

Sys.setlocale("LC_TIME", "C")

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
new_deaths <- new_cases * 0.02 + rnorm(n, 0, 2)
new_deaths <- pmax(0, new_deaths)

stringency_index <- 30 + 40 * sin(2 * pi * as.numeric(dates) / 365.25) + rnorm(n, 0, 5)
stringency_index <- pmin(100, pmax(0, stringency_index))

hosp_patients <- new_cases * 0.15 + rnorm(n, 0, 5)
hosp_patients <- pmax(0, hosp_patients)

data <- data.frame(date = dates, location = "Czech Republic", new_cases = new_cases, new_deaths = new_deaths, stringency_index = stringency_index, hosp_patients = hosp_patients, total_cases = cumsum(new_cases), reproduction_rate = 1 + 0.3 * sin(2 * pi * as.numeric(dates) / 365.25) + rnorm(n, 0, 0.1))

# Převod na časovou řadu
ts_data <- ts(data$new_cases, start = c(2020, 1), frequency = 365.25)
head(data)

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

# STL dekompozice
stl_decomp <- stl(ts_data, s.window = "periodic", t.window = 365)
plot(stl_decomp, main = "STL dekompozice časové řady nových případů")

# Klouzavé průměry různých řádů
ma_7 <- SMA(data$new_cases, n = 7)   
ma_30 <- SMA(data$new_cases, n = 30) 
ma_365 <- SMA(data$new_cases, n = 365) 

# Exponenciální vyrovnání
exp_smooth <- HoltWinters(ts_data, gamma = FALSE)

# Zajistíme stejnou délku všech vektorů
n_obs <- length(data$new_cases)
fitted_values <- fitted(exp_smooth)

# Pokud má fitted jiný začátek, vyplníme NA
exp_smooth_full <- rep(NA, n_obs)
start_idx <- length(data$new_cases) - length(fitted_values) + 1
exp_smooth_full[start_idx:n_obs] <- as.numeric(fitted_values[,1])

# Grafické srovnání
trend_data <- data.frame(date = data$date, original = data$new_cases, ma_7 = ma_7, ma_30 = ma_30, ma_365 = ma_365, exp_smooth = exp_smooth_full)

trend_plot <- trend_data %>%
  pivot_longer(cols = -date, names_to = "method", values_to = "value") %>%
  ggplot(aes(x = date, y = value, color = method)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = c("original" = "gray", "ma_7" = "blue", 
                                "ma_30" = "green", "ma_365" = "red",
                                "exp_smooth" = "purple")) +
  labs(title = "Srovnání metod vyhlazení trendu",
       x = "Datum", y = "Hodnota",
       color = "Metoda") +
  theme_minimal()

print(trend_plot)

# Analýza trendové složky
trend_component <- stl_decomp$time.series[,"trend"]
summary(trend_component)

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

# Automatické hledání SARIMA modelu
auto_sarima <- auto.arima(ts_data, seasonal = TRUE, stepwise = FALSE, 
                          approximation = FALSE, trace = TRUE)
print(summary(auto_sarima))

# Manuální testování několika SARIMA modelů
sarima_models <- list()
aic_values <- c()

# Možné kombinace parametrů
params <- expand.grid(p = 0:2, d = 0:1, q = 0:2, 
                      P = 0:2, D = 0:1, Q = 0:2)
params <- params[params$d + params$D <= 2, ]  # omezení na rozumné diferenciace

# Testování prvních 20 kombinací (pro časové důvody)
for(i in 1:min(20, nrow(params))) {
  tryCatch({
    model <- arima(ts_data, order = c(params$p[i], params$d[i], params$q[i]),
                   seasonal = list(order = c(params$P[i], params$D[i], params$Q[i]), 
                                   period = 365))
    sarima_models[[i]] <- model
    aic_values[i] <- AIC(model)
  }, error = function(e) {
    aic_values[i] <- NA
  })
}

# Nejlepší model
best_idx <- which.min(aic_values)
if(length(best_idx) > 0 && !is.na(best_idx)) {
  manual_best_sarima <- sarima_models[[best_idx]]
  print(paste("Nejlepší manuální SARIMA: ARIMA", 
              paste(manual_best_sarima$arma[c(1,6,2)], collapse = ","),
              "x", paste(manual_best_sarima$arma[c(3,7,4)], collapse = ","),
              "[", manual_best_sarima$arma[5], "]"))
  print(paste("AIC:", round(AIC(manual_best_sarima), 2)))
}

# Použijeme auto.arima výsledek jako finální
final_sarima <- auto_sarima

# Diagnostika residuí
checkresiduals(final_sarima)


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

if(grepl("poly", deparse(best_func_model$call))) {
  pred_func <- predict(best_func_model, newdata = future_data, interval = "prediction")
} else {
  pred_func <- predict(best_func_model, newdata = future_data, interval = "prediction")
}

# SARIMA predikce
pred_sarima <- forecast(final_sarima, h = horizon)

# ARIMAX predikce (potřebujeme budoucí hodnoty externích proměnných)
# Pro zjednodušení použijeme poslední hodnoty
last_external <- external_complete[nrow(external_complete), ]
future_external <- matrix(rep(last_external, horizon), nrow = horizon, byrow = TRUE)
pred_arimax <- forecast(arimax_model, xreg = future_external, h = horizon)

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