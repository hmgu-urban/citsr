set.seed(123)

T <- 1:100
E <- rep(c(0,1), each = length(T))
baseline <- 50 + 2 * T
offset <- c(rep(0, 100), rep(15, 100))
y <- baseline + offset + rnorm(200, 0, 5)

intervention_time <- 50
I <- ifelse(rep(T, 2) > intervention_time, 1, 0)

y[E == 0 & T > intervention_time] <- y[E == 0 & T > intervention_time] + 10
y[E == 1 & T > intervention_time] <- y[E == 1 & T > intervention_time] + 40

df_cits_example <- data.frame(
  y = y,
  T = rep(T,2),
  I = I,
  E = E
)

head(df_cits_example)
tail(df_cits_example)
