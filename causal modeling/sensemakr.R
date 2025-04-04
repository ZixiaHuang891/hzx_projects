# loads package
library(sensemakr)

# loads dataset
data("darfur")

# runs regression model
model <- lm(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
              pastvoted + hhsize_darfur + female + village, data = darfur)

# runs sensemakr for sensitivity analysis
sensitivity <- sensemakr(model = model, 
                         treatment = "directlyharmed",
                         benchmark_covariates = "female",
                         kd = 1:3)

# short description of results
sensitivity

# long description of results
summary(sensitivity)

# plot bias contour of point estimate
plot(sensitivity)

# plot bias contour of t-value
plot(sensitivity, sensitivity.of = "t-value")

# plot bias contour of lower limit of confidence interval
plot(sensitivity, sensitivity.of = "lwr")

# plot extreme scenario
plot(sensitivity, type = "extreme")
