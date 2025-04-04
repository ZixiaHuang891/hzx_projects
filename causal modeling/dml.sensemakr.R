# loads package
library(dml.sensemakr)
library(ggplot2)

## loads data
data("pension")

# set treatment, outcome and covariates
y <- pension$net_tfa  # net total financial assets
d <- pension$e401     # 401K eligibility
x <- model.matrix(~ -1 + age + inc  + educ+ fsize + marr + twoearn + pira + hown, 
                  data = pension)

# run DML (Partially Linear Model)
dml.401k <- dml(y, d, x, model = "plm", cf.folds = 5, cf.reps = 5)

# summary of DML results with median method (default)
summary(dml.401k)


# sensitivity analysis
sens.401k <- sensemakr(dml.401k, 
                       cf.y = 0.04, cf.d = 0.03, rho2 = 1,
                       bound_label = "Max Match")

# see results
sens.401k


# long descriptions of results
summary(sens.401k)

# contout plots
plot(sens.401k)

# heterogeneous effects

# income quartiles
g1 <- cut(x[,"inc"], quantile(x[,"inc"], c(0, 0.25,.5,.75,1), na.rm = TRUE),
          labels = c("q1", "q2", "q3", "q4"), include.lowest = T)

# compute GATE
dml.401k <- dml_gate(dml.401k, groups = g1)

# ggplot
# coefficient plot under conditional ignorability
group.names <- paste0("gate.q",1:4)
df   <- data.frame(groups = 1:4, estimate = coef(dml.401k)[group.names])
cis  <- confint(dml.401k)[group.names, ]
cis  <- setNames(as.data.frame(cis), c("lwr.ci", "upr.ci"))
df   <- cbind(df, cis)
ggplot(df, aes(x = groups, y = estimate)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.ci, ymax = upr.ci), alpha = 0.1, col = "blue", fill = "blue") +
  theme_bw() + 
  xlab("Income Groups by Quartiles") + 
  ylab("ATE")

# confidence bounds plot
bds   <- confidence_bounds(dml.401k, cf.y = 0.04, cf.d = 0.03, level = 0)
bds   <- setNames(as.data.frame(bds), c("lwr.bound", "upr.bound"))
cbds  <- confidence_bounds(dml.401k, cf.y = 0.04, cf.d = 0.03, level = .95)
cbds  <- setNames(as.data.frame(cbds), c("lwr.cbound", "upr.cbound"))
df2   <- cbind(df, bds[-1,], cbds[-1, ])
ggplot(df2, aes(x = groups, y = estimate)) + geom_line() +
  geom_ribbon(aes(ymin = lwr.bound, ymax = upr.bound),   alpha = 0.1, col = "red", fill = "red") +
  geom_ribbon(aes(ymin = lwr.cbound, ymax = upr.cbound), alpha = 0.1, col = "blue", fill = "blue") +
  theme_bw() + 
  xlab("Income Groups by Quartiles") + 
  ylab("ATE")
