

library(copula)
library(VineCopula)
library(survival)
library(doParallel)

setwd("~/Dropbox/Yohann/application - simulated data")

source("R/fcts_simulated_data3.R")

# ----------------------------------------------------------------------------
# Resultats
# ----------------------------------------------------------------------------

t1 = Sys.time()
results = make_result(v_family = c(1, 3, 4, 5),#
                       v_rate_cens = c(0, 0.1, 0.3, 0.5), #
                       v_n = c(500, 1000),
                       v_h = c(0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4), # 
                       n_repet = 100, #100
                       grid = seq(0.1, 0.9, 0.2),
                       fct_tau = function(x){return(sigmoid(-3 + 4*x))}, # kendall tau vary between 0.047 and 0.73
                       noyau = noyau_quadratique,
                       marginal = "lnorm")
Sys.time() - t1

write.csv(results, "output_2018-04-19/result_simulated_data.csv", row.names = F)


# xtable(
#   rbind(seq(0.1, 0.9, 0.2),
#       round(fct_tau(seq(0.1, 0.9, 0.2)),2),
#       round(BiCopTau2Par(family = 1, tau = fct_tau(seq(0.1, 0.9, 0.2))), 2),
#       round(BiCopTau2Par(family = 3, tau = fct_tau(seq(0.1, 0.9, 0.2))), 2),
#       round(BiCopTau2Par(family = 4, tau = fct_tau(seq(0.1, 0.9, 0.2))), 2)),
# )


