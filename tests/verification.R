library(gofBridge)

erreur <- function(expr) {
  tryCatch({ expr; FALSE }, error = function(e) TRUE)
}

set.seed(2026)
t_reg <- 0:40
t_irr <- suppressWarnings(arrondir_temps(runif(40, 0, 80)))

# Plan régulier : accroissements observés, tous les tests appliqués
Xg <- sim_HGP(t_reg, 2, 1)
Xw <- sim_HWP(t_reg, 2, 1.5)
for (m in c("Bridge", "ExtBridge")) {
  rg <- gof.process(Xg, "HGP", m, boot = 19)
  rw <- gof.process(Xw, "HWP", m)
  stopifnot(rg$methode == "direct", identical(rg$dx, diff(Xg$X)), nrow(rg$resultats) == 4,
            !anyNA(rg$resultats$Decision))
  stopifnot(rw$methode == "direct", identical(rw$dx, diff(Xw$X)), nrow(rw$resultats) == 6,
            !anyNA(rw$resultats$Decision))
}
r0 <- suppressWarnings(gof.process(Xw[-1, ], "HWP"))
stopifnot(r0$methode == "direct", r0$n_dx == 40)

# Plan irrégulier, méthode Bridge : N accroissements, pas t_N / N, X(t_N) conservé
Xg <- sim_HGP(t_irr, 2, 1)
Xw <- sim_HWP(t_irr, 2, 1.5)
N <- nrow(Xg) - 1
rg <- gof.process(Xg, "HGP", "Bridge", boot = 19)
rw <- gof.process(Xw, "HWP", "Bridge")
stopifnot(rg$pont == "gamma", rg$n_dx == N, isTRUE(all.equal(rg$pas, max(t_irr) / N)),
          isTRUE(all.equal(sum(rg$dx), Xg$X[N + 1])),
          isTRUE(all.equal(rg$parametre_pont[["a"]], as.numeric(estimer(Xg$t, Xg$X)["a"]))))
stopifnot(rw$pont == "brownien", rw$n_dx == N, isTRUE(all.equal(rw$pas, max(t_irr) / N)),
          isTRUE(all.equal(sum(rw$dx), Xw$X[N + 1])),
          isTRUE(all.equal(rw$parametre_pont[["sigma"]], as.numeric(estimer_wiener(Xw$t, Xw$X)["sigma"]))))

# Plan irrégulier, méthode Ext Bridge : observations conservées exactement
obs_conservees <- function(r, X) {
  Y <- c(0, cumsum(r$dx))
  isTRUE(all.equal(Y[match(X$t, r$t_pont)], X$X))
}
rg <- gof.process(Xg, "HGP", "ExtBridge", boot = 19)
rw <- gof.process(Xw, "HWP", "ExtBridge")
stopifnot(rg$pont == "gamma", rg$pas == 1, rg$n_dx == max(t_irr), obs_conservees(rg, Xg))
stopifnot(rw$pont == "brownien", rw$pas == 1, rw$n_dx == max(t_irr), obs_conservees(rw, Xw))

# Test non applicable : D'Agostino-Pearson avec moins de 20 accroissements
r15 <- suppressWarnings(gof.process(sim_HWP(0:15, 2, 1.5), "HWP"))
stopifnot(is.na(r15$resultats$Decision[6]), nchar(r15$resultats$remarque[6]) > 0,
          !anyNA(r15$resultats$Decision[1:5]))

# Contrôles : chaque appel doit produire une erreur
stopifnot(
  erreur(gof.process(Xg)),
  erreur(gof.process(Xg, "NHGP")),
  erreur(gof.process(Xg, c("HGP", "HWP"))),
  erreur(gof.process(Xg, "HGP", "Ext")),
  erreur(gof.process(Xg, "HGP", digits = -1)),
  erreur(gof.process(Xg, "HGP", M_max = 10.5)),
  erreur(gof.process(Xw, "HWP", alpha = 1)),
  erreur(gof.process(Xg, "HGP", boot = 10)),
  erreur(gof.process(Xg, "HGP", a_poids = -1)),
  erreur(gof.process(as.matrix(unname(Xg)), "HGP")),
  erreur(gof.process(Xg$X, "HGP")),
  erreur(gof.process(transform(Xw, X = replace(X, 3, NA)), "HWP")),
  erreur(gof.process(Xw[1, ], "HWP")),
  erreur(gof.process(Xw[c(1, 3, 2, 4:nrow(Xw)), ], "HWP")),
  erreur(gof.process(transform(Xw, X = X + 1), "HWP")),
  erreur(gof.process(sim_HWP(t_irr, 0, 1.5), "HGP")),
  erreur(gof.process(sim_HWP(0:5, 2, 1.5), "HWP")),
  erreur(gof.process(data.frame(t = c(0, 1.5, 4, 7.25, 9:20), X = 0:15), "HWP", "ExtBridge")),
  erreur(gof.process(Xw, "HWP", "ExtBridge", M_max = 50))
)
stopifnot(gof.process(data.frame(t = c(0, 1.5, 4, 7.25, 9:20), X = cumsum(c(0, rnorm(15)))), "HWP",
                 "ExtBridge", digits = 2, M_max = 10000)$pas == 0.25)

# Sous H1 : trajectoire gamma de forme faible testée comme HWP
stopifnot(all(gof.process(sim_HGP(0:60, 0.3, 1), "HWP")$resultats$Decision))
