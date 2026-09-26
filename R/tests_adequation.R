# Nombre de rééchantillonnages bootstrap par défaut et nombre minimal d'accroissements
K <- 199
n_min <- 10

# --- Contrôles communs ---
# positif = TRUE pour la loi gamma (accroissements strictement positifs)
controle_dx <- function(dx, n_min_test = n_min, positif = FALSE) {
  if (!is.numeric(dx) || !is.vector(dx)) {
    stop("'dx' doit \u00eatre un vecteur num\u00e9rique.")
  }
  if (any(!is.finite(dx))) {
    stop("'dx' ne doit contenir ni NA, ni NaN, ni valeur infinie.")
  }
  if (length(dx) < max(n_min, n_min_test)) {
    stop("'dx' doit contenir au moins ", max(n_min, n_min_test), " accroissements (", length(dx), " fournis).")
  }
  if (positif && any(dx <= 0)) {
    stop("'dx' doit \u00eatre strictement positif (", sum(dx <= 0), " valeur(s) <= 0).")
  }
  if (length(unique(dx)) == 1) {
    stop("'dx' est constant : test d'ad\u00e9quation impossible.")
  }
}

controle_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' doit \u00eatre un r\u00e9el dans ]0, 1[.")
  }
}

controle_boot <- function(boot, alpha) {
  controle_alpha(alpha)
  if (!is.numeric(boot) || length(boot) != 1 || !is.finite(boot) || boot != floor(boot) ||
      (boot + 1) * alpha < 1 - 1e-12) {
    stop("'boot' doit \u00eatre un entier tel que (boot + 1) * alpha >= 1, soit boot >= ",
         ceiling(1 / alpha - 1 - 1e-12), " pour alpha = ", alpha, ".")
  }
}

# ============================================================================
# Tests d'adéquation à la loi gamma (paramètres estimés, bootstrap paramétrique)
# ============================================================================

# --- Test de Henze-Meintanis-Ebner (HME1) ---
gof.gp.henze <- function(dx, a_poids = 1, boot = K, alpha = 0.05) {
  controle_dx(dx, positif = TRUE)
  controle_boot(boot, alpha)
  gofgamma::test.HME1(dx, a = a_poids, boot = boot, alpha = alpha)
}

# --- Test de Villaseñor-González-Estrada (VGE) ---
gof.gp.vge <- function(dx, boot = K, alpha = 0.05) {
  controle_dx(dx, positif = TRUE)
  controle_boot(boot, alpha)

  fit   <- goft::gamma_fit(dx)
  shape <- fit["shape", 1]
  scale <- fit["scale", 1]

  T_obs <- abs(as.numeric(goft::gamma_test(dx)$statistic))

  T_boot <- replicate(boot, {
    dx_star <- rgamma(length(dx), shape = shape, scale = scale)
    abs(as.numeric(goft::gamma_test(dx_star)$statistic))
  })

  if (!is.finite(T_obs) || any(!is.finite(T_boot))) {
    stop("Statistique VGE non d\u00e9finie (observ\u00e9e ou bootstrap).")
  }
  p_val <- mean(T_boot >= T_obs)

  list(statistic = T_obs, p.value = p_val, shape = shape, scale = scale,
       Decision = p_val < alpha, sig.level = alpha, boot.run = boot)
}

# --- Test de Kolmogorov-Smirnov ---
gof.gp.ks <- function(dx, boot = K, alpha = 0.05) {
  controle_dx(dx, positif = TRUE)
  controle_boot(boot, alpha)
  gofgamma::test.KS(dx, boot = boot, alpha = alpha)
}

# --- Test de Betsch-Ebner ---
gof.gp.betsch <- function(dx, boot = K, alpha = 0.05) {
  controle_dx(dx, positif = TRUE)
  controle_boot(boot, alpha)
  gofgamma::test.BE(dx, boot = boot, alpha = alpha)
}

# ============================================================================
# Tests de normalité (paramètres estimés)
# ============================================================================

resultat_test <- function(nom, statistique, p_valeur, alpha, n) {
  list(test = nom, statistic = unname(statistique), p.value = unname(p_valeur),
       Decision = unname(p_valeur) < alpha, sig.level = alpha, n = n)
}

# --- Test de Shapiro-Wilk ---
gof.wiener.shapiro <- function(dx, alpha = 0.05) {
  controle_dx(dx)
  controle_alpha(alpha)
  if (length(dx) > 5000) {
    stop("Le test de Shapiro-Wilk est limit\u00e9 \u00e0 5000 observations.")
  }
  r <- stats::shapiro.test(dx)
  resultat_test("Shapiro-Wilk", r$statistic, r$p.value, alpha, length(dx))
}

# --- Test de Lilliefors (Kolmogorov-Smirnov, paramètres estimés) ---
gof.wiener.lilliefors <- function(dx, alpha = 0.05) {
  controle_dx(dx)
  controle_alpha(alpha)
  r <- nortest::lillie.test(dx)
  resultat_test("Lilliefors (KS)", r$statistic, r$p.value, alpha, length(dx))
}

# --- Test de Jarque-Bera ---
gof.wiener.jarquebera <- function(dx, alpha = 0.05) {
  controle_dx(dx)
  controle_alpha(alpha)
  r <- tseries::jarque.bera.test(dx)
  resultat_test("Jarque-Bera", r$statistic, r$p.value, alpha, length(dx))
}

# --- Test d'Anderson-Darling ---
gof.wiener.anderson <- function(dx, alpha = 0.05) {
  controle_dx(dx)
  controle_alpha(alpha)
  r <- nortest::ad.test(dx)
  resultat_test("Anderson-Darling", r$statistic, r$p.value, alpha, length(dx))
}

# --- Test de Cramér-von Mises ---
gof.wiener.cramer <- function(dx, alpha = 0.05) {
  controle_dx(dx)
  controle_alpha(alpha)
  # cvm.test avertit quand la p-valeur sort de sa plage de calcul (p < 7.37e-10 ou p > 0.25) :
  # la valeur renvoyée est alors la borne, la décision au seuil alpha n'en dépend pas
  r <- suppressWarnings(nortest::cvm.test(dx))
  resultat_test("Cramer-von Mises", r$statistic, r$p.value, alpha, length(dx))
}

# --- Test omnibus de D'Agostino-Pearson (n >= 20) ---
gof.wiener.agostino <- function(dx, alpha = 0.05) {
  controle_dx(dx, n_min_test = 20)
  controle_alpha(alpha)
  r <- fBasics::dagoTest(dx)
  resultat_test("D'Agostino-Pearson", r@test$statistic[1], r@test$p.value[1], alpha, length(dx))
}
