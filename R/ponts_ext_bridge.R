# Pont Gamma séquentiel sur la grille fine (méthode Ext Bridge)
# Dans chaque intervalle observé [t_i, t_{i+1}] découpé en m_i pas de longueur delta :
# Y_k ~ Beta(a*delta, a*(t_{i+1} - s_{k+1})) et accroissement_k = Y_k * (x_{i+1} - x_{s_k}),
# le dernier accroissement étant le reste jusqu'à x_{i+1} (formules (3.1)-(3.2)).
# Les accroissements sont calculés directement pour éviter les zéros numériques.
simulate_gamma_grille <- function(X, t, a, digits = 0, M_max = 100) {
  if (!is.numeric(X) || !is.vector(X) || !is.numeric(t) || !is.vector(t)) {
    stop("'X' et 't' doivent \u00eatre des vecteurs num\u00e9riques.")
  }
  if (length(X) != length(t)) {
    stop("'X' et 't' doivent avoir la m\u00eame longueur (", length(X), " et ", length(t), ").")
  }
  if (length(t) == 0) {
    stop("'X' et 't' ne doivent pas \u00eatre vides.")
  }
  if (any(!is.finite(X)) || any(!is.finite(t))) {
    stop("'X' et 't' ne doivent contenir ni NA, ni NaN, ni valeur infinie.")
  }
  if (!is.numeric(a) || length(a) != 1 || !is.finite(a) || a <= 0) {
    stop("'a' doit \u00eatre un r\u00e9el strictement positif.")
  }
  if (!is.numeric(M_max) || length(M_max) != 1 || !is.finite(M_max) || M_max < 1 || M_max != floor(M_max)) {
    stop("'M_max' doit \u00eatre un entier >= 1.")
  }
  if (t[1] < 0) {
    stop("'t' doit \u00eatre positif ou nul.")
  }
  if (t[1] > 0) {
    warning("'t' ne commence pas \u00e0 0 : le point (0, 0) est ajout\u00e9.")
    t <- c(0, t)
    X <- c(0, X)
  } else if (X[1] != 0) {
    stop("La trajectoire doit partir de X(0) = 0 (X[1] = ", X[1], ").")
  }
  if (length(t) < 2) {
    stop("Il faut au moins deux observations.")
  }
  if (any(diff(t) <= 0)) {
    stop("'t' doit \u00eatre strictement croissant.")
  }
  if (any(diff(X) < 0)) {
    stop("'X' doit \u00eatre croissant : un accroissement n\u00e9gatif est incompatible avec un processus gamma.")
  }

  grille <- construire_grille(t, digits, M_max)
  pas_fine <- grille$pas_fine
  m <- grille$m
  M <- sum(m)

  dx_obs <- diff(X)
  dx_fine <- unlist(lapply(seq_along(m), function(i) {
    inc <- numeric(m[i])
    reste <- dx_obs[i]
    if (m[i] > 1) {
      for (k in 1:(m[i] - 1)) {
        y <- rbeta(1, a * pas_fine, a * (m[i] - k) * pas_fine)
        inc[k] <- y * reste
        reste <- reste - inc[k]
      }
    }
    inc[m[i]] <- reste
    inc
  }))

  Y <- c(0, cumsum(dx_fine))
  Y[grille$indices_obs] <- X

  if (length(dx_fine) != M || any(!is.finite(dx_fine)) || any(dx_fine < 0) ||
      !isTRUE(all.equal(as.numeric(tapply(dx_fine, rep(seq_along(m), m), sum)), dx_obs))) {
    stop("Sortie incoh\u00e9rente : les accroissements fins doivent \u00eatre finis, positifs et redonner les accroissements observ\u00e9s.")
  }
  if (any(dx_fine == 0)) {
    warning(sum(dx_fine == 0), " accroissement(s) reconstruit(s) exactement nul(s) : ",
            "accroissement observ\u00e9 nul ou sous-d\u00e9passement num\u00e9rique.")
  }

  list(pas_fine = pas_fine, t_fine = grille$t_fine, new_obs = Y, dx = dx_fine,
       indices_obs = grille$indices_obs)
}


# Pont brownien séquentiel sur la grille fine (méthode Ext Bridge, rapport § 4.3)
# Dans chaque intervalle observé [t_i, t_{i+1}] découpé en m_i pas de longueur delta :
# X(s_{k+1}) | X(s_k), X(t_{i+1}) ~ N(m_k, v_k) avec
#   m_k = x_{s_k} + (s_{k+1} - s_k) / (t_{i+1} - s_k) * (x_{i+1} - x_{s_k})
#   v_k = sigma^2 (s_{k+1} - s_k) (t_{i+1} - s_{k+1}) / (t_{i+1} - s_k)
# le dernier accroissement étant le reste jusqu'à x_{i+1}.
simulate_wiener_grille <- function(X, t, sigma, digits = 0, M_max = 100) {
  if (!is.numeric(X) || !is.vector(X) || !is.numeric(t) || !is.vector(t)) {
    stop("'X' et 't' doivent \u00eatre des vecteurs num\u00e9riques.")
  }
  if (length(X) != length(t)) {
    stop("'X' et 't' doivent avoir la m\u00eame longueur (", length(X), " et ", length(t), ").")
  }
  if (length(t) == 0) {
    stop("'X' et 't' ne doivent pas \u00eatre vides.")
  }
  if (any(!is.finite(X)) || any(!is.finite(t))) {
    stop("'X' et 't' ne doivent contenir ni NA, ni NaN, ni valeur infinie.")
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || !is.finite(sigma) || sigma <= 0) {
    stop("'sigma' doit \u00eatre un r\u00e9el strictement positif.")
  }
  if (!is.numeric(M_max) || length(M_max) != 1 || !is.finite(M_max) || M_max < 1 || M_max != floor(M_max)) {
    stop("'M_max' doit \u00eatre un entier >= 1.")
  }
  if (t[1] < 0) {
    stop("'t' doit \u00eatre positif ou nul.")
  }
  if (t[1] > 0) {
    warning("'t' ne commence pas \u00e0 0 : le point (0, 0) est ajout\u00e9.")
    t <- c(0, t)
    X <- c(0, X)
  } else if (X[1] != 0) {
    stop("La trajectoire doit partir de X(0) = 0 (X[1] = ", X[1], ").")
  }
  if (length(t) < 2) {
    stop("Il faut au moins deux observations.")
  }
  if (any(diff(t) <= 0)) {
    stop("'t' doit \u00eatre strictement croissant.")
  }
  
  grille <- construire_grille(t, digits, M_max)
  pas_fine <- grille$pas_fine
  m <- grille$m
  M <- sum(m)
  
  dx_obs <- diff(X)
  dx_fine <- unlist(lapply(seq_along(m), function(i) {
    inc <- numeric(m[i])
    reste <- dx_obs[i]
    if (m[i] > 1) {
      for (k in 1:(m[i] - 1)) {
        restant <- m[i] - k + 1
        moy <- reste / restant
        v <- sigma^2 * pas_fine * (restant - 1) / restant
        inc[k] <- moy + sqrt(v) * rnorm(1)
        reste <- reste - inc[k]
      }
    }
    inc[m[i]] <- reste
    inc
  }))
  
  Y <- c(0, cumsum(dx_fine))
  Y[grille$indices_obs] <- X
  
  if (length(dx_fine) != M || any(!is.finite(dx_fine)) ||
      !isTRUE(all.equal(as.numeric(tapply(dx_fine, rep(seq_along(m), m), sum)), dx_obs))) {
    stop("Sortie incoh\u00e9rente : les accroissements fins doivent \u00eatre finis et redonner les accroissements observ\u00e9s.")
  }
  
  list(pas_fine = pas_fine, t_fine = grille$t_fine, new_obs = Y, dx = dx_fine,
       indices_obs = grille$indices_obs)
}
