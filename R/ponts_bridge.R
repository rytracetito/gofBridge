# Pont gamma (Edith Grall-Maës)
gamma_bridge <- function(X, t, a) {
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
  
  n <- length(X) - 1
  t_n <- t[length(t)]
  teq <- t_n / n
  
  t_star <- seq(0, t_n, teq)
  
  Y <- numeric(length(t_star))
  Y[1] <- 0
  
  if (n > 1) {
    for (i in 2:n) {
      inter <- max(t[t <= t_star[i]])
      if (inter > t_star[i - 1]) {
        t_inf <- inter
        x_inf <- X[t == inter]
      } else {
        t_inf <- t_star[i - 1]
        x_inf <- Y[i - 1]
      }
      t_sup <- min(t[t > t_star[i]])
      x_sup <- X[t == t_sup]
      y <- rbeta(1, a * (t_star[i] - t_inf), a * (t_sup - t_star[i]))
      Y[i] <- x_inf + (x_sup - x_inf) * y
    }
  }
  Y[n + 1] <- X[n + 1]
  
  if (length(Y) != n + 1 || Y[1] != 0 || any(!is.finite(Y)) || any(diff(Y) < 0)) {
    stop("Sortie incoh\u00e9rente : new_obs doit avoir n + 1 valeurs finies, partir de 0 et \u00eatre croissant.")
  }
  if (any(diff(Y) == 0)) {
    warning(sum(diff(Y) == 0), " accroissement(s) reconstruit(s) exactement nul(s) : ",
            "accroissement observ\u00e9 nul ou sous-d\u00e9passement num\u00e9rique de rbeta.")
  }
  
  list(teq = teq, new_obs = Y, t_reg = t_star)
}

# Pont Brownien (Grall Maes) : reconstruction sur la grille régulière de pas t_N / N
# X(t) | X(t_inf) = x_inf, X(t_sup) = x_sup ~ N(x_inf + (x_sup - x_inf)(t - t_inf)/(t_sup - t_inf),
#                                              sigma^2 (t - t_inf)(t_sup - t)/(t_sup - t_inf))
wiener_bridge <- function(X, t, sigma) {
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

  n <- length(X) - 1
  t_n <- t[length(t)]
  teq <- t_n / n

  t_star <- seq(0, t_n, teq)

  Y <- numeric(length(t_star))
  Y[1] <- 0

  if (n > 1) {
    for (i in 2:n) {
      inter <- max(t[t <= t_star[i]])
      if (inter > t_star[i - 1]) {
        t_inf <- inter
        x_inf <- X[t == inter]
      } else {
        t_inf <- t_star[i - 1]
        x_inf <- Y[i - 1]
      }
      t_sup <- min(t[t > t_star[i]])
      x_sup <- X[t == t_sup]

      moy <- x_inf + (x_sup - x_inf) * (t_star[i] - t_inf) / (t_sup - t_inf)
      var_cond <- sigma^2 * (t_star[i] - t_inf) * (t_sup - t_star[i]) / (t_sup - t_inf)
      Y[i] <- rnorm(1, mean = moy, sd = sqrt(var_cond))
    }
  }
  Y[n + 1] <- X[n + 1]

  if (length(Y) != n + 1 || Y[1] != 0 || any(!is.finite(Y))) {
    stop("Sortie incoh\u00e9rente : new_obs doit avoir n + 1 valeurs finies et partir de 0.")
  }

  list(teq = teq, new_obs = Y, t_reg = t_star)
}
