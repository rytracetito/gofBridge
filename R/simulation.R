# Contrôle des instants d'observation (l'instant 0 est ajouté s'il manque)
controle_temps <- function(t_obs) {
  if (!is.numeric(t_obs) || !is.vector(t_obs) || length(t_obs) == 0) {
    stop("'t_obs' doit \u00eatre un vecteur num\u00e9rique non vide.")
  }
  if (any(!is.finite(t_obs))) {
    stop("'t_obs' ne doit contenir ni NA, ni NaN, ni valeur infinie.")
  }
  if (t_obs[1] < 0) {
    stop("'t_obs' doit \u00eatre positif ou nul.")
  }
  if (t_obs[1] > 0) {
    warning("'t_obs' ne commence pas \u00e0 0 : l'instant 0 est ajout\u00e9, avec X(0) = 0.")
    t_obs <- c(0, t_obs)
  }
  if (length(t_obs) < 2) {
    stop("'t_obs' doit contenir au moins deux instants.")
  }
  if (any(diff(t_obs) <= 0)) {
    stop("'t_obs' doit \u00eatre strictement croissant.")
  }
  as.numeric(t_obs)
}

# Contrôle d'un paramètre réel (strictement positif si positif = TRUE)
controle_parametre <- function(x, nom, positif = TRUE) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || (positif && x <= 0)) {
    stop("'", nom, "' doit \u00eatre un r\u00e9el", if (positif) " strictement positif" else "", ".")
  }
}

# Contrôle de la trajectoire : même longueur, départ en (0, 0), valeurs finies,
# croissante si croissante = TRUE (modèles gamma et gaussien inverse)
controle_trajectoire <- function(t, X, croissante = FALSE) {
  if (length(t) != length(X)) {
    stop("t et X n'ont pas la m\u00eame longueur (", length(t), " et ", length(X), ").")
  }
  if (t[1] != 0 || X[1] != 0) {
    stop("La trajectoire doit commencer en (t, X) = (0, 0).")
  }
  if (any(!is.finite(t)) || any(!is.finite(X))) {
    stop("La trajectoire contient des valeurs non finies.")
  }
  if (any(diff(t) <= 0)) {
    stop("t doit \u00eatre strictement croissant.")
  }
  if (croissante && any(diff(X) < 0)) {
    stop("La trajectoire X doit \u00eatre croissante.")
  }
}

# Processus gamma homogène (accroissements Gamma(a dt, b))
sim_HGP <- function(t_obs, a, b) {
  controle_parametre(a, "a")
  controle_parametre(b, "b")

  X <- simulate_gamma_rwa(times = t_obs, a = a, b = b, nsim = 1)
  t <- attr(X, "times")
  X <- unname(X[, 1])

  controle_trajectoire(t, X, croissante = TRUE)
  return(data.frame(t = t, X = X))
}

# Processus gamma non homogène (fonction de forme A(t) = a t^alpha)
sim_NHGP <- function(t_obs, a, b, alpha) {
  controle_parametre(a, "a")
  controle_parametre(b, "b")
  controle_parametre(alpha, "alpha")

  A_func <- function(t) a * t^alpha
  X <- simulate_gamma_rwa(times = t_obs, A = A_func, b = b, nsim = 1)
  t <- attr(X, "times")
  X <- unname(X[, 1])

  controle_trajectoire(t, X, croissante = TRUE)
  return(data.frame(t = t, X = X))
}

# Processus gaussien inverse homogène
sim_IGP <- function(t_obs, mu, lambda) {
  controle_parametre(mu, "mu")
  controle_parametre(lambda, "lambda")
  t_obs <- controle_temps(t_obs)

  n <- length(t_obs) - 1
  delta_t <- diff(t_obs)

  delta_X <- statmod::rinvgauss(n, mean = mu * delta_t, shape = lambda * delta_t^2)

  if (length(delta_X) != n || any(!is.finite(delta_X)) || any(delta_X < 0)) {
    stop("Incr\u00e9ments invalides g\u00e9n\u00e9r\u00e9s (n\u00e9gatifs ou non finis) : ",
         "v\u00e9rifie mu, lambda et t_obs")
  }

  X <- c(0, cumsum(delta_X))
  if (any(diff(X) == 0)) {
    warning(sum(diff(X) == 0), " accroissement(s) de X exactement nul(s) : valeurs inf\u00e9rieures \u00e0 la pr\u00e9cision ",
            "num\u00e9rique (sous-d\u00e9passement de rinvgauss ou absorption dans la somme cumul\u00e9e).")
  }

  controle_trajectoire(t_obs, X, croissante = TRUE)
  return(data.frame(t = t_obs, X = X))
}

# Processus de Wiener homogène X(t) = mu t + sigma B(t)
sim_HWP <- function(t_obs, mu, sigma) {
  controle_parametre(mu, "mu", positif = FALSE)
  controle_parametre(sigma, "sigma")
  t_obs <- controle_temps(t_obs)

  dt <- diff(t_obs)
  dx <- rnorm(length(dt), mean = mu * dt, sd = sigma * sqrt(dt))
  X <- c(0, cumsum(dx))

  controle_trajectoire(t_obs, X)
  return(data.frame(t = t_obs, X = X))
}

# Processus de Wiener non homogène X(t) = mu Lambda(t) + sigma B(tau(t))
# Lambda et tau : fonctions déterministes continues, nulles en 0, tau croissante
sim_NHWP <- function(t_obs, mu, sigma, Lambda, tau = function(t) t) {
  controle_parametre(mu, "mu", positif = FALSE)
  controle_parametre(sigma, "sigma")
  if (!is.function(Lambda) || !is.function(tau)) {
    stop("'Lambda' et 'tau' doivent \u00eatre des fonctions.")
  }
  t_obs <- controle_temps(t_obs)

  L <- Lambda(t_obs)
  ta <- tau(t_obs)
  if (length(L) != length(t_obs) || any(!is.finite(L))) {
    stop("Lambda(t) doit renvoyer une valeur finie pour chaque instant.")
  }
  if (length(ta) != length(t_obs) || any(!is.finite(ta))) {
    stop("tau(t) doit renvoyer une valeur finie pour chaque instant.")
  }
  if (L[1] != 0 || ta[1] != 0) {
    stop("Lambda(0) et tau(0) doivent valoir 0.")
  }
  if (any(diff(ta) < 0)) {
    stop("tau doit \u00eatre croissante.")
  }

  dx <- mu * diff(L) + sigma * sqrt(diff(ta)) * rnorm(length(t_obs) - 1)
  X <- c(0, cumsum(dx))

  controle_trajectoire(t_obs, X)
  return(data.frame(t = t_obs, X = X))
}
