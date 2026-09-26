# Méthode des moments 
estimer <- function(t, x) {
  if (!is.numeric(t) || !is.vector(t) || !is.numeric(x) || !is.vector(x)) {
    stop("'t' et 'x' doivent \u00eatre des vecteurs num\u00e9riques.")
  }
  if (length(t) != length(x)) {
    stop("'t' et 'x' doivent avoir la m\u00eame longueur (", length(t), " et ", length(x), ").")
  }
  if (length(t) < 3) {
    stop("Il faut au moins 3 observations (2 accroissements) pour estimer nu.")
  }
  if (any(!is.finite(t)) || any(!is.finite(x))) {
    stop("'t' et 'x' ne doivent contenir ni NA, ni NaN, ni valeur infinie.")
  }
  
  dt <- diff(t)
  dx <- diff(x)
  N <- length(dt)
  
  if (any(dt <= 0)) {
    stop("'t' doit \u00eatre strictement croissant.")
  }
  if (any(dx < 0)) {
    stop("'x' doit \u00eatre croissant : un accroissement n\u00e9gatif est incompatible avec un processus gamma.")
  }
  
  mu_hat <- sum(dx) / sum(dt)
  if (mu_hat <= 0) {
    stop("mu_hat = 0 : trajectoire constante, mod\u00e8le gamma inadapt\u00e9.")
  }
  
  nu_hat <- sum((dx - mu_hat * dt)^2 / dt) / N
  if (!is.finite(nu_hat)) {
    stop("nu_hat non fini : d\u00e9passement num\u00e9rique, v\u00e9rifier l'\u00e9chelle des donn\u00e9es.")
  }
  if (nu_hat <= 0) {
    stop("nu_hat = 0 : accroissements exactement proportionnels aux dur\u00e9es, estimation d\u00e9g\u00e9n\u00e9r\u00e9e.")
  }
  
  a_hat <- mu_hat^2 / nu_hat
  b_hat <- mu_hat / nu_hat
  if (!is.finite(a_hat) || !is.finite(b_hat)) {
    stop("Estimation de (a, b) non finie.")
  }
  
  c(a = a_hat, b = b_hat)
}


# Estimation des paramètres du processus de Wiener homogène
# (maximum de vraisemblance : mu = x_N / t_N, sigma2 = (1/N) somme (dx - mu dt)^2 / dt)
estimer_wiener <- function(t, x) {
  if (!is.numeric(t) || !is.vector(t) || !is.numeric(x) || !is.vector(x)) {
    stop("'t' et 'x' doivent \u00eatre des vecteurs num\u00e9riques.")
  }
  if (length(t) != length(x)) {
    stop("'t' et 'x' doivent avoir la m\u00eame longueur (", length(t), " et ", length(x), ").")
  }
  if (length(t) < 3) {
    stop("Il faut au moins 3 observations (2 accroissements) pour estimer sigma.")
  }
  if (any(!is.finite(t)) || any(!is.finite(x))) {
    stop("'t' et 'x' ne doivent contenir ni NA, ni NaN, ni valeur infinie.")
  }

  dt <- diff(t)
  dx <- diff(x)
  N <- length(dt)

  if (any(dt <= 0)) {
    stop("'t' doit \u00eatre strictement croissant.")
  }

  mu_hat <- sum(dx) / sum(dt)

  sigma2_hat <- sum((dx - mu_hat * dt)^2 / dt) / N
  if (!is.finite(sigma2_hat)) {
    stop("sigma2_hat non fini : d\u00e9passement num\u00e9rique, v\u00e9rifier l'\u00e9chelle des donn\u00e9es.")
  }
  if (sigma2_hat <= 0) {
    stop("sigma2_hat = 0 : accroissements exactement proportionnels aux dur\u00e9es, estimation d\u00e9g\u00e9n\u00e9r\u00e9e.")
  }

  c(mu = mu_hat, sigma2 = sigma2_hat, sigma = sqrt(sigma2_hat))
}
