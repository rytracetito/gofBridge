# Fonction PGCD
pgcd <- function(a, b) {
  a <- abs(a)
  b <- abs(b)
  while (b != 0) {
    temp <- b
    b <- a %% b
    a <- temp
  }
  return(a)
}

controle_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1 || !is.finite(digits) || digits < 0 || digits != floor(digits)) {
    stop("'digits' doit \u00eatre un entier positif ou nul.")
  }
}

# Arrondi des instants d'inspection à 'digits' décimales, à faire AVANT la simulation
# (ou le recueil) des données : l'instant 0 est inclus, les doublons sont supprimés
arrondir_temps <- function(t, digits = 0) {
  controle_digits(digits)
  if (!is.numeric(t) || !is.vector(t) || length(t) == 0 || any(!is.finite(t))) {
    stop("'t' doit \u00eatre un vecteur num\u00e9rique fini non vide.")
  }
  if (any(t < 0)) {
    stop("'t' doit \u00eatre positif ou nul.")
  }
  t_arr <- round(t, digits)
  t_arr <- sort(unique(c(0, t_arr)))
  n_perdus <- length(unique(t[t > 0])) - (length(t_arr) - 1)
  if (n_perdus > 0) {
    warning(n_perdus, " instant(s) confondu(s) apr\u00e8s arrondi \u00e0 ", digits, " d\u00e9cimale(s) : doublons supprim\u00e9s.")
  }
  if (length(t_arr) < 2) {
    stop("Moins de deux instants distincts apr\u00e8s arrondi : augmenter 'digits'.")
  }
  t_arr
}

# Construire grille : pas = PGCD des écarts entre instants observés (qui doivent
# être sur le réseau 10^-digits) ; tous les instants observés sont sur la grille
construire_grille <- function(t_obs, digits = 0, M_max = 1e6) {
  controle_digits(digits)
  if (!is.numeric(t_obs) || !is.vector(t_obs) || length(t_obs) < 2 || any(!is.finite(t_obs))) {
    stop("'t_obs' doit \u00eatre un vecteur num\u00e9rique fini d'au moins 2 \u00e9l\u00e9ments.")
  }
  if (any(diff(t_obs) <= 0)) {
    stop("'t_obs' doit \u00eatre strictement croissant.")
  }
  if (t_obs[1] != 0) {
    stop("'t_obs' doit commencer \u00e0 0.")
  }

  facteur <- 10^digits
  t_int <- round(t_obs * facteur)
  if (any(abs(t_obs * facteur - t_int) > 1e-6)) {
    stop("Les instants ne sont pas des multiples de 10^-", digits, " : appliquer arrondir_temps(t, digits) ",
         "avant de simuler les donn\u00e9es, ou indiquer dans 'digits' le nombre de d\u00e9cimales des instants.")
  }

  ecarts <- diff(t_int)
  delta_int <- Reduce(pgcd, ecarts)
  m <- ecarts / delta_int
  if (sum(m) > M_max) {
    stop("La grille fine comporte ", sum(m), " pas (> M_max = ", M_max, ") : ",
         "r\u00e9duire 'digits' ou augmenter 'M_max'.")
  }

  pas_fine <- delta_int / facteur
  t_fine <- (0:sum(m)) * pas_fine
  indices_obs <- c(1, 1 + cumsum(m))
  t_fine[indices_obs] <- t_obs

  list(
    t_fine      = t_fine,
    pas_fine    = pas_fine,
    m           = m,
    indices_obs = indices_obs
  )
}
