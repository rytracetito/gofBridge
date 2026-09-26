

# Simuler un processus gamma par la méthode Random Walk Approximation (RWA)
simulate_gamma_rwa <- function(times, A = NULL, a = NULL, nsim = 1, b = 1) {
  if (!is.numeric(times) || !is.vector(times) || length(times) == 0) {
    stop("'times' doit \u00eatre un vecteur num\u00e9rique non vide.")
  }
  if (any(!is.finite(times))) {
    stop("'times' ne doit contenir ni NA, ni NaN, ni valeur infinie.")
  }
  if (times[1] < 0) {
    stop("'times' doit \u00eatre positif ou nul.")
  }
  if (times[1] > 0) {
    warning("'times' ne commence pas \u00e0 0 : l'instant 0 est ajout\u00e9, avec X(0) = 0.")
    times <- c(0, times)
  }
  if (length(times) < 2) {
    stop("'times' doit contenir au moins deux points.")
  }
  if (any(diff(times) <= 0)) {
    stop("'times' doit \u00eatre strictement croissant.")
  }
  if (!is.numeric(nsim) || length(nsim) != 1 || !is.finite(nsim) || nsim < 1 || nsim != floor(nsim)) {
    stop("'nsim' doit \u00eatre un entier positif.")
  }
  if (!is.numeric(b) || length(b) != 1 || !is.finite(b) || b <= 0) {
    stop("'b' doit \u00eatre un r\u00e9el strictement positif.")
  }
  
  if (!is.null(A)) {
    if (!is.function(A)) {
      stop("'A' doit \u00eatre une fonction.")
    }
    if (!is.null(a)) {
      warning("'A' et 'a' sont fournis : 'a' est ignor\u00e9.")
    }
    A_vals <- A(times)
  } else if (!is.null(a)) {
    if (!is.numeric(a) || length(a) != 1 || !is.finite(a) || a <= 0) {
      stop("'a' doit \u00eatre un r\u00e9el strictement positif.")
    }
    A_vals <- a * times
  } else {
    stop("Il faut fournir soit une fonction 'A', soit un scalaire 'a'.")
  }
  
  if (length(A_vals) != length(times) || any(!is.finite(A_vals))) {
    stop("A(times) doit renvoyer une valeur finie pour chaque instant.")
  }
  if (any(diff(A_vals) <= 0)) {
    stop("La fonction A doit \u00eatre strictement croissante.")
  }
  
  diff_A <- diff(A_vals)
  n_intervals <- length(diff_A)
  
  increments <- matrix(rgamma(n_intervals * nsim, shape = rep(diff_A, each = nsim), rate = b),
                       nrow = n_intervals, ncol = nsim, byrow = TRUE)
  
  if (any(!is.finite(increments))) {
    stop("Accroissements simul\u00e9s non finis : v\u00e9rifier 'A', 'a' et 'b'.")
  }
  X <- rbind(rep(0, nsim), apply(increments, 2, cumsum))
  
  rownames(X) <- paste0("t", seq_along(times) - 1)
  colnames(X) <- paste0("sim", seq_len(nsim))
  attr(X, "times") <- times
  
  if (nrow(X) != length(times) || times[1] != 0 || any(X[1, ] != 0) ||
      any(!is.finite(X)) || any(diff(X) < 0)) {
    stop("Sortie incoh\u00e9rente : X et times doivent avoir la m\u00eame longueur, commencer \u00e0 0 et X doit \u00eatre croissant.")
  }
  if (any(diff(X) == 0)) {
    warning(sum(diff(X) == 0), " accroissement(s) de X exactement nul(s) : valeurs inf\u00e9rieures \u00e0 la pr\u00e9cision ",
            "num\u00e9rique (sous-d\u00e9passement de rgamma ou absorption dans la somme cumul\u00e9e), param\u00e8tre de forme trop petit.")
  }
  
  return(X)
}

