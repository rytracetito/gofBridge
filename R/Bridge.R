# Tests d'adéquation au processus gamma homogène (HGP) ou au processus de Wiener homogène (HWP)
# - plan régulier : tests sur les accroissements observés ;
# - plan irrégulier : paramètre du pont estimé sous H0 (a pour HGP, sigma pour HWP), reconstruction
#   par pont gamma ou pont brownien selon 'methode' :
#     "Bridge"    : grille régulière de pas t_N / N (Grall-Maës, 2012),
#     "ExtBridge" : grille fine contenant les observations (arguments 'digits' et 'M_max'),
#   puis tests sur les accroissements reconstruits.
# Tous les tests du processus sont appliqués aux mêmes accroissements.
# Decision = TRUE : rejet de H0 au niveau alpha.
# 'boot' et 'a_poids' ne servent que pour HGP ; 'digits' et 'M_max' que pour "ExtBridge".
Bridge <- function(X, processus, methode = "Bridge", digits = 0, M_max = 100,
                   alpha = 0.05, boot = 199, a_poids = 1) {
  if (missing(processus) || !is.character(processus) || length(processus) != 1 ||
      !(processus %in% c("HGP", "HWP"))) {
    stop("'processus' doit valoir 'HGP' (processus gamma homog\u00e8ne) ou 'HWP' (processus de Wiener homog\u00e8ne).")
  }
  if (!is.character(methode) || length(methode) != 1 || !(methode %in% c("Bridge", "ExtBridge"))) {
    stop("'methode' doit valoir 'Bridge' ou 'ExtBridge'.")
  }
  controle_digits(digits)
  if (!is.numeric(M_max) || length(M_max) != 1 || !is.finite(M_max) || M_max < 1 || M_max != floor(M_max)) {
    stop("'M_max' doit \u00eatre un entier >= 1.")
  }
  controle_alpha(alpha)
  if (processus == "HGP") {
    controle_boot(boot, alpha)
    if (!is.numeric(a_poids) || length(a_poids) != 1 || !is.finite(a_poids) || a_poids <= 0) {
      stop("'a_poids' doit \u00eatre un r\u00e9el strictement positif.")
    }
  }

  if (!(is.data.frame(X) || is.matrix(X)) || !all(c("t", "X") %in% colnames(X))) {
    stop("'X' doit \u00eatre un data.frame ou une matrice avec des colonnes nomm\u00e9es 't' et 'X'.")
  }
  if (is.data.frame(X)) {
    t <- X[["t"]]
    x <- X[["X"]]
  } else {
    t <- unname(X[, "t"])
    x <- unname(X[, "X"])
  }
  if (!is.numeric(t) || !is.numeric(x)) {
    stop("Les colonnes 't' et 'X' doivent \u00eatre num\u00e9riques.")
  }
  t <- as.numeric(t)
  x <- as.numeric(x)
  if (any(!is.finite(t)) || any(!is.finite(x))) {
    stop("Les colonnes 't' et 'X' ne doivent contenir ni NA, ni NaN, ni valeur infinie.")
  }
  if (length(t) < 2) {
    stop("Il faut au moins deux observations.")
  }
  if (any(diff(t) <= 0)) {
    stop("'t' doit \u00eatre strictement croissant.")
  }
  if (t[1] < 0) {
    stop("'t' doit \u00eatre positif ou nul.")
  }
  if (t[1] > 0) {
    warning("'t' ne commence pas \u00e0 0 : le point (0, 0) est ajout\u00e9 avant l'estimation et le pont.")
    t <- c(0, t)
    x <- c(0, x)
  } else if (x[1] != 0) {
    stop("La trajectoire doit partir de X(0) = 0 (X[1] = ", x[1], ").")
  }
  if (processus == "HGP" && any(diff(x) < 0)) {
    stop(sum(diff(x) < 0), " accroissement(s) observ\u00e9(s) n\u00e9gatif(s) : trajectoire incompatible avec ",
         "un processus gamma (croissant).")
  }

  dt <- diff(t)
  tol <- 1e-8

  if (max(abs(dt - mean(dt))) <= tol * mean(dt)) {
    methode_appliquee <- "direct"
    pont <- "aucun"
    parametre_pont <- NA_real_
    dx <- diff(x)
    pas <- mean(dt)
    t_pont <- t
  } else {
    methode_appliquee <- methode
    if (processus == "HGP") {
      pont <- "gamma"
      parametre_pont <- c(a = as.numeric(estimer(t, x)["a"]))
      if (methode == "Bridge") {
        interp <- gamma_bridge(x, t, parametre_pont[["a"]])
      } else {
        interp <- simulate_gamma_grille(x, t, parametre_pont[["a"]], digits = digits, M_max = M_max)
      }
    } else {
      pont <- "brownien"
      parametre_pont <- c(sigma = as.numeric(estimer_wiener(t, x)["sigma"]))
      if (methode == "Bridge") {
        interp <- wiener_bridge(x, t, parametre_pont[["sigma"]])
      } else {
        interp <- simulate_wiener_grille(x, t, parametre_pont[["sigma"]], digits = digits, M_max = M_max)
      }
    }
    if (methode == "Bridge") {
      dx <- diff(interp$new_obs)
      pas <- interp$teq
      t_pont <- interp$t_reg
    } else {
      dx <- interp$dx
      pas <- interp$pas_fine
      t_pont <- interp$t_fine
    }
  }

  controle_dx(dx, positif = (processus == "HGP"))

  if (processus == "HGP") {
    tests <- c("Henze", "VGE", "KS", "BET")
  } else {
    tests <- c("Shapiro", "Lilliefors", "JarqueBera", "Anderson", "Cramer", "Agostino")
  }

  lignes <- lapply(tests, function(nom) {
    r <- tryCatch(
      switch(nom,
             "Henze"      = gof.gp.henze(dx, a_poids = a_poids, boot = boot, alpha = alpha),
             "VGE"        = gof.gp.vge(dx, boot = boot, alpha = alpha),
             "KS"         = gof.gp.ks(dx, boot = boot, alpha = alpha),
             "BET"        = gof.gp.betsch(dx, boot = boot, alpha = alpha),
             "Shapiro"    = gof.wiener.shapiro(dx, alpha = alpha),
             "Lilliefors" = gof.wiener.lilliefors(dx, alpha = alpha),
             "JarqueBera" = gof.wiener.jarquebera(dx, alpha = alpha),
             "Anderson"   = gof.wiener.anderson(dx, alpha = alpha),
             "Cramer"     = gof.wiener.cramer(dx, alpha = alpha),
             "Agostino"   = gof.wiener.agostino(dx, alpha = alpha)),
      error = function(e) e)
    if (inherits(r, "error")) {
      return(data.frame(test = nom, statistique = NA_real_, valeur_critique = NA_real_,
                        p_valeur = NA_real_, Decision = NA, remarque = conditionMessage(r)))
    }
    if (nom %in% c("Henze", "KS", "BET")) {
      data.frame(test = nom, statistique = as.numeric(r$T.value), valeur_critique = as.numeric(r$cv),
                 p_valeur = NA_real_, Decision = as.logical(r$Decision), remarque = "")
    } else {
      data.frame(test = nom, statistique = as.numeric(r$statistic), valeur_critique = NA_real_,
                 p_valeur = as.numeric(r$p.value), Decision = as.logical(r$Decision), remarque = "")
    }
  })
  resultats <- do.call(rbind, lignes)
  rownames(resultats) <- NULL

  if (any(is.na(resultats$Decision))) {
    echecs <- resultats$test[is.na(resultats$Decision)]
    warning("Test(s) non appliqu\u00e9(s) : ", paste(echecs, collapse = ", "), " (voir la colonne 'remarque').")
  }

  list(resultats      = resultats,
       processus      = processus,
       methode        = methode_appliquee,
       pont           = pont,
       parametre_pont = parametre_pont,
       pas            = pas,
       n_dx           = length(dx),
       t_pont         = t_pont,
       dx             = dx)
}
