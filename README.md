# gofBridge

Package R de tests d'adéquation au **processus gamma homogène** (HGP) et au **processus de Wiener homogène** (HWP), pour des trajectoires de dégradation observées à des instants réguliers ou irréguliers.

- Plan régulier : tests sur les accroissements observés.
- Plan irrégulier : paramètre du pont estimé sous H0 (a pour HGP, sigma pour HWP), reconstruction de la trajectoire par pont gamma ou pont brownien, puis tests sur les accroissements reconstruits. Deux méthodes de reconstruction :
  - `"Bridge"` : grille régulière de pas t_N / N (Grall-Maës, 2012) ;
  - `"ExtBridge"` : grille fine contenant exactement toutes les observations.

Tests appliqués (tous sur les mêmes accroissements) :

| Processus | Tests |
|---|---|
| HGP | Henze-Meintanis-Ebner, Villaseñor-González-Estrada, Kolmogorov-Smirnov, Betsch-Ebner (paramètres estimés, bootstrap paramétrique) |
| HWP | Shapiro-Wilk, Lilliefors, Jarque-Bera, Anderson-Darling, Cramér-von Mises, D'Agostino-Pearson |

## Installation

Les dépendances (`gofgamma`, `goft`, `nortest`, `tseries`, `fBasics`, `statmod`, toutes sur CRAN) sont déclarées obligatoires dans `DESCRIPTION` : les deux méthodes ci-dessous les installent automatiquement.

```r
install.packages("remotes")
```

Depuis GitHub (dépôt privé : jeton d'accès personnel en lecture sur `gofBridge`, enregistré une fois avec `gitcreds::gitcreds_set()`) :

```r
remotes::install_github("rytracetito/gofBridge")
```

Depuis le fichier source :

```r
remotes::install_local("chemin/vers/gofBridge_0.1.0.tar.gz")
```

`install.packages("gofBridge_0.1.0.tar.gz", repos = NULL, type = "source")` n'installe pas les dépendances : l'utiliser seulement si elles sont déjà installées.

## Utilisation

```r
library(gofBridge)

set.seed(1)
t_obs <- arrondir_temps(runif(40, 0, 80))

X <- sim_HWP(t_obs, mu = 2, sigma = 1.5)
r <- Bridge(X, processus = "HWP", methode = "Bridge")
r$resultats        # une ligne par test : statistique, valeur critique ou p-valeur, décision
r$methode          # "direct" (plan régulier), "Bridge" ou "ExtBridge"
r$parametre_pont   # sigma estimé sous H0

Y <- sim_HGP(t_obs, a = 2, b = 1)
Bridge(Y, processus = "HGP", methode = "ExtBridge")$resultats
```

`Decision = TRUE` signifie le rejet de H0 au niveau `alpha`. Documentation : `?Bridge`, `?sim_HGP`, `?arrondir_temps`, `?estimer`.

## Fonctions exportées

| Fonction | Rôle |
|---|---|
| `Bridge` | Détection du plan, pont, tests d'adéquation |
| `sim_HGP`, `sim_NHGP`, `sim_IGP` | Simulation : gamma homogène, gamma non homogène, gaussien inverse |
| `sim_HWP`, `sim_NHWP` | Simulation : Wiener homogène, Wiener non homogène X(t) = mu Lambda(t) + sigma B(tau(t)) |
| `arrondir_temps` | Arrondi des instants d'inspection (méthode `"ExtBridge"`) |
| `estimer`, `estimer_wiener` | Estimation des paramètres sous H0 |

## Référence

E. Grall-Maës (2012), *Use of the Kolmogorov–Smirnov test for gamma process*, Proc. IMechE Part O: Journal of Risk and Reliability, 226(6), 624–634. DOI : [10.1177/1748006X12462522](https://doi.org/10.1177/1748006X12462522)
