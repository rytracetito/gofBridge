\name{gof.process}
\alias{gof.process}
\title{Tests d'adéquation au processus gamma homogène ou au processus de Wiener homogène}
\description{
Teste l'adéquation d'une trajectoire observée au processus gamma homogène (HGP)
ou au processus de Wiener homogène (HWP). Le plan d'inspection est détecté
automatiquement : s'il est régulier, les tests portent sur les accroissements
observés ; sinon, la trajectoire est d'abord reconstruite par pont gamma (HGP)
ou pont brownien (HWP), le paramètre du pont étant estimé sous H0.
}
\usage{
gof.process(X, processus, methode = "Bridge", digits = 0, M_max = 100,
            alpha = 0.05, boot = 199, a_poids = 1)
}
\arguments{
  \item{X}{data.frame ou matrice avec des colonnes nommées \code{t} (instants
    d'inspection, strictement croissants, positifs ou nuls) et \code{X}
    (niveaux de dégradation observés). Si \code{t[1] > 0}, le point (0, 0)
    est ajouté avec un avertissement ; si \code{t[1] = 0}, \code{X[1]} doit valoir 0.}
  \item{processus}{\code{"HGP"} (processus gamma homogène) ou \code{"HWP"}
    (processus de Wiener homogène). Obligatoire.}
  \item{methode}{méthode de reconstruction en plan irrégulier :
    \code{"Bridge"} (grille régulière de pas \eqn{t_N/N}, Grall-Maës, 2012) ou
    \code{"ExtBridge"} (grille fine contenant exactement toutes les
    observations). Ignoré en plan régulier.}
  \item{digits}{\code{"ExtBridge"} uniquement : nombre de décimales des
    instants ; les instants doivent être des multiples de \eqn{10^{-digits}}
    (voir \code{\link{arrondir_temps}}).}
  \item{M_max}{\code{"ExtBridge"} uniquement : nombre maximal de pas de la
    grille fine.}
  \item{alpha}{niveau des tests, dans \eqn{]0, 1[}.}
  \item{boot}{HGP uniquement : nombre de rééchantillonnages du bootstrap
    paramétrique, avec \code{(boot + 1) * alpha >= 1}.}
  \item{a_poids}{HGP uniquement : paramètre de pondération du test de
    Henze-Meintanis-Ebner.}
}
\details{
Le plan est considéré comme régulier si
\eqn{\max_i |\Delta t_i - \bar{\Delta t}| \le 10^{-8} \bar{\Delta t}}.

En plan irrégulier, le paramètre du pont est estimé sous H0 sur les données
observées : \eqn{\hat a = \hat\mu^2/\hat\nu} pour HGP (\code{\link{estimer}}),
\eqn{\hat\sigma} pour HWP (\code{\link{estimer_wiener}}), avec
\eqn{\hat\mu = x_N / t_N} et
\eqn{\hat\nu = \frac{1}{N}\sum_i (\Delta x_i - \hat\mu \Delta t_i)^2 / \Delta t_i}.
La trajectoire est ensuite reconstruite par pont gamma (lois bêta) ou pont
brownien (lois normales conditionnelles), sur la grille de la méthode choisie.

Tous les tests du processus sont appliqués aux mêmes accroissements :
\itemize{
  \item HGP : Henze-Meintanis-Ebner (\code{"Henze"}),
    Villaseñor-González-Estrada (\code{"VGE"}), Kolmogorov-Smirnov
    (\code{"KS"}), Betsch-Ebner (\code{"BET"}), avec paramètres estimés et
    bootstrap paramétrique. Henze, KS et BET renvoient une valeur critique
    bootstrap ; VGE renvoie la p-valeur \eqn{\frac{1}{B}\sum_b 1\{T^*_b \ge T\}}.
  \item HWP : Shapiro-Wilk, Lilliefors, Jarque-Bera, Anderson-Darling,
    Cramér-von Mises, D'Agostino-Pearson (au moins 20 accroissements).
}
Un test non applicable (par exemple D'Agostino-Pearson avec moins de 20
accroissements) renvoie \code{NA}, avec son motif dans la colonne
\code{remarque}, et un avertissement est émis ; les autres tests sont appliqués.

Avec \code{processus = "HGP"}, une trajectoire présentant un accroissement
observé négatif provoque une erreur (incompatible avec un processus gamma).
}
\value{
Une liste :
\item{resultats}{data.frame, une ligne par test : \code{test},
  \code{statistique}, \code{valeur_critique}, \code{p_valeur},
  \code{Decision} (\code{TRUE} : rejet de H0 au niveau \code{alpha}),
  \code{remarque}.}
\item{processus}{processus testé.}
\item{methode}{méthode appliquée : \code{"direct"} (plan régulier),
  \code{"Bridge"} ou \code{"ExtBridge"}.}
\item{pont}{\code{"aucun"}, \code{"gamma"} ou \code{"brownien"}.}
\item{parametre_pont}{paramètre du pont estimé sous H0 (\code{a} ou
  \code{sigma}), \code{NA} en plan régulier.}
\item{pas}{pas de la grille des accroissements testés.}
\item{n_dx}{nombre d'accroissements testés.}
\item{t_pont}{instants de la grille (instants observés en plan régulier).}
\item{dx}{accroissements testés.}
}
\references{
Grall-Maës, E. (2012). Use of the Kolmogorov-Smirnov test for gamma process.
\emph{Proceedings of the Institution of Mechanical Engineers, Part O: Journal
of Risk and Reliability}, 226(6), 624-634. \doi{10.1177/1748006X12462522}
}
\seealso{\code{\link{sim_HGP}}, \code{\link{arrondir_temps}}, \code{\link{estimer}}}
\examples{
set.seed(1)
t_obs <- arrondir_temps(runif(40, 0, 80))

X <- sim_HWP(t_obs, mu = 2, sigma = 1.5)
gof.process(X, "HWP")$resultats
gof.process(X, "HWP", methode = "ExtBridge")$resultats

Y <- sim_HGP(t_obs, a = 2, b = 1)
r <- gof.process(Y, "HGP", boot = 19)
r$resultats
r$methode
r$parametre_pont
}
