\name{simulation}
\alias{sim_HGP}
\alias{sim_NHGP}
\alias{sim_IGP}
\alias{sim_HWP}
\alias{sim_NHWP}
\title{Simulation de trajectoires de dégradation}
\description{
Simulent une trajectoire aux instants \code{t_obs} : processus gamma homogène
ou non homogène, processus gaussien inverse homogène, processus de Wiener
homogène ou non homogène.
}
\usage{
sim_HGP(t_obs, a, b)
sim_NHGP(t_obs, a, b, alpha)
sim_IGP(t_obs, mu, lambda)
sim_HWP(t_obs, mu, sigma)
sim_NHWP(t_obs, mu, sigma, Lambda, tau = function(t) t)
}
\arguments{
  \item{t_obs}{instants d'inspection, strictement croissants et positifs ou
    nuls ; l'instant 0 est ajouté (avec un avertissement) s'il manque.}
  \item{a, b}{forme par unité de temps et taux : accroissements
    \eqn{\Gamma(a \Delta t, b)} (\code{sim_HGP}) ou
    \eqn{\Gamma(a (t_{i+1}^{\alpha} - t_i^{\alpha}), b)} (\code{sim_NHGP}).}
  \item{alpha}{exposant de la fonction de forme \eqn{A(t) = a t^{\alpha}}.}
  \item{mu}{dérive (réel quelconque pour les processus de Wiener, strictement
    positif pour le processus gaussien inverse).}
  \item{lambda}{paramètre de forme du processus gaussien inverse :
    accroissements \eqn{IG(\mu \Delta t, \lambda \Delta t^2)}.}
  \item{sigma}{coefficient de diffusion, strictement positif.}
  \item{Lambda, tau}{fonctions déterministes nulles en 0, \code{tau}
    croissante : \eqn{X(t) = \mu \Lambda(t) + \sigma B(\tau(t))}.}
}
\value{data.frame avec les colonnes \code{t} et \code{X}, commençant en (0, 0).}
\examples{
set.seed(1)
t_obs <- 0:20
sim_HGP(t_obs, a = 2, b = 1)
sim_HWP(t_obs, mu = 2, sigma = 1.5)
sim_NHWP(t_obs, mu = 1, sigma = 1.5, Lambda = function(t) t^2)
}
