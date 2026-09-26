\name{estimation}
\alias{estimer}
\alias{estimer_wiener}
\title{Estimation des paramètres sous H0}
\description{
\code{estimer} : paramètres \eqn{(a, b)} du processus gamma homogène par la
méthode des moments. \code{estimer_wiener} : paramètres \eqn{(\mu, \sigma)} du
processus de Wiener homogène par maximum de vraisemblance.
}
\usage{
estimer(t, x)
estimer_wiener(t, x)
}
\arguments{
  \item{t}{instants d'inspection, strictement croissants (au moins 3).}
  \item{x}{niveaux observés (croissants pour \code{estimer}).}
}
\details{
Avec \eqn{N} accroissements :
\eqn{\hat\mu = \sum_i \Delta x_i / \sum_i \Delta t_i} et
\eqn{\hat\nu = \frac{1}{N}\sum_i (\Delta x_i - \hat\mu \Delta t_i)^2 / \Delta t_i}.
\code{estimer} renvoie \eqn{\hat a = \hat\mu^2/\hat\nu} et
\eqn{\hat b = \hat\mu/\hat\nu} ; \code{estimer_wiener} renvoie
\eqn{\hat\mu}, \eqn{\hat\sigma^2 = \hat\nu} et \eqn{\hat\sigma}.
}
\value{
\code{estimer} : vecteur nommé \code{c(a, b)}.
\code{estimer_wiener} : vecteur nommé \code{c(mu, sigma2, sigma)}.
}
\examples{
set.seed(1)
X <- sim_HGP(0:50, a = 2, b = 1)
estimer(X$t, X$X)
W <- sim_HWP(0:50, mu = 2, sigma = 1.5)
estimer_wiener(W$t, W$X)
}
