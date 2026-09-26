\name{arrondir_temps}
\alias{arrondir_temps}
\title{Arrondi des instants d'inspection}
\description{
Arrondit les instants d'inspection à \code{digits} décimales, ajoute l'instant
0 et supprime les doublons (avec un avertissement). À appliquer avant la
simulation ou le recueil des données lorsque la méthode \code{"ExtBridge"} de
\code{\link{Bridge}} est utilisée.
}
\usage{arrondir_temps(t, digits = 0)}
\arguments{
  \item{t}{instants d'inspection, positifs ou nuls.}
  \item{digits}{nombre de décimales, entier positif ou nul.}
}
\value{Vecteur strictement croissant d'instants, commençant à 0.}
\examples{
arrondir_temps(c(0.4, 2.26, 3.71, 3.74), digits = 1)
}
