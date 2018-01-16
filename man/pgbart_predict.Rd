\name{pgbart_predict}
\title{Test of particle Gibbs for Bayesian Additive Regression Trees Sampler}
\alias{pgbart_predict}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model.\cr
For numeric response \eqn{y}, we have
\eqn{y = f(x) + \epsilon}{y = f(x) + e},
where \eqn{\epsilon \sim N(0,\sigma^2)}{e ~ N(0,sigma\^2)}.\cr
For a binary response \eqn{y}, \eqn{P(Y=1 | x) = F(f(x))}, where \eqn{F}
denotes the standard normal cdf (probit link).

In both cases, \eqn{f} is the sum of many tree models.
The goal is to have very flexible inference for the uknown
function \eqn{f}.

In the spirit of \dQuote{ensemble models},
each tree is constrained by a prior to be a weak learner
so that it contributes a
small amount to the overall fit.
}
\usage{
pgbart_predict (x.test, model = "./pgbart.model")}

\arguments{
\item{x.test}{Explanatory variables for test (out of sample) data.\cr
   Should have same structure as x.train.\cr
   \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of \code{train_data}.}

\item{model}{The path to save model file. The file contains a lot of details about the trees and nodes, even information about every node. }
}
\details{
   BART is an Bayesian MCMC method.
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case
   and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
   The \code{plot} method sets mfrow to c(1,2) and makes two plots.\cr
   The first plot is the sequence of kept draws of \eqn{\sigma}{sigma}
   including the burn-in draws.  Initially these draws will decline as BART finds fit
   and then level off when the MCMC has burnt in.\cr
   The second plot has \eqn{y} on the horizontal axis and posterior intervals for
   the corresponding \eqn{f(x)} on the vertical axis.
   \code{bart} returns a list assigned class \sQuote{bart}.
   In the numeric \eqn{y} case, the list has components:

   \item{yhat.test}{
   A matrix with (ndpost/keepevery) rows and nrow(x.test) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.test.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.test.\cr
   Burn-in is dropped.
   }
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns. Only exists when \eqn{y} is not binary}

   In the binary \eqn{y} case, the returned list has the components
   yhat.test, and binaryOffset.

   Note that in the binary \eqn{y} case yhat.test are
   \eqn{f(x)} + binaryOffset.  If you want draws of the probability
   \eqn{P(Y=1 | x)} you need to apply the normal cdf (\code{pnorm})
   to these values.
}
\references{
Chipman, H., George, E., and McCulloch R. (2010)
   Bayesian Additive Regression Trees.
   \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298.

Lakshminarayanan B, Roy D, Teh Y W. Particle Gibbs for Bayesian additive regression trees[C]//Artificial Intelligence and Statistics. 2015: 553-561.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning.
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.
}

\author{
MCPRL: \email{xxx@bupt.edu.cn}.\cr
Hugh Chipman: \email{hugh.chipman@gmail.com}\cr
Robert McCulloch: \email{robert.e.mcculloch@gmail.com}\cr

}
\seealso{
\code{\link{pdpgbart}}
}
\examples{
##simulate data (example from Friedman MARS paper)
set.seed(99)
x=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
##run inference process
predictFit = pgbart_predict(x, model='./pgbart.model')
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
