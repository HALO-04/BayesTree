\name{pdpgbart}
\title{Partial Dependence Plots for PGBART}
\alias{pdpgbart}
\alias{pd2pgbart}
\description{
   Run \code{\link{pgbart_train}} at test observations constructed so that
   a  plot can be created
   displaying the effect of
   a single variable (\code{pdpgbart}) or pair of variables (\code{pd2pgbart}).
   Note that if \eqn{y} is a binary with \eqn{P(Y=1 | x) = F(f(x))}, \eqn{F} the standard
   normal cdf, then the plots are all on the \eqn{f} scale.
}
\usage{
pdpgbart(
   x.train, y.train, x.test=matrix(0.0,0,0),
   model = "./pgbart.model",
   if_trained = FALSE,train_result=NULL,
   xind=1:2, levs=NULL, levquants=c(.05,(1:9)/10,.95),
   pl=TRUE,  plquants=c(.05,.95),
   sigest=NA, sigdf=3, sigquant=.90,
   k=2.0,
   power=NA, base=.95,
   binaryOffset=0,
   ntree=200,
   ndpost=1000, nskip=100, usepg=TRUE, numparticles = 10,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=FALSE, numcut=100, printcutoffs=0,
   verbose=TRUE,...)

pd2pgbart = function (
   x.train, y.train,x.test=matrix(0.0,0,0),
   model = "./pgbart.model",
   if_trained = FALSE,train_result=NULL,
   xind=1:2,levs=NULL, levquants=c(.05,(1:9)/10,.95),
   pl=TRUE, plquants=c(.05,.95), 
   sigest=NA, sigdf=3, sigquant=.90,
   k=2.0,
   power=NA, base=.95,
   binaryOffset=0,
   ntree=200,
   ndpost=1000, nskip=100, usepg=TRUE, numparticles = 10,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=FALSE, numcut=100, printcutoffs=0,
   verbose=TRUE,...)
}
\arguments{
    \item{x.train}{
   Explanatory variables for training (in sample) data.\cr
   May be a matrix or a data frame,
   with (as usual) rows corresponding to observations and columns to variables.\cr
   If a variable is a factor in a data frame, it is replaced with dummies.
   Note that q dummies are created if q>2 and
   one dummy is created if q=2, where q is the number of levels of the factor.
   \code{makeind} is used to generate the dummies.
   \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x}
   which is a row of x.train.
   }

   \item{y.train}{
   Dependent variable for training (in sample) data.\cr
   If y is numeric a continous response model is fit (normal errors).\cr
   If y is a factor (or just has values 0 and 1) then a binary response model
   with a probit link is fit.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Should have same structure as x.train.\cr
   \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
   }
   \item{model}{
   The path to save model file. The file contains a lot of details about the trees and nodes, even information about every node.
   }
  \item{if_trained}{
    If had trained the data. if \code{TRUE}, it can plot directly without extra train time. While if \code{FALSE}, function will train data first.}

  \item{train_result}{
    \code{if_trained} is TRUE, \code{train_result} will be used to plot in direct.}

   \item{xind}{
     Integer vector indicating which variables are to be plotted.
     In \code{pdpgbart}, corresponds to the variables (columns of \code{x.train}) for which a plot is to be constructed.
     In \code{plotpdpgbart}, corresponds to the indices in list returned by \code{pdpgbart} for which plot is to be constructed.
     In \code{pd2pgbart}, the indicies of a pair of variables (columns of \code{x.train}) to plot.
   }

   \item{levs}{
     Gives the values of a variable at which the plot is to be constructed.
     Must be a list, where the \eqn{i}th component gives the values for the \eqn{i}th variable.
     In \code{pdpgbart}, it should have same length as \code{xind}. In \code{pd2pgbart}, it should have length 2.
     See also argument \code{levquants}.
   }

  \item{levquants}{
    If \code{levs} in \code{NULL}, the values of each variable used in the plot is
    set to the quantiles (in \code{x.train}) indicated by levquants.
    Must be a vector of numeric type.
  }

  \item{pl}{
    For \code{pdpgbart} and \code{pd2pgbart}, if \code{TRUE}, plot is subsequently made (by calling \code{plot.*}).
  }

  \item{plquants}{
    In the plots, beliefs about \eqn{f(x)} are indicated by plotting the
    posterior median and a lower and upper quantile.
    \code{plquants} is a double vector of length two giving the lower and upper quantiles.
  }

   \item{sigest}{
   The prior for the error variance  (\eqn{\sigma^2}{sigma\^2}) is inverted chi-squared
   (the standard conditionally conjugate prior).
   The prior is specified by choosing the degrees of freedom, a rough estimate of the
   corresponding standard deviation and a quantile to put this rough estimate at.
   If sigest=NA then the rough estimate will be the usual least squares estimator.
   Otherwise the supplied value will be used.
   Not used if y is binary.
   }

   \item{sigdf}{
   Degrees of freedom for error variance prior.
   Not used if y is binary.
   }

   \item{sigquant}{
   The quantile of the prior that the rough estimate (see sigest) is placed at.
   The closer the quantile is to 1,
   the more aggresive the fit will be as you are putting more prior weight
   on error standard deviations (\eqn{\sigma}{sigma}) less than the rough estimate.
   Not used if y is binary.
   }

   \item{k}{
   For numeric y,
   k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5.
   The response (y.train) is internally scaled to range from -.5 to .5.
   For binary y,
   k is the number of prior standard deviations \eqn{f(x)} is away from +/-3.
   In both cases, the bigger k is, the more conservative the fitting will be.
   }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

   \item{binaryOffset}{
   Used for binary \eqn{y}.\cr
   The model is \eqn{P(Y=1 | x) = F(f(x) + binaryOffset)}.\cr
   The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
   a probability other than .5.
   }

   \item{ntree}{
   The number of trees in the sum.
   }

   \item{ndpost}{
   The number of posterior draws after burn in, ndpost/keepevery will actually be returned.
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }
   
   \item{usepg}{
    Two sampling methods such as "pg" and "cpg".
    If true, sampling methods is "pg". Otherwise, sampling methods is "cgm". 
   }
   \item{numparticles}{
    The number of particles.Only when use "pg" sampling method, numparticles will work.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.\cr
   A \dQuote{draw} will consist of values of the error standard deviation (\eqn{\sigma}{sigma})
   and \eqn{f^*(x)}{f*(x)}
   at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
   the current draw of \eqn{f}.
   }

   \item{keeptrainfits}{
   If true the draws of \eqn{f(x)} for \eqn{x} = rows of x.train are returned.
   }

   \item{usequants}{
   Decision rules in the tree are of the form
   \eqn{x \le c}{x <= c}  vs. \eqn{x > c}
   for each variable corresponding to a column of x.train.
   usequants determines how the set of possible c is determined.
   If usequants is true, then the c are
   a subset of the values (xs[i]+xs[i+1])/2 where xs is
   unique sorted values obtained from
   the corresponding column of x.train.
   If usequants is false, the cutoffs are equally spaced across the
   range of values taken on by the corresponding
   column of x.train.
   }
   
   \item{numcut}{
   The number of possible values of c (see usequants).
   If a single number if given, this is used for all variables.
   Otherwise a vector with length equal to ncol(x.train) is required,
   where the \eqn{i^{th}}{i^th} element gives the number of c used for
   the \eqn{i^{th}}{i^th} variable in x.train.
   If usequants is false, numcut equally spaced cutoffs
   are used covering the range of values in the corresponding
   column of x.train.  If usequants is true, then  min(numcut, the number of unique values in the
   corresponding columns of x.train - 1) c values are used.
   }

   \item{printcutoffs}{
   The number of cutoff rules c to printed to screen before the MCMC is run.
   Give a single integer, the same value will be used for all variables.
   If 0, nothing is printed.
   }

   \item{verbose}{
   Logical, if FALSE supress printing.
   }


}
\details{
   We divide the predictor vector \eqn{x} into a subgroup of interest,
   \eqn{x_s} and the complement \eqn{x_c = x \setminus x_s}{x_c = x - x_s}.
   A prediction \eqn{f(x)} can
   then be written as \eqn{f(x_s, x_c)}. To estimate the effect of \eqn{x_s}
   on the prediction, Friedman suggests the partial dependence
   function
   \deqn{
   f_s(x_s) = \frac{1}{n}\sum_{i=1}^n f(x_s,x_{ic})
   }{
   f_s(x_s) = (1/n) \sum_{i=1}\^n f(x_s,x_{ic})
   }
   where \eqn{x_{ic}} is the \eqn{i}th observation of \eqn{x_c} in the data. Note
   that \eqn{(x_s, x_{ic})} will generally not be one of the observed data
   points. Using BART it is straightforward to then estimate and even
   obtain uncertainty bounds for \eqn{f_s(x_s)}. A draw of \eqn{f^*_s(x_s)}{f*_s(x_s)}
   from the induced BART posterior on \eqn{f_s(x_s)} is obtained by
   simply computing \eqn{f^*_s(x_s)}{f*_s(x_s)} as a byproduct of each MCMC draw
   \eqn{f^*}{f*}. The median (or average)
   of these MCMC draws \eqn{f^*_s(x_s)}{f*_s(x_s)} then yields an
   estimate of \eqn{f_s(x_s)}, and lower and upper quantiles can be used
   to obtain intervals for \eqn{f_s(x_s)}.

   In \code{pdpgbart} \eqn{x_s} consists of a single variable in \eqn{x} and in
   \code{pd2pgbart} it is a pair of variables.

   This is a computationally intensive procedure.
   For example, in \code{pdbart}, to compute the partial dependence plot
   for 5 \eqn{x_s} values, we need
   to compute \eqn{f(x_s, x_c)} for all possible \eqn{(x_s, x_{ic})} and there
   would be \eqn{5n} of these where \eqn{n} is the sample size.
   All of that computation would be done for each kept BART draw.
   For this reason running BART with \code{keepevery} larger than 1 (eg. 10)
   makes the procedure much faster.
}
\value{
  The plot methods produce the plots and don't return anything.

  \code{pdpgbart} and \code{pd2pgbart} return lists with components
  given below. The list returned by \code{pdpgbart} is assigned class
  \code{pdpgbart} and the list returned by \code{pd2pgbart} is assigned
  class \code{pd2pgbart}.

  \item{fd}{
    A matrix whose \eqn{(i, j)} value is the \eqn{i}th
    draw of \eqn{f_s(x_s)} for the \eqn{j}th value of \eqn{x_s}.
    \dQuote{fd} is for \dQuote{function draws}.

    For \code{pdpgbart} \code{fd} is actually a list whose
    \eqn{k}th component is the matrix described above
    corresponding to the \eqn{k}th variable chosen by argument xind.\
    The number of columns in each matrix will equal the number of values
    given in the corresponding component of argument \code{levs} (or number of values in \code{levquants}).

    For \code{pd2pgbart}, \code{fd} is a single matrix.
    The columns correspond to all possible pairs of values for the pair
    of variables indicated by \code{xind}.
    That is, all possible \eqn{(x_i, x_j)} where \eqn{x_i} is a value in
    the levs component corresponding to the first \eqn{x} and
    \eqn{x_j} is a value in the levs components corresponding to the second one.
    The first \eqn{x} changes first.
  }

  \item{levs}{
    The list of levels used, each component corresponding to a variable.
    If argument \code{levs} was supplied it is unchanged.
    Otherwise, the levels in \code{levs} are as constructed using argument \code{levquants}.
  }

  \item{xlbs}{
    A vector of character strings which are the plotting labels used for the variables.
  }
   
   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{sigma}{post burn in draws of sigma, length = ndpost/keepevery.}
   \item{first.sigma}{burn-in draws of sigma.}
   \item{varcount}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   \item{sigest}{
   The rough error standard deviation (\eqn{\sigma}{sigma}) used in the prior.
   }
 
   \item{y}{
   The input dependent vector of values for the dependent variable.\cr
   This is used in plot.bart.
   }

  The remaining components returned in the list are the same as in the value of \code{\link{pgbart_train}}.
  They are simply passed on from the BART run used to create the partial dependence plot.

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
\examples{
\dontrun{
## simulate data
f <- function(x) { return(0.5 * x[,1] + 2 * x[,2] * x[,3]) }

sigma <- 0.2
n     <- 100

set.seed(27)
x <- matrix(2 * runif(n * 3) -1, ncol = 3);
colnames(x) <- c('rob', 'hugh', 'ed')

Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

## first two plot regions are for pdbart, third for pd2bart
par(mfrow = c(1, 3))

## pdbart: one dimensional partial dependence plot
set.seed(99)
xind=1:ncol(x.train)
pdb1 <-
 pdpgbart(
    x, y,x.test=matrix(0.0,0,0),
    model = "./pgbart.model",
    if_trained = FALSE,train_result=NULL,
    xind=xind,levs=NULL, levquants=c(.05,(1:9)/10,.95),
    pl=TRUE, plquants=c(.05,.95), 
    sigest=NA, sigdf=3, sigquant=.90,
    k=2.0,
    power=NA, base=.95,
    binaryOffset=0,
    ntree=200,
    ndpost=100, nskip=100, usepg=TRUE, numparticles = 10,
    printevery=100, keepevery=1, keeptrainfits=TRUE,
    usequants=FALSE, numcut=100, printcutoffs=0,
    verbose=TRUE
  )

## pd2bart: two dimensional partial dependence plot
set.seed(99)
pdb2 <-
 pdpgbart(
    x, y,x.test=matrix(0.0,0,0),
    model = "./pgbart.model",
    if_trained = FALSE,train_result=NULL,
    xind=1:2,levs=NULL, levquants=c(.05,(1:9)/10,.95),
    pl=TRUE, plquants=c(.05,.95), 
    sigest=NA, sigdf=3, sigquant=.90,
    k=2.0,
    power=NA, base=.95,
    binaryOffset=0,
    ntree=200,
    ndpost=100, nskip=100, usepg=TRUE, numparticles = 10,
    printevery=100, keepevery=1, keeptrainfits=TRUE,
    usequants=FALSE, numcut=100, printcutoffs=0,
    verbose=TRUE
  )
plot(pdb2)

## compare BART fit to linear model and truth = Ey
lmFit  <- lm(y ~., data.frame(x, y))
fitmat <- cbind(y, Ey, lmFit$fitted, pdb1$yhat.train.mean)
colnames(fitmat) <- c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))
}}
