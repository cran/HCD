\name{HCD}
\alias{HCD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
hierarchical community detection with recursive spectral methods
}
\description{
Hierarchical community by recursive spectral partitioning. It includes
the splitting methods of spectral clustering and sign splitting, as well
stopping rules for fixed stopping, non-backtracking matrix checking and
edge cross-validation.
}
\usage{
HCD(A, method = "SS", stopping = "NB", reg = FALSE, n.min = 25, D = NULL,notree=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
adjacency matrix. Can be standard R matrix or dsCMatrix (or other type in package Matrix)
}
  \item{method}{
splitting method. "SS" (default) for sign splitting, "SC" for spectral clustering
}
  \item{stopping}{
stopping rule. "NB" (default) for non-backtracking matrix spectrum,
"ECV" for edge cross-validation, "Fix"for fixed D layers of partitioning
(needs D value)
}
  \item{reg}{
logic value on whether regularization is needed. By default it is FALSE.Set it to be TRUE will add reguarlization, which help the performance on sparse networks, but it will make the computation slower.
}
  \item{n.min}{
integer number. The algorithm will stop splitting if the current size is <= 2*n.min.
}
  \item{D}{
the number of layers to partition, if stopping=="Fix".
}
  \item{notree}{
logical value on whether the tree and the corresponding similarity will
be computed. If TRUE (default), will not produce the data.tree object or the community similarity matrix. Only the cluster label and the tree path strings will be returned. This typically makes the runing faster.
}
}
\details{
For stopping rules, ECV is nonparametric rank evaluation by
cross-validation, a more generally applicable approach without assuming
SBM or its variants. ECV is also applicable for weighted networks.So it is
believed to be more robust than NB but less effective if the true model
is close to BTSBM. However, the ECV is computationally much more
intensive.

Notice that the algorithm does not reply on the assumption of the BTSBM.
But
the estimated probabiilty matrix from the output is based on the
BTSBM.
}
\value{
  A list of the following objects:
  \item{labels}{detected community labels of nodes}
  \item{ncl}{number of clusters from the algorithm}
  \item{cluster.tree}{a data.tree object for the binary tree between
    communities}
  \item{P}{estimated connection probability matrix between n nodes,
    according to BTSBM}
  \item{node.bin.sim.mat}{binary string similarity between nodes}
  \item{comm.bin.sim.mat}{binary string similarity between communities}
  \item{tree.path}{a list of strings to describe the path from root to
    each community along the tree}
}
\references{
Tianxi Li, Lihua Lei, Sharmodeep Bhattacharyya, Purnamrita Sarkar, Peter
Bickel, and Elizaveta Levina.  Hierarchical community detection by
recursive partitioning. arXiv:1810.01509
}
\author{
  Tianxi Li, Lihua Lei, Sharmodeep Bhattacharyya, Purnamrita Sarkar, Peter Bickel, and Elizaveta Levina.\cr

Maintainer: Tianxi Li <tianxili@virginia.edu>
}
\examples{
dt <- BTSBM(n=1600,d=4,a.seq=0.2^seq(0,4),lambda=50)
A <- dt$A.list[[1]]
# you can try various versions of the algorithm as below: the Fix is fastest and ECV is slowest.
system.time(HCD.result <- HCD(A,method="SC",stopping="Fix",D=4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ HCD}% use one of  RShowDoc("KEYWORDS")
\keyword{ BTSBM }% __ONLY ONE__ keyword per line
