\name{BTSBM}
\alias{BTSBM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates networks from binary tree stochastic block model
}
\description{
Generates networks from binary tree stochastic block model, with
provided sequence of connection probability along the tree
}
\usage{
BTSBM(n, d, a.seq, lambda, alpha = NULL, N = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
 number of nodes in the network
}
  \item{d}{
number of layers until leaves (excluding the root)
}
  \item{a.seq}{
the connection probability sequence along the tree, a_r, see details in the paper
}
  \item{lambda}{
average node degree, only used when alpha is not provided
}
\item{alpha}{
  the common scaling of the a_r sequence. So at the end, essentially the a_r sequence is a.seq*alpha
}
  \item{N}{
the number of networks to generate from the same model
}
}

\value{
  A list of objections of
  \item{A.list}{the generated network adjacency matrices}
  \item{B}{the connection probability matrix between K communities,
    where K = 2^d}
  \item{label}{the vector of community labels for n nodes}
  \item{P}{the connection probability matrix between the n nodes. It is
    the expectation of adjacency matrices, except on the diagonal}
  \item{comm.sim.mat}{the binary string similarity matrix between
    communities}
  \item{node.sim.mat}{the binary string similarity matrix between nodes}
  }
\references{
Li, T., Lei, L., Bhattacharyya, S., Van den Berge, K., Sarkar, P., Bickel, P.J. and Levina, E., 2022. Hierarchical community detection by recursive partitioning. Journal of the American Statistical Association, 117(538), pp.951-968.
}
\author{
  Tianxi Li, Lihua Lei, Sharmodeep Bhattacharyya, Purnamrita Sarkar, Peter Bickel, and Elizaveta Levina.\cr

Maintainer: Tianxi Li <tianxili@umn.edu>

}
\examples{
dt <- BTSBM(n=1600,d=4,a.seq=0.2^seq(0,4),lambda=50)
A <- dt$A.list[[1]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BTSBM }