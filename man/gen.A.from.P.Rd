\name{gen.A.from.P}
\alias{gen.A.from.P}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
generates a network from the given connection probability
}
\description{
  Generates an adjacency matrix from a given probability matrix,
  according independent Bernoulli -- the so-called inhomogeneous
  Erdos-Renyi model. It is used to generate new networks from a given model.
}
\usage{
gen.A.from.P(P, undirected = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
connection probability between nodes
}
  \item{undirected}{
logic value. FALSE (default) if the network is undirected, so the
adjacency matrix will be symmetric with only upper diagonal entries
being generated as independent Bernoulli.
}
}
\value{
An adjacency matrix
%% ...
}
\references{
Li, T., Lei, L., Bhattacharyya, S., Van den Berge, K., Sarkar, P., Bickel, P.J. and Levina, E., 2022. Hierarchical community detection by recursive partitioning. Journal of the American Statistical Association, 117(538), pp.951-968.
}
\author{
  Tianxi Li, Lihua Lei, Sharmodeep Bhattacharyya, Purnamrita Sarkar, Peter Bickel, and Elizaveta Levina.\cr

Maintainer: Tianxi Li <tianxili@umn.edu>}

