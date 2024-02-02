\name{HCDplot}
\alias{HCDplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot the result of hierarchical community detection
}
\description{
Generate dendrogram of the HCD result. 
}
\usage{
HCDplot(hcd,mode="community",labels=NULL,main=NULL,label.cex=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hcd}{
The result of an HCD call.
}
  \item{mode}{
plotting community hierarchy or node hierarchy. The default value is "community", indicating plotting hierarchy between communities. Alternatively, the plot is for all nodes, which is not recommended because usually there are too many of them.
}
  \item{labels}{
the labels of the each leaf of the tree. By default, the community/node index is used. The user can also specify another sequence of characters.
}
  \item{main}{
title of the plot.
}
  \item{label.cex}{
size of the leaf label in the plot. When plotting node hierarchy, typically there are too many nodes so the labels will seriously overlap. Use a smaller size (say, label.cex=0.3) may help.
}
}
\value{No return value, called for visualization.}
\references{
Li, T., Lei, L., Bhattacharyya, S., Van den Berge, K., Sarkar, P., Bickel, P.J. and Levina, E., 2022. Hierarchical community detection by recursive partitioning. Journal of the American Statistical Association, 117(538), pp.951-968.
}
\author{
  Tianxi Li, Lihua Lei, Sharmodeep Bhattacharyya, Purnamrita Sarkar, Peter Bickel, and Elizaveta Levina.\cr

Maintainer: Tianxi Li <tianxili@umn.edu>
}
\examples{
dt <- BTSBM(n=80,d=4,a.seq=0.2^seq(0,4),lambda=20)
A <- dt$A.list[[1]]
system.time(HCD.result <- HCD(A,method="SC",stopping="Fix",D=4,notree=FALSE,n.min=5))

HCDplot(HCD.result,mode="community",main="Community Tree")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ HCD}% use one of  RShowDoc("KEYWORDS")
