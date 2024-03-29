% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{match_gset_nms}
\alias{match_gset_nms}
\title{Matches gset_nms to prnGSPA.}
\usage{
match_gset_nms(gset_nms = NULL)
}
\arguments{
\item{gset_nms}{Character string or vector containing the shorthanded name(s),
full file path(s), or both, to gene sets for enrichment analysis. For
species among \code{"human", "mouse", "rat"}, the default of
\code{c("go_sets", "c2_msig", "kinsub")} will utilize terms from gene
ontology (\code{GO}), molecular signatures (\code{MSig}) and
kinase-substrate network (\code{PSP Kinase-Substrate}). Custom \code{GO},
\code{MSig} and other data bases at given species are also supported. See
also: \code{\link{prepGO}} for the preparation of custom \code{GO};
\code{\link{prepMSig}} for the preparation of custom \code{MSig}. For other
custom data bases, follow the same format of list as \code{GO} or
\code{MSig}.}
}
\description{
Not currently used.
}
