% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psmtable.R
\name{cleanPSMRT_by_sets}
\alias{cleanPSMRT_by_sets}
\title{Clean up PSM retention times for a set of TMTset[i]LCMSinj[j]}
\usage{
cleanPSMRT_by_sets(filelist, lfq_mbr = TRUE, max_n_apexes = 2L)
}
\arguments{
\item{filelist}{A file names of PSM tables TMTSet1_LCMSinj1_PSM_N.txt etc.}

\item{lfq_mbr}{Logical; perform LFQ MBR nor not.}

\item{max_n_apexes}{The maximum number of apexes for consideration.}
}
\description{
No clean-ups across TMTset[i]LCMSinj[j].
}
