#' Helper: switches among ion types for calculating MS2 masses.
#'
#' @param aas2 A sequence of amino-acid residues with \emph{masses}. Residues
#'   are in names and masses in values (note that argument \code{aas}
#'   corresponds to residues without masses).
#' @inheritParams calcpep_a1_t1_nl1
#' @inheritParams calc_ms2ions
ms2ions_by_type <- function (aas2, ntmass, ctmass, type_ms2ions, digits) {
  
  switch(type_ms2ions, 
         by = byions(ntmass, ctmass, aas2, digits), 
         cz = czions(ntmod, ctmod, aas2, digits), 
         ax = axions(ntmod, ctmod, aas2, digits), 
         stop("Unknown type.", call. = FALSE))
}



#' Masses of singly-charged b- and y-ions.
#'
#' @rdname bions_base
#' @seealso \link{add_complement_ions}
byions <- function (ntmass, ctmass, aas2, digits = 4L) {
  bs <- bions_base(aas2, ntmass, digits)
  ys <- yions_base(aas2, ctmass, digits)
  c(bs, ys)
}


#' Masses of singly-charged c- and z-ions.
#'
#' @rdname bions_base
czions <- function (ntmass, ctmass, aas2, digits = 4L) {
  cs <- cions_base(aas2, ntmass, digits)
  zs <- zions_base(aas2, ctmass, digits)
  c(cs, zs)
}


#' Masses of singly-charged a- and x-ions.
#'
#' @rdname bions_base
axions <- function (ntmass, ctmass, aas2, digits = 4L) {
  as <- aions_base(aas2, ntmass, digits)
  xs <- xions_base(aas2, ctmass, digits)
  c(as, xs)
}


#' B-ions.
#'
#' For (1) "amods- tmod- vnl- fnl-", (2) "amods- tmod+ vnl- fnl-".
#'
#' @param aas2 A sequence of amino-acid residues with \emph{masses}. Residues
#'   are in names and masses in values.
#'
#'   The masses reflects fixed/variable modifications, and/or fixed/variable
#'   neutral losses.
#'   
#' @param digits Integer; the number of decimal places to be used.
#' @param tmass The mass of a fixed or variable N-term or C-term modification.
#'
#' @importFrom stringr str_split
#' @examples
#' \donttest{
#' ## (1) "amods- tmod- vnl- fnl-"
#' # (Fixed N-term mods; also for no N-term mod)
#'
#' fixedmods = c("TMT6plex (N-term)", "TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods = c("Oxidation (M)", "Deamidated (N)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods, add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' aa_masses = aa_masses_all[[1]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549 # - electron
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#'
#' if (is_empty(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147 # + (H) + (H+)
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MAKEMASSPECFUN"
#' aas <- str_split(aa_seq, "", simplify = TRUE)
#' aas2 <- aa_masses[aas]
#'
#' b <- bions_base(aas2, ntmass)
#' y <- yions_base(aas2, ctmass)
#'
#'
#' ## (2) "amods- tmod+ vnl- fnl-"
#' # (2a, N-term)
#' fixedmods = c("TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods = c("TMT6plex (N-term)", "Acetyl (Protein N-term)", "Oxidation (M)",
#'             "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods, add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' aa_masses = aa_masses_all[[3]]
#'
#' # (Fixed or variable C-term mods +/- makes no difference on b-ions;
#' # and vice versa for y-ions)
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#'
#' if (is_empty(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' aa_seq <- "MAKEMASSPECFUN"
#' aas <- str_split(aa_seq, "", simplify = TRUE)
#' aas2 <- aa_masses[aas]
#'
#' b <- bions_base(aas2, ntmass)
#' y <- yions_base(aas2, ctmass)
#'
#'
#' # (2b, C-term)
#' fixedmods = c("TMT6plex (K)", "Carbamidomethyl (C)")
#' varmods = c("TMT6plex (N-term)", "Amidated (Protein C-term)", "Oxidation (M)",
#'             "Deamidated (N)", "Gln->pyro-Glu (N-term = Q)")
#'
#' mod_indexes <- seq_along(c(fixedmods, varmods)) %>%
#'   as.hexmode() %>%
#'   `names<-`(c(fixedmods, varmods))
#'
#' aa_masses_all <- calc_aamasses(fixedmods, varmods, add_varmasses = FALSE,
#'                                add_nlmasses = FALSE)
#'
#' # `TMT6plex (N-term)`; `Amidated (Protein C-term)`
#' aa_masses = aa_masses_all[[7]]
#'
#' ntmod <- attr(aa_masses, "ntmod", exact = TRUE)
#' ctmod <- attr(aa_masses, "ctmod", exact = TRUE)
#'
#' if (is_empty(ntmod)) {
#'   ntmass <- aa_masses["N-term"] - 0.000549
#' } else {
#'   ntmass <- aa_masses[names(ntmod)] - 0.000549
#' }
#'
#' if (is_empty(ctmod)) {
#'   ctmass <- aa_masses["C-term"] + 2.01510147
#' } else {
#'   ctmass <- aa_masses[names(ctmod)] + 2.01510147
#' }
#'
#' b <- bions_base(aas2, ntmass)
#' y <- yions_base(aas2, ctmass)
#'
#' }
bions_base <- function (aas2, tmass, digits = 4L) {
  c(tmass, aas2) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Y-ions.
#' 
#' @rdname bions_base
yions_base <- function (aas2, tmass, digits = 4L) {
  
  # (1) OH (C-term), + H (neutralizes the N-term on a fragment) + H+
  # (2) Other C-term (other than OH) + H + H+: X + 1.007825 + 1.00727647
  
  ions <- c(tmass, rev(aas2)) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' B2-ions.
#' 
#' @rdname bions_base
b2ions_base <- function (aas2, tmass, digits = 4L) {
  ions <- bions_base(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' B*-ions.
#' 
#' @rdname bions_base
bstarions <- function (aas2, tmass, digits = 4L) {
  # -NH3:17.026549
  c(tmass - 17.026549, aas2) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' B*2-ions.
#' 
#' @rdname bions_base
bstar2ions <- function (aas2, tmass, digits = 4L) {
  ions <- bstarions(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' B0-ions.
#' 
#' @rdname bions_base
b0ions <- function (aas2, tmass, digits = 4L) {
  # -H2O 18.010565
  c(tmass - 18.010565, aas2) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' B02-ions.
#' 
#' @rdname bions_base
b02ions <- function (aas2, tmass, digits = 4L) {
  ions <- b0ions(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' Y2-ions.
#' 
#' @rdname bions_base
y2ions <- function (aas2, tmass, digits = 4L) {
  ions <- yions_base(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' Y*-ions.
#' 
#' @rdname bions_base
ystarions <- function (aas2, tmass, digits = 4L) {
  ions <- c(tmass - 17.026549, rev(aas2)) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Y*2-ions.
#' 
#' @rdname bions_base
ystar2ions <- function (aas2, tmass, digits = 4L) {
  ions <- ystarions(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' Y0-ions.
#' 
#' @rdname bions_base
y0ions <- function (aas2, tmass, digits = 4L) {
  ions <- c(tmass - 18.010565, rev(aas2)) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Y02-ions.
#' 
#' @rdname bions_base
y02ions <- function (aas2, tmass, digits = 4L) {
  ions <- y0ions(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' C-ions.
#' 
#' @rdname bions_base
cions_base <- function (aas2, tmass, digits = 4L) {
  c(tmass + 17.026549, aas2) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' C2-ions.
#' 
#' @rdname bions_base
c2ions <- function (aas2, tmass, digits = 4L) {
  ions <- cions_base(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' Z-ions.
#' 
#' @rdname bions_base
zions_base <- function (aas2, tmass, digits = 4L) {
  ions <- c(tmass - 17.026549, rev(aas2)) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' Z2-ions.
#' 
#' @rdname bions_base
z2ions <- function (aas2, tmass, digits = 4L) {
  ions <- zions_base(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' A-ions.
#' 
#' @rdname bions_base
aions_base <- function (aas2, tmass, digits = 4L) {
  c(tmass - 27.9949146, aas2) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' A2-ions.
#' 
#' @rdname bions_base
a2ions <- function (aas2, tmass, digits = 4L) {
  ions <- aions_base(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' A*-ions.
#' 
#' @rdname bions_base
astarions <- function (aas2, tmass, digits = 4L) {
  # -CO -NH3 = -(27.9949146 + 17.026549)
  c(tmass - 45.0214636, aas2) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' A*2-ions.
#' 
#' @rdname bions_base
astar2ions <- function (aas2, tmass, digits = 4L) {
  ions <- astarions(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' A0-ions.
#' 
#' @rdname bions_base
a0ions <- function (aas2, tmass, digits = 4L) {
  # -CO -H2O = -(27.9949146 + 18.010565)
  c(tmass - 46.0054796, aas2) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' A02-ions.
#' 
#' @rdname bions_base
a02ions <- function (aas2, tmass, digits = 4L) {
  ions <- a0ions(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


#' X-ions.
#' 
#' @rdname bions_base
xions <- function (aas2, tmass, digits = 4L) {
  # +CO -H2 = 27.9949146 - 2*1.007825
  ions <- c(tmass + 25.9792646, rev(aas2)) %>% 
    cumsum() %>% 
    `[`(-1) %>% 
    round(digits = digits)
}


#' X2-ions.
#' 
#' @rdname bions_base
x2ions <- function (aas2, tmass, digits = 4L) {
  ions <- xions(aas2, tmass, digits)
  (ions + 1.00727647)/2
}


