#' Entrez IDs of human molecular signatures.
#'
#' A dataset containing human entrez IDs by the gene sets of molecular signatures.
#'
#' @format A list at a length of 5501 
#' \describe{
#'   \item{hs_...}{human entries}
#' }
#' @source \url{https://www.gsea-msigdb.org/gsea/index.jsp}
"c2_msig_hs"


#' Entrez IDs of mouse molecular signatures.
#'
#' A dataset containing mouse entrez IDs by the gene sets of molecular signatures.
#'
#' @format A list at a length of 5501 
#' \describe{
#'   \item{mm_...}{mouse entries}
#' }
#' @source \url{https://www.gsea-msigdb.org/gsea/index.jsp}
"c2_msig_mm"


#' Entrez IDs of rat molecular signatures.
#'
#' A dataset containing rat entrez IDs by the gene sets of molecular signatures.
#'
#' @format A list at a length of 5501 
#' \describe{
#'   \item{rn_...}{rat entries}
#' }
#' @source \url{https://www.gsea-msigdb.org/gsea/index.jsp}
"c2_msig_rn"


#' Entrez IDs of human gene-ontology annotations.
#'
#' A dataset containing human entrez IDs by the gene sets of GO.
#'
#' @format A list at a length of 18320 
#' \describe{
#'   \item{hs_...}{human entries}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
"go_sets_hs"


#' Entrez IDs of mouse gene-ontology annotations.
#'
#' A dataset containing mouse entrez IDs by the gene sets of GO.
#'
#' @format A list at a length of 18364 
#' \describe{
#'   \item{mm_...}{mouse entries}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
"go_sets_mm"


#' Entrez IDs of rat gene-ontology annotations.
#'
#' A dataset containing rat entrez IDs by the gene sets of GO.
#'
#' @format A list at a length of 18726 
#' \describe{
#'   \item{rn_...}{rat entries}
#' }
#' @source \url{http://current.geneontology.org/products/pages/downloads.html}
"go_sets_rn"


#' Entrez IDs of human KEGG annotations.
#'
#' A dataset containing human entrez IDs by the gene sets of KEGG pathways.
#'
#' @format A list at a length of 229 
#' \describe{
#'   \item{hsa...}{human entries}
#' }
#' @source \url{https://www.genome.jp/kegg/pathway.html}
"kegg_sets_hs"


#' Entrez IDs of mouse KEGG annotations.
#'
#' A dataset containing mouse entrez IDs by the gene sets of KEGG pathways.
#'
#' @format A list at a length of 225 
#' \describe{
#'   \item{mmu...}{mouse entries}
#' }
#' @source \url{https://www.genome.jp/kegg/pathway.html}
"kegg_sets_mm"


#' Entrez IDs of rat KEGG annotations.
#'
#' A dataset containing rat entrez IDs by the gene sets of KEGG pathways.
#'
#' @format A list at a length of 225 
#' \describe{
#'   \item{mmu...}{rat entries}
#' }
#' @source \url{https://www.genome.jp/kegg/pathway.html}
"kegg_sets_rn"


#' Lookups of human or mouse kinases.
#'
#' A dataset containing the \code{refseq}, \code{uniprot}, \code{gene names} of
#' human or mouse kinases.
#'
#' @format A data frame with 2050 rows and 9 variables:
#' \describe{
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{gene}{gene names}
#'   \item{kin_attr}{the attribute being a kinase or not}
#'   \item{kin_class}{the class of a kinase; for example: TK, tyrosine kinase}
#'   \item{kin_order}{the number order of a kinase class by the kinase tree}
#'   \item{uniprot_acc}{UniProt accession number}
#'   \item{uniprot_id}{UniProt ID}
#'   \item{prot_desc}{protein description}
#'   \item{entrez}{Entrez ID}
#' }
#' @source \url{http://kinase.com/human/kinome/phylogeny.html}
"kinase_lookup"


#' BioMart-Ensembl databases for human proteins.
#' 
#' A S4 object from \link[biomaRt]{useMart} at "dataset = hsapiens_gene_ensembl". 
#' @source \url{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}
"mart_hs"


#' BioMart-Ensembl for mouse proteins.
#' 
#' A S4 object from \link[biomaRt]{useMart} at "dataset = mmusculus_gene_ensembl". 
#' @source \url{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}
"mart_mm"


#' BioMart-Ensembl for rat proteins.
#' 
#' A S4 object from \link[biomaRt]{useMart} at "dataset = rnorvegicus_gene_ensembl". 
#' @source \url{https://bioconductor.org/packages/release/bioc/html/biomaRt.html}
"mart_rn"


#' Lookups of cRAP proteins. 
#'
#' A dataset containing the IDs of cRAP proteins. 
#'
#' @format A data frame with 195 rows and 10 variables:
#' \describe{
#'   \item{uniprot_acc}{UniProt accession number}
#'   \item{uniprot_id}{UniProt ID}
#'   \item{status}{the review status according to UniProt}
#'   \item{prot_desc}{protein description}
#'   \item{gene}{gene names}
#'   \item{organism}{the organism attribute according to UniProt}
#'   \item{length}{the number of amino acid residues under a proposed protein}
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{entrez}{Entrez ID (not currently used)}
#'   \item{fasta_name}{Fasta name according to \code{read_fasta(pattern = ...)}}
#' }
#' @source \url{https://www.thegpm.org/crap/}
"prn_annot_crap"


#' Lookups among RefSeq accessions, Entrez IDs and gene names for human
#' proteins.
#'
#' A dataset containing human \code{refseq} accessions, \code{entrez} IDs and
#' \code{gene} names.
#'
#' @format A data frame with 275081 rows and 4 variables:
#' \describe{
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{gene}{Gene name}
#'   \item{entrez}{Entrez ID}
#'   \item{species}{Species name}
#' }
"refseq_entrez_hs"


#' Lookups among RefSeq accessions, Entrez IDs and gene names for mouse
#' proteins.
#'
#' A dataset containing mouse \code{refseq} accessions, \code{entrez} IDs and
#' \code{gene} names.
#'
#' @format A data frame with 209384 rows and 4 variables: 
#' \describe{
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{gene}{Gene name}
#'   \item{entrez}{Entrez ID}
#'   \item{species}{Species name}
#' }
"refseq_entrez_mm"


#' Lookups among RefSeq accessions, Entrez IDs and gene names for rat
#' proteins.
#'
#' A dataset containing mouse \code{refseq} accessions, \code{entrez} IDs and
#' \code{gene} names.
#'
#' @format A data frame with 158959 rows and 4 variables:
#' \describe{
#'   \item{refseq_acc}{RefSeq accession number}
#'   \item{gene}{Gene name}
#'   \item{entrez}{Entrez ID}
#'   \item{species}{Species name}
#' }
"refseq_entrez_rn"


#' Lookups between UniProt accessions and Entrez IDs for human proteins.
#'
#' A dataset containing human \code{entrez} IDs and \code{UniProt} accessions.
#'
#' @format A data frame with 34383 rows and 4 variables:
#' \describe{
#'   \item{uniprot_acc}{UniProt accession number}
#'   \item{gene}{Gene name}
#'   \item{entrez}{Entrez ID}
#'   \item{species}{Species name}
#' }
"uniprot_entrez_hs"


#' Lookups between UniProt accessions and Entrez IDs for mouse proteins.
#'
#' A dataset containing mouse \code{entrez} IDs and \code{UniProt} accessions.
#'
#' @format A data frame with 35818 rows and 4 variables:
#' \describe{
#'   \item{uniprot_acc}{UniProt accession number}
#'   \item{gene}{Gene name}
#'   \item{entrez}{Entrez ID}
#'   \item{species}{Species name}
#' }
"uniprot_entrez_mm"


#' Lookups between UniProt accessions and Entrez IDs for rat proteins.
#'
#' A dataset containing rat \code{entrez} IDs and \code{UniProt} accessions.
#'
#' @format A data frame with 23097 rows and 4 variables:
#' \describe{
#'   \item{uniprot_acc}{UniProt accession number}
#'   \item{gene}{Gene name}
#'   \item{entrez}{Entrez ID}
#'   \item{species}{Species name}
#' }
"uniprot_entrez_rn"


#' A list of UniProt taxon IDs and Species names.
#'
#' The dataset is based on UniProt.ws::availableUniprotSpecies().
#'
#' @format A data frame with 10818 rows and 2 variables:
#' \describe{
#'   \item{taxid}{UniProt taxon ID}
#'   \item{organism}{UniProt Species name}
#' }
"uniprot_species"


#' A list of MaxQuant modifications.
#'
#' The dataset is based on "modificatins.xml" in MaxQuant.
#'
#' @format A data frame with 549 rows and 4 variables:
#' \describe{
#'   \item{title}{Modification}
#'   \item{position}{Position of a modification}
#'   \item{composition}{Composition of elements under a modification}
#'   \item{mass}{Monoisotopic mass of a modification}
#' }
"mq_mods"


#' Masses of amino acid residues.
#'
#' The dataset is based on Mascot.
#'
#' @format A data frame with 26 rows and 6 variables:
#' \describe{
#'   \item{one_letter}{A one-letter representation of amino acid residues}
#'   \item{three_letter}{A three-letter representation of amino acid residues}
#'   \item{fullname}{The full name of amino acid residues}
#'   \item{monoisotopic_da}{Monoisotopic mass of a amino acid residue in Dalton}
#'   \item{average_da}{Average mass of a amino acid residue in Dalton}
#'   \item{composition}{Elemental composition of a amino acid residue}
#' }
"aa_residues"

