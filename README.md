  - [Introduction to proteoQ](#introduction-to-proteoq)
  - [Installation](#installation)
  - [1 Data normalization](#data-normalization)
      - [1.1 Set up experiment for Mascot
        workflow](#set-up-experiment-for-mascot-workflow)
          - [1.1.1 Prepare fasta databases](#prepare-fasta-databases)
          - [1.1.2 Prepare Mascot PSM data](#prepare-mascot-psm-data)
          - [1.1.3 Prepare metadata](#prepare-metadata)
          - [1.1.4 Load experiment](#load-experiment)
      - [1.2 Summarise Mascot PSMs](#summarise-mascot-psms)
          - [1.2.1 Process PSMs](#process-psms)
          - [1.2.2 Summarize PSMs to
            peptides](#summarize-psms-to-peptides)
          - [1.2.3 Summarize peptides to
            proteins](#summarize-peptides-to-proteins)
      - [1.3 Renormalize data against column
        subsets](#renormalize-data-against-column-subsets)
      - [1.4 Renormalize data against row
        subsets](#renormalize-data-against-row-subsets)
      - [1.5 Summarize MaxQuant results](#summarize-maxquant-results)
      - [1.6 Summarize Spectrum Mill
        results](#summarize-spectrum-mill-results)
      - [1.7 Workflow scripts](#workflow-scripts)
  - [2 Basic informatics](#basic-informatics)
      - [2.1 MDS and PCA plots](#mds-and-pca-plots)
      - [2.2 Correlation plots](#correlation-plots)
      - [2.3 Heat maps](#heat-maps)
      - [2.4 Significance tests and volcano plot
        visualization](#significance-tests-and-volcano-plot-visualization)
      - [2.5 Gene sets under volcano
        plots](#gene-sets-under-volcano-plots)
      - [2.6 Trend Analysis](#trend-analysis)
      - [2.7 NMF Analysis](#nmf-analysis)
      - [2.8 STRING Analysis](#string-analysis)
      - [2.9 Missing value imputation](#missing-value-imputation)
  - [3 Labs](#labs)
      - [3.1 Reference choices](#reference-choices)
          - [3.1.1 References on data
            scaling](#references-on-data-scaling)
          - [3.1.2 References on data CV](#references-on-data-cv)
      - [3.2 Data subsets](#data-subsets)
      - [3.3 Random effects](#random-effects)
          - [3.3.1 Single random effect](#single-random-effect)
          - [3.3.2 Multiple random effects](#multiple-random-effects)
  - [4 Column keys](#column-keys)
      - [4.1 Mascot](#mascot)
          - [4.1.1 PSMs](#psms)
          - [4.1.2 Peptides](#peptides)
          - [4.1.3 Proteins](#proteins)
      - [4.2 MaxQuant](#maxquant)
          - [4.2.1 PSMs](#psms-1)
          - [4.2.2 Peptides](#peptides-1)
          - [4.2.3 Proteins](#proteins-1)
  - [References](#references)

<style>
p.comment {
background-color: #e5f5f9;
padding: 10px;
border: 1px solid black;
margin-left: 0px;
border-radius: 5px;
}

</style>

## Introduction to proteoQ

Chemical labeling using tandem mass tag
([TMT](https://en.wikipedia.org/wiki/Tandem_mass_tag)) has been commonly
applied in mass spectrometry (MS)-based quantification of proteins and
peptides. The `proteoQ` tool is designed for automated and reproducible
analysis of proteomics data. It interacts with an `Excel` spread sheet
for dynamic sample selections, aesthetics controls and statistical
modelings. It further integrates the utilities of data filtration and
ordering into functions at the users' interface. The arrangements allow
users to put *ad hoc* manipulation of data behind the scene and instead
apply metadata to openly address biological questions using various
informatic tools. In addition, the entire workflow is documented and can
be conveniently reproduced upon revisiting.

The tool currently processes the peptide spectrum matches (PSM) tables
from [Mascot](https://http://www.matrixscience.com/),
[MaxQuant](https://www.maxquant.org/) and [Spectrum
Mill](https://www.agilent.com/en/products/software-informatics/masshunter-suite/masshunter-for-life-science-research/spectrum-mill)
searches, for 6-, 10- or 11-plex TMT experiments using Thermo's orbitrap
mass analyzers. Peptide and protein results are then produced with
users' selection of parameters in data filtration, alignment and
normalization. The package further offers a suite of tools and
functionalities in statistics, informatics and data visualization by
creating 'wrappers' around published R routines.

<p class="comment">

Click
<strong>[here](https://htmlpreview.github.io/?https://github.com/qzhang503/proteoQ/blob/master/README.html)</strong>
to render a html version of the README.

</p>

