# Trans-QTL pipeline

Main script: [pipeline.R](pipeline.R)

The main goal of this analysis is to compute genome-wide trans scores in a target cohort which will be used to learn the “core genes”. This procedure consists of the following points:

1.  Compute the regional genotypic scores in the target cohort using weights from published eQTL and pQTL studies.
2.  Classify the scores into cis, cis-x and trans.
3.  Construct genome-wide trans scores.
4.  Perform downstream analysis to identify “core genes”.

## Details of implementation and code

### GENOSCORES

We use [GENOSCORES](https://genoscores.cphs.mvm.ed.ac.uk) application to 
compute regional genotypic scores in a target cohort using weights from published GWAS, and adjust them for LD.
 
You will need to install GENOSCORES R package on a computing cluster or powerful Linux box. The detailed steps are provided in https://genoscores.cphs.mvm.ed.ac.uk/quickstart

#### Overview:

1.  Register on our website which will allow using the API: https://genoscores.cphs.mvm.ed.ac.uk/register
2.  Install Genoscores R package in a usual way with `R CMD INSTALL genoscores_xxx.tar.gz`
3.  Install all the dependencies (see https://genoscores.cphs.mvm.ed.ac.uk/quickstart#install-the-dependencies)
4.  Initialise the package and generate API token to access our database (https://genoscores.cphs.mvm.ed.ac.uk/quickstart#initialise-the-package)

#### Compute regional genotypic scores

Compute the regional genotypic scores for the target cohort:

1.  Specify imputed target genotypes in PLINK format. If the genotypes are not in build hg38 specify a liftOver chain file to lift the positions to hg38. GENOSCORES will do the lifting.
2.  Specify LD reference panel in PLINK. We have the most recent release of 1000 Genomes panel which I can share with you if needed.
3.  Specify that weights from eQTLGen study (eQTLs) or Somalogic pQTLs are to be used to compute the scores.
4.  Run the analysis (preferably parallelising by chromosome).

See details in [pipeline.R](pipeline.R).

We suggest running this analysis per chromosome as this can be very computationally intensive.

Use the included bash script `analysis.cli.sh` which will invoke per chromosome analysis.

#### Classify scores and aggregate them into genome-wide trans scores

Having computed the regional scores, we now need to combine them into genome-wide trans
scores and (optionally) remove any HLA scores. This can be done with R script (gw.trans.scores.R) attached.

The script does the following:

1.  Combine scores (if computed per chromosome).
2.  Classify the scores into cis, cis-x and trans, and annotate by genes.
3.  Exclude any scores in HLA (optional).
4.  For each gene sum up all the regional trans scores into a single score.
5.  Also include any cis or cis-x scores for the gene which will be used later.
6.  Process score metadata and write everything into .Rdata.gz files
