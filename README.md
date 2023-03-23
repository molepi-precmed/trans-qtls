# Genomewide aggregated trans- effects (GATE) computational pipeline

This pipeline is used to compute genome-wide trans scores in a target cohort and learn the “core genes”. 
The code in this repository was used in the analysis reported in:


*Iakovliev, A., McGurnaghan, S., Hayward, C., Colombo, M., Lipschutz, D., Spiliopoulou, A., Colhoun, H.M., McKeigue, P.M., 2023. Trans-eQTL effects on risk of type 1 diabetes: a test of the sparse effector (omnigenic) hypothesis of complex trait genetics. American Journal of Human Genetics (in press).*


The procedure consists of the following:

1.  Compute the regional genotypic scores in the target cohort using weights from published eQTL, pQTL and other studies available at [https://genoscores.cphs.mvm.ed.ac.uk/studies](https://genoscores.cphs.mvm.ed.ac.uk/studies).
2.  Classify the scores into cis, cis-x and trans depending on the distance from the corresponding transcript.
3.  Construct genome-wide trans scores by aggregating regional trans scores.
4.  Perform downstream analysis to identify “core genes”.

Main script: [pipeline.R](pipeline.R)

## Details

### GENOSCORES

- We use [GENOSCORES](https://genoscores.cphs.mvm.ed.ac.uk) application and R package (distributed separately) to compute regional genotypic scores in a target cohort using weights from published GWAS, adjusting for LD.
 
- You will need to install GENOSCORES R package on a computing cluster or powerful Linux box. The installation process is described at [https://genoscores.cphs.mvm.ed.ac.uk/quickstart](https://genoscores.cphs.mvm.ed.ac.uk/quickstart)

#### Installation overview:

1.  Register on our website to start using the API: [https://genoscores.cphs.mvm.ed.ac.uk/register](https://genoscores.cphs.mvm.ed.ac.uk/register).
2.  Install package dependencies as listed in [https://genoscores.cphs.mvm.ed.ac.uk/quickstart#install-the-dependencies](https://genoscores.cphs.mvm.ed.ac.uk/quickstart#install-the-dependencies) .
3. Install GENOSCORES R package with `R CMD INSTALL genoscores_xxx.tar.gz`.
4.  Initialise the package and generate API token to access our database [https://genoscores.cphs.mvm.ed.ac.uk/quickstart#initialise-the-package](https://genoscores.cphs.mvm.ed.ac.uk/quickstart#initialise-the-package).

#### Overview of the procedure to compute regional genotypic scores

See [pipeline.R](pipeline.R), [README](genoscores/README.md), [example.analysis.R](genoscores/example.analysis.R) and [analysis.cli.sh](genoscores/example.analysis.cli.sh) for details.

1. Specify imputed target genotypes in PLINK format. If the genotypes are not in build hg38 specify a liftOver chain file to lift the positions to hg38. GENOSCORES will do the lifting.
2. Specify LD reference panel in PLINK. We have the most recent release of 1000 Genomes panel which can be shared upon request.
3. Specify studies from which the weights are to be used in the analysis (eQTLGen, Somalogic pQTL, etc.).
4.  Run the analysis parallelising by chromosome (this is facilitated by [analysis.cli.sh](genoscores/example.analysis.cli.sh)).

#### Classify scores and aggregate them into genome-wide trans scores

Having computed the regional scores, we now need to combine them into genome-wide trans scores. This can be done using [gw.trans.scores.R](gw.trans.scores.R). 
This script can be used as a template to write R code customised for your specific analysis.

The steps already implemented are:

1. Load the scores and the metadata.
2. Classify the scores into cis, cis-x and trans, annotate eQTL by genes. Annotation of pQTL can only be performed on Genoscores server. Contact us if you need pQTL annotation, or annotate pQTLs manually using [uniprot](https://www.uniprot.org) table. We provide [annotate.pqtls.R](annotate.pqtls.R) as an example.
3. Optionally exclude any scores in HLA region.
4. Optionally compute genome-wide scores for Immunce cells study.
5. For each gene sum up all the regional trans scores into a single genome-wide score.
6. Save processed scores into specified output directory.

### Compute associations between the scores and a phenotype of interest

The scripts in this section were written for T1D case/control study but can be adapted to other autoimmune disorders.

#### Create R object with genes implicated in disease

We provide [diabetesgenes.R](diabetesgenes.R) script as an example of how we created an R object containing all known T1D and T2D genes that we use to annotate findings from trans-eQTL analysis.

To help with annotation we also include a precomputed [allgenes.RData](helper-data/allgenes.RData) object which contains all gene symbols with annotation from `ensembl`.

#### Prepare scores and metadata for the association analysis

Two options are available in [qtl.metadata.R](qtl.metadata.R):

- Load and prepare eQTL scores (`analysis=eQTL`)
- Load and prepare pQTL scores (`analysis=pQTL`)

#### Prepare phenotype for the association analysis

The script [phenotype.R](phenotype.R) can be used to write code to load and process phenotype data. This will be unique for each analysis.

The only requirements are:

1. Phenotype data should be in the data.table format.
2. No N/A are allowed in the phenotype data.table.
3. IDs have to be matched between scores and phenotypes.

#### Compute associations

The script [associations.R](association.R) can be used to compute associations between the scores and the phenotype of interest. It also has some useful code to post-process and annotate associations.

### Next steps

The above code defines a framework sufficient to identify the "core genes" in the specified target cohort. However, we are planning to add some helper scripts to outline various downstream analyses which can be performed with the computed associations.

## License

This code was developed by Paul McKeigue and Andri Iakovliev and is licensed under [GPL-3 license](https://www.gnu.org/licenses/gpl-3.0.txt)
