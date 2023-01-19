## read individual-level phenotype dataset
phenotype <- fread(pheno.file)

phenotype <- phenotype[!is.na(PHENO)]
setkey(phenotype, id)
