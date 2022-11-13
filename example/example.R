library(devtools)
library(tidyverse)
library(GSA)

data("1000g")

# generate allele count table for each gene and each sample
gene.ac <- count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "allele")

# generate the table of allele carrier status for each gene and each sample
gene.carrier <- count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "carrier")

# here just pick some genes and use as an example as a gene set input
gene.list <- gene.ac$GENE[10:14]

# count the total number of alleles in the gene set for each sample
geneset.ac <- count_vars_in_gene_set_per_sample(table = gene.ac, gene.set = gene.list, gene.set.name = "genesetA" )

# count the total number of genes carrying at least one variant in the gene set for each sample
geneset.carrier <- count_vars_in_gene_set_per_sample(table = gene.carrier, gene.set = gene.list, gene.set.name = "genesetA" )

# the input ped must be assigned with correct class/type for each variable. e.g.
# * the outcome/dependent variable must be a binary variable, either 0 or 1.
# * according to the natural property of some independent variables, they may need to be assigned with a proper class/type. e.g. usually the variable representing sex should be assigned as a factor.
#
ped.with.geneset.count <- ped %>%
  as_tibble() %>%
  left_join(geneset.ac, by = c("IND_ID" = "sampleID")) %>%
  replace_na(list(genesetA = 0)) %>%
  mutate(DISEASE = DISEASE - 1,
         SEX = factor(SEX))


# outcome variable: it must be a binary variable. The value should be either 0 or 1.
# interested.variables : all variables here must have the same class/type. e.g. genesetA and AGE are both numeric types in our example here. If you add SEX to interested.variables, you will get an error. Because SEX is factor.
glm.results.tidy <- glm_for_multiple_tests(ped = ped.with.geneset.count, outcome = "DISEASE", interested.variables = c("genesetA", "AGE"), covariates= c("QT","SEX"), tidy.output = TRUE, conf.int = TRUE, exponentiate = F)

glm.results.tidy

glm.results <- glm_for_multiple_tests(ped = ped.with.geneset.count, outcome = "DISEASE", interested.variables = c("genesetA", "AGE"), covariates= c("QT","SEX"))

glm.results
