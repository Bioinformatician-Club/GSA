library(devtools)
library(tidyverse)
check()

data("1000g")

# function 1 based on maidou version
gene.ac <- count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "allele")

gene.carrier <- count_vars_per_gene_per_sample(table = variants, gene.col = 8, gt.col = 12, type = "carrier")

# function 2 based on xiaoyi version
gene.list <- gene.ac$GENE[10:14]

geneset.ac <- count_vars_in_gene_set_per_sample(table = gene.ac, gene.set = gene.list, gene.set.name = "genesetA" )

geneset.carrier <- count_vars_in_gene_set_per_sample(table = gene.carrier, gene.set = gene.list, gene.set.name = "genesetA" )

# function 3
# the input ped must be set with correct class/type for each variable. e.g. outcome/y/dependent variable must be a binary variable. Some of independent variables need to be set in proper type, for instance, sex variable normally should be set to factor.

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
