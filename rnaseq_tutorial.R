##### This is test 

dpois(x = 3, lambda = 5)

0:12

dpois( x= 0:12, lambda = 5)

barplot(dpois(0:12, 5), names.arg = 0:12, col = "green")

genotype = c("AA","AO","BB","AO","OO","AO","AA","BO","BO",
             "AO","BB","AO","BO","AB","OO","AB","BB","AO","AO")

table(genotype)

genotypeF = factor(genotype)
levels(genotypeF)

table(genotypeF)