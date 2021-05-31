source('AMAT.R')


library(phyloseq) # See https://joey711.github.io/phyloseq/install.html for installation
library(OMiAT)    # See https://github.com/hk1785/OMiAT for installation


## A toy dataset

data(MiData)
otu.tab <- otu_table(MiData)
n<-nrow(otu.tab)

y.con <- rnorm(n)
y.bin <- as.numeric(runif(n) < 0.5)

x1 <- sample_data(MiData)[[3]]
x2 <- sample_data(MiData)[[4]]
cov <- cbind(x1, x2)

# AMAT with continuous outcome

set.seed(1)
test<-AMAT(y=y.con,Z=otu.tab,X=cov,model="c",total = NULL)

test[[1]]  ## AMAT's p value
test[[2]]  ## Testing subset: a subset of the column names of otu.tab

# AMAT with binary outcome

set.seed(1)
test<-AMAT(y=y.bin,Z=otu.tab,X=cov,model="b",total = NULL)

# Test involving a subset of the entire community

total<-apply(otu.tab,1,sum)  # Total reads in the entire community for each sample.
sub.otu.tab<-otu.tab[,1:50]  # Without loss of generality assume that we want to conduct a test that involves the first fifty columns of the OTU table.

set.seed(1)
test<-AMAT(y=y.con,Z=sub.otu.tab,X=cov,model="c",total = total)


# References

## Hyunwook Koh (2020). OMiAT: Optimal Microbiome-based Association Test (OMiAT). R package version 6.0. https://github.com/hk1785/OMiAT

## phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. Paul J. McMurdie and Susan Holmes (2013) PLoS ONE 8(4):e61217.