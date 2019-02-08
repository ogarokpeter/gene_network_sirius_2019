library(RankAggreg)
filename1 <- '/home/user/Sirius/gene_network_sirius_2019/RankAggregation/data/clr_exps_10_ranked.txt'
filename2 <- '/home/user/Sirius/gene_network_sirius_2019/RankAggregation/data/aracne_exps_10_ranked.txt'

d <- read.table(filename1)
d <- as.matrix(d)
print(d)
CES <- RankAggreg(d, 3, method="GA", maxIter=10, verbose=False)
res <- CES$top.list
print(res)
