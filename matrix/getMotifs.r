# your input (example):

# T	G	C	A
# 0.2927	0.1933	0.2472	0.2668
# 0.0554	0.9208	0.0233	0.0005
# 0.2316	0.7674	0.0005	0.0005
# 0.0005	0.0005	0.9985	0.0005
# 0.9861	0.0036	0.0098	0.0005
# 0.0005	0.9985	0.0005	0.0005
# 0.9208	0.0741	0.0046	0.0005
# 0.8037	0.0461	0.0451	0.1052
# 0.5498	0.0813	0.2585	0.1104

# inverse the letter if you want to get the reverse complement

library(Cairo)
library(motifStack)

# a <- read.table("ER7_new_matrix.txt",sep="\t",header=TRUE)
# Cairo(600, 400,type="png", bg="white",file="cistrome_degenerated_logo.png")
# b <- apply(a, 2, rev) #add this line if you want to get the reverse complement
# pfm <- t(as.matrix(b))
# plotMotifLogo(pfm)
# dev.off()

a <- read.table("ARF2_DR_matrix.txt",sep="\t",header=TRUE)
a<-data.frame(a[,4],a[,3],a[,2],a[,1])
names(a)<-c("T","G","C","A")
Cairo(600, 400,type="png", bg="white",file="ARF2_DR_matrix.png")
# CairoSVG("complete_ER7_from_sequences_frequences.svg",12,6, bg="white")
#a <- apply(a, 2, rev) 
pfm <- t(as.matrix(a))
plotMotifLogo(pfm)
dev.off()