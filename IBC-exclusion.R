args = commandArgs(trailingOnly=TRUE)
ibcFile <- args[1]# "F.ibc"

cat(paste0("ibcFile = ",ibcFile,".\n"))

## Read PCs from 455k samples
ibc <- read.table(ibcFile,h=T,stringsAsFactors = F,comment.char = "!")
out <- which(abs(ibc[,"Fhat3"]-mean(ibc[,"Fhat3"]))>3*sd(ibc[,"Fhat3"]))
cat(paste0("Exclusion of ",length(out)," individuals based on inbreeding coefficients.\n"))

png("Histogram_of_F.png",width=1000,height=1000,res=200)
hist(ibc[-out,"Fhat3"],nclass=50,xlim=0.025*c(-1,1),
     xlab="Inbreeding coefficient",main="Histogram of F")
dev.off()

eur_with_wgs_F_qced <- ibc[-out,"IID"]
cat(paste0("Finally ",length(eur_with_wgs_F_qced)," are kept.\n"))
write.table(cbind(eur_with_wgs_F_qced,eur_with_wgs_F_qced),"UKB.EUR_qced_with_wgs_no_F_outliers.id",quote=F,row.names=F,col.names=F)
