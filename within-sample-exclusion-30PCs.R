args = commandArgs(trailingOnly=TRUE)
pcaFile <- args[1]# "PCA_approx30PCs.eigenvec"
wgsFile <- args[2]# "PCA_approx30PCs.eigenvec"

cat(paste0("pcaFile = ",pcaFile,".\n"))

## Read PCs from 455k samples
pcs <- read.table(pcaFile,h=T,stringsAsFactors = F,comment.char = "!")[,-1]
npc <- ncol(pcs)-1
cat(paste0("Exclusion based on ",npc," PCs.\n"))

OUT <- matrix(0,nrow=nrow(pcs),ncol=npc)
for(i in 1:npc){
  x <- c(scale(pcs[,paste0("PC",i)]))
  OUT[which(abs(x)>4),i] <- 1
}
l   <- rowSums(OUT)
table(l==0)
out <- which(l>0)

png("Histogram_of_within_sample_PCs.png",width=4000,height=3000,res=200)
op <- par(mfrow=c(3,4))
for(i in 1:12){
  s <- sd(pcs[,paste0("PC",i)])
  hist(pcs[,paste0("PC",i)]/s,nclass=100,main=paste0("Within-sample PC",i))
  hist(pcs[-out,paste0("PC",i)]/s,nclass=100,col="pink",add=T)
  abline(v=range(pcs[-out,paste0("PC",i)])/s,col="red",lty=2)
  legend("topleft",legend="Kept For analysis",fill="pink",box.lty=0)
}
par(op)
dev.off()

wgs <- read.table(wgsFile,colClasses = "numeric")
cat(paste0("Found ",nrow(wgs)," individuals with WGS data.\n"))

idsInWGS <- wgs[,2]
eur_with_wgs_qced <- intersect(pcs[-out,"IID"],idsInWGS)
cat(paste0("Finally ",length(eur_with_wgs_qced)," European ancestry individuals with WGS data are kept.\n"))
write.table(cbind(eur_with_wgs_qced,eur_with_wgs_qced),"UKB.EUR_qced_with_wgs.id",quote=F,row.names=F,col.names=F)




 
# # table(rowSums(abs(scale(ukb[as.character(eurid),]))>3)>1)
# # 
# # bp  <- read.table("9280_12505_UKBiobank_SCZ_whereYouLive_15032018.tab",h=T,stringsAsFactors=FALSE,sep="\t")
# # bp  <- bp[,c("f.eid","f.129.0.0","f.130.0.0","f.20075.0.0","f.20074.0.0")]
# # colnames(bp)      <- c("IID","BirthNorthCoord","BirthEastCoord","HomeNorthCoord","HomeEastCoord")
# # bp <- na.omit(bp)
# # bp <- bp[which(apply(bp[,-1],1,min)>0),]
# 
# # load("idsInWGS.RData")
# # idsInWGS <- as.numeric(idsInWGS)
# # idsInWES <- read.table("ukb23158_c1_b0_v1.fam",stringsAsFactors = F)[,2]
# # 
# # s <- intersect(eurid,bp$IID)
# # s <- intersect(s,idsInWES)
# # s <- intersect(s,idsInWGS)
# # 
# # rel <- read.table("ukb_inferredRelationship.txt")
# # colnames(rel) <- c("IID1","IID2","IBD0","IBD1","IBD2","YOB1","YOB2","SEX1","SEX2","REL")
# # 
# # r  <- rel[which(rel$IID1%in%s & rel$IID2%in%s),]
# # i  <- table(c(r$IID1,r$IID2))
# # sx <- s
# # n  <- 0
# # while(max(i)>1){
# #   n <- n + 1
# #   idOut <- names(which.max(i))
# #   sx <- sx[-which(sx%in%idOut)]
# #   r  <- r[which(r$IID1%in%sx & r$IID2%in%sx),]
# #   i  <- table(c(r$IID1,r$IID2))
# # }
# # 
# # r  <- r[which(r$IID1%in%sx & r$IID2%in%sx),]
# # i  <- table(c(r$IID1,r$IID2))
# # s  <- sample(c(1,2),size=nrow(r),replace = T)
# # o  <- sapply(1:nrow(r), function(k) r[k,paste0("IID",s[k])])
# # sx <- sx[-which(sx%in%o)]
# # 
# # if(length(which(r$IID1%in%sx & r$IID2%in%sx))==0){
# #   cat(paste0("\tFinito [n=",length(sx)," samples selected].\n"))
# # }
# 
# 
# 
# 
