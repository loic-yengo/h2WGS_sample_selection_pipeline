args = commandArgs(trailingOnly=TRUE)
refFile <- args[1]# "1kg10PCs.eigenvec"
ukbFile <- args[2]# "ukbForPCAproj10PCs.proj.eigenvec"
popFile <- args[3]# "1kg_pop2504.txt"

cat(paste0("refFile = ",refFile,".\n"))
cat(paste0("ukbFile = ",ukbFile,".\n"))
cat(paste0("popFile = ",popFile,".\n"))

## Load PCS
ref <- read.table(refFile,stringsAsFactors=FALSE)
ukb <- read.table(ukbFile,stringsAsFactors=FALSE)

## get parameters
spop  <- read.table(popFile,header=TRUE,stringsAsFactors=FALSE)
spop  <- spop[which(spop[,"FID"]%in%ref[,2]),]; rownames(spop) <- spop[,"FID"]
Pops  <- unique(spop[,"SuperPOP_CODE"])
Cols  <- c("pink","lightblue","lightgreen","grey","khaki3")
names(Cols) <- Pops

nPCS  <- ncol(ref)-2
for(k in 1:nPCS){
  ref[,2+k] <- ref[,2+k] 
  ukb[,2+k] <- ukb[,2+k] 
}

rownames(ref) <- ref[,2]
rownames(ukb) <- ukb[,2]

ref <- as.matrix(ref[,-(1:2)])
ukb <- as.matrix(ukb[,-(1:2)])

n <- nrow(ref)
N <- nrow(ukb)

meanEUR   <- apply(ref[which(spop[rownames(ref),"SuperPOP_CODE"]=="EUR"),],2,mean)
sdEUR     <- apply(ref[which(spop[rownames(ref),"SuperPOP_CODE"]=="EUR"),],2,sd)

nSD       <- 3 # Wainschtein et al. (2022)
criterion <- sapply(1:nrow(ukb),function(i) sum(abs(ukb[i,]-meanEUR)<nSD*sdEUR))
eurid     <- as.numeric(rownames(ukb[which(criterion==nPCS),]))

## Draw PC plots
png("PC-plots.png",width=3000,height=2000,res=200)
op <- par(mfrow=c(2,3))
for(i in 1:3){
  for(j in (i+1):4){
    plot(ref[,i],ref[,j],col=Cols[spop[,"SuperPOP_CODE"]],pch=19,
         xlab=paste0("Eig",i),ylab=paste0("Eig",j),axes=FALSE)
    axis(1);axis(2);abline(h=0,v=0,col="coral1")
    l <- which(criterion==nPCS)
    points(ukb[l,i],ukb[l,j],col="darkred",pch=19,cex=0.5)
    legend("bottomright",legend=c(names(Cols),"UKB-EUR"),col=c(Cols,"darkred"),
           box.lty=0,pch=19)
    
  }
}
par(op)
dev.off()

eurid <- rownames(ukb[l,])
write.table(cbind(eurid,eurid),"UKB.EUR.id",quote=F,row.names=F,col.names=F)
cat(paste0("Output file: UKB.EUR.id (N=",length(eurid),").\n"))

# ## Read PCs from 455k samples
# pcs <- read.table("PCA_approx30PCs.eigenvec",h=T,stringsAsFactors = F,comment.char = "!")[,-1]
# OUT <- matrix(0,nrow=nrow(pcs),ncol=30)
# for(i in 1:30){
#   x <- c(scale(pcs[,paste0("PC",i)]))
#   OUT[which(abs(x)>4),i] <- 1
# }
# l   <- rowSums(OUT)
# table(l==0)
# 
# out <- which(l>0)#which(rowSums(abs(scale(pcs[,paste0("PC",1:30)]))>3)>1)
# 
# png("Histogram_of_within_sample_PCs.png",width=4000,height=3000,res=200)
# op <- par(mfrow=c(3,4))
# for(i in 1:12){
#   s <- sd(pcs[,paste0("PC",i)])
#   hist(pcs[,paste0("PC",i)]/s,nclass=100,main=paste0("Within-sample PC",i))
#   hist(pcs[-out,paste0("PC",i)]/s,nclass=100,col="pink",add=T)
#   abline(v=range(pcs[-out,paste0("PC",i)])/s,col="red",lty=2)
#   legend("topleft",legend="Kept For analysis",fill="pink",box.lty=0)
# }
# par(op)
# dev.off()
# 
# load("idsInWGS.RData")
# eur_with_wgs_qced <- intersect(pcs[-out,"IID"],idsInWGS)
# write.table(cbind(eur_with_wgs_qced,eur_with_wgs_qced),"UKB.EUR_qced_with_wgs.id",quote=F,row.names=F,col.names=F)
# 
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
