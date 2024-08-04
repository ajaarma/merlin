# extra awk script om laatste line eraf te halen
#
skip<-NULL
for(chr in 1:n.chr){
add<-paste(" awk '!/^Swap/' output",chr," > output",chr,sep="")
# sed -i '/^Swap/d' out*
skip<-rbind(skip,add)
}
skip<-as.vector(skip)

write.table(skip,file="skip",
quote=F,sep="\t",row.names=F,col.names=F)


# haal alle output[Chr] files naar outputLOD dir op Windows
# inlezen output files

setwd("//readynasgenet1/Doof/Erik/Linux/cvc_knaepen/output_lod")

n.chr = 22

LOD.scores<-vector("list",n.chr)

for (i in 1:n.chr){
    #temp<-read.table(paste("tstoutput",i,sep=""),header=T,sep="")
    #temp <- read.table(paste("output",i,"_out.txt",sep=""),header=T)
    temp <- read.table(paste("output",i,sep=""),header=T,sep="")
    temp2<-temp[,c(1:2)]
    names(temp2)<-c("map","LOD")
    
LOD.scores[[i]]<-temp2
i}

#Adding up of LOD score
for (i in 1:n.chr)
{
  a = sum(LOD.scores[[i]][2])
  cat(i,"\t",a,"\n",file="LOD_ADD.txt",append=TRUE)

}


# alles afprinten

pdf("graph.pdf")
par(mfrow=c(2,3),oma=c(10,0,10,0),mar=c(3,2,4,2))
for (i in 1:n.chr)
{
attach(LOD.scores[[i]])
plot(map,LOD,type="l",ylim=c(-15,5),col="blue",lwd=1.5,mgp=c(2,0.5,0),
main=paste("Chr ",i,sep=""),xlab="",ylab="",cex.axis=0.8,las=1)
abline(h=0,col="grey")
abline(h=3,col="red",lwd=1.5)
abline(h=-2,col="red",lwd=1.5)
detach(LOD.scores[[i]])
}
dev.off()

# highest and lowest LOD scores
max.lod<-rep(NA,n.chr)
chrom<-seq(1,n.chr)

for (i in 1:n.chr){
max.lod[i]<-max(LOD.scores[[i]]$LOD)
}
lod.peaks<-cbind(chrom,max.lod)

# zoom in on linked regions
# positieve LOD op  chr8
m<-13
peak.region2<-LOD.scores[[m]][LOD.scores[[m]]$LOD>2,]
peak.lowlim2<-range(peak.region2$map)[1]
peak.uplim2<-range(peak.region2$map)[2]

peak.snps2<-mapfile[[m]][mapfile[[m]]$POSITION>=peak.lowlim2 & mapfile[[m]]$POSITION<=peak.uplim2,-1]

temp<-peak.snps2[,c(2,1)]
names(temp)<-c("SNP","Position")
peak.region2$map<-round(peak.region2$map,1)

write.table(peak.region2,"clipboard",row.names=F,quote=F, sep="\t")
#write.table(peak.snps2,"clipboard",row.names=F,quote=F, sep="\t")





#######################################################################
#
# Printing the Ouput of Sum of LOD scores obtained in previous step
#
#
#
#######################################################################
#BRANCH_1 = "/home/aakumar/WORK/LINK_ANAL/inputfiles_BRANCH_1/alloutput"

n.chr = 22
LOD.scores_1 <-vector("list",n.chr)
LOD.scores_2 <- vector("list",n.chr)
LOD.sum12 <- vector("list",n.chr)
#pdf("graph_LOD_RUN6_V1.pdf")
#pdf("graph_LOD_DBIE.pdf")
#pdf("graph_LOD_VASS.pdf")
#system(" sed -i '/^Swap/d' /home/aakumar/WORK/LINK_ANAL/DVD_BRS/DBEC_inputfiles/alloutput/output*") 
pdf("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/DBIE_inputfiles/graph_LOD_KIV_45cM_95Pent.pdf")
par(mfrow=c(2,3),oma=c(10,0,10,0),mar=c(3,2,4,2))
for (i in 1:n.chr){
    print(i)
    #temp_1 <- read.table(paste("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/DBEC_inputfiles/alloutput/output",i,sep=""),header=T,sep="")
    temp_1 <- read.table(paste("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/DBIE_inputfiles/alloutput/output",i,sep=""),header=T,sep="")
    #temp_2 <- read.table(paste("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/DBIE_inputfiles/alloutput/output",i,sep=""),header=T,sep="")
    temp12<-temp_1[,c(1:2)]
    #temp22 <- temp_2[,c(1:2)]
    names(temp12)<-c("map","LOD")
    #names(temp22) <- c("map","LOD")
    LOD.scores_1[[i]]<-temp12
    #LOD.scores_2[[i]]<-temp22
    pos = as.numeric(temp12[,1])
    #lod.sum <- as.numeric(temp12[,2])+ as.numeric(temp22[,2])
    lod.sum <- as.numeric(temp12[,2])#+ as.numeric(temp22[,2])
    plot(pos,lod.sum,type="l",col="blue",lwd=1.5,main=paste("Chr ",i,sep=""),xlab="",ylab="",cex.axis=0.8,las=1,ylim=c(-10,5))
    abline(h=0,col="grey")
    abline(h=3,col="red",lwd=1.5)
    abline(h=-2,col="red",lwd=1.5)
    LOD.sum12[[i]] = data.frame(pos,lod.sum)
    #cat("Chromsome: ",i,"\n",file="LOD_sum_V2.txt",append=TRUE)
    #for(j in 1:length(lod.sum))
    #{
    #   cat("Chromsome: ",i,"\n",file="LOD_sum_V2.txt",append=TRUE)
    #   cat(pos[j],"\t",lod.sum[j],"\n",file="LOD_sum_V2.txt",append=TRUE)
    #}
    #cat("\n",file="LOD_sum_V2.txt",append=TRUE)
     
    }
dev.off()

n.chr="X"
#LOD.scores_1 <- vector("list",n.chr)
pdf("graph_LOD_X.pdf")
for(i in n.chr){
   print(i)
   temp_1 <- read.table("/home/aakumar/WORK/LINK_ANAL/Linkage_data_ilse/inputfilesX/input1/tempoutput",header=T,sep="")
   pos = as.numeric(temp_1[,1])
   lod.score = as.numeric(temp_1[,2])
   plot(pos,lod.score,type="l",col="blue",lwd=1.5,main=paste("Chr ",i,sep=""),xlab="",ylab="",cex.axis=0.8,las=1,ylim=c(-10,5))
   abline(h=0,col="grey")
   abline(h=2,col="red",lwd=1.5)
   abline(h=-2,col="red",lwd=1.5)
   abline(v=33.705,col="lightgray",lwd=1.5,xlab="33.705")
   abline(v=93.705,col="lightgray",lwd=1.5,xlab="93.705")
   }

dev.off()

###################
## DEBIE Family
#
#
####################

n.chr="X"
#LOD.scores_1 <- vector("list",n.chr)
pdf("graph_LOD_X.pdf")
for(i in n.chr){
   print(i)
   #temp_1 <- read.table("/home/aakumar/WORK/LINK_ANAL/Linkage_data_ilse/inputfilesX/input1/tempoutput",header=T,sep="")
   temp_1 <- read.table("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/DBIE_inputfiles/alloutput/output",header=T,sep="")
   pos = as.numeric(temp_1[,1])
   lod.score = as.numeric(temp_1[,2])
   plot(pos,lod.score,type="l",col="blue",lwd=1.5,main=paste("Chr ",i,sep=""),xlab="",ylab="",cex.axis=0.8,las=1,ylim=c(-10,5))
   abline(h=0,col="grey")
   abline(h=2,col="red",lwd=1.5)
   abline(h=-2,col="red",lwd=1.5)
   abline(v=33.705,col="lightgray",lwd=1.5,xlab="33.705")
   abline(v=93.705,col="lightgray",lwd=1.5,xlab="93.705")
   }

dev.off()





