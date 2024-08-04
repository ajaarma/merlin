# Formatteren linkage data naar  MEGA2 files
##################################################
# Names file
# Pedigree (pre-makeped) file
# Map file

# dit deel zelfde voor Simwalk of MERLIN


# read in full data
##############################

# komen uit Excell, met rijen = merker en kolommen = ID
system("mkdir inputfiles",wait=TRUE)

Data<-read.table("../Full_Data_Table.txt",header=T,sep='\t',dec=",",na.strings="NC",stringsAsFactors=F)
root <- Data[,1:10]
som.cols<-grep("^KIV",names(Data))
data.som<-as.data.frame(cbind(root,Data[,som.cols]))
Data = data.som
#colnames(Data) =  gsub("95BY","000095",colnames(Data))

snp.chr<-table(Data[,4])[-c(23:26)]
snp.chr<-as.data.frame(cbind(as.numeric(snp.chr),as.numeric(row.names(snp.chr))))
names(snp.chr)<-c("n.snp","Chr")

# X, Y en XYi en MT merkers eruit
# enkel naam, chromos en positie overhouden, plus GTs
autosom<-Data[Data$Chr!="X"&Data$Chr!="Y"&Data$Chr!="XY"&Data$Chr!="MT"&Data$Chr!=0,-c(1,3,6:10)]
autosom$Chr<-as.numeric(autosom$Chr)
autosom = autosom[substr(autosom$Name,1,2)=="rs",]
o<-order(autosom$Chr,autosom$Name)
autosom<-autosom[o,]

map.0<-autosom[,c(2,3,1)]
#names(map.0)<-c("CHROMOSOME","physical.p","MARKER")
names(map.0)<-c("CHROMOSOME","physical.p","MARKER")

n.ID<-18
n.col<-n.ID+3

autosom.num<-autosom[,c(1:3)]
for(i in 4:n.col){
  temp<-as.numeric(as.factor(autosom[,i]))-1
  autosom.num<-cbind(autosom.num,temp)
}

autosom.num$sum.GT<-apply(autosom.num[,c(4:n.col)],1,sum,na.rm=T)
autosom.num$inform<-(autosom.num$sum.GT>1&autosom.num$sum.GT<(2*n.ID-1))
table(autosom.num$inform)

##########################
# FALSE   TRUE 
# 32394 215821

inform.set<-autosom.num[autosom.num$inform==T,1]
map.inform<-map.0[map.0$MARKER%in%inform.set,]

o<-order(map.inform$CHROMOSOME,map.inform$physical.p)
map.inform<-map.inform[o,]

n.chr<-22
mapfile.0<-vector("list",n.chr)

for(i in 1:n.chr){
  temp<-map.inform[map.inform$CHROMOSOME==i,]
  temp$POSITION<-temp$physical.p
  mapfile.0[[i]]<-temp[,c(1,4,3)]
}

interval<-200000

# fix number of segments for each chromosome1, by intermarker distance
# boundaries per segment 
chr.length<-rep(NA,n.chr)
n.markerseg<-rep(NA,n.chr)
boundaries<-vector("list",n.chr)
for (i in 1:n.chr){

  chr.length[i]<-max(mapfile.0[[i]]$POSITION)
  n.markerseg[i]<-round(chr.length[i]/interval,0)
  boundaries[[i]]<-seq(1,chr.length[i],interval)
}

marker.seg<-vector("list",n.chr)
marker.seg[[i]]<-vector("list",n.markerseg[i])

selected.marker<-c(0,0,"noname")

for (i in 1:n.chr){
  for (j in 1:n.markerseg[i]){

    # upper and bottom boundary from segments, for i chromosoom, j intervals

    lowlim<-boundaries[[i]][j]
    uplim<-boundaries[[i]][j+1]

    marker.seg[[i]][[j]]<-mapfile.0[[i]][mapfile.0[[i]]$POSITION>lowlim & mapfile.0[[i]]$POSITION<uplim,]


    # within each segment, select first merker
    r<-1
    # some are leeg
    temp<-marker.seg[[i]][[j]][r,]

    if(is.na(temp$MARKER)==FALSE)
     {
        selected.marker<-rbind(selected.marker,temp)
     }
     else
     {
       cat(i,"\t",j,"\n")
     }
  }
}

selected<-selected.marker[-1,]

mapfile<-vector("list",n.chr)
for(i in 1:n.chr){
  temp<-selected[selected$CHROMOSOME==i,]
  temp$POSITION<-round(as.numeric(temp$POSITION)/1000000,digits=4)
  mapfile[[i]]<-temp[,c(1,2,3)]
  #write.table(mapfile[[i]],paste("inputfiles/map.",i,sep=""),quote=F,row.names=F,sep="\t")
}

snp.chr<-table(selected[,1])
snp.chr<-as.data.frame(cbind(as.numeric(snp.chr),as.numeric(row.names(snp.chr))))
names(snp.chr)<-c("n.snp","Chr")


# sorteren op chrosomoom en NAAM (NIET : position)

o<-order(autosom$Chr,autosom$Name)
autosom<-autosom[o,]




# Names file - met chromosoom om later op te splitsen
Names.1<-autosom[,c(1,2)]
n.rows<-dim(Names.1)[1]
Type<-as.character(rep("M",length=n.rows))
Names.2<-data.frame(Type,Names.1,stringsAsFactors=F)

Names.3<-Names.2[Names.2$Name%in%selected$MARKER,]

trait.row<-c("A","TRAIT")
n.chr<-22
namesfile<-vector("list",n.chr)

for (i in 1:n.chr){

  temp<-Names.3[Names.3$Chr==i,c(1,2)]
  namesfile[[i]]<-rbind(trait.row,temp)

write.table(namesfile[[i]],paste("inputfiles/names.",i,sep=""),
quote=F,row.names=F,col.names=F,sep="\t")
}


# names file for each chromosome
##################################"


# Map file : chromos, physical position en merkernaam
###################################################
# Pedigree (pedin) file
###############################

# Eerste 6 kolommen van pre-file
# prefile moet nog getrimt ! bevat kolom teveel en stalen teveel!!!

#raw.Pre<-read.table("pedfile_branch_1_affected_only.txt",header=F,sep="\t",stringsAsFactors=F)
raw.Pre<-read.table("pedfile_branch_1.txt",header=F,sep="\t",stringsAsFactors=F)

raw.Pre$V2<-noquote(raw.Pre$V2)

# staal 10 mislukt - weghalen uit Prefile

#overbodig<-c("10")
#Pre.6<-raw.Pre[raw.Pre$V2%in%overbodig==F,]

Pre.6 = raw.Pre[-1,-c(7:8)]

names(Pre.6)<-c("Pedigree","ID","Father","Mother","Sex","TRAIT.A")
Pre.6$Pedigree<-1

# blijft over : 36 stalen, waarvan 23 getypeerd en 13 niet

# sorteer Pre.file op subjectID
o<-order(Pre.6$ID)
Pre.6<-Pre.6[o,]

# voeg genotypes toe

gt.1<-autosom[,-3]
gt.1<-autosom[autosom$Name%in%Names.3$Name,-c(3,4:14,21)]

# voeg lege kolom bij voor ongetypeerde familieleden 
# (die wel in de ped-file zitten)

# 37 rijen in pre-file, 24 stalen getypeerd
# voeg 13 lege kolommen toe van ongetypeerde founders

nr.snps<-dim(gt.1)[1]
KIV_000025<-matrix(NA,nrow=nr.snps,ncol=1)
KIV_000028<-matrix(NA,nrow=nr.snps,ncol=1)

#HOP.39<-matrix(NA,nrow=nr.snps,ncol=1)
#HOP.40<-matrix(NA,nrow=nr.snps,ncol=1)
#HOP.41<-matrix(NA,nrow=nr.snps,ncol=1)
#HOP.42<-matrix(NA,nrow=nr.snps,ncol=1)
#HOP.43<-matrix(NA,nrow=nr.snps,ncol=1)

empty<-as.data.frame(cbind(KIV_000025,KIV_000028))

names(empty)<-c("KIV_000025","KIV_000028")

gt.2<-data.frame(gt.1,empty,stringsAsFactors=F)


# moet getransposeerd
library(reshape2)

temp<-melt(gt.2,id.var=c("Name","Chr"))

gt.melt<-melt(gt.2,id.var=c("Name","Chr"))
gt.melt$ID<-noquote(substr(gt.melt$variable,5,10))

gt.melt$value<-as.factor(gt.melt$value)

# recodeer genotype naar 0.1.2
gt.melt$value<-as.numeric(gt.melt$value)-1

### dim(gt.melt)
##  9183955       5
###
#################################################""
# sorteer op ID, binnen chromosoom
# dus zelfde sibj orde als in Pre.6
#library("data.table")
#oo<-order(gt.melt$Chr,gt.melt$ID,gt.melt$Name)
#gt.melt<-gt.melt[oo,]
#gt.melt.table <- data.table(gt.melt)
#oo = setkey(gt.melt.table,Chr, ID)
#oo <- oo[,.SD[order(Name)],by=list(Chr,ID)]
#gt.melt <- data.frame(oo)

# hiertussen : op NA zetten van foute GT's (mistyping analyse) (post_mistyping.R)


# genotypes opsplitsen per chromosoon
# per chromosoom naar breed formaat
n.chr<-22
genotypes<-vector("list",n.chr)

for (i in 1:n.chr){
temp<-gt.melt[gt.melt$Chr==i,-c(2,3)]
genotypes[[i]]<-dcast(temp,ID~Name)
}

#system("mkdir inputfiles",wait=FALSE)
system("mkdir inputfiles/input_segments",wait=TRUE)
system("mkdir inputfiles/input_segments/rename",wait=TRUE)
system("mkdir inputfiles/input_segments/mega2_lod",wait=TRUE)
system("mkdir inputfiles/input_segments/mega2_mistyp",wait=TRUE)

save(gt.melt,file="inputfiles/gt.melt")

save(snp.chr,file="inputfiles/snp_chr_df")
save(genotypes,file="inputfiles/genotypes_df")
save(namesfile,file="inputfiles/namesfile_df")

# met deze genotypes-list worden gesegmenteerde genotype-files gemaakt (zie segmentalte.R)
# GOTO segmentatie
