library(reshape)

system(" mkdir /home/aakumar/WORK/LINK_ANAL/DVD_BRS/VASS_inputfiles")

Data <-read.table("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/DATA/Full_Data_Table_22_indiv.txt",header=T,sep='\t',dec=",",na.strings="NC",stringsAsFactors=F,fill=T)
root <-Data[,1:10]

# VASS
vass.cols <-grep("^VASS",names(Data))
data.vass <-as.data.frame(cbind(root,Data[,vass.cols]))
Data = data.vass

#snp.chr = table(Data$Chr)[-c(23:26)]
snp.chr = table(Data$Chr)[-c(1,24:27)]
snp.chr <-as.data.frame(cbind(as.numeric(snp.chr),as.numeric(row.names(snp.chr))))
autosom.0<-Data[Data$Chr!="X"&Data$Chr!="Y"&Data$Chr!="XY"& Data$Chr!="MT"&Data$Chr!=0,-c(1,3,6:10)]
autosom.0$Chr<-as.numeric(autosom.0$Chr)

#autosom = autosom.0[substr(autosom.0$Name,1,2)=="rs",]
autosom = autosom.0[grep("rs",autosom.0$Name),]
o<-order(autosom$Chr,autosom$Name)
autosom<-autosom[o,]

map.0<-autosom[,c(2,3,1)]
names(map.0)<-c("CHROMOSOME","physical.p","MARKER")

#n.ID<-18
#n.ID : It is the number of Patients having genotype data.
n.ID<-8  
n.col<-n.ID+3

# make GT numeric and add up
autosom.num<-autosom[,c(1:3)]
for(i in 4:n.col){
  temp<-as.numeric(as.factor(autosom[,i]))-1
  autosom.num<-cbind(autosom.num,temp)
}

autosom.num$sum.GT<-apply(autosom.num[,c(4:n.col)],1,sum,na.rm=T)
autosom.num$inform<-(autosom.num$sum.GT>1&autosom.num$sum.GT<(2*n.ID-1))
table(autosom.num$inform)
#FALSE   TRUE 
# 103955 174226
#FALSE   TRUE 
# 39645 208570 

inform.set<-autosom.num[autosom.num$inform==T,1]
map.inform<-map.0[map.0$MARKER%in%inform.set,]

o<-order(map.inform$CHROMOSOME,map.inform$physical.p)
map.inform<-map.inform[o,]

n.chr<-22
mapfile.0<-vector("list",n.chr)

for(i in 1:n.chr){
  temp<-map.inform[map.inform$CHROMOSOME==i,]
  colnames(temp)[2] <- "POSITION"
  #temp$POSITION<-temp$physical.p
  #mapfile.0[[i]]<-temp[,c(1,4,3)]
  mapfile.0[[i]]<-temp
}



#interval<-300000
#interval<-200000
interval<-450000
#interval <-50000

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

# segment the markerlist per chromosoom
# using boundaries

marker.seg<-vector("list",n.chr)
marker.seg[[i]]<-vector("list",n.markerseg[i])

selected.marker<-c(0,0,"noname")
#### Selected #######
##  5922

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

# remove first row
selected<-selected.marker[-1,]

mapfile<-vector("list",n.chr)
for(i in 1:n.chr){
  temp<-selected[selected$CHROMOSOME==i,]
  temp$POSITION<-round(as.numeric(temp$POSITION)/1000000,digits=4)
  mapfile[[i]]<-temp[,c(1,3,2)]
  #write.table(mapfile[[i]],paste("VASS_inputfiles/map.",i,sep=""), quote=F,row.names=F,sep="\t")
  write.table(mapfile[[i]],paste("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/VASS_inputfiles/map.",i,sep=""), quote=F,row.names=F,sep="\t")
}

# per Chr calculate how many SNps are left
snp.chr<-table(selected[,1])
snp.chr<-as.data.frame(cbind(as.numeric(snp.chr),as.numeric(row.names(snp.chr))))
names(snp.chr)<-c("n.snp","Chr")


##################################
# names file for each chromosome
##################################"

# Names file - met chromosoom om later op te splitsen
Names.1<-autosom[,c(1,2)]
n.rows<-dim(Names.1)[1]
Type<-as.character(rep("M",length=n.rows))
Names.2<-data.frame(Type,Names.1,stringsAsFactors=F)


# only selected merkers
Names.3<-Names.2[Names.2$Name%in%selected$MARKER,]

trait.row<-c("A","TRAIT")
n.chr<-22
namesfile<-vector("list",n.chr)

for (i in 1:n.chr){
  
  temp<-Names.3[Names.3$Chr==i,c(1,2)]
  namesfile[[i]]<-rbind(trait.row,temp)
  
  write.table(namesfile[[i]],paste("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/VASS_inputfiles/names.",i,sep=""), quote=F,row.names=F,col.names=F,sep="\t")
}

###############################
# Pedigree (pedin) file
###############################

# Eerste 6 kolommen van pre-file
# prefile moet nog getrimt ! bevat kolom teveel en stalen teveel!!!
#raw.pre<-read.table("pedfile_branch_1_affected_only.txt",header=F,sep="\t")
raw.pre<-read.table("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/PED/VASS_Ped_file.csv",header=F,sep="\t")
#pre.6<-raw.pre[,-7]
#pre.6 = raw.pre[-1,-c(7:8)]
pre.6 = raw.pre[-1,-7]

names(pre.6)<-c("Pedigree","ID","Father","Mother","Sex","TRAIT.A")

pre.6$Pedigree<-as.numeric(pre.6$Pedigree)

# sorteer pre.file op subjectID
o<-order(pre.6$ID)
pre.6<-pre.6[o,]

# remove 2 samples with unknownGT
##pre.6<-pre.6[pre.6$TRAIT.A!=0,] #Done Just for the affected only analysis on 26.03.2014

# voeg genotypes toe
# eerste kolommen verwijderen van stalen met onbekend	 GT
# en kolom met position

#gt.1<-autosom[autosom$Name%in%Names.3$Name,-c(3,10:11)]
#gt.1<-autosom[autosom$Name%in%Names.3$Name,-c(3,4:14,21)]
gt.1<-autosom[autosom$Name%in%Names.3$Name,-3]

# add empty coilumn ffor untyped family members
# (that are present in ped-file )

# voeg founder (vader) toe

n.snp<-dim(gt.1)[1]
#KIV_000025<-matrix(NA,nrow=n.snp,ncol=1)
#KIV_000028<-matrix(NA,nrow=n.snp,ncol=1)
VASS_11 <-matrix(NA,nrow=n.snp,ncol=1)
VASS_12 <-matrix(NA,nrow=n.snp,ncol=1)
VASS_21 <-matrix(NA,nrow=n.snp,ncol=1)
VASS_25 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_29 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_31 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_33 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_35 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_36 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_37 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_38 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_39 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_41 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_42 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_210 <-matrix(NA,nrow=n.snp,ncol=1)
VASS_211 <-matrix(NA,nrow=n.snp,ncol=1) #
VASS_212 <-matrix(NA,nrow=n.snp,ncol=1)
VASS_310 <-matrix(NA,nrow=n.snp,ncol=1)



#empty<-as.data.frame(cbind(KIV_000025,KIV_000028))
empty<-as.data.frame(cbind(VASS_11, VASS_12, VASS_21, VASS_25, VASS_29,VASS_31,VASS_33,VASS_35,VASS_36,VASS_37,VASS_38,VASS_39,VASS_41,VASS_42,VASS_210,VASS_211,VASS_212,VASS_310))
#names(empty)<-c("KIV_000025","KIV_000028")
#names(empty)<-c("VASS_11","VASS_12","VASS_21","VASS_24","VASS_26","VASS_27","VASS_210","VASS_212","VASS_311")
names(empty)<-c("VASS_11", "VASS_12","VASS_21","VASS_25","VASS_29","VASS_31","VASS_33","VASS_35","VASS_36","VASS_37", "VASS_38","VASS_39","VASS_41","VASS_42","VASS_210","VASS_211", "VASS_212","VASS_310")
gt.2<-data.frame(gt.1,empty,stringsAsFactors=F)


# moet getransposeerd
library(reshape2)

#temp<-melt.data.frame(gt.2,id.var=c("Name","Chr"))

#gt.melt<-melt.data.frame(gt.2,id.var=c("Name","Chr"))
gt.melt <- melt(gt.2,id.var=c("Name","Chr"))
#gt.melt$ID<-(substr(gt.melt$variable,5,10)) #as.numeric
#gt.melt$ID<- unlist(strsplit(gt.melt$variable,".GType")) #as.numeric
#tmp = apply(gt.melt,1,func(x) { a = strsplit(paste("\"",x,"\"",sep=""),".GType"); return(a)}
tmp = apply(as.matrix(gt.melt$variable),1,function(x){ a1 = unlist(strsplit(x,".GType")); return(a1) })
gt.melt$ID <- tmp

a = names(table(gt.melt$ID))
#gt.melt.1 = gt.melt[gt.melt$ID%in%a[c(1,2,12:19)],] 
gt.melt.1 = gt.melt[gt.melt$ID%in%a,]
gt.melt = gt.melt.1


# recodeer genotype naar 0.1.2
gt.melt$value<-as.numeric(as.factor(gt.melt$value))-1

# genotypes opsplitsen per chromosoon
# per chromosoom naar breed formaat
n.chr<-22
genotypes<-vector("list",n.chr)


for (i in 1:n.chr){
temp<-gt.melt[gt.melt$Chr==i,-c(2,3)]
genotypes[[i]]<-cast(temp,ID~Name)
}



# MERLIN : verdergaan met input files per chromosoom

# voor 1 chromosoom (ongesegmenteerd): 
#########################################
# voor elk chromosoom, genotypes omzetten in 2 allele-kolommen

n.chr<-22
snpnames<-vector("list",n.chr)

for (i in 1:n.chr){
snpnames[[i]]<-names(genotypes[[i]])[-1]
}

# creer 2 kolommen voor elk genotype
# elke kolom = allel
# buitenste loop : chromosomen (elementen van de lijst)
# binnenste loops : eerst opvullen met allelen,
# dan kolom headers
alleles<-vector("list",n.chr)
n.snp<-rep(NA,n.chr)

for (i in 1:n.chr){
	n.snp[i]<-dim(genotypes[[i]])[2]-1
	alleles[[i]]<-matrix(NA,nrow=nrow(genotypes[[i]]),ncol=n.snp[i]*2)

    	for (m in 1:n.snp[i]){
    		j<-m*2
    		k<-m*2-1
    		l<-m+1
    		alleles[[i]][,k]<-genotypes[[i]][,l]+1-(genotypes[[i]][,l]>0)
    		alleles[[i]][,j]<-genotypes[[i]][,l]+1-(genotypes[[i]][,l]>1)
    	}

	alleles[[i]]<-as.data.frame(alleles[[i]])    
# kolom names
	for (m in 1:n.snp[i]){
    		j<-m*2
    		k<-m*2-1
    		names(alleles[[i]])[k]<-paste(snpnames[[i]][m],".1",sep="")
    		names(alleles[[i]])[j]<-paste(snpnames[[i]][m],".2",sep="")
    	}
}


# mergen van genotypes met eerste zes kolommen
# rijen staan in zelfde volgorde
# dus cbind en geen merge

prefile<-vector("list",n.chr)

for (i  in 1:n.chr){
	prefile[[i]]<-data.frame(pre.6,alleles[[i]])
	write.table(prefile[[i]],file=paste("/home/aakumar/WORK/LINK_ANAL/DVD_BRS/VASS_inputfiles/pedin.",i,sep=""),sep="\t", quote=F,row.names=F,col.names=F,na="0")
}

system("cp -r /home/aakumar/WORK/LINK_ANAL/DVD_BRS/VASS_inputfiles  /home/aakumar/WORK/LINK_ANAL/DVD_BRS/ORIG_VASS_inputfiles")
system("cp /home/aakumar/WORK/LINK_ANAL/DVD_BRS/PAR/* /home/aakumar/WORK/LINK_ANAL/DVD_BRS/VASS_inputfiles/")
