
library(snpStats)
library(ggplot2)

# Code requires manuel input and saving of dataframes below#

##################################Input files###################################

chr01=readRDS("ind_m45_g50_1_5mb_res.RData")
chr02=readRDS("ind_m45_g50_2_5mb_res.RData")
chr03=readRDS("ind_m45_g50_3_5mb_res.RData")
chr04=readRDS("ind_m45_g50_4_5mb_res.RData")
chr05=readRDS("ind_m45_g50_5_5mb_res.RData")
chr06=readRDS("ind_m45_g50_6_5mb_res.RData")
chr07=readRDS("ind_m45_g50_7_5mb_res.RData")
chr08=readRDS("ind_m45_g50_8_5mb_res.RData")
chr09=readRDS("ind_m45_g50_9_5mb_res.RData")
chr10=readRDS("ind_m45_g50_10_5mb_res.RData")
chr11=readRDS("ind_m45_g50_11_5mb_res.RData")
chr12=readRDS("ind_m45_g50_12_5mb_res.RData")
chr13=readRDS("ind_m45_g50_13_5mb_res.RData")
chr14=readRDS("ind_m45_g50_14_5mb_res.RData")
chr15=readRDS("ind_m45_g50_15_5mb_res.RData")
chr16=readRDS("ind_m45_g50_16_5mb_res.RData")
chr17=readRDS("ind_m45_g50_17_5mb_res.RData")
chr18=readRDS("ind_m45_g50_18_5mb_res.RData")
chr19=readRDS("ind_m45_g50_19_5mb_res.RData")
chr20=readRDS("ind_m45_g50_20_5mb_res.RData")
chr21=readRDS("ind_m45_g50_21_5mb_res.RData")
chr22=readRDS("ind_m45_g50_22_5mb_res.RData")

genelist=read.table("BiomartProteinCodingGeneList_GRCh37_160117.bed")

nrdims=sum(length(chr01[1,1,]),length(chr02[1,1,]),length(chr03[1,1,]),length(chr04[1,1,])
           ,length(chr05[1,1,]),length(chr06[1,1,]),length(chr07[1,1,]),length(chr08[1,1,])
           ,length(chr09[1,1,]),length(chr10[1,1,]),length(chr11[1,1,]),length(chr12[1,1,])
           ,length(chr13[1,1,]),length(chr14[1,1,]),length(chr15[1,1,]),length(chr16[1,1,])
           ,length(chr17[1,1,]),length(chr18[1,1,]),length(chr19[1,1,]),length(chr20[1,1,])
           ,length(chr21[1,1,]),length(chr22[1,1,]))

chr=array(data=c(chr01,chr02,chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10
                 ,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21
                 ,chr22)
          ,dim=c(length(chr01[,1,1]),length(chr01[1,,1]),nrdims))


################################Input choices###################################

chri=chr22
chr.nr=22

plink.gen=read.plink("individuals_m45_g50_22.bed",
                     "individuals_m45_g50_22.bim",
                     "individuals_m45_g50_22.fam")

total.gen.length.from=min(plink.gen$map$position)
total.gen.length.to=max(plink.gen$map$position)
block.length=5000000
min.snps.in.block=100

genlist_chr1=subset(genelist, genelist$V1=="chr22")

###################


plink.snps=plink.gen$map$position
snp.order=sort(plink.snps)

blocks=split(snp.order, cut(plink.gen$map$position, seq(total.gen.length.from, 
                                                        total.gen.length.to, 
                                                        block.length)))
bsize=sapply(blocks, length)
blocks=blocks[bsize>=min.snps.in.block]
block.nr=length(blocks)
b.names=names(blocks)


startsin=split(genlist_chr1$V4, cut(genlist_chr1$V2, seq(total.gen.length.from, 
                                                         total.gen.length.to, 
                                                         block.length)))

endsin=split(genlist_chr1$V4, cut(genlist_chr1$V3, seq(total.gen.length.from, 
                                                       total.gen.length.to, 
                                                       block.length)))

allblocks=array(data=NA,dim=c(length(startsin),5,1))
allblocks[,1,1]=names(startsin)
for(i in 1:length(startsin)){
  allblocks[i,2,1]=total.gen.length.from+(block.length*(i-1))
  allblocks[i,3,1]=as.numeric(allblocks[i,2,1])+block.length
}

for (i in 1:length(allblocks[,1,1])){
  for (j in 1:length(startsin)){
    if (allblocks[i,1,1]==names(startsin[j])){
      allblocks[i,4,1]=toString(union(unlist(endsin[j]),unlist(startsin[j])))
    }
  }
}

for (i in 1:length(allblocks[,1,1])){
  for (j in 1:length(bsize)){
    if (allblocks[i,1,1]==names(bsize[j])){
      allblocks[i,5,1]=noquote(bsize[[j]])
    }
  }
}

chr1res=array(data=NA,dim=c(block.nr,16,1))
chr1res=as.data.frame(chr1res)
colnames(chr1res)=c("chr","b.nr","posrough","pos.start","pos.end","SNP.nr","meanSFA","meanSFApos",
                    "meanCor","7500","7500Cor","5500","5500Cor","3500","3500Cor","genes")
chr1res[,1]=chr.nr
chr1res[,2]=seq(1,block.nr,1)
chr1res[,3]=noquote(b.names)

for (i in 1:nrow(chr1res)){
  for (j in 1:nrow(allblocks)){
    if (chr1res[i,3]==allblocks[j,1,1]){
      chr1res[i,4]=allblocks[j,2,1]
      chr1res[i,5]=allblocks[j,3,1]
      chr1res[i,16]=allblocks[j,4,1]
      chr1res[i,6]=allblocks[j,5,1]
    }
  }
}

for (i in 1:nrow(chr1res)){
  chr1res[i,7]=mean(chri[seq(1,41,1),5,i])
}


pctl=function(vector,value){
  out=((sum(value>vector)+(sum(value==vector)/2))/length(vector))*100
  return(out)
}

chr.w.ranks=array(data=NA,dim=c(41,nrow(chr1res),1))
for (i in 1:41){
  for (j in 1:nrow(chr1res)){
    chr.w.ranks[i,j,1]=pctl(chr[i,5,],chri[i,5,j])
  }
}

for (i in 1:nrow(chr1res)){
  chr1res[i,8]=mean(chr.w.ranks[,i,1])
}

for (i in 1:nrow(chr1res)){
  chr1res[i,9]=mean(chri[seq(1,41,1),6,i])
}

for (i in 1:nrow(chr1res)){
  chr1res[i,10]=pctl(chr[1,5,],chri[1,5,i])
}

for (i in 1:nrow(chr1res)){
  chr1res[i,11]=chri[1,6,i]
}

for (i in 1:nrow(chr1res)){
  chr1res[i,12]=pctl(chr[21,5,],chri[21,5,i])
}

for (i in 1:nrow(chr1res)){
  chr1res[i,13]=chri[21,6,i]
}

for (i in 1:nrow(chr1res)){
  chr1res[i,14]=pctl(chr[41,5,],chri[41,5,i])
}

for (i in 1:nrow(chr1res)){
  chr1res[i,15]=chri[41,6,i]
}



chrres1=chr1res
chrres2=chr1res
chrres3=chr1res
chrres4=chr1res
chrres5=chr1res
chrres6=chr1res
chrres7=chr1res
chrres8=chr1res
chrres9=chr1res
chrres10=chr1res
chrres11=chr1res
chrres12=chr1res
chrres13=chr1res
chrres14=chr1res
chrres15=chr1res
chrres16=chr1res
chrres17=chr1res
chrres18=chr1res
chrres19=chr1res
chrres20=chr1res
chrres21=chr1res
chrres22=chr1res

chr_5mb_res=rbind(chrres1,chrres2,chrres3,chrres4,chrres5,chrres6,chrres7,chrres8,
                  chrres9,chrres10,chrres11,chrres12,chrres13,chrres14,chrres15,chrres16,
                  chrres17,chrres18,chrres19,chrres20,chrres21,chrres22)
nrow(chr_5mb_res)

write.table(chr_5mb_res, file="5mb_table_TBCor_res.tsv", quote=FALSE, sep="\t",row.names=FALSE)

chr_5mb_res_round=chr_5mb_res

chr_5mb_res_round[,7]=round(chr_5mb_res_round[,7], digits=0)
chr_5mb_res_round[,8]=round(chr_5mb_res_round[,8], digits=0)
#chr_5mb_res_round[,9]=round(chr_5mb_res_round[,9], digits=0)
chr_5mb_res_round[,10]=round(chr_5mb_res_round[,10], digits=0)
#chr_5mb_res_round[,11]=round(chr_5mb_res_round[,11], digits=0)
chr_5mb_res_round[,12]=round(chr_5mb_res_round[,12], digits=0)
#chr_5mb_res_round[,13]=round(chr_5mb_res_round[,13], digits=0)
chr_5mb_res_round[,14]=round(chr_5mb_res_round[,14], digits=0)
#chr_5mb_res_round[,15]=round(chr_5mb_res_round[,15], digits=0)

write.table(chr_5mb_res_round, file="5mb_table_TBCor_res_round.tsv", quote=FALSE, sep="\t",row.names=FALSE)


#### SNP-Correlation Regression Analysis ####

mc=as.numeric(chr_5mb_res$meanCor)
snp=as.numeric(chr_5mb_res$SNP.nr)

mcsnp=data.frame(mc,snp)
head(mcsnp)

plot(m1)
abline(h=0)

m1=lm(mc~snp)
summary(m1)
m1$coefficients[1]
theme_set(theme_bw())
p=ggplot(data=mcsnp,aes(x=snp,y=mc))
p=p+geom_point(colour="royalblue")
p=p+geom_abline(slope=m1$coefficients[2],intercept=m1$coefficients[1],
                colour="indianred")
#p=p+geom_line(data=m1,colour="darkred",linetype="dashed")
#p=p+geom_vline(xintercept=3500,linetype="dotted",colour="darkred")
p=p+labs(y="Mean correlation of block",x="Number of SNPs in block",title="Regression analysis of correlation and number of SNPs",subtitle="5 Mb blocks")
p=p+theme(axis.title.x=element_text(margin=margin(t=10),size=25),
          axis.title.y=element_text(margin=margin(r=10),size=25),
          plot.title=element_text(margin=margin(r=10,0,0,0),size=40,hjust=0.5),
          plot.subtitle=element_text(hjust=0.5,size=20),
          axis.text.x=element_text(size=15),
          plot.margin=margin(t=10,r=10,b=10,l=10),
          axis.text.y=element_text(size=15))
p=p+scale_x_continuous(breaks=seq(0,5000,500),limits=c(0,5000))
p=p+scale_y_continuous(breaks=seq(0,0.2,0.05),limits=c(0,0.2))

print(p)

### Histogram of SNPs per block ### 

chr_5mb_res=read.table("5mb_table_TBCor_res.tsv",sep="\t",header=TRUE)

ggplot(chr_5mb_res,aes(x=as.numeric(SNP.nr)))+
  geom_histogram(binwidth=100,boundary=0,color="black",fill="royalblue")+
  labs(y="Number of blocks",x="Number of SNPs",title="Distribution of blocks by number of SNPs",
       subtitle="EUR 5 Mb blocks")+
  #coord_cartesian(xlim=c(100,5000),ylim=c(0,60))+
  theme(axis.title.x=element_text(margin=margin(t=10),size=25),
        axis.title.y=element_text(margin=margin(r=10),size=25),
        plot.title=element_text(margin=margin(r=10,0,0,5),size=35,hjust=0.5),
        plot.subtitle=element_text(hjust=0.5,size=25),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15))
  #scale_x_continuous(breaks=seq(0,5000,500),limits=c(-25,5025))+
  #scale_y_continuous(breaks=seq(0,60,5),limits=c(0,60))
range(as.numeric(chr_5mb_res$SNP.nr))


### histogram of number of genes per block ###

genenr5=array(data=NA,dim=c(nrow(chr_5mb_res),2,1))

for(i in 1:nrow(chr_5mb_res)){
  genenr5[i,1,1]=chr_5mb_res[i,3]
  genenr5[i,2,1]=lengths(gregexpr(",",chr_5mb_res[i,16]))+1
}

genenr5=as.data.frame(genenr5)
colnames(genenr5)=c("block","gene.nr")

genenr5$gene.nr=as.numeric(genenr5$gene.nr)
range(as.numeric(genenr5$gene.nr))

ggplot(genenr5,aes(x=as.numeric(gene.nr)))+
  geom_histogram(binwidth=10,boundary=0,color="black",fill="royalblue")+
  labs(y="Number of blocks",x="Number of genes",title="Distribution of blocks by number of genes",
       subtitle="EUR 5 Mb blocks")+
  #coord_cartesian(xlim=c(100,5000),ylim=c(0,60))+
  theme(axis.title.x=element_text(margin=margin(t=10),size=25),
        axis.title.y=element_text(margin=margin(r=10),size=25),
        plot.title=element_text(margin=margin(r=10,0,0,5),size=35,hjust=0.5),
        plot.subtitle=element_text(hjust=0.5,size=25),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15))+
scale_x_continuous(breaks=seq(0,200,20),limits=c(0,200))+
scale_y_continuous(breaks=seq(0,150,25),limits=c(0,150))










