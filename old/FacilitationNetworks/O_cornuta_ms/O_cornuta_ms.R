###################################################
### chunk number 1: 
###################################################
#line 78 "O_cornuta_ms.Rnw"
  library(gdata)
library(bipartite)
library(xtable)
                                        #Use the response statistic = (I-O)/(I+O)
rs=function(x){if ((x[1]-x[2])==0){0}else{(x[1]-x[2])/(sum(x))}}

                                        #Presence absence function
bin.sum=function(x){x[x!=0]<-1;sum(x)}

                                        #Binomial distibution test for the significance of each response statistic
rs.binom=function(x,alpha=0.05){
  total=sum(x)
  p=sum(dbinom(max(x):total,total,0.5))+sum(dbinom(0:min(x),total,0.5))
  if (p<=alpha){x=x}else if (p>alpha){x=c(0,0)}
  return(x)
}


                                        #Calculate response statistic
rii=function(x,rep,site,alpha=0.05,decreasing=TRUE,b.test=TRUE,label='n',sep='',zero.na=TRUE){
  
  y <- array(NA,c(length(unique(rep)),ncol(x)))
  y <- data.frame(y)
  
  for (i in seq(along=unique(rep))){
    z=x[rep==unique(rep)[i],]
    z=z[order(site[rep==unique(rep)[i]],decreasing=decreasing),]
    for (j in seq(along=(1:ncol(z)))){
      
      if (b.test == TRUE){z[,j]=rs.binom(z[,j],alpha=(alpha/ncol(z)))}else{}
      
    }
    y[i,]=apply(z,2,rs)
    rownames(y)[i] <- paste(label,unique(rep)[i],sep=sep)
  }
  
  colnames(y) <- colnames(x)
  if (zero.na == TRUE){y[is.na(y)] <- 0}
  
  return(y)	
}



###################################################
### chunk number 2: 
###################################################
#line 128 "O_cornuta_ms.Rnw"
  Oc.com <- read.csv('data/OnobrychisCushions.csv')
Oc.com <- Oc.com[,7:ncol(Oc.com)]
Oc.env <- read.csv('data/OnobrychisCushions.csv')[,1:2]

attach(Oc.env)
names(Oc.env)

                                        #Calculate response statistic (rii)
Oc.rs <- rii(Oc.com,rep,microsite,b.test=FALSE,label='Oc') #Without Test of significance
Oc.rs.. <- rii(Oc.com,rep,microsite,b.test=TRUE,label='Oc') #Test of significance (..)
Oc.order <- Oc.rs
Oc.order.. <- Oc.rs..

detach(Oc.env)

                                        #Import and sort the surface area values
surface <- read.csv('data/OnobrychisCushions.csv')[Oc.env$microsite=='on cushion',5]
Oc.surf <- surface[order(apply(Oc.order,1,bin.sum),decreasing=TRUE)]
Oc.surf.. <- surface[order(apply(Oc.order..,1,bin.sum),decreasing=TRUE)]
names(Oc.surf) <- rownames(Oc.rs)
names(Oc.surf..) <- rownames(Oc.rs..)

                                        #Reorder the matrices by their rii
                                        #Without test
Oc.binary <- Oc.rs
Oc.binary[Oc.binary != 0] <- 1
Oc.order <- Oc.rs
Oc.order <- Oc.order[order(apply(Oc.order,1,bin.sum),decreasing=TRUE),]
Oc.order <- Oc.order[,order(apply(Oc.order,2,bin.sum),decreasing=TRUE)]
Oc.order[rownames(Oc.order) == 'Oc1',order(colnames(Oc.order))] - Oc.rs[rownames(Oc.rs) == 'Oc1',order(colnames(Oc.rs))]
                                        #With test
Oc.binary.. <- Oc.rs..
Oc.binary..[Oc.binary.. != 0] <- 1
Oc.order.. <- Oc.rs..
Oc.order.. <- Oc.order..[order(apply(Oc.order..,1,bin.sum),decreasing=TRUE),]
Oc.order.. <- Oc.order..[,order(apply(Oc.order..,2,bin.sum),decreasing=TRUE)]
Oc.order..[rownames(Oc.order..) == 'Oc1',order(colnames(Oc.order..))] - Oc.rs..[rownames(Oc.rs..) == 'Oc1',order(colnames(Oc.rs..))]

##Cushion phenotype groupings based on the numbers of individuals (another dataset suggests a strong pattern of correlation between the phenotypes Dense, Loose and Intermediate and the number of individuals). The number of individuals also correlates strongly with the number of flowers produced (more beneficiary individuals means fewer flowers.)

#####SORT THE PHENOTYPE ACCORDING TO THE SORTING OF THE CUSHIONS!!!!!

                                        #Without test
pheno <- read.xls('data/Affectation phenotype.xls')
Oc.pheno <- as.character(sapply(rownames(Oc.order),function(x) unlist(strsplit(x,split='c'))[2]))
marker <- unique(Oc.pheno)
for (i in 1:length(marker)){
  Oc.pheno[Oc.pheno == marker[i]] <- as.character(pheno[,2])[pheno[,1] == marker[i]]
}
                                        #With test
pheno.. <- read.xls('data/Affectation phenotype.xls')
Oc.pheno.. <- as.character(sapply(rownames(Oc.order..),function(x) unlist(strsplit(x,split='c'))[2]))
marker.. <- unique(Oc.pheno..)
for (i in 1:length(marker..)){
  Oc.pheno..[Oc.pheno.. == marker..[i]] <- as.character(pheno..[,2])[pheno..[,1] == marker..[i]]
}




###################################################
### chunk number 3: 
###################################################
#line 190 "O_cornuta_ms.Rnw"
  par(mfrow=c(1,2))
                                        #Histograms
hist(Oc.order[Oc.order != 0],xlab='RII',main='All RII')
hist(Oc.order..[Oc.order.. != 0],xlab='RII',main='Significant RII')



###################################################
### chunk number 4: 
###################################################
#line 198 "O_cornuta_ms.Rnw"
                                        #Without test
Oc.bin <- as.matrix(Oc.order)
Oc.bin[Oc.bin != 0] <- 1
Oc.temp <- nestedtemp(abs(as.matrix(Oc.bin)))
nestedchecker(Oc.order)
Oc.nest <- oecosimu(as.matrix(Oc.order),nestedchecker,'quasiswap',nsimul=1000,alternative='less',burnin=0)
                                        #With test
Oc.bin.. <- as.matrix(Oc.order..)
Oc.bin..[Oc.bin.. != 0] <- 1
Oc.temp.. <- nestedtemp(abs(as.matrix(Oc.bin..)))
nestedchecker(Oc.order..)
Oc.nest.. <- oecosimu(as.matrix(Oc.order),nestedchecker,'quasiswap',nsimul=1000,alternative='less',burnin=0)

Oc.nest
Oc.nest..



###################################################
### chunk number 5: 
###################################################
#line 218 "O_cornuta_ms.Rnw"
  par(mfrow=c(2,2))
                                        #Nestedness Plots
plot(Oc.temp,kind='temp',col=grey(seq(1,0,length=25)),main='All RII')
plot(Oc.temp..,kind='temp',col=grey(seq(1,0,length=25)),main='Significant RII')
                                        #Bipartite Graph
plotweb(abs(Oc.order),method='normal',text.rot=90,empty=FALSE,col.low=as.numeric(factor(Oc.pheno)))
plotweb(abs(Oc.order..),method='normal',text.rot=90,empty=FALSE,col.low=as.numeric(factor(Oc.pheno..)))



###################################################
### chunk number 6: 
###################################################
#line 230 "O_cornuta_ms.Rnw"
  par(mfrow=c(1,2))
                                        #Without test
Oc.edge <- apply(Oc.order,2,function(x) length(x[x!=0]))
Oc.abun <- apply(Oc.com,2,sum)
Oc.edge <- Oc.edge[order(names(Oc.edge))]
Oc.abun <- Oc.abun[order(names(Oc.abun))]
Oc.edge.cush <- apply(Oc.order,1,function(x) length(x[x!=0]))
Oc.edge.cush <- Oc.edge.cush[order(names(Oc.edge.cush))]
names(Oc.surf) <- rownames(Oc.rs)
Oc.surf <- Oc.surf[order(names(Oc.surf))]

                                        #With test
Oc.edge.. <- apply(Oc.order..,2,function(x) length(x[x!=0]))
Oc.abun.. <- apply(Oc.com,2,sum)
Oc.edge.. <- Oc.edge..[order(names(Oc.edge..))]
Oc.abun.. <- Oc.abun..[order(names(Oc.abun..))]
Oc.edge.cush.. <- apply(Oc.order..,1,function(x) length(x[x!=0]))
Oc.edge.cush.. <- Oc.edge.cush..[order(names(Oc.edge.cush..))]
names(Oc.surf..) <- rownames(Oc.rs..)
Oc.surf <- Oc.surf[order(names(Oc.surf))]




###################################################
### chunk number 7: 
###################################################
#line 255 "O_cornuta_ms.Rnw"
  par(mfrow=c(2,2))
                                        #Without test

plot(Oc.edge~Oc.abun,pch=19,xlab='Species Total Abundance (in and out)',ylab='Number of Connections to Cushions',main='All RII')
abline(lm(Oc.edge~Oc.abun))
plot(Oc.edge.cush~Oc.surf,pch=19,xlab='Cushion Surface Area',ylab='Number of Connections to Species',main='All RII')
abline(lm(Oc.edge.cush~Oc.surf))
                                        #With test

plot(Oc.edge..~Oc.abun..,pch=19,xlab='Species Total Abundance (in and out)',ylab='Number of Connections to Cushions',main='Significant RII')
abline(lm(Oc.edge..~Oc.abun..))
plot(Oc.edge.cush..~Oc.surf..,pch=19,xlab='Cushion Surface Area',ylab='Number of Connections to Species',main='Significant RII')
abline(lm(Oc.edge.cush..~Oc.surf..))



###################################################
### chunk number 8: 
###################################################
#line 272 "O_cornuta_ms.Rnw"
  ##Test of difference in number of connections
  tapply(Oc.edge.cush,Oc.pheno,sum)
tapply(Oc.edge.cush..,Oc.pheno..,sum)
xtable(summary(aov(Oc.edge.cush~Oc.pheno)))
xtable(summary(aov(Oc.edge.cush..~Oc.pheno..)))



###################################################
### chunk number 9: 
###################################################
#line 281 "O_cornuta_ms.Rnw"
#Connections by phenotype
par(mfrow=c(1,2))
barplot(tapply(Oc.edge.cush,Oc.pheno,mean),xlab='Mean Number of Interactions')
barplot(tapply(Oc.edge.cush..,Oc.pheno..,mean),xlab='Mean Number of Interactions')



