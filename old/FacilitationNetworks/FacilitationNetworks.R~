#Modeling facilitation networks

#The following data are from Richard Michalet's sites in Lebenon. Datasets consist of species counts of plant species from within and next to cushions of two species (Astragalus Zachlensis and Onobrychis cornuta), with the sampling area standardized. There is also environmental data from sites in the original .xls file. 

library(bipartite)

#Use the response statistic = (I-O)/(I+O)
rs=function(x){if ((x[1]-x[2])==0){0}else{(x[1]-x[2])/(sum(x))}}

#Presence absence function
bin.sum=function(x){x[x!=0]<-1;sum(x)}

#response significance test
#returns the expected value given the significance test
#rs.binom=function(x,alpha=0.05){
#	s=x[1]
#	total=sum(x)
#	if (x[1]>=x[2]){
#		p=sum(dbinom(s:total,total,0.5))
#		}
#	else {
#		p=1-sum(dbinom(s:total,total,0.5))
#		}
#	if (p<=alpha){x=x}else{x=c(0,0)}
#	return(x)
#	}

rs.binom=function(x,alpha=0.05){
	total=sum(x)
	p=sum(dbinom(max(x):total,total,0.5))+sum(dbinom(0:min(x),total,0.5))
	if (p<=alpha){x=x}else if (p>alpha){x=c(0,0)}
	return(x)
	}


#Calculate response statistic
rii=function(x,rep,site,alpha=0.05,decreasing=TRUE,test=TRUE){
	
	y=array(NA,c(length(unique(rep)),ncol(x)))

for (i in seq(along=unique(rep))){
	z=x[rep==unique(rep)[i],]
	z=z[order(site[rep==unique(rep)[i]],decreasing=decreasing),]
for (j in seq(along=(1:ncol(z)))){

  if (test == TRUE){z[,j]=rs.binom(z[,j],alpha=(alpha/ncol(z)))
                  }else{z[,j] = z[,j]}
        
	}
	y[i,]=apply(z,2,rs)
	}
y[is.na(y)]<-0
return(y)	
	}

#A. Zachlensis

Az.com=read.csv('/Users/Aeolus/Documents/Active_Projects/FacilitationNetworks/AstragalusCushions.csv')
Az.com=Az.com[,7:ncol(Az.com)]

Az.env=read.csv('/Users/Aeolus/Documents/Active_Projects/FacilitationNetworks/AstragalusCushions.csv')[,1:2]
attach(Az.env)
names(Az.env)

#Calculate response statistics
Az.rs=rii(Az.com,rep,microsite)

rownames(Az.rs)=paste('Az',seq(along=unique(rep)),sep='')
colnames(Az.rs)=colnames(Az.com)
detach(Az.env)

Az.binary=Az.rs
Az.binary[Az.binary!=0]<-1
Az.order=Az.rs
Az.order=Az.order[order(apply(Az.order,1,bin.sum),decreasing=TRUE),]
Az.order=Az.order[,order(apply(Az.order,2,bin.sum),decreasing=TRUE)]

Az.order[rownames(Az.order)=='Az21',order(colnames(Az.order))]-Az.rs[rownames(Az.rs)=='Az21',order(colnames(Az.rs))]

#Graph
plotweb(abs(Az.rs))
plotweb(abs(Az.order),method='normal')

#Distribution patterns
hist(Az.rs)
edge.num=apply(Az.binary,1,sum)
hist(edge.num)

edge.num=apply(Az.binary,2,sum)
edge.hist=hist(edge.num)

plot(log(edge.hist$counts+1)~log(edge.hist$mids))
abline(lm(log(edge.hist$counts+1)~log(edge.hist$mids)))

#Is there a bias in the sampling (i.e., size of cushion) that could affect the network structures?


#Nestedness
visweb(Az.binary,type='nested')

#O. cornuta

Oc.com=read.csv('/Users/Aeolus/Documents/Active_Projects/FacilitationNetworks/OnobrychisCushions.csv')
Oc.com=Oc.com[,7:ncol(Oc.com)]

Oc.env=read.csv('/Users/Aeolus/Documents/Active_Projects/FacilitationNetworks/OnobrychisCushions.csv')[,1:2]
attach(Oc.env)
names(Oc.env)

#Calculate response statistic
Oc.rs=rii(Oc.com,rep,microsite)
Oc.order=Oc.rs
rownames(Oc.rs)=paste('Oc',seq(along=unique(rep)),sep='')
colnames(Oc.rs)=colnames(Oc.com)
detach(Oc.env)

surface=read.csv('/Users/Aeolus/Documents/Active_Projects/FacilitationNetworks/OnobrychisCushions.csv')[Oc.env$microsite=='on cushion',5]
surface=surface[order(apply(Oc.order,1,bin.sum),decreasing=TRUE)]

Oc.binary=Oc.rs
Oc.binary[Oc.binary!=0]<-1
Oc.order=Oc.rs
Oc.order=Oc.order[order(apply(Oc.order,1,bin.sum),decreasing=TRUE),]
Oc.order=Oc.order[,order(apply(Oc.order,2,bin.sum),decreasing=TRUE)]

Oc.order[rownames(Oc.order)=='Oc1',order(colnames(Oc.order))]-Oc.rs[rownames(Oc.rs)=='Oc1',order(colnames(Oc.rs))]


#Graph
plotweb(abs(Oc.rs))
plotweb(abs(Oc.order),method='normal')


#Distribution patterns
hist(Oc.rs[Oc.rs!=0])
edge.num=apply(Oc.binary,1,sum)
hist(edge.num)

edge.num=apply(Oc.binary,2,sum)
edge.hist=hist(edge.num)

plot(log(edge.hist$counts+1)~log(edge.hist$mids))
abline(lm(log(edge.hist$counts+1)~log(edge.hist$mids)))

#Is there a bias in the sampling (i.e., size of cushion) that could affect the network structures?
edge.num=apply(Oc.binary,1,sum)
summary(edge.surf<-lm(edge.num~surface))
summary(edge.polysurf<-lm(edge.num~poly(surface,2)))
plot(surface,edge.num)
abline(edge.surf)
lines(spline(surface,predict(edge.polysurf)),lty=2)

AIC(edge.surf);AIC(edge.polysurf)

summary(log.e.s<-lm(log(edge.num+1)~log(surface)))
plot(log(surface),log(edge.num+1))
abline(log.e.s)
hist(residuals(log.e.s))

#Nestedness
visweb(Oc.binary,type='nested')

heatmap(Az.order)
#Compare
plotweb(abs(Az.order),method='normal')
quartz()
plotweb(abs(Oc.order),method='normal')

visweb(Az.binary,type='nested')
quartz()
visweb(Oc.binary,type='nested')


#Separate Out Positive and Negative Interactions

#A. Zachalensis
Az.fac=Az.order
Az.fac[Az.fac<0]<-0
Az.fac=Az.fac[order(apply(Az.fac,1,bin.sum),decreasing=TRUE),order(apply(Az.fac,2,bin.sum),decreasing=TRUE)]
Az.com=Az.order
Az.com[Az.com>0]<-0
Az.com=Az.com[order(apply(Az.com,1,bin.sum),decreasing=TRUE),order(apply(Az.com,2,bin.sum),decreasing=TRUE)]
Az.f.bin=Az.fac
Az.f.bin[Az.f.bin!=0]<-1

Az.c.bin=Az.com
Az.c.bin[Az.c.bin!=0]<-1

quartz('A. Zachalensis Beneficiaries')
plotweb(abs(Az.fac),method='normal')
quartz('A. Zachalensis Excluded')
plotweb(abs(Az.com),method='normal')

quartz('A. Zachalensis Beneficiaries')
visweb(Az.f.bin)
quartz('A. Zachalensis Excluded')
visweb(Az.c.bin)

#O. cornuta
Oc.fac=Oc.order
Oc.fac[Oc.fac<0]<-0
Oc.fac=Oc.fac[order(apply(Oc.fac,1,bin.sum),decreasing=TRUE),order(apply(Oc.fac,2,bin.sum),decreasing=TRUE)]
Oc.com=Oc.order
Oc.com[Oc.com>0]<-0
Oc.com=Oc.com[order(apply(Oc.com,1,bin.sum),decreasing=TRUE),order(apply(Oc.com,2,bin.sum),decreasing=TRUE)]
Oc.f.bin=Oc.fac
Oc.f.bin[Oc.f.bin!=0]<-1

Oc.c.bin=Oc.com
Oc.c.bin[Oc.c.bin!=0]<-1

quartz('O. cornuta Beneficiaries')
plotweb(abs(Oc.fac),method='normal')
quartz('O. cornuta Excluded')
plotweb(abs(Oc.com),method='normal')

quartz('O. cornuta Beneficiaries')
visweb(Oc.f.bin)
quartz('O. cornuta Excluded')
visweb(Oc.c.bin)




###GEUM###
#Geum low site
GLow=read.csv('/Users/Aeolus/Documents/Active_Projects/FacilitationNetworks/GeumData/GeumLow.csv')
head(GLow)

GL.com=GLow[,7:ncol(GLow)]
head(GL.com)

attach(GLow)

#Calculate response statistic
levels(Cushion.number)=c('IN','OUT')
GL.rs=rii(GL.com,Number,Cushion.number,alpha=0.2,decreasing=FALSE)

rownames(GL.rs)=paste('GL',seq(along=unique(Number)),sep='')
colnames(GL.rs)=colnames(GL.com)
detach(GLow)

GL.surf=seq(along=unique(GLow$Number))
for (i in seq(along=unique(GLow$Number))){
	GL.surf[i]=GLow$Surface[GLow$Number==unique(GLow$Number)[i]][1]
	}
GL.surf=GL.surf[order(apply(GL.rs,1,bin.sum),decreasing=TRUE)]

GL.rs=GL.rs[order(apply(GL.rs,1,bin.sum),decreasing=TRUE),]
GL.rs=GL.rs[,order(apply(GL.rs,2,bin.sum),decreasing=TRUE)]


plotweb(abs(GL.rs),method='normal')

#Geum High site

GHigh=read.csv('/Users/Aeolus/Documents/Active_Projects/FacilitationNetworks/GeumData/GeumHigh.csv')
head(GHigh)

GH.com=GHigh[,7:ncol(GHigh)]
head(GH.com)

attach(GHigh)

#Calculate response statistic
levels(Cushion.number)=c('IN','OUT')
GH.rs=rii(GH.com,Number,Cushion.number,alpha=0.1,decreasing=FALSE)

rownames(GH.rs)=paste('GH',seq(along=unique(Number)),sep='')
colnames(GH.rs)=colnames(GH.com)
detach(GHigh)

GH.surf=seq(along=unique(GHigh$Number))
for (i in seq(along=unique(GHigh$Number))){
	GH.surf[i]=GHigh$Surface[GHigh$Number==unique(GHigh$Number)[i]][1]
	}
GH.surf=GH.surf[order(apply(GH.rs,1,bin.sum),decreasing=TRUE)]

GH.rs=GH.rs[order(apply(GH.rs,1,bin.sum),decreasing=TRUE),]
GH.rs=GH.rs[,order(apply(GH.rs,2,bin.sum),decreasing=TRUE)]

plotweb(abs(GH.rs),method='normal')

#Node Stats
GL.ew=apply(GL.rs,1,bin.sum)
GL.LMfit=lm(GL.ew~GL.surf)
plot(GL.ew~GL.surf)
abline(GL.LMfit)
summary(GL.LMfit)

GH.ew=apply(GH.rs,1,bin.sum)
GH.LMfit=lm(GH.ew~GH.surf)
plot(GH.ew~GH.surf)
abline(GH.LMfit)
summary(GH.LMfit)

#ANCOVA
macsite=rep('L',length=length(GL.surf))
ew=GL.ew
surf=GL.surf
x1=data.frame(macsite,ew,surf)
macsite=rep('H',length=length(GH.surf))
ew=GH.ew
surf=GH.surf
x2=data.frame(macsite,ew,surf)
x=data.frame(rbind(x1,x2))

summary(ancova.fit<-aov(log(ew+1)~log(surf)*macsite,data=x))
plot(log(ew+1)~log(surf),data=x,col=as.numeric(macsite))
abline(lm(log(GL.ew+1)~log(GL.surf)),col=1)
abline(lm(log(GH.ew+1)~log(GH.surf)),col=2)


#Nestedness

GL.fac=GL.rs
GL.fac[GL.fac<0]<-0

GH.fac=GH.rs
GH.fac[GH.fac<0]<-0

#Modularity
quartz(title='Low')
plotweb(GL.fac,method='normal')
quartz(title='High')
plotweb(GH.fac,method='normal')
