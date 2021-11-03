
# Orbicella annularis structural complexity model for rates of multidecadal framework erosion (1998 - 2018)  
# parameterised using the field-data from Long Cay (Roff et al 2015 Global Change Biology, https://doi.org/10.1111/gcb.13006)
# g.roff@uq.edu.au [26th February 2019]


rm(list=ls()) # clear history
set.seed(101) # set seed


###### Loop model for 1000 simulations ----

df = NULL
for (k in 1:1000)
    {
     cat("model loop",k,"\n")


### # 1. Create cross section of O.annularis parameterised from field data ----

## Colony width derived from 95 colonies of Orbicella annularis measured at Long Cay in 2000 (Roff et al 2015 GCB)
cDiameter <- rnorm(1, mean = 78.85, sd = 30.98) # sample colony width from distribution
cDiameter[cDiameter <30] <- 30 # set minimum colony size at 20cm based on field data
cHeight <- (cDiameter*0.56) + 45.5 # convert colony width to height using equation from Roff et al (2015) GCB (Figure S2)

nrams<-round(cDiameter/(9.08 + 4.73)) # nrams is number of ramets within a colony (colony diameter divided by average ramet diameter and ramet spacing)

## Ramet diameters derived from 30 colonies of Orbicella annularis measured at Long Cay in 2000 (Roff et al 2015 GCB)
nDiameterR<-rnorm(round(nrams), mean = 9.08, sd = 3.15) # for the number of ramets in colony, randomly sample ramet diameter from distribution
nDiameterR[nDiameterR <2] <-2 # set minimum size at 2cm based on field data

## Ramet spacing derived from 50 colonies of Orbicella annularis measured at Long Cay in 1998
nSpacing <- rnorm(round(nrams), mean = 4.73, sd = 0.83) # for the number of ramets in colony, randomly sample ramet spacing from distribution
nSpacing[nSpacing <0.4] <- 0.4 # set minimum spacing at 0.4cm based on field data

nDiameter<-(nDiameterR+nSpacing)

# Ramet heights derived from 30 colonies of Orbicella annularis measured at Long Cay in 2000 (Roff et al 2015 GCB)
nHeight<-rnorm(nrams, mean = 6.77, sd = 2.7)  # for the number of ramets in colony, randomly sample ramet height from distribution
nHeight[nHeight <2.5] <-2.5 # set minimum size at 2.5cm

coral_list<-cbind(nDiameter,nHeight) # make matrix at colony scale of ramet diameters/spacing and ramet heights

# calculate rugosity of live coral (calculated as topographic "tooth-step" profile [see Figure 1, Figure S1])
topo<-sum(nDiameter)+(sum(nHeight)*2)+(cHeight*2) # sum the ramet diameters + spacings, colony height (*2)
rugosity<-topo/cDiameter # calculate rugosity as topography / colony diameter (Ri)

### 2. Determine  live/dead ramets within colonies prior to the 1998 bleaching event 

# calculate the number of live and dead ramets prior to the 1998 bleaching event
presurv <- (rnorm(1, mean = 0.977, sd = 0.086)) # average proportion of dead ramets determined from surveys in 1998 mean SD
presurv[presurv>1] <-1 # remove negative survival from rnorm distribution 
prelive<-round(presurv*nrams) # number of live ramets pre-1998
predead<-round(nrams)-prelive # number of dead ramets pre-1998

### 3. Simulate mortality following the 1998 bleaching event 

surv <- rnorm(1, mean = 0.141, sd = 0.11) # average survival determined from surveys in 2018 mean SD
surv[surv<0] <-0 # remove negative survival from rnorm distribution 
postrams<-round(surv*nrams) # number of surviving ramets post-1998
deadrams<-nrams-postrams # number of dead ramets post-1998

### 4. Define erosion and growth

# linear external bioerosion based on U-th measurements from Roff et al 2015 GCB
erosion <- rnorm(nrams, mean = -0.110329383, sd = 0.160802775) # calculate erosion per ramet based on random sampling of mean and SD from U-th [static = 2.2]
erosion[erosion>0] <-0 # remove positive erosion measurements stemming from mean SD
# linear extension (growth) based on CT scan measuerments of ramets 2006-2011 from Roff et al 2015 GCB
growth<- rnorm(nrams, mean = 0.686, sd = 0.127) # calculate growth per ramet based on random sampling of mean and SD from CT scans [static = 13.8]


### 5. Calculate multidecadal timesequence----

# create matrix at colony level and randomly sample growth and erosion from distributions in each timestep.

# pre-1998 assumes erosion of existing dead ramets (predead) and growth of living ramets (prelive) 
difference95<-nHeight+as.numeric(c(sample(erosion, predead, replace=TRUE),sample(growth, prelive, replace=TRUE)))
difference96<-difference95+as.numeric(c(sample(erosion, predead, replace=TRUE),sample(growth, prelive, replace=TRUE)))
difference97<-difference96+as.numeric(c(sample(erosion, predead, replace=TRUE),sample(growth, prelive, replace=TRUE)))
difference98<-difference97+as.numeric(c(sample(erosion, predead, replace=TRUE),sample(growth, prelive, replace=TRUE)))

# post-1998 mortality assumse growth of surviving ramets (postrams) and erosion of dead ramets (deadrams)
difference99<-difference98+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))
difference00<-difference99+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))
difference01<-difference00+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE))) 
difference02<-difference01+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference03<-difference02+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference04<-difference03+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference05<-difference04+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference06<-difference05+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference07<-difference06+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference08<-difference07+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference09<-difference08+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference10<-difference09+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference11<-difference10+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference12<-difference11+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference13<-difference12+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference14<-difference13+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference15<-difference14+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference16<-difference15+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference17<-difference16+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  
difference18<-difference17+as.numeric(c(sample(erosion, deadrams, replace=TRUE),sample(growth, postrams, replace=TRUE)))  

coral_list<-cbind(coral_list,difference95,difference96,difference97,difference98,difference99,difference00,difference01,difference02,difference03,difference04,difference05,difference06,difference07,difference08,difference09,
                  difference10,difference11,difference12,difference13,difference14,difference15,difference16,difference17,difference18)

coral_list<-as.data.frame(coral_list)


### 6. Calculate annual change rugosity at a colony scale (Ri)----

rugosity1995<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference95)*2))/cDiameter
rugosity1996<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference96)*2))/cDiameter
rugosity1997<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference97)*2))/cDiameter
rugosity1998<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference98)*2))/cDiameter
rugosity1999<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference98)*2))/cDiameter
rugosity2000<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference99)*2))/cDiameter
rugosity2001<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference00)*2))/cDiameter
rugosity2002<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference01)*2))/cDiameter
rugosity2003<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference02)*2))/cDiameter
rugosity2004<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference03)*2))/cDiameter
rugosity2005<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference04)*2))/cDiameter
rugosity2006<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference05)*2))/cDiameter
rugosity2007<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference06)*2))/cDiameter
rugosity2008<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference07)*2))/cDiameter
rugosity2009<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference08)*2))/cDiameter
rugosity2010<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference09)*2))/cDiameter
rugosity2011<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference10)*2))/cDiameter
rugosity2012<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference11)*2))/cDiameter
rugosity2013<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference12)*2))/cDiameter
rugosity2014<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference13)*2))/cDiameter
rugosity2015<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference14)*2))/cDiameter
rugosity2016<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference15)*2))/cDiameter
rugosity2017<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference16)*2))/cDiameter
rugosity2018<-((cHeight*2)+sum(coral_list$nDiameter)+(sum(coral_list$difference17)*2))/cDiameter



### 7. populate for output & end loop ----


df = rbind(df, data.frame(rugosity1995,rugosity1996,rugosity1997,rugosity1998,rugosity1999,rugosity2000,rugosity2001,rugosity2002,rugosity2003,
                          rugosity2004,rugosity2005,rugosity2006,rugosity2007,rugosity2008,rugosity2009,rugosity2010,rugosity2011,rugosity2012,
                          rugosity2013,rugosity2014,rugosity2015,rugosity2016,rugosity2017,rugosity2018,nrams,cDiameter,surv,prelive,predead,postrams,deadrams))


}



### 8. Plot annual changes in structural complexity for 1000 simulated colonies ----


dfmean<-colMeans(df)
dfmean<-as.data.frame(t(dfmean))
TimePoints=1996:2018
matplot(t(df[, c(2:24)]), lwd = 0.5,type="l", pch=1, col="#b3b3cc", xlab="Year",ylab="Rugosity index",xaxt = "n") # plots continuous
axis(1, at = 1:23, las=2, labels = paste(TimePoints))
lines(rowMeans(t(df[, c(2:24)])), col = "red", lwd = 2) 
abline(v=3)
  
mean(df$rugosity1998)
mean(df$rugosity2018)
sum(df$delta < 0)


### 9. Plot change in structural complexity with percent bleaching and colony size ----

df$delta<-(df$rugosity2018-df$rugosity1998)

par(mfrow=c(1,2))

cor.test(df$delta,df$cDiameter)
fm1 <- cor.test(df$delta,df$prop)
plot(df$delta~df$cDiameter, xlab="Colony width (cm)",  xlim=c(25,200), ylab=expression(paste(Delta,sep=" ", Ri)))
abline(h=0)
abline(lm(df$delta~df$cDiameter), col="red",lwd = 3)
text(120,0.8,format(round(fm1$estimate,2),nsmall=2), cex=1.4)


df$prop<-(df$surv*100)
cor.test(df$delta,df$prop)
fm2 <- cor.test(df$delta,df$prop)
plot(df$delta~df$prop, xlab="Percent survival",xlim=c(0,50), ylab=expression(paste(Delta,sep=" ", Ri)))
abline(h=0)
abline(lm(df$delta~df$prop), col="red",  lwd = 3)
text(12,0.8,format(round(fm2$estimate,2),nsmall=2), cex=1.4)


par(mfrow=c(1,1))




### 10. cross section model of ramets, ramet spacings and ramet heights
x<-c(0,0,1,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11)
y<-c(0,1,1,1,0.5,0.5,1,1,0.5,0.5,1,1,0.5,0.5,1,1,0.5,0.5,1,1,0.5,0.5,1,1,0)
plot(x,y, ylim=c(0,1.5), pch=NA, xaxt='n', yaxt='n',ann=FALSE, asp=5)
polygon(x,y, col="#ffeb99")

title("Structural model of Orbicella annularis", cex=1.5)

text(1.5,0.45,"nSpacing", cex=0.8)
text(3.5,0.45,"nSpacing", cex=0.8)
text(5.5,0.45,"nSpacing", cex=0.8)
text(7.5,0.45,"nSpacing", cex=0.8)
text(9.5,0.45,"nSpacing", cex=0.8)

text(0.5,1.05,"nDiameterR", cex=0.8)
text(2.5,1.05,"nDiameterR", cex=0.8)
text(4.5,1.05,"nDiameterR", cex=0.8)
text(6.5,1.05,"nDiameterR", cex=0.8)
text(8.5,1.05,"nDiameterR", cex=0.8)
text(10.5,1.05,"nDiameterR", cex=0.8)

text(1.25,0.75,"nHeight", srt=270, cex=0.8) 
text(3.25,0.75,"nHeight", srt=270, cex=0.8) 
text(5.25,0.75,"nHeight", srt=270, cex=0.8) 
text(7.25,0.75,"nHeight", srt=270, cex=0.8) 
text(9.25,0.75,"nHeight", srt=270, cex=0.8) 

text(5.5,0.05,"cWidth") 

text(2,1.5,"nHeight = ramet height") 
text(2.5,1.44,"nDiameterR = ramet diameter") 
text(2.2,1.37,"nSpacing = ramet spacing") 
text(1.9,1.3,"cWidth = colony width") 

###### 

