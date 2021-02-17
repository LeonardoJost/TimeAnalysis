### Analysis of treatment using random sampling of participants
#     Copyright (C) 2021  Leonardo Jost
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(lme4)
library(optimx)
source("functions/helpers.R")
#center variable(vector) within groups
center=function(var,group) {
  return(var-tapply(var,group,mean,na.rm=T)[group])
}
#random simulations
randomSims=function(myDataTest,myDataControl,randomSeed,numSims=1000,noTime=TRUE,time=TRUE,timeCov=TRUE){
  set.seed(randomSeed)
  dataOfSims=data.frame(matrix(ncol=8,nrow=numSims))
  names(dataOfSims)=c("randomSample","limit","pNoTime","pTime","coefNoTime","coefTime","pTimeCov","coefTimeCov")
  #generate random simulations
  i=1
  while(i <=numSims) {
    #41 participants in total, select either 20 or 21 by random for treatment and control group
    randomSample=sample(levels(as.factor(myDataTest$ID)))
    limit=19+sample(2,1)
    randomSampleTreatment=randomSample[1:limit]
    randomSampleControl=randomSample[(limit+1):length(randomSample)]
    #get data of random participants
    treatmentGroup=myDataTest[which(myDataTest$ID %in% randomSampleTreatment),]
    treatmentGroup$group="treatment"
    controlGroup=myDataControl[which(myDataControl$ID %in% randomSampleControl),]
    controlGroup$group="control"
    randomGroupData=getReactionTimeDataset(rbind(treatmentGroup,controlGroup))
    randomGroupData$timeCovariate=center(randomGroupData$time,randomGroupData$block)
    repeat{ #this is actually only run once but simpler to cancel the rest of the code in this block
      #analysis without time
      if(noTime){
        mNoTime=lmer(reactionTime~deg*block*group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
        if(!is.null(mNoTime@optinfo$conv$lme4$messages)){
          print(paste(i,": convergence issues for notime model",sep=" "))
          break
        }
        mNoTime2=lmer(reactionTime~deg*block*group+deg*correctSide+MRexperience-block:group+(deg+block|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
        if(!is.null(mNoTime2@optinfo$conv$lme4$messages)){
          print(paste(i,": convergence issues for notime model",sep=" "))
          break
        }
        pNoTime=anova(mNoTime,mNoTime2)[[8]][2]
        coefNoTime=coef(summary(mNoTime))[9]
      } else {
        pNoTime=0
        coefNoTime=0
      }
      #analysis with time
      if(time){
        mTime=lmer(reactionTime~deg*time*block*group+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
        if(!is.null(mTime@optinfo$conv$lme4$messages)){
          print(paste(i,": convergence issues for time model",sep=" "))
          break
        }
        mTime2=lmer(reactionTime~deg*time*block*group+deg*correctSide+MRexperience-block:group+(deg+time|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
        if(!is.null(mTime2@optinfo$conv$lme4$messages)){
          print(paste(i,": convergence issues for time model2",sep=" "))
          break
        }
        pTime=anova(mTime,mTime2)[[8]][2]
        coefTime=coef(summary(mTime))[13]
      } else {
        pTime=0
        coefTime=0
      }
      #analysis with time as covariate
      if(timeCov){
        mTimeCov=lmer(reactionTime~deg*block*group+block*group*timeCovariate+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
        if(!is.null(mTimeCov@optinfo$conv$lme4$messages)){
          print(paste(i,": convergence issues for timeCov model",sep=" "))
          break
        }
        mTimeCov2=lmer(reactionTime~deg*block*group+block*group*timeCovariate+deg*correctSide+MRexperience-block:group+(deg+block|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
        if(!is.null(mTimeCov2@optinfo$conv$lme4$messages)){
          print(paste(i,": convergence issues for timeCov model",sep=" "))
          break
        }
        pTimeCov=anova(mTimeCov,mTimeCov2)[[8]][2]
        coefTimeCov=coef(summary(mTimeCov))[10]
      } else {
        pTimeCov=0
        coefTimeCov=0
      }
      #save data
      dataOfSims$randomSample[i]=list(randomSample)
      dataOfSims$limit[i]=limit
      dataOfSims$pNoTime[i]=pNoTime
      dataOfSims$pTime[i]=pTime
      dataOfSims$pTimeCov[i]=pTimeCov
      dataOfSims$coefNoTime[i]=coefNoTime
      dataOfSims$coefTime[i]=coefTime
      dataOfSims$coefTimeCov[i]=coefTimeCov
      print(paste(i,":",round(pNoTime,3),",",round(pTimeCov,3),",",round(pTime,3),",",round(coefNoTime),",",round(coefTime),",",round(coefTimeCov),sep=" "))
      #only increase loop variable if values are ok
      i=i+1
      break
    }
  }
  return(dataOfSims)
}

#plot selected indices
plotSelectedIndex=function(dataOfSims,myDataTest,myDataControl,index,name){
  #get data from index
  randomSample=unlist(dataOfSims$randomSample[index])
  limit=dataOfSims$limit[index]
  randomSampleTreatment=randomSample[1:limit]
  randomSampleControl=randomSample[(limit+1):length(randomSample)]
  #get data of random participants
  #get data of random participants
  treatmentGroup=myDataTest[which(myDataTest$ID %in% randomSampleTreatment),]
  treatmentGroup$group="treatment"
  controlGroup=myDataControl[which(myDataControl$ID %in% randomSampleControl),]
  controlGroup$group="control"
  randomGroupData=rbind(treatmentGroup,controlGroup)
  randomGroupData$cond=paste(randomGroupData$group,ifelse(randomGroupData$block=="main1","pretest","posttest"),sep="*")
  randomGroupData$condLinetype=randomGroupData$group
  randomGroupData$condColor=randomGroupData$group
  #plot graph
  generateTableAndGraphsForCondition(randomGroupData,name,FALSE,TRUE,legendProp=list(color="Group",linetypes="Group",shape="Group"))
  #print values
  print(dataOfSims[index,])
}

#generate random seeds
#sample(100000)[1]
#96235
#50793
dataOfSims=randomSims(myDataTest,myDataControl,96235)
#use control group for both
dataOfSimsControl=randomSims(myDataControl,myDataControl,50793)
sum(dataOfSimsControl$pNoTime<0.05)
sum(dataOfSimsControl$pTime<0.05)
sum(dataOfSimsControl$pTimeCov<0.05)
#difference between analyses
sum(dataOfSimsControl$pNoTime>0.05 & dataOfSimsControl$pTime>0.05)
sum(dataOfSimsControl$pNoTime>0.05 & dataOfSimsControl$pTime<=0.05)
sum(dataOfSimsControl$pNoTime<=0.05 & dataOfSimsControl$pTime>0.05)
sum(dataOfSimsControl$pNoTime<=0.05 & dataOfSimsControl$pTime<=0.05)
#save(dataOfSims,file="functions\\time as fixed Effect\\dataOfSims.RData")
#group by significance and type of analysis
sum(dataOfSims$pNoTime>0.05 & dataOfSims$pTime>0.05)
sum(dataOfSims$pNoTime>0.05 & dataOfSims$pTime<=0.05)
sum(dataOfSims$pNoTime<=0.05 & dataOfSims$pTime>0.05)
sum(dataOfSims$pNoTime<=0.05 & dataOfSims$pTime<=0.05)

sum(dataOfSims$pNoTime<=0.05 & dataOfSims$pTimeCov<=0.05)
sum(dataOfSims$pNoTime>0.05 & dataOfSims$pTimeCov>0.05)
#difference between time as covariate and not
sum(dataOfSims$pNoTime<=0.05 & dataOfSims$pTimeCov>0.05)
sum(dataOfSims$pNoTime>0.05 & dataOfSims$pTimeCov<=0.05)
sum(dataOfSimsControl$pNoTime<=0.05 & dataOfSimsControl$pTimeCov>0.05)
sum(dataOfSimsControl$pNoTime>0.05 & dataOfSimsControl$pTimeCov<=0.05)
#coefficients in wrong direction
dataOfSims[which(dataOfSims$coefNoTime>0),]
dataOfSims[which(dataOfSims$coefTime>0),]
dataOfSims[which(dataOfSims$coefTimeCov>0),]

#plot dataset for maximum p-values
#maximum p value of Time analysis while no time significant
maxPTime=which(dataOfSims$pTime==max(dataOfSims$pTime[which(dataOfSims$pNoTime<0.05)]))
plotSelectedIndex(dataOfSims,myDataTest,myDataControl,maxPTime,"randomTreatmentMaxPTime")

#maximum p value of no time analysis while time significant
maxPNoTime=which(dataOfSims$pNoTime==max(dataOfSims$pNoTime[which(dataOfSims$pTime<0.05)]))
plotSelectedIndex(dataOfSims,myDataTest,myDataControl,maxPNoTime,"randomTreatmentMaxPNoTime")

#maximum of both (sum of both)
maxPBoth=which(dataOfSims$pNoTime+dataOfSims$pTime==max(dataOfSims$pNoTime+dataOfSims$pTime))
plotSelectedIndex(dataOfSims,myDataTest,myDataControl,maxPBoth,"randomTreatmentMaxPBoth")

#plot data for no treatment condition
#minimum p value of time analysis while no time not significant
minPTimeControl=which(dataOfSimsControl$pTime==min(dataOfSimsControl$pTime[which(dataOfSimsControl$pNoTime>=0.05)]))
plotSelectedIndex(dataOfSimsControl,myDataControl,myDataControl,minPTimeControl,"randomControlMinPTime")

#minimum p value of no time analysis while time not significant
minPNoTimeControl=which(dataOfSimsControl$pNoTime==min(dataOfSimsControl$pNoTime[which(dataOfSimsControl$pTime>=0.05)]))
plotSelectedIndex(dataOfSimsControl,myDataControl,myDataControl,minPNoTimeControl,"randomControlMinPNoTime")

combineImages(c("figs/MR/Timed/randomTreatmentMaxPTimeLinePlotByCondTime.png",
                "figs/MR/Timed/randomTreatmentMaxPNoTimeLinePlotByCondTime.png",
                "figs/MR/Timed/randomTreatmentMaxPBothLinePlotByCondTime.png"),
              1,3,"figs/MR/Timed/randomTreatmentMaxP.png")

combineImages(c("figs/MR/Timed/randomControlMinPTimeLinePlotByCondTime.png",
                "figs/MR/Timed/randomControlMinPNoTimeLinePlotByCondTime.png"),
              1,2,"figs/MR/Timed/randomControlMinP.png")


#plot p-values of no treatment simulation to check type1 error rate
library(ggplot2)
dataOfSimsControl$pTimeSorted=sort(dataOfSimsControl$pTime)
dataOfSimsControl$pNoTimeSorted=sort(dataOfSimsControl$pNoTime)
dataOfSimsControl$pTimeCovSorted=sort(dataOfSimsControl$pTimeCov)
ggplot(dataOfSimsControl,aes(x=c(0:999))) + 
  geom_line(aes(y=pTimeSorted,color="1")) +
  geom_line(aes(y=pNoTimeSorted,color="2")) +
  geom_line(aes(y=pTimeCovSorted,color="3")) +
  geom_line(aes(y=c(0:999)/999,color="4")) +
  scale_color_manual(
    values=c("1"="red","2"="blue","3"="green","4"="black"), 
    labels = c("time", "no time", "time as covariate", "expected"),
    name="analysis type") +
  labs(x="number of simulations",y="magnitude of p-values") +
  theme_classic() + theme(legend.position = c(0.2,0.8))
