### Analysis of treatment using random sampling of participants
#     Copyright (C) 2019  Leonardo Jost
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

getReactionTimeDataset=function(myData){
  datasetForLMM=myData
  #scaling
  datasetForLMM$deg=datasetForLMM$deg/100
  datasetForLMM$endTime=datasetForLMM$endTime/30 #30 minutes (time is already in minutes)
  #prepare dataset
  dataset.noOutlier=datasetForLMM[which(!datasetForLMM$outlier),]
  dataset.rt=dataset.noOutlier[which(dataset.noOutlier$typeOutlier=="hit"),]
  dataset.rt$deg=dataset.rt$deg-mean(dataset.rt$deg) #center degree
  #normalizing time and centering degree are necessary to analyze main effects of partial interaction (block*group) when higher-
  #order interactions are present (deg*block*group+time*block*group). Main effects are calculated for value 0
  #0 of degree: average effect due to centering (this is "standard" main effect of removing higher order interaction)
  #0 of time: difference between blocks
  return(dataset.rt)
}

center=function(var,group) {
  return(var-tapply(var,group,mean,na.rm=T)[group])
}

set.seed(783881)
numSims=100
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
  controlGroup=myDataControl[which(myDataControl$ID %in% randomSampleControl),]
  randomGroupData=getReactionTimeDataset(rbind(treatmentGroup,controlGroup))
  randomGroupData$time=center(randomGroupData$endTime,randomGroupData$block)
  #analysis without time
  mNoTime=lmer(reactionTime~deg*block*group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
  mNoTime2=lmer(reactionTime~deg*block*group+deg*correctSide+MRexperience-block:group+(deg+block|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
  pNoTime=anova(mNoTime,mNoTime2)[[8]][2]
  coefNoTime=coef(summary(mNoTime))[9]
  if(!isSingular(mNoTime) & !isSingular(mNoTime2)){ #retry in case of singular fit
    #analysis with time
    mTime=lmer(reactionTime~deg*endTime*block*group+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
    mTime2=lmer(reactionTime~deg*endTime*block*group+deg*correctSide+MRexperience-block:group+(deg+endTime|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
    pTime=anova(mTime,mTime2)[[8]][2]
    coefTime=coef(summary(mTime))[13]
    #analysis with time as covariate
    mTimeCov=lmer(reactionTime~deg*block*group+block*group*time+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
    mTimeCov2=lmer(reactionTime~deg*block*group+block*group*time+deg*correctSide+MRexperience-block:group+(deg+block|ID)+(1|modelNumber),data=randomGroupData,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
    pTimeCov=anova(mTimeCov,mTimeCov2)[[8]][2]
    coefTimeCov=coef(summary(mTimeCov))[10]
    if(!isSingular(mTime) & !isSingular(mTime2) & !isSingular(mTimeCov) & !isSingular(mTimeCov2)){ #retry in case of singular fit      
      #save data
      dataOfSims$randomSample[i]=list(randomSample)
      dataOfSims$limit[i]=limit
      dataOfSims$pNoTime[i]=pNoTime
      dataOfSims$pTime[i]=pTime
      dataOfSims$pTimeCov[i]=pTimeCov
      dataOfSims$coefNoTime[i]=coefNoTime
      dataOfSims$coefTime[i]=coefTime
      dataOfSims$coefTimeCov[i]=coefTimeCov
      print(paste(i,":",round(pNoTime,3),",",round(pTimeCov,3),",",round(pTime,3),",",round(coefNoTime),",",round(coefTime),sep=" "))
      #only increase loop variable if values are ok
      i=i+1
    } else {
      print(paste(i,": singular fit for time model",sep=" "))
    }
  } else {
    print(paste(i,": singular fit for noTime model",sep=" "))
  }
}
save(dataOfSims,file="functions\\time as fixed Effect\\dataOfSims.RData")
#group by significance and type of analysis
sum(dataOfSims$pNoTime>0.05 & dataOfSims$pTime>0.05)
sum(dataOfSims$pNoTime>0.05 & dataOfSims$pTime<=0.05)
sum(dataOfSims$pNoTime<=0.05 & dataOfSims$pTime>0.05)
sum(dataOfSims$pNoTime<=0.05 & dataOfSims$pTime<=0.05)

sum(dataOfSims$pNoTime<=0.05 & dataOfSims$pTimeCov<=0.05)
sum(dataOfSims$pNoTime>0.05 & dataOfSims$pTimeCov>0.05)
#coefficients in wrong direction?
dataOfSims[which(dataOfSims$coefNoTime>0),]
dataOfSims[which(dataOfSims$coefTime>0),]
#plot selected indices
plotSelectedIndex=function(index,name){
  #get data from index
  randomSample=unlist(dataOfSims$randomSample[index])
  limit=dataOfSims$limit[index]
  randomSampleTreatment=randomSample[1:limit]
  randomSampleControl=randomSample[(limit+1):length(randomSample)]
  #get data of random participants
  treatmentGroup=myDataTest[which(myDataTest$ID %in% randomSampleTreatment),]
  controlGroup=myDataControl[which(myDataControl$ID %in% randomSampleControl),]
  randomGroupData=getReactionTimeDataset(rbind(treatmentGroup,controlGroup))
  randomGroupData$cond=paste(randomGroupData$group,ifelse(randomGroupData$block=="main1","pretest","posttest"),sep="*")
  randomGroupData$condForLineTypes=randomGroupData$group
  #plot graph
  generateTableAndGraphsForCondition(randomGroupData,name,FALSE,TRUE,"Group",TRUE)
  #print values
  print(dataOfSims[index,])
}
#maximum p value of Time analysis while no time significant
maxPTime=which(dataOfSims$pTime==max(dataOfSims$pTime[which(dataOfSims$pNoTime<0.05)]))
plotSelectedIndex(maxPTime,"randommaxPTime")

#maximum p value of no time analysis while time significant
maxPNoTime=which(dataOfSims$pNoTime==max(dataOfSims$pNoTime[which(dataOfSims$pTime<0.05)]))
plotSelectedIndex(maxPNoTime,"randommaxPNoTime")

#maximum of both (sum of both)
maxPBoth=which(dataOfSims$pNoTime+dataOfSims$pTime==max(dataOfSims$pNoTime+dataOfSims$pTime))
plotSelectedIndex(maxPBoth,"randommaxPBoth")
