### Analysis of all blocks without treatment
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

#load dataset
dataset.rt=datasetA1
#get random slopes
mBaseTime=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg*time*block|ID)+(deg|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))


#load dataset
dataset.rt=datasetA1
## analysis without time
mBase=lmer(reactionTime~deg*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase.summary=modelSummary(mBase,0)
mBlock=lmer(reactionTime~deg*block-block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBase,mBlock)

## analysis with time
mBaseTime=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTime.summary=modelSummary(mBaseTime,0)
#stepwise remove nonsignificant effects
mTime2=lmer(reactionTime~time*block+deg*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime2.summary=modelSummary(mTime2,0)
#split time*block
mTime3=lmer(reactionTime~deg*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime3.summary=modelSummary(mTime3,0)
#split deg*block
mTime4=lmer(reactionTime~block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime4.summary=modelSummary(mTime4,0)
#remove block
mTime5=lmer(reactionTime~deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime5.summary=modelSummary(mTime5,0)
#all effects significant

#nonsignificant effects
#block
mTimeBlock=lmer(reactionTime~block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeBlock.summary=modelSummary(mTimeBlock,0)
#deg*block
mTimeDegXBlock=lmer(reactionTime~deg*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXBlock.summary=modelSummary(mTimeDegXBlock,0)
#time*block
mTimeTimeXBlock=lmer(reactionTime~time*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeTimeXBlock.summary=modelSummary(mTimeTimeXBlock,0)


###center time
#move all to 0
dataset.rt$time2=dataset.rt$time
dataset.rt$time2[which(dataset.rt$block=="main3")]=dataset.rt$time[which(dataset.rt$block=="main3")]-20/30
dataset.rt$time2[which(dataset.rt$block=="main2")]=dataset.rt$time[which(dataset.rt$block=="main2")]-10/30
center=function(var,group) {
  return(var-tapply(var,group,mean,na.rm=T)[group])
}
dataset.rt$time=center(dataset.rt$time,dataset.rt$block)

mBaseTimeCovariate=lmer(reactionTime~deg*block+time+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
summary(mBaseTimeCovariate)
mBaseTimeCovariate2=lmer(reactionTime~deg*block+time2+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
summary(mBaseTimeCovariate2)

m1=lmer(reactionTime~block+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
m2=lmer(reactionTime~block+time+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
summary(m1)
summary(m2)
dataset.rt$cond=ifelse(dataset.rt$block=="main3","20-30min",ifelse(dataset.rt$block=="main2","10-20min","0-10min"))
generateTableAndGraphsForCondition(dataset.rt,"block2",TRUE,TRUE,legendProp=list(color="Block",linetypes="Block",shape="Block"))

