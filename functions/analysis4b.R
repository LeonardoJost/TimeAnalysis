### Analysis of a treatment and between-groups differences in pretest performance (better group as control)
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
dataset.rt=datasetA4b
## analysis without time
mBase=lmer(reactionTime~deg*block*group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase.summary=modelSummary(mBase,0)
#significant
#block*group
mBlockXGroup=lmer(reactionTime~deg*block*group-block:group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBase,mBlockXGroup)
#split deg*block*group
mBase2=lmer(reactionTime~block*group+deg*group+deg*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase2.summary=modelSummary(mBase2,0)
#all effects significant

## analysis with time
mBaseTime=lmer(reactionTime~deg*time*block*group+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTime.summary=modelSummary(mBaseTime,0)
#stepwise remove nonsignificant effects
mTime2=lmer(reactionTime~time*block*group+deg*block*group+deg*time*group+deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime2.summary=modelSummary(mTime2,0)
#split deg*group*time
mTime3=lmer(reactionTime~time*block*group+deg*block*group+deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime3.summary=modelSummary(mTime3,0)
#split deg*time*block
mTime4=lmer(reactionTime~time*block*group+deg*block*group+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime4.summary=modelSummary(mTime4,0)
#all effects significant

#nonsignificant effects
#deg*block*time
mTimeDegXBlockXTime=lmer(reactionTime~deg*block*time+time*block*group+deg*block*group+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXBlockXTime.summary=modelSummary(mTimeDegXBlockXTime,0)

#main effects of block*group
mTime5=lmer(reactionTime~time*block*group+deg*block*group+deg*time+deg*correctSide+MRexperience-block:group+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mTime4,mTime5)

#significant


## separate groups
#control group
dataset.rt=datasetA4b[which(datasetA4b$group=="control"),]

mBaseTimeControl=lmer(reactionTime~time*block+deg*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeControl.summary=modelSummary(mBaseTimeControl,0)
#split deg*block
mBaseTimeControl2=lmer(reactionTime~time*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeControl2.summary=modelSummary(mBaseTimeControl2,0)
#main effect of block
mBaseTimeControl3=lmer(reactionTime~time*block-block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBaseTimeControl3,mBaseTimeControl2)
#block non significant

#treatment group
dataset.rt=datasetA4b[which(datasetA4b$group=="treatment"),]

mBaseTimeTreatment=lmer(reactionTime~time*block+deg*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeTreatment.summary=modelSummary(mBaseTimeTreatment,0)
#split deg*block
mBaseTimeTreatment2=lmer(reactionTime~time*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeTreatment2.summary=modelSummary(mBaseTimeTreatment2,0)
#time*block non significant
#effect of deg*block
mBaseTimeTreatment3=lmer(reactionTime~deg*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeTreatment3.summary=modelSummary(mBaseTimeTreatment3,0)
#main effect of block
mBaseTimeTreatment4=lmer(reactionTime~block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeTreatment4.summary=modelSummary(mBaseTimeTreatment4,0)
#block significant