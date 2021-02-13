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

## analysis with time
mBaseTime=lmer(reactionTime~deg*time*block*group+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTime.summary=modelSummary(mBaseTime,0)
#stepwise remove nonsignificant effects
mTime2=lmer(reactionTime~time*block*group+deg*block*group+deg*time*group+deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime2.summary=modelSummary(mTime2,0)
#split deg*block*time
mTime3=lmer(reactionTime~time*block*group+deg*block*group+deg*time*group+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime3.summary=modelSummary(mTime3,0)
#split deg*time*group
mTime4=lmer(reactionTime~time*block*group+deg*block*group+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime4.summary=modelSummary(mTime4,0)
#split time*block*group
mTime5=lmer(reactionTime~time*group+time*block+deg*block*group+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime5.summary=modelSummary(mTime5,0)
#split time*block
mTime6=lmer(reactionTime~time*group+deg*block*group+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime6.summary=modelSummary(mTime6,0)
#split time*group
mTime7=lmer(reactionTime~deg*block*group+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime7.summary=modelSummary(mTime7,0)
#all effects significant

#main effects of block*group
mTime8=lmer(reactionTime~deg*block*group+deg*time+deg*correctSide+MRexperience-block:group+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mTime7,mTime8)
#significant


## separate groups
#control group
dataset.rt=datasetA4b[which(datasetA4b$group=="control"),]

#without time
mBaseControl=lmer(reactionTime~deg*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseControl.summary=modelSummary(mBaseControl,0)
#split deg*block
mBaseControl2=lmer(reactionTime~block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseControl2.summary=modelSummary(mBaseControl2,0)
#block significant

mBaseTimeControl=lmer(reactionTime~deg*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeControl.summary=modelSummary(mBaseTimeControl,0)
#main effect of block
mBaseTimeControl2=lmer(reactionTime~deg*block+deg*time+deg*correctSide+MRexperience-block+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBaseTimeControl,mBaseTimeControl2)
#block non significant

#treatment group
dataset.rt=datasetA4b[which(datasetA4b$group=="treatment"),]

mBaseTimeTreatment=lmer(reactionTime~deg*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeTreatment.summary=modelSummary(mBaseTimeTreatment,0)
#split deg*time
mBaseTimeTreatment2=lmer(reactionTime~deg*block+time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeTreatment2.summary=modelSummary(mBaseTimeTreatment2,0)
#main effect of block
mBaseTimeTreatment3=lmer(reactionTime~deg*block+time+deg*correctSide+MRexperience-block+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBaseTimeTreatment2,mBaseTimeTreatment3)
#block significant