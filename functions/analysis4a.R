### Analysis of a treatment and between-groups differences in pretest performance (better group as treatment)
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
dataset.rt=datasetA4a
## analysis without time
mBase=lmer(reactionTime~deg*block*group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase.summary=modelSummary(mBase,0)
#split deg*block*group
mBase2=lmer(reactionTime~block*group+deg*group+deg*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase2.summary=modelSummary(mBase2,0)
#split block*group
mBase3=lmer(reactionTime~deg*group+deg*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase3.summary=modelSummary(mBase3,0)
#all effects significant

## analysis with time
mBaseTime=lmer(reactionTime~deg*time*block*group+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTime.summary=modelSummary(mBaseTime,0)
#stepwise remove nonsignificant effects
mTime2=lmer(reactionTime~time*block*group+deg*block*group+deg*time*group+deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime2.summary=modelSummary(mTime2,0)
#split deg*time*block
mTime3=lmer(reactionTime~time*block*group+deg*block*group+deg*time*group+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime3.summary=modelSummary(mTime3,0)
#split deg*group*block
mTime4=lmer(reactionTime~time*block*group+deg*block+deg*time*group+
              deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime4.summary=modelSummary(mTime4,0)
#split deg*block
mTime5=lmer(reactionTime~time*block*group+deg*time*group+
              deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime5.summary=modelSummary(mTime5,0)
#split deg*time*group
mTime6=lmer(reactionTime~time*block*group+deg*group+deg*time+
              deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime6.summary=modelSummary(mTime6,0)
#all effects significant

#nonsignificant effects
#deg*group*time
mTimeDegXGroupXTime=lmer(reactionTime~time*block*group+deg*time*group+
              deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXGroupXTime.summary=modelSummary(mTimeDegXGroupXTime,0)
#deg*block*time
mTimeDegXBlockXTime=lmer(reactionTime~time*block*group+deg*time*block+
                           deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXBlockXTime.summary=modelSummary(mTimeDegXBlockXTime,0)
#deg*block*group
mTimeDegXBlockXGroup=lmer(reactionTime~time*block*group+deg*group*block+
                           deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXBlockXGroup.summary=modelSummary(mTimeDegXBlockXGroup,0)

#main effects of two-way-interactions
#block*group
mTime7=lmer(reactionTime~time*block*group+deg*group+deg*time-block:group+
              deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mTime7,mTime6)

## separate groups
#control group
dataset.rt=datasetA4a[which(datasetA4a$group=="control"),]

mBaseTimeControl=lmer(reactionTime~time*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeControl.summary=modelSummary(mBaseTimeControl,0)
#main effect of block
mBaseTimeControl2=lmer(reactionTime~time*block-block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBaseTimeControl2,mBaseTimeControl)
#block non significant

#treatment group
dataset.rt=datasetA4a[which(datasetA4a$group=="treatment"),]

mBaseTimeTreatment=lmer(reactionTime~time*block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeTreatment.summary=modelSummary(mBaseTimeTreatment,0)
#effect of deg*block
mBaseTimeTreatment2=lmer(reactionTime~time*block-block+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBaseTimeTreatment2,mBaseTimeTreatment)
#block significant