### Analysis of no treatment and between-groups differences in pretest performance
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

#load dataset
dataset.rt=datasetA5a
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
#remove block*group
mBase3=lmer(reactionTime~deg*group+deg*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase3.summary=modelSummary(mBase3,0)


## analysis with time
mBaseTime=lmer(reactionTime~deg*endTime*block*group+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTime.summary=modelSummary(mBaseTime,0)
#stepwise remove nonsignificant effects
mTime2=lmer(reactionTime~endTime*block*group+deg*block*group+deg*endTime*group+deg*endTime*block+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime2.summary=modelSummary(mTime2,0)
#split time*group*deg
mTime3=lmer(reactionTime~endTime*block*group+deg*block*group+deg*endTime*block+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime3.summary=modelSummary(mTime3,0)
#split time*block*deg
mTime4=lmer(reactionTime~endTime*block*group+deg*block*group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime4.summary=modelSummary(mTime4,0)

#nonsignificant effects
#deg*group*time
mTimeDegXGroupXTime=lmer(reactionTime~deg*endTime*group+
                           endTime*block*group+deg*block*group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXGroupXTime.summary=modelSummary(mTimeDegXGroupXTime,0)
#deg*block*time
mTimeDegXBlockXTime=lmer(reactionTime~deg*endTime*block+
                           endTime*block*group+deg*block*group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXBlockXTime.summary=modelSummary(mTimeDegXBlockXTime,0)

#main effects of two-way-interactions
#block*group
mTime5=lmer(reactionTime~endTime*block*group+deg*block*group-block:group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mTime4,mTime5)


## separate groups
#control group
dataset.rt=datasetA5a[which(datasetA5a$group=="control"),]

mBaseTimeControl=lmer(reactionTime~endTime*block+deg*block+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeControl.summary=modelSummary(mBaseTimeControl,0)
#split block*time
mBaseTimeControl2=lmer(reactionTime~deg*block+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeControl2.summary=modelSummary(mBaseTimeControl2,0)
#split deg*block
mBaseTimeControl3=lmer(reactionTime~block+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeControl3.summary=modelSummary(mBaseTimeControl3,0)
#block non significant
#main effect of block*time
mBaseTimeControl4=lmer(reactionTime~endTime*block+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeControl4.summary=modelSummary(mBaseTimeControl4,0)

#treatment group
dataset.rt=datasetA5a[which(datasetA5a$group=="treatment"),]

mBaseTimeTreatment=lmer(reactionTime~endTime*block+deg*block+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTimeTreatment.summary=modelSummary(mBaseTimeTreatment,0)
#main effect of block
mBaseTimeTreatment2=lmer(reactionTime~endTime*block+deg*block-block+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBaseTimeTreatment,mBaseTimeTreatment2)
