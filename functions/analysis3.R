### Analysis of a treatment and comparable pretest performance between groups
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
dataset.rt=datasetA3
## analysis without time
mBase=lmer(reactionTime~deg*block*group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase.summary=modelSummary(mBase,0)
#all effects significant
#main effect of block*group
mBlockXGroup=lmer(reactionTime~deg*block*group-block:group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBase,mBlockXGroup)

#main effects of two-way-interactions
mBase2=lmer(reactionTime~block*group+deg*group+deg*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase2.summary=modelSummary(mBase2,0)
#stepwise remove nonsignificant effects
#remove group*deg
mBase3=lmer(reactionTime~block*group+deg*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase3.summary=modelSummary(mBase3,0)
#all effects significant

## analysis with time
mBaseTime=lmer(reactionTime~deg*endTime*block*group+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTime.summary=modelSummary(mBaseTime,0)
#stepwise remove nonsignificant effects
mTime2=lmer(reactionTime~endTime*block*group+deg*block*group+deg*endTime*group+deg*endTime*block+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime2.summary=modelSummary(mTime2,0)
#split time*block*group
mTime3=lmer(reactionTime~deg*block*group+deg*endTime*group+deg*endTime*block+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime3.summary=modelSummary(mTime3,0)
#split deg*group*time
mTime4=lmer(reactionTime~deg*block*group+endTime*group+deg*endTime*block+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime4.summary=modelSummary(mTime4,0)
#split group*time
mTime5=lmer(reactionTime~deg*block*group+deg*endTime*block+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime5.summary=modelSummary(mTime5,0)
#split deg*block*time
mTime6=lmer(reactionTime~deg*block*group+deg*endTime+block*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime6.summary=modelSummary(mTime6,0)
#split block*time
mTime7=lmer(reactionTime~deg*block*group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime7.summary=modelSummary(mTime7,0)
#all effects significant

#nonsignificant effects
#deg*block*time
mTimeDegXBlockXTime=lmer(reactionTime~deg*block*endTime+deg*block*group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXBlockXTime.summary=modelSummary(mTimeDegXBlockXTime,0)
#block*group*time
mTimeBlockXGroupXTime=lmer(reactionTime~block*group*endTime+deg*block*group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeBlockXGroupXTime.summary=modelSummary(mTimeBlockXGroupXTime,0)
#deg*group*time
mTimeDegXGroupXTime=lmer(reactionTime~deg*group*endTime+deg*block*group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTimeDegXGroupXTime.summary=modelSummary(mTimeDegXGroupXTime,0)

#main effects of two-way-interactions
#block*group
mTime8=lmer(reactionTime~deg*block*group-block:group+deg*endTime+deg*correctSide+MRexperience+(deg+endTime|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mTime7,mTime8)
