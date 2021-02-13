### Analysis of a treatment and comparable pretest performance between groups
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
dataset.rt=datasetA3
## analysis without time
mBase=lmer(reactionTime~deg*block*group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBase.summary=modelSummary(mBase,0)
#main effect of block*group
mBlockXGroup=lmer(reactionTime~deg*block*group-block:group+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mBase,mBlockXGroup)

## analysis with time
mBaseTime=lmer(reactionTime~deg*time*block*group+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mBaseTime.summary=modelSummary(mBaseTime,0)
#stepwise remove nonsignificant effects
mTime2=lmer(reactionTime~time*block*group+deg*block*group+deg*time*group+deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime2.summary=modelSummary(mTime2,0)
#split time*block*group
mTime3=lmer(reactionTime~deg*block*group+deg*time*group+deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime3.summary=modelSummary(mTime3,0)
#split deg*group*time
mTime4=lmer(reactionTime~deg*block*group+time*group+deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime4.summary=modelSummary(mTime4,0)
#split group*time
mTime5=lmer(reactionTime~deg*block*group+deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime5.summary=modelSummary(mTime5,0)
#split deg*block*time
mTime6=lmer(reactionTime~deg*block*group+deg*time+block*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime6.summary=modelSummary(mTime6,0)
#split block*time
mTime7=lmer(reactionTime~deg*block*group+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
mTime7.summary=modelSummary(mTime7,0)
#all effects significant

#main effects of two-way-interactions
#block*group
mTime8=lmer(reactionTime~deg*block*group-block:group+deg*time+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(mTime7,mTime8)
