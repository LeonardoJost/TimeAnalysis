### random slopes for full dataset
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
#models are reduced if the complex model produces a singular fit or the less complex model is not significantly worse (at alpha=0.2)
m0=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg*block|ID)+(deg*block|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#singular fit
m1=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg*block||ID)+(deg*block||modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#singular fit
m2=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg+block||ID)+(deg+block||modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#singular fit
m3=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg+block||ID)+(deg||modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#good
m4=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg+block||ID)+(block||modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#singular fit
m4=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg+block|ID)+(deg|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#singular fit
m5=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg+block|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m3,m5)
#m5 is better
#for time analysis use also time?
m6=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg+block+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
#singular fit
m7=lmer(reactionTime~deg*time*block+deg*correctSide+MRexperience+(deg+time|ID)+(1|modelNumber),data=dataset.rt,REML=FALSE,control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
anova(m5,m7)
#m7 is better