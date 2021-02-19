### create all datasets for analysis and plot some graphs for these
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

### functions
source("functions/helpers.R")
source("functions/generateGraphsAndTables.R", encoding="utf-8")
#get the dataset for analysis with rescaled variables
getReactionTimeDataset=function(myData){
  #scaling
  myData$time=myData$time/30 #30 minutes (time is already in minutes)
  #center degree
  myData$deg=scale(myData$deg)
  #myData$block=sapply(myData$block,function(i) contr.sum(3)[i,])
  myData$correctSide=sapply(as.factor(myData$correctSide),function(i) contr.sum(2)[i,])
  myData$MRexperience=sapply(as.factor(myData$MRexperience),function(i) contr.sum(2)[i,])
  #normalizing time and centering degree are necessary to analyze main effects of partial interaction (block*group) when higher-
  #order interactions are present (deg*block*group+time*block*group). Main effects are calculated for value 0
  #0 of degree: average effect due to centering (this is "standard" main effect of removing higher order interaction)
  #0 of time: difference between blocks
  return(myData)
}
##create output directories, if they don't exist (outputs warnings otherwise)
dir.create("figs")
dir.create("output")
dir.create("figs/MR")
dir.create("figs/MR/allData")
dir.create("figs/MR/meanData")
dir.create("figs/MR/Timed/")
dir.create("figs/MR/accData/")


### main script
#load full dataset
myData=read.csv(file="data\\dataset.csv",sep=";")
#rename variables, convert ID to factor
myData$time=myData$endTime
myData$ID=as.factor(myData$ID)
#reaction time dataset
#use only correct answers and remove outliers
myData=myData[which(!myData$outlier),]
myData=myData[which(myData$typeOutlier=="hit"),]
#split dataset into three blocks by time of 10 minutes
myData$block=toChar(myData$block)
myData$block=ifelse(myData$time>10*60*1000,ifelse(myData$time>20*60*1000,"main3","main2"),"main1")
#normalize time to minutes
myData$time=myData$time/60000
#calculate mean improvement from block 1 to 2 for each participant
myData$meanImprovement=0
for(id in levels(as.factor(myData$ID))){
  myData$meanImprovement[which(myData$ID==id)]=
    mean(myData$reactionTime[which(myData$ID==id & myData$block=="main1")])-
    mean(myData$reactionTime[which(myData$ID==id & myData$block=="main2")])
}
#calculate average performance in first block for each participant
myData$meanPerformance=0
for(id in levels(as.factor(myData$ID))){
  myData$meanPerformance[which(myData$ID==id)]=mean(myData$reactionTime[which(myData$ID==id & myData$block=="main1")])
}

### create datasets for each analysis
# analysis 1: all blocks
datasetA1=getReactionTimeDataset(myData)

#generate some plots for this dataset
myData$cond=ifelse(myData$block=="main3","20-30min",ifelse(myData$block=="main2","10-20min","0-10min"))
generateTableAndGraphsForCondition(myData,"block",legendProp=list(color="Block",linetypes="Block",shape="Block"))

#normalize time 0 to end of first block
myData$time=myData$time-10

#analysis 2: blocks 1 and 3 but move block 3 forward 10 minutes
myData13=myData[which(myData$block!="main2"),]
myData13$time[which(myData13$block=="main3")]=myData13$time[which(myData13$block=="main3")]-10
datasetA2=getReactionTimeDataset(myData13)

### analysis 3-5
#preparation: create datasets for treatment and control group
#treatment dataset: blocks 1 and 3
myDataTest=myData13
myDataTest$group="treatment"
#rename main3 to main2 to compare with control
myDataTest$block[which(myDataTest$block=="main3")]="main2"
#control group: blocks 1 and 2
myDataControl=myData[which(myData$block!="main3"),]
myDataControl$group="control"
#create new ids for control (simulate different subjects in both groups)
levels(myDataControl$ID)=paste(levels(myDataControl$ID),"control",sep="")
#combina datasets into one
myDataTestControl=rbind(myDataTest,myDataControl)

#analysis 3: use values from all participants
datasetA3=getReactionTimeDataset(myDataTestControl)

#analysis 4a: select by performance and group (better group as treatment)
myDataTestControl2=myDataTestControl[which((myDataTestControl$group=="treatment" & myDataTestControl$meanPerformance<median(unique(myDataTestControl$meanPerformance))) |
                                             (myDataTestControl$group=="control" & myDataTestControl$meanPerformance>=median(unique(myDataTestControl$meanPerformance)))
),]

datasetA4a=getReactionTimeDataset(myDataTestControl2)
#plot block*group interaction over time
myDataTestControl2$cond=paste(myDataTestControl2$group,ifelse(myDataTestControl2$block=="main1","pretest","posttest"),sep="*")
myDataTestControl2$condLinetype=myDataTestControl2$group
myDataTestControl2$condColor=myDataTestControl2$group
generateTableAndGraphsForCondition(myDataTestControl2,"blockXgroup4a",legendProp=list(color="Group",linetypes="Group",shape="Group"))

#analysis 4b: select by performance and group (better group as control)
myDataTestControl3=myDataTestControl[which((myDataTestControl$group=="treatment" & myDataTestControl$meanPerformance>=median(unique(myDataTestControl$meanPerformance))) |
                                             (myDataTestControl$group=="control" & myDataTestControl$meanPerformance<median(unique(myDataTestControl$meanPerformance)))
),]

datasetA4b=getReactionTimeDataset(myDataTestControl3)
#plot block*group interaction over time
myDataTestControl3$cond=paste(myDataTestControl3$group,ifelse(myDataTestControl3$block=="main1","pretest","posttest"),sep="*")
myDataTestControl3$condLinetype=myDataTestControl3$group
myDataTestControl3$condColor=myDataTestControl3$group
generateTableAndGraphsForCondition(myDataTestControl3,"blockXgroup4b",legendProp=list(color="Group",linetypes="Group",shape="Group"))

#analysis 5a: select by improvement (better learners in control group)
myDataTestControl4=myDataTestControl[which((myDataTestControl$group=="treatment" & myDataTestControl$meanImprovement<median(unique(myDataTestControl$meanImprovement))) |
                                             (myDataTestControl$group=="control" & myDataTestControl$meanImprovement>=median(unique(myDataTestControl$meanImprovement)))
),]

datasetA5a=getReactionTimeDataset(myDataTestControl4)
#plot block*group interaction over time
myDataTestControl4$cond=paste(myDataTestControl4$group,ifelse(myDataTestControl4$block=="main1","pretest","posttest"),sep="*")
myDataTestControl4$condLinetype=myDataTestControl4$group
myDataTestControl4$condColor=myDataTestControl4$group
generateTableAndGraphsForCondition(myDataTestControl4,"blockXgroup5a",legendProp=list(color="Group",linetypes="Group",shape="Group"))

#analysis 5a: select by improvement, but without treatment
myDataTestControl5=myDataTestControl[which(myDataTestControl$group=="control"),]
myDataTestControl5$group=ifelse(myDataTestControl5$meanImprovement>=median(myDataTestControl5$meanImprovement),"treatment","control")

datasetA5b=getReactionTimeDataset(myDataTestControl5)
#plot block*group interaction over time
myDataTestControl5$cond=paste(myDataTestControl5$group,ifelse(myDataTestControl5$block=="main1","pretest","posttest"),sep="*")
myDataTestControl5$condLinetype=myDataTestControl5$group
myDataTestControl5$condColor=myDataTestControl5$group
generateTableAndGraphsForCondition(myDataTestControl5,"blockXgroup5b",legendProp=list(color="Group",linetypes="Group",shape="Group"))


### random simulations
#get one set for control and treatment data each from all participants
#participants are later randomly selected from groups
#note that ids are same in this case, otherwise datasets are as before
myDataTest=myData13
myDataTest$group="treatment"
#rename main3 to main2 to compare with control
myDataTest$block[which(myDataTest$block=="main3")]="main2"
#control group with only block 1 and 2
myDataControl=myData[which(myData$block!="main3"),]
myDataControl$group="control"