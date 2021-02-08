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
getReactionTimeDataset=function(myData, acc=FALSE){
  datasetForLMM=myData
  #scaling
  datasetForLMM$deg=datasetForLMM$deg/100
  datasetForLMM$endTime=datasetForLMM$endTime/30 #30 minutes (time is already in minutes)
  #prepare dataset
  dataset.noOutlier=datasetForLMM[which(!datasetForLMM$outlier),]
  dataset.rt=dataset.noOutlier[which(dataset.noOutlier$typeOutlier=="hit"),]
  dataset.rt$deg=dataset.rt$deg-mean(dataset.rt$deg) #center degree
  dataset.noOutlier$deg=dataset.noOutlier$deg-mean(dataset.noOutlier$deg) #center degree
  #normalizing time and centering degree are necessary to analyze main effects of partial interaction (block*group) when higher-
  #order interactions are present (deg*block*group+time*block*group). Main effects are calculated for value 0
  #0 of degree: average effect due to centering (this is "standard" main effect of removing higher order interaction)
  #0 of time: difference between blocks
  if(acc)
    return(dataset.noOutlier)
  else
    return(dataset.rt)
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
#split block by time of 10 minutes
myData$block=toChar(myData$block)
myData$block=ifelse(myData$endTime>10*60*1000,ifelse(myData$endTime>20*60*1000,"main3","main2"),"main1")
#normalize time to minutes
myData$endTime=myData$endTime/60000 

### analysis 1
#all blocks
datasetA1=getReactionTimeDataset(myData)

#generate some plots of the separation by blocks
myData$cond=ifelse(myData$block=="main3","20-30min",ifelse(myData$block=="main2","10-20min","0-10min"))
generateTableAndGraphsForCondition(myData,"block",TRUE,TRUE,legendProp=list(color="Block",linetypes="Block",shape="Block"))

#normalize time 0 to end of first block
myData$endTime=myData$endTime-10

### analysis 2
#blocks 1 and 3 but move block 3 in time to block 2
myData13=myData[which(myData$block!="main2"),]
myData13$endTime[which(myData13$block=="main3")]=myData13$endTime[which(myData13$block=="main3")]-10
datasetA2=getReactionTimeDataset(myData13)

### analysis 3-5
#add duplicated control group
myDataTest=myData13
myDataTest$group="treatment"
#rename main3 to main2 to compare with control
myDataTest$block[which(myDataTest$block=="main3")]="main2"
#control group with only block 1 and 2
myDataControl=myData[which(myData$block!="main3"),]
myDataControl$group="control"
#create new ids for control (simulate different subjects in both groups)
levels(myDataControl$ID)=paste("idControl",sample.int(length(levels(myDataControl$ID))),sep="")
myDataTestControl=rbind(myDataTest,myDataControl)

### analysis 3
datasetA3=getReactionTimeDataset(myDataTestControl)

### analysis 4a
#calculate average performance in first block
myDataTestControl$meanPerformance=0
for(id in levels(as.factor(myDataTestControl$ID))){
  myDataTestControl$meanPerformance[which(myDataTestControl$ID==id)]=mean(myDataTestControl$reactionTime[which(myDataTestControl$ID==id & myDataTestControl$block=="main1")])
}
#select by performance and group (better group as treatment)
myDataTestControl2=myDataTestControl[which((myDataTestControl$group=="treatment" & myDataTestControl$meanPerformance<median(myDataTestControl$meanPerformance)) |
                                             (myDataTestControl$group=="control" & myDataTestControl$meanPerformance>=median(myDataTestControl$meanPerformance))
),]

datasetA4a=getReactionTimeDataset(myDataTestControl2)
#plot block*group interaction over time
myDataTestControl2$cond=paste(myDataTestControl2$group,ifelse(myDataTestControl2$block=="main1","pretest","posttest"),sep="*")
myDataTestControl2$condLinetype=myDataTestControl2$group
myDataTestControl2$condColor=myDataTestControl2$group
generateTableAndGraphsForCondition(myDataTestControl2,"blockXgroup4a",FALSE,TRUE,legendProp=list(color="Group",linetypes="Group",shape="Group"))

### analysis 4b
#select by performance and group (better group as control)
myDataTestControl3=myDataTestControl[which((myDataTestControl$group=="treatment" & myDataTestControl$meanPerformance>=median(myDataTestControl$meanPerformance)) |
                                             (myDataTestControl$group=="control" & myDataTestControl$meanPerformance<median(myDataTestControl$meanPerformance))
),]

datasetA4b=getReactionTimeDataset(myDataTestControl3)
#plot block*group interaction over time
myDataTestControl3$cond=paste(myDataTestControl3$group,ifelse(myDataTestControl3$block=="main1","pretest","posttest"),sep="*")
myDataTestControl3$condLinetype=myDataTestControl3$group
myDataTestControl3$condColor=myDataTestControl3$group
generateTableAndGraphsForCondition(myDataTestControl3,"blockXgroup4b",FALSE,TRUE,legendProp=list(color="Group",linetypes="Group",shape="Group"))

### analysis 5a
#compare block 1 and 2/3 without actual treatment but separated by performance
myDataTestControl4=myData
#calculate average performance in first block
myDataTestControl4$meanPerformance=0
for(id in levels(as.factor(myDataTestControl4$ID))){
  myDataTestControl4$meanPerformance[which(myDataTestControl4$ID==id)]=mean(myDataTestControl4$reactionTime[which(myDataTestControl4$ID==id & myDataTestControl4$block=="main1")])
}
myDataTestControl4$group=ifelse(myDataTestControl4$meanPerformance<median(myDataTestControl4$meanPerformance),"treatment","control")
#merge block 2 and 3
myDataTestControl4$block[which(myDataTestControl4$block=="main3")]="main2"

datasetA5a=getReactionTimeDataset(myDataTestControl4)
#plot block*group interaction over time
myDataTestControl4$cond=paste(myDataTestControl4$group,ifelse(myDataTestControl4$block=="main1","pretest","posttest"),sep="*")
myDataTestControl4$condLinetype=myDataTestControl4$group
myDataTestControl4$condColor=myDataTestControl4$group
generateTableAndGraphsForCondition(myDataTestControl4,"blockXgroup5a",FALSE,TRUE,legendProp=list(color="Group",linetypes="Group",shape="Group"))

### analysis 5b
#separate groups not by performance but by improvement
myDataTestControl5=myDataTestControl4
for(id in levels(as.factor(myDataTestControl4$ID))){
  myDataTestControl5$meanPerformance[which(myDataTestControl5$ID==id)]=
    mean(myDataTestControl5$reactionTime[which(myDataTestControl5$ID==id & myDataTestControl5$block=="main1")])-
    mean(myDataTestControl5$reactionTime[which(myDataTestControl5$ID==id & myDataTestControl5$block=="main2")])
}
myDataTestControl5$group=ifelse(myDataTestControl5$meanPerformance<median(myDataTestControl5$meanPerformance),"control","treatment")

datasetA5b=getReactionTimeDataset(myDataTestControl5)
#plot block*group interaction over time
myDataTestControl5$cond=paste(myDataTestControl5$group,ifelse(myDataTestControl5$block=="main1","pretest","posttest"),sep="*")
myDataTestControl5$condLinetype=myDataTestControl5$group
myDataTestControl5$condColor=myDataTestControl5$group
generateTableAndGraphsForCondition(myDataTestControl5,"blockXgroup5b",FALSE,TRUE,legendProp=list(color="Group",linetypes="Group",shape="Group"))

### random simulations
#get control and treatment data from all participants
myDataTest=myData13
myDataTest$group="treatment"
#rename main3 to main2 to compare with control
myDataTest$block[which(myDataTest$block=="main3")]="main2"
#control group with only block 1 and 2
myDataControl=myData[which(myData$block!="main3"),]
myDataControl$group="control"