### generate graph and table output
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

source("functions/helpers.R")

#generate table and graphs for cond from MRData dataset
generateTableAndGraphsForCondition=function(MRData,conditionString,degreeGraphs=TRUE,timeGraphs=TRUE,legendProp=list(),ylab="Reaction time(ms)"){
  #calculate means by angle and interesting condition (important for plotting accuracy)
  #careful with outliers
  library(plyr)
  if(is.null(MRData$cond2))
    MRData$cond2=1
  if(is.null(MRData$condLinetype))
    MRData$condLinetype=MRData$cond
  if(is.null(MRData$condShape))
    MRData$condShape=MRData$cond
  if(is.null(MRData$condColor))
    MRData$condColor=MRData$cond
  #averages for each participant and condition
  MRDataMeansByIDDegcond=ddply(MRData,
                               .(ID,deg,cond,condLinetype,condShape,condColor),
                               summarize,
                               reactionTime1=weighted.mean(reactionTime,(typeOutlier=="hit")*cond2),
                               reactionTime2=max(c(0,weighted.mean(reactionTime,(typeOutlier=="hit")*(1-cond2))),na.rm=T),
                               reactionTime=reactionTime1-reactionTime2,
                               hits1=sum((typeOutlier=="hit")*cond2),
                               misses1=sum((typeOutlier=="incorrect")*cond2),
                               acc1=hits1/(hits1+misses1),
                               hits2=sum((typeOutlier=="hit")*(1-cond2)),
                               misses2=sum((typeOutlier=="incorrect")*(1-cond2)),
                               acc2=max(c(0,hits2/(hits2+misses2)),na.rm=T),
                               acc=acc1-acc2,
                               outliers=sum(outlier))
  #average and sd over all participants
  MRDataMeansByDegcond=ddply(MRDataMeansByIDDegcond,
                             .(deg,cond),
                             summarize,
                             reactionTimeSd=sd(reactionTime,na.rm=T),
                             reactionTime=mean(reactionTime,na.rm=T),
                             accSd=sd(acc),
                             acc=mean(acc))
  #format digits
  MRDataMeansByDegcond=format(MRDataMeansByDegcond, digits=4)
  #write table
  write.table(MRDataMeansByDegcond,file=paste("output\\MRDataMeansByDeg",conditionString,".csv",sep=""),sep=";", col.names=NA)
  
  #generate plots
  if (degreeGraphs) {
    #all data
    #generateGraphs(MRData,paste("MR/allData/",conditionString,sep=""),legendProp,ylab)
    #average by participants
    generateGraphs(MRDataMeansByIDDegcond,paste("MR/meanData/",conditionString,sep=""),legendProp,ylab)
    #accuracy is always only for averages
    generateAccGraphs(MRDataMeansByIDDegcond,paste("MR/accData/",conditionString,sep=""),legendProp)
  }
  if (timeGraphs) {
    #plot line graphs for changes over time
    generateLineGraphsByTime(MRData[which(MRData$typeOutlier=="hit"),],paste("MR/Timed/",conditionString,sep=""),legendProp,ylab)
  }
}

#generate reaction time graphs
generateGraphs=function(dataset,title,legendProp=list(),ylab="Reaction time(ms)") {
  if(is.null(legendProp$pos))
     legendProp$pos=c(0.8,0.2)
  library(ggplot2)
  #plot data as line graph (mean Data by degree and condition)
  ggplot(dataset,aes(y=reactionTime,x=deg,fill=condLinetype, shape=condShape,color=condColor)) + 
    stat_summary(na.rm=TRUE, fun=mean, geom="line",aes(linetype=condLinetype)) +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_cl_normal,geom="errorbar",position = "dodge", width=5) +
    scale_x_continuous(breaks=c(0:4)*45)+
    labs(x="degrees(°)",y=ylab,color=legendProp$color,linetype=legendProp$linetype,shape=legendProp$shape) + 
    guides(fill=FALSE) + 
    theme_classic() + theme(legend.position = legendProp$pos) + 
    scale_colour_discrete(drop=TRUE,limits = levels(dataset$condColor)) + 
    scale_linetype_discrete(drop=TRUE,limits = levels(dataset$condLinetype)) + 
    scale_shape_discrete(drop=TRUE,limits = levels(dataset$condShape))
  ggsave(paste("figs/",title,"LinePlotByCondDegree.png",sep=""))
}

#generate lins graphs by time
generateLineGraphsByTime=function(dataset,title,legendProp=list(),ylab="Reaction time(ms)") {
  if(is.null(legendProp$pos))
    legendProp$pos=c(0.9,0.9)
  library(ggplot2)
  #plot data as line graph (mean Data by degree and condition)
  ggplot(dataset,aes(y=reactionTime,x=time, color=condColor, linetype=condLinetype)) + 
    geom_smooth(aes(fill=condColor, group=cond)) +
    labs(x="time(min)",y=ylab,color=legendProp$color,fill=legendProp$color,linetype=legendProp$linetype,shape=legendProp$shape) + 
    theme_classic() + theme(legend.position = legendProp$pos) + 
    scale_colour_discrete(drop=TRUE,limits = levels(dataset$condColor)) + 
    scale_linetype_discrete(drop=TRUE,limits = levels(dataset$condLinetype)) + 
    scale_fill_discrete(drop=TRUE,limits = levels(dataset$condColor))
  ggsave(paste("figs/",title,"LinePlotByCondTime.png",sep=""))
  #plot again with linear smoothing
  ggplot(dataset,aes(y=reactionTime,x=time, color=condColor, linetype=condLinetype)) + 
    geom_smooth(aes(fill=condColor, group=cond),method="lm") +
    labs(x="time(min)",y=ylab,color=legendProp$color,fill=legendProp$color,linetype=legendProp$linetype,shape=legendProp$shape) + 
    theme_classic() + theme(legend.position = legendProp$pos) + 
    scale_colour_discrete(drop=TRUE,limits = levels(dataset$condColor)) + 
    scale_linetype_discrete(drop=TRUE,limits = levels(dataset$condLinetype)) + 
    scale_fill_discrete(drop=TRUE,limits = levels(dataset$condColor))
  ggsave(paste("figs/",title,"LinePlotByCondTimeLinear.png",sep=""))
}

#generate graphs for accuracy
generateAccGraphs=function(dataset,title,legendProp=list()) {
  if(is.null(legendProp$pos))
    legendProp$pos=c(0.8,0.8)
  library(ggplot2)
  #plot data as line graph (mean Data by degree and condition)
  ggplot(dataset,aes(y=acc,x=deg, fill=condLinetype, shape=condShape,color=condColor)) + 
    stat_summary(na.rm=TRUE, fun=mean, geom="line",aes(linetype=condLinetype)) +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_cl_normal,geom="errorbar",position = "dodge", width=5) +
    scale_x_continuous(breaks=c(0:4)*45)+
    labs(x="degrees(°)",y="Proportion of correct answers",color=legendProp$color,linetype=legendProp$linetype,shape=legendProp$shape) + 
    guides(fill=FALSE) + 
    theme_classic() + theme(legend.position = legendProp$pos) + 
    scale_colour_discrete(drop=TRUE,limits = levels(dataset$condColor)) + 
    scale_linetype_discrete(drop=TRUE,limits = levels(dataset$condLinetype)) + 
    scale_shape_discrete(drop=TRUE,limits = levels(dataset$condShape))
  ggsave(paste("figs/",title,"LinePlotByCondDegree.png",sep=""))
}

#combine multiple images into one
combineImages=function(imagesList,rows,columns,outputFile,outputWidth=1028){
  library(png) #for reading in PNGs
  # setup plot
  par(mai=rep(0,4)) # no margins
  
  # layout the plots into a matrix
  layout(matrix(1:(rows*columns), ncol=columns, byrow=TRUE))
  
  # do the plotting
  for(i in 1:(rows*columns)) {
    img =readPNG(imagesList[i])
    plot(NA,xlim=0:1,ylim=0:1,bty="n",axes=0,xaxs = 'i',yaxs='i')
    rasterImage(img,0,0,1,1)
  }
  
  # write to new image
  dev.print(png,outputFile,width=outputWidth*columns,height=outputWidth*rows)
  dev.off()
  gc()
}
