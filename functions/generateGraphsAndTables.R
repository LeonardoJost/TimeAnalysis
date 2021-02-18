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
generateTableAndGraphsForCondition=function(MRData,conditionString,legendProp=list(),ylab="Reaction time(ms)"){
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
  generateLineGraphsByTime(MRData[which(MRData$typeOutlier=="hit"),],paste("MR/Timed/",conditionString,sep=""),legendProp,ylab)
}

#generate lins graphs by time
generateLineGraphsByTime=function(dataset,title,legendProp=list(),ylab="Reaction time(ms)") {
  if(is.null(legendProp$pos))
    legendProp$pos=c(0.8,0.9)
  library(ggplot2)
  #plot data as line graph (mean Data by degree and condition)
  ggplot(dataset,aes(y=reactionTime,x=time, color=condColor, linetype=condLinetype)) + 
    geom_smooth(aes(fill=condColor, group=cond)) +
    labs(x="time(min)",y=ylab,color=legendProp$color,fill=legendProp$color,linetype=legendProp$linetype,shape=legendProp$shape) + 
    theme_classic(base_size=8) + theme(legend.position = legendProp$pos) + 
    scale_colour_discrete(drop=TRUE,limits = levels(dataset$condColor)) + 
    scale_linetype_discrete(drop=TRUE,limits = levels(dataset$condLinetype)) + 
    scale_fill_discrete(drop=TRUE,limits = levels(dataset$condColor))
  ggsave(paste("figs/",title,"LinePlotByCondTime.tiff",sep=""),height=8,width=8,unit="cm",dpi=1000)
}

#combine multiple images into one
combineImages=function(imagesList,rows,columns,outputFile){
  library(magick)
  initImage=image_read(imagesList[1])
  imageRows=rep(initImage,columns)
  imageColumns=rep(initImage,rows)
  counter=1
  for(i in 1:rows){
    for(j in 1:columns){
      imageRows[j]=image_read(imagesList[counter])
      counter=counter+1
    }
    imageRow=image_append(imageRows)
    imageColumns[i]=imageRow
  }
  image=image_append(imageColumns,stack=TRUE)
  image_write(image, path = outputFile)
  gc()
}
