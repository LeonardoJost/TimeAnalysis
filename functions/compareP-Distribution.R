### Analysis of random slopes on power and type1 error rates
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

#this script is used to inspect the distribution of p-values to check
#power and type 1 error rates for different analyses
#especially for comparison of random slopes (because of ongoing discussion?)


#plot p-values of no treatment simulation to check type1 error rate
library(ggplot2)
dataOfSimsControl$pTimeSorted=sort(dataOfSimsControl$pTime)
dataOfSimsControl$pNoTimeSorted=sort(dataOfSimsControl$pNoTime)
dataOfSimsControl$pTimeCovSorted=sort(dataOfSimsControl$pTimeCov)
ggplot(dataOfSimsControl,aes(x=c(0:999))) + 
  geom_line(aes(y=pTimeSorted,color="1")) +
  geom_line(aes(y=pNoTimeSorted,color="2")) +
  geom_line(aes(y=pTimeCovSorted,color="3")) +
  geom_line(aes(y=c(0:999)/999,color="4")) +
  scale_color_manual(
    values=c("1"="red","2"="blue","3"="green","4"="black"), 
    labels = c("time", "no time", "time as covariate", "expected"),
    name="analysis type") +
  labs(x="number of simulations",y="magnitude of p-values") +
  theme_classic() + theme(legend.position = c(0.2,0.8))