
from __future__ import division

numSlaves = 16
#instance price dollars/hour
slaveCost = 1.4 #r3.4xlarge = 1.4
bigMachine = 6.82 #i2.8xlarge = 6.82
#There's a master
clusterPrice = (numSlaves+1)*slaveCost

"""
Initial times for reference
adam = 199
ups3 = 7*60
dwns3 = 5*60
reorder = 692
sort = (19*60) + 16
var = (35*60) + 49
"""

#times in minutes
adam = 199
ups3 = 7*60
dwns3 = 5*60
reorder = 692
sort = (19*60) + 16
var = (35*60) + 49

time = adam + ups3 + dwns3 + reorder + sort + var
cost = (adam + ups3)*clusterPrice + (dwns3 + reorder + sort + var)*bigMachine

print("time",time)

#for reference, the all-gatk pipeline
gatTime = (148*60)
gatCost = gatTime*bigMachine
print("gatTime", gatTime)

timeCompare = time/gatTime
costCompare = cost/gatCost

print("time ratio comb/gatk", "%.2f" % timeCompare)
print("cost ratio comb/gatk", "%.2f" % costCompare)