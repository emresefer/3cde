#Evaluation utilities
import os
import sys
import math
import numpy as np
from munkres import Munkres
from copy import deepcopy


def getMI(part1,part2,nodecount):
    """returns mi score
    Args:
       part1:
       part2:
       nodecount:
    Returns:
       miscore:
    """
    pdist1 = {index: float(len(part1[index]))/nodecount for index in xrange(len(part1))}
    pdist2 = {index: float(len(part2[index]))/nodecount for index in xrange(len(part2))}
    miscore = 0.0
    for ind1 in xrange(len(part1)):
         for ind2 in xrange(len(part2)):
             p = float(len(set(part1[ind1]).intersection(set(part2[ind2]))))/nodecount
             if p != 0:
                miscore += p*(math.log(p/(pdist1[ind1]*pdist2[ind2]),2.0))
    return miscore

def getVI(part1,part2,nodecount):
    """returns vi score
    Args:
       part1:
       part2:
       nodecount:
    Returns:
       miscore:
    """
    pdist1 = {index: float(len(part1[index]))/nodecount for index in xrange(len(part1))}
    pdist2 = {index: float(len(part2[index]))/nodecount for index in xrange(len(part2))}
    viscore = getEntropy(pdist1) + getEntropy(pdist2) - 2.0*getMI(part1,part2,nodecount)
    if abs(viscore) < 0.001:
       return 0.0  
    return viscore
    
def getEntropy(probdist):
    """gets entropy
    Args:
       probdist:
    Returns:
       entsum:
    """
    entsum = sum([-1.0*probdist[item]*math.log(probdist[item],2.0) for item in probdist.keys() if probdist[item] != 0])
    return entsum

def addEmptyClusters(domains,allnodes):
    """adds empty clusters
    """        
    part1 = sorted([range(start,end+1) for start,end in domains], key=lambda item: item[0])
    curin,part2 = 1,[]
    for clust in part1:
        start,end = min(clust),max(clust)
        if curin <= start:
           part2.append(range(curin,start))
        curin = end+1
    if curin <= len(allnodes):
       part2.append(range(curin,len(allnodes)+1))
    if TESTMODE:   
       for clust1 in part1:
           for clust2 in part2:
               assert len(set(clust1) & set(clust2)) == 0
    part1.extend(part2)
    part1 = [part for part in part1 if len(part)>0]
    if TESTMODE:  
       for ind1 in xrange(len(part1)):
           range1 = set(part1[ind1]) 
           for ind2 in xrange(ind1+1,len(part1)):
               range2 = set(part1[ind2])
               assert len(range1 & range2) == 0
       assert len(set(node for part in part1 for node in part) ^ set(allnodes)) == 0        
    return part1

TESTMODE = False
def matchScore(truecomp2domain,truecomp2scale,estcomp2domain,estcomp2scale,freqmat,sctype):
    """bipartite matching based score
    Args:
       truecomp2domain:
       truecomp2scale:
       estcomp2domain:
       estcomp2scale:
       freqmat:
       sctype:
    Returns:
       score:
       align: alignment tuple
    """
    assert sctype in ["error","vi"]
    compcount = len(truecomp2domain.keys()) 
    scoremat = np.zeros((compcount,compcount),dtype=np.float)
    allnodes = range(1,np.shape(freqmat)[0]+1)
    for truecomp in xrange(compcount):
        truedomains = truecomp2domain[truecomp]
        if sctype == "vi":
           part1 = addEmptyClusters(truedomains,allnodes)
        elif sctype == "error":
           truemat = np.zeros(np.shape(freqmat),dtype=np.float)
           for start,end in truedomains:
               truemat[start-1:end,start-1:end] += truecomp2scale[truecomp]
        for estcomp in xrange(compcount):
            estdomains = estcomp2domain[estcomp]
            if sctype == "vi":
               part2 = addEmptyClusters(estdomains,allnodes)
               scoremat[truecomp,estcomp] = getVI(part1,part2,len(allnodes))
            elif sctype == "error":             
               estmat = np.zeros(np.shape(freqmat),dtype=np.float)
               for start,end in estdomains:
                   estmat[start-1:end,start-1:end] += estcomp2scale[estcomp]
               scoremat[truecomp,estcomp] = np.sum(np.fabs(truemat-estmat))
    m = Munkres()
    scoremat2 = deepcopy(scoremat)  
    indexes = m.compute(scoremat2)
    score = sum([scoremat[row,column] for row, column in indexes])
    if sctype == "vi":
       avgscore = score/(float(compcount)*math.log(len(allnodes),2.0))
    elif sctype == "error":
       avgscore = score/(float(compcount)*(len(allnodes)**2))
    if TESTMODE:                              
       if sctype == "vi":
          for ind1 in xrange(np.shape(scoremat)[0]):
              for ind2 in xrange(np.shape(scoremat)[1]):
                  assert scoremat[ind1,ind2] >= 0 and scoremat[ind1,ind2] <= math.log(len(allnodes),2.0)
       elif sctype == "error":
          for ind1 in xrange(np.shape(scoremat)[0]):
              for ind2 in xrange(np.shape(scoremat)[1]):
                  assert scoremat[ind1,ind2] >= 0      
    return avgscore,indexes
