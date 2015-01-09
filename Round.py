#Rounding related methods
import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import operator
import os
import sys
import itertools
import EmbedUtilities
from copy import deepcopy
sys.path.append("Tests")
from MainTest import MainTest

def getDom2Dom(domains):
    """get domain to domain interaction but start position include
    """    
    dom2dom = {(domin1,domin2): False for domin1 in xrange(len(domains))  for domin2 in xrange(len(domains))}
    for domin1 in xrange(len(domains)):
        dom1 = domains[domin1]
        for domin2 in xrange(len(domains)):
            dom2 = domains[domin2]
            if domin1 == domin2:
               continue 
            if dom2[0] <= dom1[0] and dom2[1] >= dom1[0]:
               dom2dom[(domin1,domin2)] = True
    return dom2dom


def partBack(comp2scale,scales):
    """partition back scales
    Args:
       comp2scale:
       scales:
    Returns:
       comp2sclist: comp2 scale list
    """
    sortedscales = sorted(scales,reverse=True)
    comp2sclist = {comp:[] for comp in comp2scale.keys()}
    tcomp2scale = {comp:comp2scale[comp] for comp in comp2scale.keys()}
    for comp in tcomp2scale.keys():
        for scale in sortedscales:
            if tcomp2scale[comp] / scale >= 1:
               tcomp2scale[comp] -= scale
               comp2sclist[comp].append(scale)
    for comp in comp2sclist.keys():
        assert sum(comp2sclist[comp]) == comp2scale[comp]
    return comp2sclist


def runIterRoundPartial(iterparams,fqsum,method,k):
    """rounds iter partial rounding
    Args:
       iterparams:
       fqsum:
       method:
       k:
    Returns:
       curobjval:
       cursol:
    """
    [compcount,freqmat,domains,scales,A,b,comp2scale,dom2dom,comp2sclist,xdict,var2index,node2dom] = iterparams
    comp2R = {comp:set() for comp in xrange(compcount)}
    allkeys = list(xdict.keys())
    random.shuffle(allkeys)
    for (domin,comp) in allkeys:
        if random.random() < method(xdict[(domin,comp)]):
           comp2R[comp].add(domin) 
    comp2mark = {comp:set() for comp in xrange(compcount)}
    for (domin,comp) in allkeys:
        for domin2 in comp2R[comp]:
            if dom2dom[(domin,domin2)]:
               comp2mark[comp].add(domin)
               break
    comp2R = {comp: comp2R[comp]-comp2mark[comp] for comp in comp2mark.keys()}
    dom2comp = {}
    for comp in comp2R.keys():
        for domin in comp2R[comp]:
            if not dom2comp.has_key(domin):
               dom2comp[domin] = set() 
            dom2comp[domin].add(comp)
    #real constraint part
    for domin in dom2comp.keys():
        if len(dom2comp[domin]) > k:
           comps = list(dom2comp[domin])
           random.shuffle(comps)
           for comp in comps[k:]:
               comp2R[comp].remove(domin)
    binxdict = {(domin,comp,scale): 0.0 for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales}
    for comp in comp2R.keys():
        for domin in comp2R[comp]:
            for scale in comp2sclist[comp]:
                binxdict[(domin,comp,scale)] = 1
    binxdict2 = {(domin,comp,comp2scale[comp]): 0.0 for domin in xrange(len(domains)) for comp in xrange(compcount)}
    for comp in comp2R.keys():
        for domin in comp2R[comp]:
            binxdict2[(domin,comp,comp2scale[comp])] = 1
    curobjval = EmbedUtilities.estFracObjectiveMat(A,b,fqsum,var2index,binxdict2)
    
    for comp in comp2R.keys():
        alldoms = list(comp2R[comp])
        for dom in alldoms:
            assert dom not in comp2mark[comp]
    return curobjval,comp2R 


def RoundPipage(xdict,ydict,objval,yset):
    """does pipage rounding
    Args:
       xdict:
       ydict:
       objval:
       yset:
    Returns:
       sol:
    """
    sol = {}
    print objval
    print "sol: ",len(xdict.keys()), len(ydict.keys()),len(yset)
    xvals = set(xdict.values())
    yvals = set(ydict.values())
    for domin in xrange(len(domains)):
        for comp in xrange(compcount):
            mysum = sum([xdict[(domin,comp,scale)] for scale in scales if xdict.has_key((domin,comp,scale))])
            assert abs(mysum-1) < 0.01
    print xvals
    print yvals 
    print "done"
    #fqsum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    #print fqsum
    #print objval-fqsum
    return sol


def runIterRound(iterparams,fqsum,method,kernmat):
    """rounds iter rounding
    Args:
       iterparams:
       fqsum:
       method:
    Returns:
       curobjval:
       cursol:
    """
    [compcount,freqmat,domains,predomains,scales,A,b,comp2scale,dom2dom,comp2sclist,xdict,var2index,node2dom,algopar2] = iterparams
    comp2R = {comp:set() for comp in xrange(compcount)}
    allkeys = list(xdict.keys())
    random.shuffle(allkeys)
    for (domin,comp) in allkeys:
        if random.random() < method(xdict[(domin,comp)]):
           comp2R[comp].add(domin) 
    comp2mark = {comp:set() for comp in xrange(compcount)}
    for (domin,comp) in allkeys:
        for domin2 in comp2R[comp]:
            if dom2dom[(domin,domin2)]:
               comp2mark[comp].add(domin)
               break
    comp2R = {comp: comp2R[comp]-comp2mark[comp] for comp in comp2mark.keys()}
    binxdict = {(domin,comp,scale): 0.0 for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales}
    for comp in comp2R.keys():
        for domin in comp2R[comp]:
            for scale in comp2sclist[comp]:
                binxdict[(domin,comp,scale)] = 1
    binxdict2 = {(domin,comp,comp2scale[comp]): 0.0 for domin in xrange(len(domains)) for comp in xrange(compcount)}
    for comp in comp2R.keys():
        for domin in comp2R[comp]:
            binxdict2[(domin,comp,comp2scale[comp])] = 1
    curobjval = EmbedUtilities.estFracObjectiveMat(A,b,fqsum,var2index,binxdict2)
    if False:
       curobjval2 = EmbedUtilities.estFracObjective(binxdict,freqmat,node2dom,scales,compcount,kernmat)
       curobjval3 = EmbedUtilities.getRatio(freqmat,comp2R,comp2scale,domains,kernmat)[1]
       if algopar2.has_key("lambda"):
          priorobj = EmbedUtilities.addPriorObj(comp2R,domains,predomains,algopar2["lambda"])
          curobjval2 += priorobj
          curobjval3 += priorobj
       assert abs(curobjval-curobjval2) < curobjval*1e-09 and abs(curobjval-curobjval3) < curobjval*1e-09
       for comp in comp2R.keys():
           alldoms = list(comp2R[comp])
           for dom in alldoms:
               assert dom not in comp2mark[comp]
    return curobjval,comp2R  


def roundCRPartial(xdict,sideparams,comp2scale,k,count):
    """round based on contention resolution partial case
    Args:
       xdict:
       sideparams:
       comp2scale:
       scales:
       k:
       count: round count
    Returns:
       comp2dominds:
    """
    [compcount,freqmat,interdom,domains,scales,predomains,prob,algopar] = sideparams
    for domin in xrange(len(domains)):
        assert sum([xdict[(domin,comp)] for comp in xrange(compcount) if xdict.has_key((domin,comp))])-k < 0.001
    comp2sclist = partBack(comp2scale,scales)
    var2index = EmbedUtilities.mapVarsSpec(domains,compcount,comp2scale)[0]
    A,b = EmbedUtilities.genCoefs(freqmat,compcount,domains,interdom,var2index)
    if prob == "map":
       EmbedUtilities.addPriorCoef(b,algopar["lambda"],predomains,compcount,scales,domains,var2index,"min")  
    fqsum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    dom2dom = getDom2Dom(domains)
    iterparams = [compcount,freqmat,domains,scales,A,b,comp2scale,dom2dom,comp2sclist,xdict,var2index,node2dom] 
    minobjval,bestsol = fqsum, None
    for param in ["exp","normal"]:
        if param == "exp":
           method = lambda val: 1.0-math.exp(-1.0*val)
        elif param == "normal":
           method = lambda val: val
        for index in xrange(count):
            curobjval, cursol = runIterRoundPartial(iterparams,fqsum,method,k)
            for domin in xrange(len(domains)):    
                assert len(set(comp for comp in cursol.keys() if domin in cursol[comp])) <= k
            if curobjval < minobjval:
               minobjval = curobjval
               bestsol = deepcopy(cursol)
    return bestsol
    

def runIterRoundfrac(iterparams,fqsum,method,kernmat):
    """rounds iter rounding
    Args:
       iterparams:
       fqsum:
       method:
    Returns:
       curobjval:
       cursol:
    """
    [compcount,freqmat,domains,predomains,A,b,comp2scale,dom2dom,xdict,var2index,node2dom,algopar2] = iterparams
    comp2R = {comp:set() for comp in xrange(compcount)}
    allkeys = list(xdict.keys())
    random.shuffle(allkeys)
    for (domin,comp) in allkeys:
        if random.random() < method(xdict[(domin,comp)]):
           comp2R[comp].add(domin) 
    comp2mark = {comp:set() for comp in xrange(compcount)}
    for (domin,comp) in allkeys:
        for domin2 in comp2R[comp]:
            if dom2dom[(domin,domin2)]:
               comp2mark[comp].add(domin)
               break
    comp2R = {comp: comp2R[comp]-comp2mark[comp] for comp in comp2mark.keys()}
    binxdict2 = {(domin,comp,comp2scale[comp]): 0.0 for domin in xrange(len(domains)) for comp in xrange(compcount)}
    for comp in comp2R.keys():
        for domin in comp2R[comp]:
            binxdict2[(domin,comp,comp2scale[comp])] = 1
    curobjval = EmbedUtilities.estFracObjectiveMat(A,b,fqsum,var2index,binxdict2)
    if False:
       curobjval2 = EmbedUtilities.getRatio(freqmat,comp2R,comp2scale,domains,kernmat)[1]
       if algopar2.has_key("lambda"):
          priorobj = EmbedUtilities.addPriorObj(comp2R,domains,predomains,algopar2["lambda"])
          curobjval2 += priorobj
       assert abs(curobjval-curobjval2) < curobjval*1e-10
       for comp in comp2R.keys():
           alldoms = list(comp2R[comp])
           for dom in alldoms:
               assert dom not in comp2mark[comp]
    return curobjval,comp2R  


def roundCRfrac(xdict,sideparams,comp2scale,algopar,count):
    """round based on contention resolution 1, 1/e for fractional case
    Args:
       xdict:
       sideparams:
       comp2scale:
       scales:
       algopar:
       count: round count
    Returns:
       comp2dominds:
    """
    if algopar.has_key("k"):
       return roundCRPartial(xdict,sideparams,comp2scale,algopar["k"],count) 
    [compcount,freqmat,interdom,domains,predomains,prob,algopar2,kernmat] = sideparams
    var2index = EmbedUtilities.mapVarsSpec(domains,compcount,comp2scale)[0]
    A,b = EmbedUtilities.genCoefs(freqmat,compcount,domains,interdom,var2index,kernmat)
    if prob == "map":
       EmbedUtilities.addPriorCoef(b,algopar["lambda"],predomains,compcount,domains,var2index,"min") 
    if kernmat == None:  
       fqsum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    else:
       fqsum = sum([(freqmat[ind1,ind2]**2)*kernmat[ind1,ind2] for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])    
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    dom2dom = getDom2Dom(domains)
    iterparams = [compcount,freqmat,domains,predomains,A,b,comp2scale,dom2dom,xdict,var2index,node2dom,algopar] 
    minobjval,bestsol = max(fqsum,10000000000.0), None
    for param in ["exp","normal"]:
        if param == "exp":
           method = lambda val: 1.0-math.exp(-1.0*val)
        elif param == "normal":
           method = lambda val: val
        for index in xrange(count):
            curobjval, cursol = runIterRoundfrac(iterparams,fqsum,method,kernmat)
            if curobjval < minobjval:
               minobjval = curobjval
               bestsol = deepcopy(cursol)
    assert bestsol != None                   
    return bestsol,minobjval


def roundCR(xdict,sideparams,comp2scale,algopar,count):
    """round based on contention resolution 1, 1/e
    Args:
       xdict:
       sideparams:
       comp2scale:
       scales:
       algopar:
       count: round count
    Returns:
       comp2dominds:
    """
    if algopar.has_key("k"):
       return roundCRPartial(xdict,sideparams,comp2scale,algopar["k"],count) 
    [compcount,freqmat,interdom,domains,scales,predomains,prob,algopar2,kernmat] = sideparams
    comp2sclist = partBack(comp2scale,scales)
    var2index = EmbedUtilities.mapVarsSpec(domains,compcount,comp2scale)[0]
    A,b = EmbedUtilities.genCoefs(freqmat,compcount,domains,interdom,var2index,kernmat)
    if prob == "map":
       EmbedUtilities.addPriorCoef(b,algopar["lambda"],predomains,compcount,domains,var2index,"min") 
    if kernmat == None:  
       fqsum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    else:
       fqsum = sum([(freqmat[ind1,ind2]**2)*kernmat[ind1,ind2] for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])    
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    dom2dom = getDom2Dom(domains)
    iterparams = [compcount,freqmat,domains,predomains,scales,A,b,comp2scale,dom2dom,comp2sclist,xdict,var2index,node2dom,algopar] 
    minobjval,bestsol = max(fqsum,10000000000.0), None
    for param in ["exp","normal"]:
        if param == "exp":
           method = lambda val: 1.0-math.exp(-1.0*val)
        elif param == "normal":
           method = lambda val: val
        for index in xrange(count):
            curobjval, cursol = runIterRound(iterparams,fqsum,method,kernmat)
            if curobjval < minobjval:
               minobjval = curobjval
               bestsol = deepcopy(cursol)
    assert bestsol != None           
    return bestsol,minobjval

   
class Round():

    @staticmethod
    def greedy3Rand(xdict,sideparams):
        """greedy rounding iterative loop
        Args:
           xdict:
           sideparams:
        Returns:
           comp2dominds: domain indices of each component
           comp2scale: component to scale map
        """
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        comp2dominds = {comp: set() for comp in xrange(compcount)}
        comp2scale = {comp: None for comp in xrange(compcount)}
        scales,compcount,freqmat,interdom,domains = sideparams
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        sorted_x = sorted(xdict.iteritems(), key=operator.itemgetter(1),reverse=True)
        cursol = {(domin,comp,scale): 0 for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales}
        curobjval = EmbedUtilities.estFracObjective(cursol,freqmat,node2dom,scales,compcount)
        for (domin,comp,scale),val in sorted_x:
            cursol[(domin,comp,scale)] = 1
            if violateCons(cursol,compcount,domains):
               cursol[(domin,comp,scale)] = 0
               continue
            estobjval = EmbedUtilities.estFracObjective(cursol,freqmat,node2dom,scales,compcount)
            if estobjval <= curobjval:
               curobjval = estobjval
            else:
               cursol[(domin,comp,scale)] = 0
        assert violateCons(cursol,compcount,domains)
        for (domin,comp,scale) in cursol.keys():
            comp2dominds[comp].add(domin)
            comp2scale[comp].add(scale)
        return comp2dominds,comp2scale
      

    @staticmethod
    def greedy2Rand(xdict,sideparams):
        """greey rounding 2
        Args:
           xdict:
           sideparams:
        Returns:
           comp2dominds: domain indices of each component
           comp2scale: component to scale map
        """
        scales,compcount,freqmat,interdom,domains = sideparams
        comp2dominds = {comp: set() for comp in xrange(compcount)}
        comp2scale = {comp: None for comp in xrange(compcount)}
        #divxdict = {(domin,comp,scale): xdict[(domin,comp,scale)]/ydict[(comp,scale)] for domin,comp,scale in xdict.keys()}
        divxdict = {(domin,comp,scale): xdict[(domin,comp,scale)] for domin,comp,scale in xdict.keys()}
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
         
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        soldict = {(domin,comp,scale): 0.0 for (domin,comp,scale) in xdict.keys()}  
        curobjval = EmbedUtilities.estFracObjective(soldict,freqmat,node2dom,scales,compcount)
        allcomps = range(compcount)
        random.shuffle(allcomps)
        for comp in allcomps:
            print "before {0} {1}".format(comp,curobjval)
            bestscaleobj = curobjval
            bestscaledict = {}
            bestscale = None
            bestdominds = None
            for scale in scales:
                scaledict = deepcopy(soldict)
                
                x = {domin: divxdict[(domin,comp,scale)] for domin in xrange(len(domains))}
                curdoms = set()
                remdomset = set(range(len(domains)))
                imprflag = True
                while imprflag and len(remdomset) > 0:
                   imprflag = False
                   remdoms = list(remdomset) 
                   random.shuffle(remdoms)
                   p = random.random()
                   found, prob = remdoms[0], 0.0
                   divsum = sum(x.values())
                   #if divsum == 0:
                   #   break 
                   normprob = {domin: x[domin]/divsum for domin in remdomset}
                   for dom in remdoms:
                       prob += normprob[dom]
                       if prob >= p:
                          found = dom
                          imprflag = True
                          break
                   if imprflag:
                      delset = set(domin for domin in remdomset if intersect(domains[domin],domains[found]))
                      remdomset -= delset
                      for deldom in delset:
                          del x[deldom]
                      curdoms.add(found)
            
                #x = {domin: divxdict[(domin,comp,scale)] for domin in xrange(len(domains))}     
                #sorted_x = sorted(x.iteritems(), key=operator.itemgetter(1),reverse=True)
                #curdoms = set()
                #while len(sorted_x) > 0:
                #   found,foundval = sorted_x[0]
                #   intflag = False
                #   for curdom in curdoms:
                #       if intersect(domains[curdom],domains[found]):
                #          intflag = True
                #          break
                #   del sorted_x[0]   
                #   if not intflag:
                #      curdoms.add(found)

                      
                for curdomin in curdoms:
                    scaledict[(curdomin,comp,scale)] = 1.0
                scaleobjval = EmbedUtilities.estFracObjective(scaledict,freqmat,node2dom,scales,compcount)     
                if scaleobjval <= bestscaleobj:
                   bestscaleobj = scaleobjval
                   bestscaledict = deepcopy(scaledict)
                   bestscale = scale
                   bestdominds = set(curdoms)
            if bestscale != None:
               soldict = deepcopy(bestscaledict)
               curobjval = bestscaleobj
               comp2scale[comp] = bestscale
               comp2dominds[comp] = set(bestdominds)
            else:
               break
        print "lastsol ",EmbedUtilities.estFracObjective(soldict,freqmat,node2dom,scales,compcount)  
        emdict = {(domin,comp,scale): 0.0 for (domin,comp,scale) in xdict.keys()}  
        print "emptysol ",EmbedUtilities.estFracObjective(emdict,freqmat,node2dom,scales,compcount)

        print np.shape(freqmat)[0]*np.shape(freqmat)[1]
        print math.sqrt(EmbedUtilities.estFracObjective(emdict,freqmat,node2dom,scales,compcount) / float(np.shape(freqmat)[0]*np.shape(freqmat)[1]))
        print "max: ", np.amax(freqmat)
        print "min: ", np.amin(freqmat)
        print "info"
        for comp in comp2scale.keys():
            print comp,comp2scale[comp]
        print allcomps    
        exit(1)
    
        for comp in xrange(compcount):
            for domin,comp2,scale in newxdict.keys():
                if comp2 == comp:
                   print domin,comp,scale,newxdict[(domin,comp,scale)] 
            break
    
    @staticmethod
    def roundRand(xdict,ydict,sideparams):
        """randomized rounding
        Args:
           xdict:
           ydict:
           sideparams:
        Returns:
           comp2dominds: domain indices of each component
           comp2scale: component to scale map
        """
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        compcount,scales,domains = sideparams
        comp2dominds = {comp: set() for comp in xrange(compcount)}
        comp2scale = {comp: None for comp in xrange(compcount)}
        #comp scale assignment by randomized round
        for comp in xrange(compcount):
            p = random.random()
            tot = 0.0
            curscales = list(scales)
            random.shuffle(curscales)
            selscale = curscales[0]
            for scale in curscales:
                tot += ydict[(comp,scale)]
                if tot >= p:
                   selscale = scale 
                   break 
            comp2scale[comp] = selscale
    
        domcom2val = {(domin,comp): sum([xdict[(dom,tcomp,scale)] for dom,tcomp,scale in xdict.keys() if tcomp==comp and dom==domin]+[0.0]) for domin in xrange(len(domains)) for comp in xrange(compcount)}
        com2dom2val = {comp:{} for comp in xrange(compcount)}
        for domin,comp in domcom2val.keys():
            com2dom2val[comp][domin] = domcom2val[(domin,comp)]
        #comp domain assignment by modified randomized rounding
        for comp in xrange(compcount):
            domset = set()
            remdomset = set(com2dom2val[comp].keys())
            imprflag = True
            while imprflag and len(remdomset) > 0:
               imprflag = False
               curdoms = list(remdomset) 
               random.shuffle(curdoms)
               p = random.random()
               found, prob = curdoms[0], 0.0
               divsum = float(sum(com2dom2val[comp][remdom] for remdom in remdomset))
               normprob = {dom: com2dom2val[comp][dom]/divsum for dom in remdomset}
               for dom in curdoms:
                   prob += normprob[dom]
                   if prob >= p:
                      found = dom
                      imprflag = True
                      break
               if imprflag:
                  delset = set(domin for domin in remdomset if intersect(domains[domin],domains[found]))
                  remdomset -= delset
                  domset.add(found)
            comp2dominds[comp] = set(domset)
        return comp2dominds, comp2scale

    @staticmethod 
    def greedyRand(xdict,ydict,sideparams):
        """randomized rounding
        Args:
           xdict:
           ydict:
           sideparams:
        Returns:
           compdict: domains of each component
           comp2s: component to scale map
        """
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        compcount,scales,domains = sideparams
        comp2dominds = {comp: set() for comp in xrange(compcount)}
        comp2scale = {comp: None for comp in xrange(compcount)}
    
        #greedy comp scale assignment
        #maxvals = {comp: 0.0 for comp in xrange(compcount)}
        #for comp,scale in ydict.keys():
        #    if ydict[(comp,scale)] >= maxvals[comp]:
        #       maxvals[comp] = ydict[(comp,scale)]
        #       comp2scale[comp] = scale
                                 
        domcom2val = {(domin,comp): sum([xdict[(dom,tcomp,scale)] for dom,tcomp,scale in xdict.keys() if tcomp==comp and dom==domin]+[0.0]) for domin in xrange(len(domains)) for comp in xrange(compcount)}
        com2dom2val = {comp:{} for comp in xrange(compcount)}
        for domin,comp in domcom2val.keys():
            com2dom2val[comp][domin] = domcom2val[(domin,comp)]
        #comp domain assignment by modified randomized rounding
        for comp in xrange(compcount):
            domset = set()
            remdomset = set(com2dom2val[comp].keys())
            imprflag = True
            while imprflag and len(remdomset) > 0:
               imprflag = False
               curdoms = list(remdomset) 
               random.shuffle(curdoms)
               p = random.random()
               found, prob = curdoms[0], 0.0
               divsum = float(sum(com2dom2val[comp][remdom] for remdom in remdomset))
               normprob = {dom: com2dom2val[comp][dom]/divsum for dom in remdomset}
               for dom in curdoms:
                   prob += normprob[dom]
                   if prob >= p:
                     found = dom
                     imprflag = True
                     break
               if imprflag:
                  delset = set(domin for domin in remdomset if intersect(domains[domin],domains[found]))
                  remdomset -= delset
                  domset.add(found)
            comp2dominds[comp] = set(domset)
        return comp2dominds, comp2scale

