#BiSubmodular Maximization Tests
import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import os
import sys
import itertools
import operator
import EmbedUtilities
import Round
from copy import deepcopy
from MainTest import MainTest


class TestfdBiSub(MainTest):

    @staticmethod
    def testObjMethods(sideparams,comp2scale,fqsum):
        """tests matrix based and other objective method
        Args:
           sideparams:
           comp2scale:
        Returns:
           bool: true/false
        """
        [compcount,freqmat,interdom,domains,scales,predomains,prob,algopar,kernmat] = sideparams
        allscales = range(1,2*max(scales)+1)
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        var2index,index = {},0
        for domin in xrange(len(domains)):
            for comp in xrange(compcount):
                for scale in allscales:
                    var2index[(domin,comp,scale)] = index
                    index += 1
        A,b = EmbedUtilities.genCoefs(freqmat,compcount,domains,interdom,var2index,kernmat)
        if prob == "map":
           EmbedUtilities.addPriorCoef(b,algopar["lambda"],predomains,compcount,domains,var2index,"min")  
        for iter in xrange(5):
            newxdict = {(domin,comp,scale): 0.0 for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in allscales}
            alldoms = list(range(len(domains)))
            comp2dominds = {comp:[] for comp in xrange(compcount)}
            for comp in xrange(compcount):
                random.shuffle(alldoms)
                for domin in alldoms:
                    if random.random() < 0.2:
                       newxdict[(domin,comp,comp2scale[comp])] = 1.0
                       comp2dominds[comp].append(domin)       
            obj1 = EmbedUtilities.estFracObjectiveMat(A,b,fqsum,var2index,newxdict)
            obj2 = EmbedUtilities.estFracObjective(newxdict,freqmat,node2dom,allscales,compcount,kernmat)
            if prob == "map":
               priorobj = EmbedUtilities.addPriorObj(comp2dominds,domains,predomains,algopar["lambda"]) 
               obj2 += priorobj
            assert abs(obj1 - obj2) < obj1*1e-08
        return True        
        
    @staticmethod
    def testLPData(sideparams,objstr):
        """tests lp data, script string code etc for the first step
        Args:
          sideparams:
          objstr:
        Returns:
          bool: true or false
        """
        if len(sideparams) == 9:
           [compcount,freqmat,interdom,domains,scales,predomains,prob,algopar,kernmat] = sideparams
        elif len(sideparams) == 8:
           [compcount,freqmat,interdom,domains,predomains,prob,algopar,kernmat] = sideparams  
        #Compares lp coefs in string with the expected ones
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        dom2index = {domains[index]:index for index in xrange(len(domains))}
        var2index, varcount, index2var = {}, 0, {}
        for dom,comp in list(itertools.product(domains,range(compcount))):
            var2index[(dom2index[dom],comp)] = varcount
            index2var[varcount] = (dom2index[dom],comp)
            varcount += 1
        coefmat = np.zeros((varcount,varcount),dtype=np.float)
        impstr = objstr.split("[")[1].split("]")[0]
        for part in impstr.split(" + "):
            splitted = part.split()
            assert len(splitted) == 4 and splitted[2] == "*"
            domin1,comp1 = [int(item) for item in splitted[1].replace("x","").split("_")]
            domin2,comp2 = [int(item) for item in splitted[3].replace("x","").split("_")]
            ind1 = var2index[(domin1,comp1)]
            ind2 = var2index[(domin2,comp2)]
            coefmat[ind1,ind2] += float(splitted[0])/4
            coefmat[ind2,ind1] += float(splitted[0])/4
        assert np.allclose(coefmat.transpose(), coefmat)  
        return True

    @staticmethod
    def getXdict(domains,compcount,scales,comp2scale,comp2dominds):
        """
        """       
        xdict = {(domin,comp,scale):0 for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales}
        for comp in comp2scale.keys():
            for domin in comp2dominds[comp]:
                xdict[(domin,comp,comp2scale[comp])] = 1
        return xdict

    @staticmethod
    def testSDPOutput(sideparams,curcomp2scale,comp2dominds):  
        """tests SDP output including objective val
        Args:
           sideparams:
           curcomp2scale:
           comp2dominds:
        Returns:
           bool: true/false
        """
        [compcount,freqmat,interdom,domains,scales,predomains,prob,algopar,kernmat] = sideparams
        assert len(set(curcomp2scale.keys()) ^ set(range(compcount))) == 0
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        scalelist = range(1,2*max(scales)+1)
        txdict = TestfdBiSub.getXdict(domains,compcount,scalelist,curcomp2scale,comp2dominds)
        step2objval = EmbedUtilities.estFracObjective(txdict,freqmat,node2dom,scalelist,compcount)
        ratio,tobjval = EmbedUtilities.getRatio(freqmat,comp2dominds,curcomp2scale,domains)
        assert abs(tobjval - step2objval) < step2objval*1e-08
        return True
          
    @staticmethod
    def testLPOutput(mainxdict,sideparams,cliques,fracobjval,comp2scale,comp2dominds,algopar):
        """tests LP output including objective val
        Args:
           mainxdict:
           sideparams:
           cliques:
           fracobjval:
           comp2scale,comp2dominds:
           algopar:
        Returns:
           bool: true/false
        """
        if len(sideparams) == 9:
           [compcount,freqmat,interdom,domains,scales,predomains,prob,algopar2,kernmat] = sideparams
        elif len(sideparams) == 8:
           [compcount,freqmat,interdom,domains,predomains,prob,algopar2,kernmat] = sideparams  
        if algopar.has_key("k"):
           for domin in xrange(len(domains)):    
               assert sum([mainxdict[(domin,comp)] for comp in xrange(compcount) if mainxdict.has_key((domin,comp))]) - algopar["k"] < 0.001
               assert len(set(comp for comp in comp2dominds.keys() if domin in comp2dominds[comp])) <= algopar["k"]
        assert TestfdBiSub.testIntersect(comp2dominds,domains)
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        if kernmat == None:
           fqsum = sum([freqmat[in1,in2]**2 for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1])])
        else:
           fqsum = sum([(freqmat[ind1,ind2]**2)*kernmat[ind1,ind2] for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])      
        for index in xrange(compcount):
            for clique in cliques:
                sumval = sum([mainxdict[(domin,comp)] for domin,comp in mainxdict.keys() if comp == index and domains[domin] in clique])
                assert sumval <= 1.001
        if len(sideparams) == 9:        
           scalelist = range(1,2*max(scales)+1)
           txdict = TestfdBiSub.getXdict(domains,compcount,scalelist,comp2scale,comp2dominds)
           step1objval = EmbedUtilities.estFracObjective(txdict,freqmat,node2dom,scalelist,compcount,kernmat)
           if prob == "map":
              step1objval += EmbedUtilities.addPriorFracObj(txdict,domains,predomains,algopar["lambda"])
        ratio,tobjval = EmbedUtilities.getRatio(freqmat,comp2dominds,comp2scale,domains,kernmat)
        if prob == "map":
           tobjval += EmbedUtilities.addPriorObj(comp2dominds,domains,predomains,algopar["lambda"])
        if len(sideparams) == 9:   
           print "afters:"   
           print step1objval,tobjval,fracobjval   
           assert abs(tobjval - step1objval) < step1objval*1e-08
           if prob != "map":
              assert abs(ratio*fqsum - tobjval) < step1objval*1e-08 and tobjval-fracobjval > -0.1
        elif len(sideparams) == 8:
           print "afters:"   
           print tobjval,fracobjval   
           if prob != "map":
              assert abs(ratio*fqsum - tobjval) < step1objval*1e-08 and tobjval-fracobjval > -0.1                 
        return True

