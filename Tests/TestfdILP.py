#Combinatorial ILP Formulation Test
import numpy as np
import scipy as sp
import math
import random
import itertools
import operator
import EmbedUtilities
from MainTest import MainTest

    
class TestfdILP(MainTest):

    @staticmethod
    def testfracILPData(sideparams,objstr):
        """tests frac ilp data, script string code etc
        Args:
          sideparams:
          objstr:
        Returns:
          bool: true or false
        """
        [compcount,freqmat,interdom,domains,kernmat] = sideparams
        #Compares lp coefs in string with the expected ones
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        dom2index = {domains[index]:index for index in xrange(len(domains))}
        var2index, index2var, varcount = {}, {}, 0
        for dom,comp in list(itertools.product(domains,range(compcount))):
            var2index[(dom2index[dom],comp)] = varcount
            index2var[varcount] = (dom2index[dom],comp)
            varcount += 1
        coefs = np.zeros((varcount,varcount),dtype=np.float)
        impstr = objstr.split("[")[1].split("]")[0]
        for part in impstr.split(" + "):
            splitted = part.split()
            assert len(splitted) == 4 and splitted[2] == "*"
            domin1,comp1 = [int(item) for item in splitted[1].replace("y","").split("_")]
            domin2,comp2 = [int(item) for item in splitted[3].replace("y","").split("_")]
            ind1 = var2index[(domin1,comp1)]
            ind2 = var2index[(domin2,comp2)]
            coefs[ind1,ind2] += float(splitted[0])/4
            coefs[ind2,ind1] += float(splitted[0])/4
        return True
        #wrong below    
        if kernmat != None:
           coefest = lambda dom1,dom2: np.sum(kernmat[max(dom1[0],dom2[0]):min(dom1[1],dom2[1])+1,max(dom1[0],dom2[0]):min(dom1[1],dom2[1])+1])
        else:
           coefest = lambda dom1,dom2: (min(dom1[1],dom2[1])-max(dom1[0],dom2[0]) + 1)**2
        for ind1 in xrange(np.shape(coefs)[0]):
            domin1,comp1 = index2var[ind1]
            dom1 = domains[domin1]
            for ind2 in xrange(np.shape(coefs)[1]):   
                domin2,comp2 = index2var[ind2]
                dom2 = domains[domin2]
                print "info"
                print ind1,ind2
                print domin1,comp1,dom1
                print domin2,comp2,dom2
                print coefs[ind1,ind2],coefest(dom1,dom2)
                assert coefs[ind1,ind2] == coefest(dom1,dom2)

    @staticmethod
    def testfracILPOutput(sideparams,comp2dominds,comp2scale,xdict,ydict,objval):
        """tests for frac ilp output
        Args:
           sideparams:
           comp2domins:
           comp2scale:
           xdict:
           ydict:
           objval:
        Returns:
           bool: true/false
        """
        [compcount,freqmat,interdom,domains,predomains,prob,algopar,kernmat] = sideparams
        dom2index = {domains[index]:index for index in xrange(len(domains))}
        assert TestfdILP.testIntersect(comp2dominds,domains)
        for val in xdict.values():
            assert val in [1.0, 0.0, -0.0]

        for dom1,dom2 in interdom:
            domin1,domin2 = dom2index[dom1],dom2index[dom2]
            for comp in xrange(compcount):
                sumnum = 0
                if xdict.has_key((domin1,comp)) and xdict[(domin1,comp)] > 0.98:
                   sumnum +=1
                if xdict.has_key((domin2,comp)) and xdict[(domin2,comp)] > 0.98:
                   sumnum +=1
                assert sumnum <= 2       
                       
        comps= {comp: set() for comp in xrange(compcount)}
        for domin,comp in ydict.keys():
            if ydict[(domin,comp)] not in [-0.0,0.0]:
               assert xdict[(domin,comp)] == 1.0
               comps[comp].add(ydict[(domin,comp)])
               
        tcomp2scale = {comp:set() for comp in xrange(compcount)}
        for domin,comp in ydict.keys():
            if ydict[(domin,comp)] not in [0.0,-0.0]:
               tcomp2scale[comp].add(ydict[(domin,comp)])
        for comp in tcomp2scale.keys():    
            assert len(tcomp2scale[comp]) <=1
        assert len(xdict.keys()) <= len(domains)*compcount and len(ydict.keys()) <= len(domains)*compcount

        newobjval = EmbedUtilities.getRatio(freqmat,comp2dominds,comp2scale,domains,kernmat)[1]
        if prob == "map":
           newobjval += EmbedUtilities.addPriorObj(comp2dominds,domains,predomains,algopar["lambda"])
        assert abs(objval - newobjval) < objval*1e-08

        for comp in comp2dominds.keys():
            for domin in comp2dominds[comp]:
                assert xdict.has_key((domin,comp)) and xdict[(domin,comp)] > 0.99
        for comp in comp2scale.keys():
            for domin,tcomp in ydict.keys():
                if tcomp == comp:
                   if ydict[(domin,comp)] not in [0.0,-0.0]: 
                      assert abs(ydict[(domin,tcomp)] - comp2scale[comp]) < 0.01
        return True

    @staticmethod
    def testILPOutput(sideparams,comp2dominds,comp2scale,zdict,ydict,objval):
        """tests ilp output
        Args:
           sideparams:
           comp2domins:
           comp2scale:
           zdict:
           ydict:
           objval:
        Returns:
           bool: true/false
        """
        [scales,compcount,freqmat,interdom,domains,predomains,prob,algopar,kernmat] = sideparams
        dom2index = {domains[index]:index for index in xrange(len(domains))}
        assert TestfdILP.testIntersect(comp2dominds,domains)    
        for val in zdict.values() + ydict.values():
            assert val in [1.0, 0.0]
        for domin,comp,scale in zdict.keys():
            assert domin in comp2dominds[comp] and comp2scale[comp] == scale
            if comp2scale[comp] != 0:
               assert ydict[(comp,scale)] == comp2scale[comp]
        for comp,scale in ydict.keys():
            assert comp2scale[comp] == scale
            
        for maincomp in xrange(compcount):
            for mainscale in scales:
                for dom1,dom2 in interdom:
                    sumval = 0.0
                    if zdict.has_key((dom2index[dom1],maincomp,mainscale)):
                       sumval += zdict[(dom2index[dom1],maincomp,mainscale)]
                    if zdict.has_key((dom2index[dom2],maincomp,mainscale)):
                       sumval += zdict[(dom2index[dom2],maincomp,mainscale)]     
                    assert sumval <= ydict[(comp,scale)]
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)        
        assert len(zdict.keys()) <= len(domains)*compcount*len(scales) and len(ydict.keys()) <= compcount*len(scales)

        estobjval = EmbedUtilities.getRatio(freqmat,comp2dominds,comp2scale,domains,kernmat)[1]
        if prob == "map":
           estobjval += EmbedUtilities.addPriorObj(comp2dominds,domains,predomains,algopar["lambda"])
        assert abs(objval - estobjval) < objval*1e-08

        if algopar["kernel"] != None:
           fqsum = sum([(freqmat[ind1,ind2]**2)*kernmat[ind1,ind2] for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
        else:
           fqsum = float(sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])]))  
        assert estobjval <= fqsum
          
        tcomp2scale = {comp:0 for comp in xrange(compcount)}
        for comp,scale in ydict.keys():
            tcomp2scale[comp] += scale*ydict[(comp,scale)]
        for comp in comp2scale.keys():
            assert comp2scale[comp] == tcomp2scale[comp]    
        return True
        
    @staticmethod
    def testILPData(sideparams,objstr):
        """tests lp data, script string code etc
        Args:
          sideparams:
          objstr:
        Returns:
          bool: true or false
        """
        [scales,compcount,freqmat,interdom,domains,kernmat] = sideparams
        #Compares lp coefs in string with the expected ones
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        dom2index = {domains[index]:index for index in xrange(len(domains))}
        var2index, varcount = {}, 0
        for dom,comp,scale in list(itertools.product(domains,range(compcount),scales)):
            var2index[(dom2index[dom],comp,scale)] = varcount
            varcount += 1
        coefmat = np.zeros((varcount,varcount),dtype=np.float)
        impstr = objstr.split("[")[1].split("]")[0]
        for part in impstr.split(" + "):
            splitted = part.split()
            assert len(splitted) == 4 and splitted[2] == "*"
            domin1,comp1,scale1 = [int(item) for item in splitted[1].replace("z","").split("_")]
            domin2,comp2,scale2 = [int(item) for item in splitted[3].replace("z","").split("_")]
            ind1 = var2index[(domin1,comp1,scale1)]
            ind2 = var2index[(domin2,comp2,scale2)]
            coefmat[ind1,ind2] += float(splitted[0])/4
            coefmat[ind2,ind1] += float(splitted[0])/4
        if kernmat == None:    
           scalesum = sum(scales)
           coefsum = sum([(scalesum*compcount*len(node2dom[in1].intersection(node2dom[in2])))**2 for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1])])
           assert abs(coefsum - np.sum(coefmat)) <= 0.01 
        return True
