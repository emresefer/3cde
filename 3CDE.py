#BiSubmodular Maximization Formulation of Frequency Deconvolution
import networkx as nx
import numpy as np
import math
import random
import os
import sys
import time
import string
sys.path.append("lib")
import itertools
import operator
import EmbedUtilities
import Bls_Sdp
import Round
from copy import deepcopy
sys.path.append("Tests")
from TestfdBiSub import TestfdBiSub


def genBoundStr(domains,compcount):
    """generates lp relaxation bound constraints
    Args:
       domains: set of domains
       compcount: number of components 
    Returns:
       boundstr:
    """
    boundstr = " Bounds\n "
    boundstr += "\n"+ " ".join([" 0 <= x{0}_{1} <= 1\n".format(domin,comp) for domin in xrange(len(domains)) for comp in xrange(compcount)])
    return boundstr

def genConsStr(domains,compcount,cliques,algopar):
    """generates lp relaxation constraints
    Args:
       domains: set of domains
       compcount: number of components
       cliques:
       algopar:
    Returns:
       consstr:
    """
    dom2index = {domains[index]: index for index in xrange(len(domains))}
    consstr = " Subject To \n"
    for clique in cliques:
        consstr += " ".join([" + ".join([" x{0}_{1} ".format(dom2index[dom],comp) for dom in clique]) + " <= 1 \n " for comp in xrange(compcount)])
    if algopar.has_key("k"):
       consstr += "".join([" + ".join([" x{0}_{1} ".format(dom2index[dom],comp) for comp in xrange(compcount)]) + " <= {0}\n".format(algopar["k"]) for dom in domains]) 
    return consstr



def genObjStr(freqmat,compcount,domains,interdom,comp2sc,predomains,prob,algopar,kernmat):
    """generates objective function string
    Args:
       freqmat: frequency matrix
       compcount: number of components
       domains: all domains
       interdom: intersecting domains
       comp2sc: current scale assignment
       predomains:
       prob:
       algopar:
       kernmat: only used in kernel version
    Returns:
       objstr: 
    """
    if algopar["kernel"] != None:
       coefest = lambda dom1,dom2: np.sum(kernmat[max(dom1[0],dom2[0]):min(dom1[1],dom2[1])+1,max(dom1[0],dom2[0]):min(dom1[1],dom2[1])+1])
    elif algopar["obj"] == "frobenius":
       coefest = lambda dom1,dom2: (min(dom1[1],dom2[1])-max(dom1[0],dom2[0]) + 1)**2
    dom2index = {domains[index]:index for index in xrange(len(domains))}
    objstr = " Minimize\n obj: " 
    singles, quads = [], []
    for dom1,dom2 in interdom:
        domin1 = dom2index[dom1]
        domin2 = dom2index[dom2]
        qcoef = coefest(dom1,dom2)
        quads.extend([" %.16f x{0}_{1} * x{2}_{3} ".format(domin1,comp1,domin2,comp2) %(4*qcoef*comp2sc[comp1]*comp2sc[comp2]) for comp1 in xrange(compcount) for comp2 in xrange(compcount)])
    for dom in domains:
        domin = dom2index[dom]
        qcoef = coefest(dom,dom)
        quads.extend([" %.16f x{0}_{1} * x{0}_{2} ".format(domin,comp1,comp2) %(2*qcoef*comp2sc[comp1]*comp2sc[comp2]) for comp1 in xrange(compcount) for comp2 in xrange(compcount)])
        if algopar["kernel"] == None:
           fsum = np.sum(freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1])
        else:
           fsum = np.sum(np.multiply(kernmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1],freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1]))   
        singles.extend([" - %.16f x{0}_{1} ".format(domin,comp) %(2*fsum*comp2sc[comp]) for comp in xrange(compcount)]) 
    if prob == "map":
       singles.extend([" - {0} x{1}_{2} ".format(algopar["lambda"],dom2index[dom],comp) for dom in predomains for comp in xrange(compcount)])
    if len(singles) > 0:
       objstr += " ".join(singles)
    if len(quads) > 0:
       objstr += " + [ " + " + ".join(quads) + " ] "                   
    return objstr


def convertCplexOut(cplexoutpath):
    """reads output solution and returns xdict,ydict,objval 
    Args:
       cplexoutpath: Cormin output file
    Returns:
       xdict: 
       objval:
    """
    retvalues,objval = EmbedUtilities.readCplexOut(cplexoutpath,specific = ["x"])
    xdict = {}
    for key in retvalues.keys():
        if key.startswith("x"):
           xdict[tuple(int(part) for part in key.replace("x","").split("_"))] = retvalues[key]
    return xdict,objval


def getComp2scale(compcount,scales,binsol,index2var):
    """return comp2scale
    Args:
       compcount:
       scales:
       binsol:
       index2var:
    Returns:
       comp2scale:
    """        
    comp2scale = {comp:1 for comp in xrange(compcount)}
    for ind in xrange(len(binsol)):
        if binsol[ind] == 1:
           tdom,comp,scale = index2var[ind]
           comp2scale[comp] += scale
    return comp2scale


def getX(compcount,scales,binsol,index2var):
    """
    """
    xdict = {(0,comp,scale):-1 for comp in xrange(compcount) for scale in scales}
    for ind in xrange(len(binsol)):
        if binsol[ind] == 1:
           tdom,comp,scale = index2var[ind]
           xdict[(0,comp,scale)] = 1
    return xdict


def runScaleOpt(comp2dominds,sideparams,roundnum,PATHS):
    """runs second step of scale optimization
    Args:
       comp2dominds:
       sideparams:
       roundnum:
       PATHS:
       coefmat:
    Returns:
       minval:
       comp2scale:
    """
    compcount,freqmat,interdom,domains,scales,predomains,prob,algopar,kernmat = sideparams
    var2index,index2var = EmbedUtilities.mapVars([1],compcount,scales) #=EmbedUtilities.mapVars2(comp2dominds,compcount,scales)
    A,b = genSdpCoef(freqmat,scales,compcount,domains,interdom,comp2dominds,var2index,kernmat)
    Amat = np.dot(A.transpose(),A)
    side = -1.0*np.dot(b.transpose(),A)
    numsum = np.sum([item**2 for item in b])
    P = np.append(Amat,[side.transpose()],0)
    P = np.column_stack([P, np.append(side,numsum)])
    runfolder = "desdpbisub"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    os.makedirs(runfolder)
    Xmat,objval = Bls_Sdp.runSDPRelax(P,runfolder,PATHS)
    minval,bestsol = sys.maxint, None
    for citer in xrange(roundnum):
        binsol = Bls_Sdp.gaussRound(Xmat,objval,P)
        comp2scale = getComp2scale(compcount,scales,binsol,index2var)
        xdict = getX(compcount,scales,binsol,index2var)
        score = EmbedUtilities.estFracObjectiveMat(Amat,side,numsum,var2index,xdict)
        if TESTMODE:
           score3 = EmbedUtilities.getRatio(freqmat,comp2dominds,comp2scale,domains,kernmat)[1]
           print score,score3,abs(score-score3)/float(score)
           assert abs(score-score3) <= score*1e-06
        #score3 = EmbedUtilities.getRatio(freqmat,comp2dominds,comp2scale,domains,kernmat)[1]
        #assert abs(score-score3) <= score*1e-08
        if score < minval:
           minval = score
           bestsol = deepcopy(comp2scale)

    if TESTMODE:
       newmat = np.dot(P,Xmat)
       assert abs(sum([newmat[ind1,ind1] for ind1 in xrange(np.shape(newmat)[0])]) - objval) < objval*1e-07
       for eig in (np.linalg.eigh(Xmat))[0]:
           assert eig >= 0.0
       assert np.allclose(P.transpose(), P)    
    return minval,bestsol

    
def genSdpCoef(freqmat,scales,compcount,domains,interdom,comp2dominds,var2index,kernmat):
    """generates sdp matrix coefs
    Args:
       freqmat: frequency matrix
       scales: set of scales
       compcount: number of components
       domains: all domains
       interdom: intersecting domains
       comp2dominds: first step solution
       var2index:
       kernmat:
    Returns:
       coefmat:
       b:
       fsum:
    """
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    dom2index = {domains[index]:index for index in xrange(len(domains))}
    vcount = np.shape(freqmat)[0]
    scalesum = sum(scales)
    A = np.zeros((vcount**2,compcount*len(scales)),dtype=np.float)
    b = np.array([0.0] * vcount**2)
    for ind1 in xrange(np.shape(freqmat)[0]):
        for ind2 in xrange(np.shape(freqmat)[1]):
            b[(ind1*vcount)+ind2] = freqmat[ind1,ind2]
            for comp in comp2dominds.keys():
                doms = set(dom for dom in comp2dominds[comp] if dom in node2dom[ind1] and dom in node2dom[ind2])
                assert len(doms) <= 1
                if len(doms) == 1:
                   for scale in scales:
                       A[(ind1*vcount)+ind2,var2index[(0,comp,scale)]] = 0.5*scale
                   b[(ind1*vcount)+ind2] -= (scalesum*0.5)+1
    if kernmat != None:
       for row in xrange(np.shape(A)[0]):
           ind1,ind2 = row/vcount, row%vcount 
           for col in xrange(np.shape(A)[1]):
               A[row,col] *= math.sqrt(kernmat[ind1,ind2])
           b[row] *= math.sqrt(kernmat[ind1,ind2])
    return A,b


def getScales(unitscale,freqmat,compcount):
    """get scales
    Args:
       unitscale:
       freqmat:
       compcount:
    Returns:
       scales:
    """
    if unitscale:
       return [1]
    return [2**index for index in xrange(0,int(math.ceil(math.log(np.amax(freqmat)+1,2))))]

 
ROUNDCOUNT = 25 
TESTMODE = False
def runBisub(freqmat,domains,compcount,outprefix,scalemode,prob,algopar,PATHS):
    """runs bisubmodular optimization with rounding
    Args:
       freqmat: frequency matrix
       domains: list of domains(list of set of nodes)
       compcount: number of components
       outprefix: folder+prefix
       scalemode,prob,algopar:
       PATHS:
    Returns:
       comp2doms:
       comp2scale:
    """
    unitscale = False
    if scalemode == "binary":
       unitscale = True
    stime = time.time()
    predomains = set(domains)
    if prob == "map":
       alldomains = EmbedUtilities.getPosDomains(freqmat,domains)
       sentdoms = EmbedUtilities.preprocDomains(alldomains,domains)
       domains = list(sentdoms)     
    interdom = EmbedUtilities.getInterDomain(domains)
    cliques = EmbedUtilities.findDomainCliqueDecomp(domains,interdom)
    scales = getScales(unitscale,freqmat,compcount)
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    if algopar["kernel"] != None:
       coefmat = EmbedUtilities.genKernelCoefMat(freqmat,algopar) 
       fqsum = sum([(freqmat[ind1,ind2]**2)*coefmat[ind1,ind2] for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    else:
       coefmat = None      
       fqsum = float(sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])]))
    comp2scale,comp2dominds = {comp:1 for comp in xrange(compcount)}, {}
    sideparams = [compcount,freqmat,interdom,domains,scales,predomains,prob,algopar,coefmat]
    objlist = [fqsum]
    curobjval = fqsum
    while True:
       #first step
       objstr = genObjStr(freqmat,compcount,domains,interdom,comp2scale,predomains,prob,algopar,coefmat)
       consstr =  genConsStr(domains,compcount,cliques,algopar)
       boundstr = genBoundStr(domains,compcount)
       outmethod = globals()["convertCplexOut"]
       runfolder = "delpbisub"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
       xdict, objval = EmbedUtilities.runCplexCode(consstr,objstr,boundstr,"",runfolder,outmethod)
       objval += fqsum
       if unitscale:
          curcomp2dominds,step1objval = Round.roundCR(xdict,sideparams,comp2scale,algopar,ROUNDCOUNT)
          curobjval, comp2dominds = step1objval,dict(curcomp2dominds)
          objlist.append(curobjval)   
       else:
          curcomp2dominds,step1objval = Round.roundCR(xdict,sideparams,comp2scale,algopar,ROUNDCOUNT)
       if TESTMODE:
          tempobjval = EmbedUtilities.getRatio(freqmat,curcomp2dominds,comp2scale,domains,coefmat)[1]
          if prob == "map":
             priorobj = EmbedUtilities.addPriorObj(curcomp2dominds,domains,predomains,algopar["lambda"]) 
             tempobjval += priorobj
          assert abs(tempobjval-step1objval) < step1objval*1e-08
          assert TestfdBiSub.testLPData(sideparams,objstr)
          dommap = lambda x: comp2dominds if x else curcomp2dominds
          assert TestfdBiSub.testLPOutput(xdict,sideparams,cliques,objval,comp2scale,dommap(unitscale),algopar)
       if not unitscale and curobjval - step1objval >= 2:
          print "step1 ",step1objval 
          curobjval,comp2dominds = step1objval, deepcopy(curcomp2dominds)
          objlist.append(step1objval)
       else:
          break
       #second step
       step2objval,curcomp2scale = runScaleOpt(comp2dominds,sideparams,1000,PATHS)
       if prob == "map":
          priorobj = EmbedUtilities.addPriorObj(comp2dominds,domains,predomains,algopar["lambda"])
          step2objval += priorobj
       if TESTMODE:
          tempobjval = EmbedUtilities.getRatio(freqmat,comp2dominds,curcomp2scale,domains,coefmat)[1]
          if prob == "map":
             tempobjval += priorobj
          assert abs(step2objval - tempobjval) < step2objval*1e-09
          assert TestfdBiSub.testSDPOutput(sideparams,curcomp2scale,comp2dominds)   
       if curobjval - step2objval >= 2:
          print "step2 ",step2objval 
          curobjval,comp2scale = step2objval, deepcopy(curcomp2scale)
          objlist.append(step2objval)
       else:
          break
    etime = time.time()
    print "iterated objectives: ",objlist
    metadata = {"objval": curobjval,"ratio": float(curobjval)/fqsum, "fqsum":fqsum, "time":etime-stime, "objlist":objlist}
    
    if TESTMODE:
       newobjval = EmbedUtilities.getRatio(freqmat,comp2dominds,comp2scale,domains,coefmat)[1]
       if prob == "map":
          newobjval += EmbedUtilities.addPriorObj(comp2dominds,domains,predomains,algopar["lambda"]) 
       assert abs(curobjval - newobjval) < 0.1
       assert TestfdBiSub.testObjMethods(sideparams,comp2scale,fqsum)
       sideparams = [scales,compcount,freqmat,interdom,domains]
       assert TestfdBiSub.testData(cliques,sideparams)
       assert TestfdBiSub.testImprovement(sideparams,comp2dominds,comp2scale,coefmat)
    
    EmbedUtilities.makeDeconOutput(outprefix,domains,comp2dominds,comp2scale,metadata)   

    
def makeParser():
    """
    """
    parser = argparse.ArgumentParser(description='Process 3CDE parameters')
    parser.add_argument('-f', dest='freqfilename', type=str, action='store', default='freq.gz', help='Ensemble Interaction Matrix(default: freq.gz)')
    parser.add_argument('-d', dest='domainfile', type=str, action='store', default=None, help='Prior Possible Domains(default: None)')
    parser.add_argument('-c', dest='compcount', type=int, action='store', default=5, help='number of classes(default: 5)')
    parser.add_argument('-o', dest='outprefix', type=str, action='store', default='deconout', help='output prefix(default: deconout)')
    parser.add_argument('-s', dest='scalemode', type=str, action='store', default='uniform', help='scale type uniform or integer(default: uniform)')
    parser.add_argument('-p', dest='prob', type=str, action='store', default='map', help='Estimation type, ML or MAP prob(default: map)')
    parser.add_argument('-l', dest='lambdaval', type=float, action='store', default=1.0, help='lambda for prior domain effects(default: 1.0)')
    parser.add_argument('-kern', dest='kernel', type=str, action='store', default=None, help='Kernel type, exp or gauss. Leave blank for None(Default: None)')
    parser.add_argument('-coef', dest='kerncoef', type=float, action='store', default=None, help='Kernel coefficient such as 0.2. Leave blank for None(Default: None)')
    parser.add_argument('-sdpt', dest='sdptdir', type=str, action='store', default='SDPT3-4.0', help='SDPT Path(Default: SDPT3-4.0)')
    parser.add_argument('-mat', dest='matlabdir', type=str, action='store', default='matlab', help='MATLAB Path(Default: matlab)')
    return parser

def checkParams(prob,scalemode,algopar):
    """
    """
    if prob not in ["map","ml"] or algopar["obj"] not in ["frobenius"] or scalemode not in ["uniform","binary","integer"] or algopar["kernel"] not in ["exp","gauss",None]:
       print "Parameter mismatch, Please run with -h for help!!"
       exit(1)
    if algopar["kernel"] in ["exp","gauss"] and algopar["kerncoef"] in [None,'None']:
       print "Parameter mismatch, Please run with -h for help!!"
       exit(1)     
       
if  __name__ =='__main__':
    """runs
    """
    import argparse
    parser = makeParser()
    args = parser.parse_args(sys.argv[1:])
    algopar = {'kernel': args.kernel, 'obj': 'frobenius', 'kerncoef': args.kerncoef,'lambda': args.lambdaval}
    PATHS = {'SDPTDIR': args.sdptdir, 'MATLABPATH': args.matlabdir}
    globals().update(vars(args))
    checkParams(prob,scalemode,algopar)
    scalemode = {"uniform":"binary","binary":"binary","integer":"integer"}[scalemode]
    freqmat,nodenames,in2pos = EmbedUtilities.readFreqFile(freqfilename)
    if domainfile == None:
       domains = EmbedUtilities.selectSubsetDomain(freqmat)
    else:   
       tdomains = EmbedUtilities.readDomainFile(domainfile)
       domains = [(start-1,end-1) for start, end in tdomains]
    runBisub(freqmat,domains,compcount,outprefix,scalemode,prob,algopar,PATHS)   
