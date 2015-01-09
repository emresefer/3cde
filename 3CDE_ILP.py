#Combinatorial ILP Formulation for Optimal Solution
import numpy as np
import os
import sys
import time
import random
import string
sys.path.append("lib")
import EmbedUtilities
sys.path.append("Tests")
from TestfdILP import TestfdILP


def genBoundStr(domains,compcount,scales):
    """generates ilp bound constraints
    Args:
       domains: set of domains
       compcount: number of components
       scales: set of scales
    Returns:
       boundstr:
    """
    boundstr = " Bounds\n "
    boundstr += "\n"+ " ".join([" 0 <= z{0}_{1}_{2} <= 1\n".format(domin,comp,scale) for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales])
    boundstr += " ".join([" 0 <= y{0}_{1} <= 1\n".format(comp,scale) for comp in xrange(compcount) for scale in scales])
    return boundstr


def genConsStr(domains,interdom,scales,compcount):
    """generates ilp constraints
    Args:
       domains: set of domains
       interdom: interacting domain pairs
       scales: set of scales
       compcount: number of components
    Returns:
       objstr:
    """ 
    dom2index = {domains[index]: index for index in xrange(len(domains))}
    consstr = " Subject To \n"
    for dom1,dom2 in interdom:
        consstr += " ".join(["z{0}_{1}_{2} + z{3}_{1}_{2} - y{1}_{2} <= 0\n ".format(dom2index[dom1],comp,scale,dom2index[dom2]) for comp in xrange(compcount) for scale in scales])
    remdoms = set(domains) - set(item for part in interdom for item in part)
    for remdom in remdoms:
        consstr += " ".join(["z{0}_{1}_{2} - y{1}_{2} <= 0\n ".format(dom2index[remdom],comp,scale) for comp in xrange(compcount) for scale in scales])    
    consstr += " ".join([" + ".join([ "y{0}_{1}".format(comp,scale) for scale in scales]) + " <= 1 \n"  for comp in xrange(compcount)])
    return consstr


def genVarStr(domains,compcount,scales):
    """generates ilp bound constraints
    Args:
       domains: set of domains
       compcount: number of components
       scales: set of scales
    Returns:
       varstr:
    """
    varstr = ""
    varstr += "\n Binary \n"
    varstr +=  " \n ".join([" z{0}_{1}_{2} ".format(domin,comp,scale) for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales])
    varstr += "\n" + "\n".join([" y{0}_{1} ".format(comp,scale) for comp in xrange(compcount) for scale in scales])
    return varstr


def genObjStr(freqmat,scales,compcount,domains,interdom,predomains,prob,algopar,kernmat):
    """generates objective function string
    Args:
       freqmat: frequency matrix
       scales: set of scales
       compcount: number of components
       domains: all domains
       interdom: intersecting domains
       predomains:
       prob,algopar:
       kernmat:
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
    pairs = [(comp,scale) for comp in xrange(compcount) for scale in scales]
    for dom1,dom2 in interdom:
        domin1 = dom2index[dom1]
        domin2 = dom2index[dom2]
        qcoef = coefest(dom1,dom2)
        #quads.extend([" {0} z{1}_{2}_{3} * z{4}_{5}_{6} ".format(4*qcoef*scale1*scale2,domin1,comp1,scale1,domin2,comp2,scale2) for comp1,scale1 in pairs for comp2,scale2 in pairs])
        quads.extend([" %.16f z{0}_{1}_{2} * z{3}_{4}_{5} ".format(domin1,comp1,scale1,domin2,comp2,scale2) %(4*qcoef*scale1*scale2) for comp1,scale1 in pairs for comp2,scale2 in pairs])
    for dom in domains:
        domin = dom2index[dom]
        qcoef = coefest(dom,dom)
        #quads.extend([" {0} z{1}_{2}_{3} * z{1}_{4}_{5} ".format(2*qcoef*scale1*scale2,domin,comp1,scale1,comp2,scale2) for comp1,scale1 in pairs for comp2,scale2 in pairs])
        quads.extend([" %.16f z{0}_{1}_{2} * z{0}_{3}_{4} ".format(domin,comp1,scale1,comp2,scale2) %(2*qcoef*scale1*scale2) for comp1,scale1 in pairs for comp2,scale2 in pairs])
        if algopar["kernel"] == None:
           fsum = np.sum(freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1])
        else:
           fsum = np.sum(np.multiply(kernmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1],freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1]))  
        #singles.extend([" - {0} z{1}_{2}_{3} ".format(2*scale*fsum,domin,comp,scale) for comp,scale in pairs])
        singles.extend([" - %.16f z{0}_{1}_{2} ".format(domin,comp,scale) %(2*scale*fsum) for comp,scale in pairs]) 
    if prob == "map":
       singles.extend([" - {0} z{1}_{2}_{3} ".format(algopar["lambda"],dom2index[dom],comp,scale) for dom in predomains for comp,scale in pairs])    
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
    retvalues,objval = EmbedUtilities.readCplexOut(cplexoutpath,specific = ["z","y"])
    zdict,ydict = {}, {}
    for key in retvalues.keys():
        if key.startswith("z"):
           zdict[tuple(int(part) for part in key.replace("z","").split("_"))] = retvalues[key]
        elif key.startswith("y"):
           ydict[tuple(int(part) for part in key.replace("y","").split("_"))] = retvalues[key]   
    return zdict,ydict,objval


def addMissing(cliques,domains):
    """adds cliques for unseen nodes if exists
    Args:
       cliques:
       domains:
    Returns:
       cliques:
    """    
    clidomains = set(dom for cli in cliques for dom in cli)
    difset = set(domains) - clidomains
    if len(difset) == 0:
       return cliques
    for difdom in difset:
        cliques.append([difdom])      
    return cliques 

def getILPScales(unitscale,freqmat,compcount):
    """
    """
    if unitscale:
       return [1]
    return range(1,max(2,int(np.amax(freqmat))+1-compcount+1))

def toOutFormat(zdict,ydict,compcount):
    """
    Args:
       zdict:
       ydict:
       compcount:
    Returns:
       comp2dominds:
       comp2scale:
    """    
    comp2dominds = {comp:set() for comp in xrange(compcount)}
    for domin,comp,scale in zdict.keys():
        comp2dominds[comp].add(domin)
    comp2scale = {comp:0 for comp in xrange(compcount)}
    for comp,scale in ydict.keys():
        comp2scale[comp] = scale
    return comp2dominds, comp2scale


TESTMODE = False
def runILP(freqmat,domains,compcount,outprefix,scalemode,prob,algopar):
    """runs ilp formulation
    Args:
       freqmat: frequency matrix
       domains: list of domains(list of set of nodes)
       compcount: number of components
       outprefix: folder+prefix
       scalemode:
       prob:
       algopar:
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
    scales = getILPScales(unitscale,freqmat,compcount)
    if algopar["kernel"] != None:
       coefmat = EmbedUtilities.genKernelCoefMat(freqmat,algopar) 
       fqsum = sum([(freqmat[ind1,ind2]**2)*coefmat[ind1,ind2] for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    else:
       coefmat = None      
       fqsum = float(sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])]))  
    objstr = genObjStr(freqmat,scales,compcount,domains,interdom,predomains,prob,algopar,coefmat)
    consstr = genConsStr(domains,interdom,scales,compcount)
    boundstr = genBoundStr(domains,compcount,scales)
    varstr =  genVarStr(domains,compcount,scales)
    outmethod = globals()["convertCplexOut"]
    runfolder = "deconilp"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    zdict,ydict, objval = EmbedUtilities.runCplexCode(consstr,objstr,boundstr,varstr,runfolder,outmethod)
    etime = time.time()
    comp2dominds, comp2scale = toOutFormat(zdict,ydict,compcount)
    objval += fqsum 
    metadata = {"objval": objval,"ratio": float(objval)/fqsum, "fqsum":fqsum, "time":etime-stime}

    if TESTMODE:
       newobjval = EmbedUtilities.getRatio(freqmat,comp2dominds,comp2scale,domains,coefmat)[1]
       print objval
       print newobjval
       if prob == "map":
          priorobj = EmbedUtilities.addPriorObj(comp2dominds,domains,predomains,algopar["lambda"])  
          newobjval += priorobj
          print priorobj
       print newobjval     
       assert abs(objval - newobjval) < 0.1
       sideparams = [scales,compcount,freqmat,interdom,domains,predomains,prob,algopar,coefmat]
       assert TestfdILP.testILPOutput(sideparams,comp2dominds,comp2scale,zdict,ydict,objval)
       sideparams = [scales,compcount,freqmat,interdom,domains]
       cliques = EmbedUtilities.findDomainCliqueDecomp(domains,interdom)
       assert TestfdILP.testData(cliques,sideparams)
       assert TestfdILP.testImprovement(sideparams,comp2dominds,comp2scale)
       sideparams = [scales,compcount,freqmat,interdom,domains,coefmat]
       assert TestfdILP.testILPData(sideparams,objstr)

    EmbedUtilities.makeDeconOutput(outprefix,domains,comp2dominds,comp2scale,metadata)   

    
   
def makeParser():
    """
    """
    parser = argparse.ArgumentParser(description='Process 3CDE parameters')
    parser.add_argument('-f', dest='freqfilename', type=str, action='store', default='freq.gz', help='Ensemble Interaction Matrix(default: freq.gz)')
    parser.add_argument('-d', dest='domainfile', type=str, action='store', default=None, help='Prior Possible Domains(default: None)')
    parser.add_argument('-c', dest='compcount', type=int, action='store', default=5, help='number of classes(default: 5)')
    parser.add_argument('-o', dest='outprefix', type=str, action='store', default='deconilpout', help='output prefix(default: deconilpout)')
    parser.add_argument('-s', dest='scalemode', type=str, action='store', default='binary', help='scale type, binary or integer(default: binary)')
    parser.add_argument('-p', dest='prob', type=str, action='store', default='map', help='Estimation type, ML or MAP prob(default: map)')
    parser.add_argument('-l', dest='lambdaval', type=float, action='store', default=1.0, help='lambda for prior domain effects(default: 1.0)')
    parser.add_argument('-kern', dest='kernel', type=str, action='store', default=None, help='Kernel type, exp or gauss. Leave blank for None(Default: None)')
    parser.add_argument('-coef', dest='kerncoef', type=float, action='store', default=None, help='Kernel coefficient such as 0.2. Leave blank for None(Default: None)')
    return parser

def checkParams(prob,scalemode,algopar):
    """
    """
    if prob not in ["map","ml"] or algopar["obj"] not in ["frobenius"] or scalemode not in ["binary","integer","uniform"] or algopar["kernel"] not in ["exp","gauss",None]:
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
    globals().update(vars(args))
    checkParams(prob,scalemode,algopar)
    scalemode = {"uniform":"binary","binary":"binary","integer":"integer"}[scalemode]
    freqmat,nodenames,in2pos = EmbedUtilities.readFreqFile(freqfilename)
    if domainfile == None:
       domains = EmbedUtilities.selectSubsetDomain(freqmat)
    else:   
       tdomains = EmbedUtilities.readDomainFile(domainfile)
       domains = [(start-1,end-1) for start, end in tdomains]
    runILP(freqmat,domains,compcount,outprefix,scalemode,prob,algopar)         
