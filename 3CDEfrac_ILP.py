#Combinatorial ILP Formulation for Optimal Solution
import numpy as np
import os
import sys
import time
import random
import string
import itertools
sys.path.append("lib")
import EmbedUtilities
sys.path.append("Tests")
from TestfdILP import TestfdILP


def genBoundStr(domains,compcount,maxscale):
    """generates ilp bound constraints
    Args:
       domains: set of domains
       compcount: number of components
       maxscale: set of scales
    Returns:
       boundstr:
    """
    boundstr = " Bounds\n "
    boundstr += "\n"+ " ".join([" 0 <= y{0}_{1} <= {2}\n".format(domin,comp,maxscale) for domin in xrange(len(domains)) for comp in xrange(compcount)])
    boundstr += " ".join([" 0 <= x{0}_{1} <= 1\n".format(domin,comp) for domin in xrange(len(domains)) for comp in xrange(compcount)])
    return boundstr


def genConsStr(domains,interdom,maxscale,compcount):
    """generates ilp constraints
    Args:
       domains: set of domains
       interdom: interacting domain pairs
       maxscale: maximum scale
       compcount: number of components
    Returns:
       objstr:
    """ 
    dom2index = {domains[index]: index for index in xrange(len(domains))}
    consstr = " Subject To \n"
    for dom1,dom2 in interdom:
        consstr += " ".join(["x{0}_{1} + x{2}_{1} <= 1\n ".format(dom2index[dom1],comp,dom2index[dom2]) for comp in xrange(compcount)])
    consstr += " ".join(["y{0}_{1} - {2} x{0}_{1} <= 0\n ".format(dom2index[dom],comp,maxscale) for dom in domains for comp in xrange(compcount)])
    for dom1,dom2 in itertools.combinations(domains,2):   
        if (dom1,dom2) not in interdom and (dom2,dom1) not in interdom:
           domin1,domin2 = dom2index[dom1], dom2index[dom2] 
           consstr +=" ".join(["y{0}_{2} - y{1}_{2} + {4} x{0}_{2} + {4} x{1}_{2} <= {3}\n ".format(domin1,domin2,comp,2*maxscale,maxscale) for comp in xrange(compcount)])  
           consstr +=" ".join(["y{0}_{2} - y{1}_{2} - {4} x{0}_{2} - {4} x{1}_{2} >= {3}\n ".format(domin1,domin2,comp,-2*maxscale,maxscale) for comp in xrange(compcount)])
           consstr +=" ".join(["y{1}_{2} - y{0}_{2} + {4} x{0}_{2} + {4} x{1}_{2} <= {3}\n ".format(domin1,domin2,comp,2*maxscale,maxscale) for comp in xrange(compcount)])  
           consstr +=" ".join(["y{1}_{2} - y{0}_{2} - {4} x{0}_{2} - {4} x{1}_{2} >= {3}\n ".format(domin1,domin2,comp,-2*maxscale,maxscale) for comp in xrange(compcount)])
    return consstr


def genVarStr(domains,compcount):
    """generates ilp bound constraints
    Args:
       domains: set of domains
       compcount: number of components
    Returns:
       varstr:
    """
    varstr = ""
    varstr += "\n Binary \n"
    varstr +=  " \n ".join([" x{0}_{1} ".format(domin,comp) for domin in xrange(len(domains)) for comp in xrange(compcount)])
    varstr += "\n"
    return varstr


def genObjStr(freqmat,compcount,domains,interdom,predomains,prob,algopar,kernmat):
    """generates objective function string
    Args:
       freqmat: frequency matrix
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
    for dom1,dom2 in interdom:
        domin1 = dom2index[dom1]
        domin2 = dom2index[dom2]
        qcoef = coefest(dom1,dom2)
        quads.extend([" %.16f y{0}_{1} * y{2}_{3} ".format(domin1,comp1,domin2,comp2) %(4*qcoef) for comp1 in xrange(compcount) for comp2 in xrange(compcount)])     
    for dom in domains:
        domin = dom2index[dom]
        qcoef = coefest(dom,dom)
        quads.extend([" %.16f y{0}_{1} * y{0}_{2} ".format(domin,comp1,comp2) %(2*qcoef) for comp1 in xrange(compcount) for comp2 in xrange(compcount)])
        if algopar["kernel"] == None:
           fsum = np.sum(freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1])
        else:
           fsum = np.sum(np.multiply(kernmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1],freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1]))
        singles.extend([" - %.16f y{0}_{1} ".format(domin,comp) %(2*fsum) for comp in xrange(compcount)]) 
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
    retvalues,objval = EmbedUtilities.readCplexOut(cplexoutpath,specific = ["x","y"])
    xdict,ydict = {}, {}
    for key in retvalues.keys():
        if key.startswith("x"):
           xdict[tuple(int(part) for part in key.replace("x","").split("_"))] = retvalues[key]
        elif key.startswith("y"):
           ydict[tuple(int(part) for part in key.replace("y","").split("_"))] = retvalues[key]   
    return xdict,ydict,objval


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


def toOutFormat(xdict,ydict,compcount):
    """
    Args:
       xdict:
       ydict:
       compcount:
    Returns:
       comp2dominds:
       comp2scale:
    """    
    comp2dominds = {comp:set() for comp in xrange(compcount)}
    for domin,comp in xdict.keys():
        if xdict[(domin,comp)] not in [0.0,-0.0]:
           comp2dominds[comp].add(domin)
    comp2scale = {comp:set() for comp in xrange(compcount)}
    for domin,comp in ydict.keys():
        if ydict[(domin,comp)] not in [0.0,-0.0]:
           comp2scale[comp].add(ydict[(domin,comp)])
    print comp2scale    
    for comp in comp2scale.keys():
        #assert len(comp2scale[comp]) <=1
        if len(comp2scale[comp]) == 0:
           comp2scale[comp] = 0.0
        elif len(comp2scale[comp]) != 1:
           comp2scale[comp] = sorted(list(comp2scale[comp]))[-1]   
        else:
           comp2scale[comp] = list(comp2scale[comp])[0]
    return comp2dominds, comp2scale


TESTMODE = False
def runILPfrac(freqmat,domains,compcount,outprefix,prob,algopar):
    """runs ilp formulation
    Args:
       freqmat: frequency matrix
       domains: list of domains(list of set of nodes)
       compcount: number of components
       outprefix: folder+prefix
       prob:
       algopar:
    Returns:
       comp2doms:
       comp2scale:
    """
    stime = time.time()
    predomains = set(domains)
    if prob == "map":
       alldomains = EmbedUtilities.getPosDomains(freqmat,domains)
       sentdoms = EmbedUtilities.preprocDomains(alldomains,domains)
       domains = list(sentdoms)
    if algopar["kernel"] != None:
       coefmat = EmbedUtilities.genKernelCoefMat(freqmat,algopar) 
       fqsum = sum([(freqmat[ind1,ind2]**2)*coefmat[ind1,ind2] for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    else:
       coefmat = None      
       fqsum = float(sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])]))   
    interdom = EmbedUtilities.getInterDomain(domains)
    maxscale = max(2,np.amax(freqmat)+1-compcount+1)
    objstr = genObjStr(freqmat,compcount,domains,interdom,predomains,prob,algopar,coefmat)
    consstr = genConsStr(domains,interdom,maxscale,compcount)
    boundstr = genBoundStr(domains,compcount,maxscale)
    varstr =  genVarStr(domains,compcount)
    outmethod = globals()["convertCplexOut"]
    runfolder = "deconilpmixed"+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    xdict,ydict,objval = EmbedUtilities.runCplexCode(consstr,objstr,boundstr,varstr,runfolder,outmethod)
    comp2dominds, comp2scale = toOutFormat(xdict,ydict,compcount)
    objval += fqsum       
    etime = time.time()
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
       sideparams = [compcount,freqmat,interdom,domains,predomains,prob,algopar,coefmat]
       assert TestfdILP.testfracILPOutput(sideparams,comp2dominds,comp2scale,xdict,ydict,objval)
       sideparams = [None,compcount,freqmat,interdom,domains]
       cliques = EmbedUtilities.findDomainCliqueDecomp(domains,interdom)
       assert TestfdILP.testData(cliques,sideparams)
       assert TestfdILP.testImprovement(sideparams,comp2dominds,comp2scale)
       sideparams = [compcount,freqmat,interdom,domains,coefmat]
       assert TestfdILP.testfracILPData(sideparams,objstr)
       
    EmbedUtilities.makeDeconOutput(outprefix,domains,comp2dominds,comp2scale,metadata)


def makeParser():
    """
    """
    parser = argparse.ArgumentParser(description='Process 3CDE parameters')
    parser.add_argument('-f', dest='freqfilename', type=str, action='store', default='freq.gz', help='Ensemble Interaction Matrix(default: freq.gz)')
    parser.add_argument('-d', dest='domainfile', type=str, action='store', default=None, help='Prior Possible Domains(default: None)')
    parser.add_argument('-c', dest='compcount', type=int, action='store', default=5, help='number of classes(default: 5)')
    parser.add_argument('-o', dest='outprefix', type=str, action='store', default='deconilpout', help='output prefix(default: deconilpout)')
    parser.add_argument('-p', dest='prob', type=str, action='store', default='map', help='Estimation type, ML or MAP prob(default: map)')
    parser.add_argument('-l', dest='lambdaval', type=float, action='store', default=1.0, help='lambda for prior domain effects(default: 1.0)')
    parser.add_argument('-kern', dest='kernel', type=str, action='store', default=None, help='Kernel type, exp or gauss. Leave blank for None(Default: None)')
    parser.add_argument('-coef', dest='kerncoef', type=float, action='store', default=None, help='Kernel coefficient such as 0.2. Leave blank for None(Default: None)')
    return parser

def checkParams(prob,algopar):
    """
    """
    if prob not in ["map","ml"] or algopar["obj"] not in ["frobenius"] or algopar["kernel"] not in ["exp","gauss",None]:
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
    checkParams(prob,algopar)
    freqmat,nodenames,in2pos = EmbedUtilities.readFreqFile(freqfilename)
    if domainfile == None:
       domains = EmbedUtilities.selectSubsetDomain(freqmat)
    else:   
       tdomains = EmbedUtilities.readDomainFile(domainfile)
       domains = [(start-1,end-1) for start, end in tdomains]
    runILPfrac(freqmat,domains,compcount,outprefix,prob,algopar)
