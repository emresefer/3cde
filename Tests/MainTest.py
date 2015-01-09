#Main Test Class
import EmbedUtilities
import networkx as nx
import numpy as np
from copy import deepcopy

class MainTest():

    @staticmethod
    def testIntersect2(comp2domains):
        """tests overlap
        """
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        for comp in comp2domains.keys():
            doms = list(comp2domains[comp])
            for ind1 in xrange(len(doms)):
                dom1 = doms[ind1]
                for ind2 in xrange(ind1+1,len(doms)):
                    dom2 = doms[ind2]
                    assert not intersect(dom1,dom2)
        return True
    
    @staticmethod
    def testIntersect(comp2dominds,domains):
        """tests overlap
        """
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        for comp in comp2dominds.keys():
            alldoms = list(comp2dominds[comp])
            for ind1 in xrange(len(alldoms)):
                domin1 = alldoms[ind1]
                for ind2 in xrange(ind1+1,len(alldoms)):
                    domin2 = alldoms[ind2]
                    assert not intersect(domains[domin1],domains[domin2])
        return True
    
    @staticmethod
    def testImprovement(sideparams,comp2dominds,comp2scale,kernmat=None):
        """tests improvement in the solution
        Args:
           sideparams:
           comp2dominds:
           comp2scale:
           kernmat:
        Returns:
           bool:
        """
        [scales,compcount,freqmat,interdom,domains] = sideparams
        firstmat = np.zeros(np.shape(freqmat),dtype=np.float)
        secondmat = np.zeros(np.shape(freqmat),dtype=np.float)
        for comp in comp2dominds.keys():
            for domin in comp2dominds[comp]:
                start,end = domains[domin]
                firstmat[start:end+1,start:end+1] += 1
                secondmat[start:end+1,start:end+1] += comp2scale[comp]
        difmat1 = freqmat - firstmat        
        difmat2 = freqmat - secondmat
        if kernmat == None:
           firstdif = sum([difmat1[ind1,ind2]**2 for ind1 in xrange(np.shape(firstmat)[0]) for ind2 in xrange(np.shape(firstmat)[1])])
           freqsum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
           seconddif = sum([difmat2[ind1,ind2]**2 for ind1 in xrange(np.shape(secondmat)[0]) for ind2 in xrange(np.shape(secondmat)[1])])
        else:
           firstdif = sum([(difmat1[ind1,ind2]**2)*kernmat[ind1,ind2] for ind1 in xrange(np.shape(firstmat)[0]) for ind2 in xrange(np.shape(firstmat)[1])])
           freqsum = sum([(freqmat[ind1,ind2]**2)*kernmat[ind1,ind2] for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
           seconddif = sum([(difmat2[ind1,ind2]**2)*kernmat[ind1,ind2] for ind1 in xrange(np.shape(secondmat)[0]) for ind2 in xrange(np.shape(secondmat)[1])])    
        assert firstdif <= freqsum and seconddif <= firstdif
        return True
    
    @staticmethod  
    def testData(cliques,sideparams):
        """tests data
        Args:
          cliques:
          sideparams:
        Returns:
          bool: true or false
        """
        [scales,compcount,freqmat,interdom,domains] = sideparams
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        storeG = EmbedUtilities.inter2Graph(interdom,domains)
        tinterG = deepcopy(storeG)
        for clique in cliques:
            for ind1 in xrange(len(clique)):
                node1 = clique[ind1]
                for ind2 in xrange(ind1+1,len(clique)):
                    node2 = clique[ind2]
                    if tinterG.has_edge(node1,node2):
                       tinterG.remove_edge(node1,node2)
        assert tinterG.number_of_edges() == 0
        tinterG = deepcopy(storeG)
        for clique in cliques:
            for ind1 in xrange(len(clique)):
                node1 = clique[ind1]
                for ind2 in xrange(ind1+1,len(clique)):
                    node2 = clique[ind2]
                    assert tinterG.has_edge(node1,node2)
        assert len(set([dom for clique in cliques for dom in clique]) ^ set(domains)) == 0
        for dom1,dom2 in interdom:
            assert intersect(dom1,dom2)
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        for node in node2dom.keys():
            for domin in node2dom[node]:
                assert domains[domin][0] <= node and node <= domains[domin][1]    
        return True
    
