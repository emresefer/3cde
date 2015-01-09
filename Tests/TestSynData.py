#Synthetic Data Generator Test
import numpy as np
import sys
sys.path.append("lib")
import EmbedUtilities
import itertools
import myutilities as myutil


class TestSynData:

    @staticmethod
    def testSynData(freqmat,comp2domain,comp2scale):
        """tests the synthetic generated data
        Args:
          freqmat:
          comp2domain:
          comp2scale:
        Returns:
          bool: true fale
        """
        maxfreq = np.amax(freqmat)
        maxindices = set((in1,in2) for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1]) if freqmat[in1,in2] == maxfreq)
        #freq matrix sum check
        newfreqmat = np.zeros(np.shape(freqmat),dtype=np.float)
        for comp in comp2domain.keys():
            doms = comp2domain[comp]
            for s,e in doms:
                newfreqmat[s-1:e,s-1:e] += comp2scale[comp]
        assert np.allclose(freqmat,newfreqmat)
        #each motif can only contribute as its scale
        maxtotcount = 0
        for in1,in2 in maxindices:
            totcount = 0
            for comp in comp2domain.keys():
                count = 0
                for dom in comp2domain[comp]:
                    if in1 >= dom[0] and in2 >= dom[0] and in1 <= dom[1] and in2 <= dom[1]:
                       count += 1
                       totcount += comp2scale[comp]
                assert count <= 1
            if totcount > maxtotcount:
               maxtotcount = totcount
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        #generated domains must not intersect 
        for comp in comp2domain.keys():
            curdoms = list(comp2domain[comp])
            for dom1,dom2 in list(itertools.combinations(curdoms,2)):
                assert not intersect(dom1,dom2)
        assert np.amax(freqmat) - sum(comp2scale.values()) < 0.001
        return True
     
    @staticmethod
    def testReadWrite(freqfile,freqmat,nodecount,truefile,comp2domain,comp2scale):
        """
        Args:
          freqfile:
          freqmat:
          nodecount:
          truefile:
          comp2domain:
          rcomp2scale:
        Returns:
          bool:
        """
        freqmat2, allnodes,in2pos = EmbedUtilities.readFreqFile(freqfile)
        assert len(allnodes) == nodecount
        for ind1 in xrange(np.shape(freqmat)[0]):
            for ind2 in xrange(np.shape(freqmat)[1]):
                assert abs(freqmat[ind1,ind2] - freqmat2[ind1,ind2]) <= 0.01
        if truefile != None: 
           rcomp2domain,rcomp2scale = EmbedUtilities.readDeconOut(truefile)
           for comp in rcomp2scale.keys():
               assert abs(rcomp2scale[comp] - comp2scale[comp]) < 0.0001
           for comp in rcomp2domain.keys():
              assert len(set(rcomp2domain[comp]) ^ set(comp2domain[comp])) == 0
        return True


    

    

    
