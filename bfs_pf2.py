import pylon as pn
import networkx as nx
from scipy.sparse import hstack, vstack,csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, splu
#import distsys as ds
from distsys import *
#from pylon.util import Named, Serializable
from iterator import *

class bfs_pf(Traverser):
    """Implement backward and forward sweep method
    """
    def __init__(self,case):
        Traverser.__init__(self,case)
        self.max_iter=200
        self.tolerance=1e-6
        self.verbose=False
        self.oldE=np.array([b.E for b in self.case.buses])
        self.presentE=np.array([b.E for b in self.case.buses])
 
    def updateCurrent(self):
        """ Performing backward sweep level by level. First, please write the steps needed in this sweep.
        """
        max_level=max(self.case.bus_level.values())

        # start with the maximum level. decrease step bu step.
        k=max_level # the lecel counter
        while k>0 :
            buses=[bus for bus in self.getBusesAtLevel(k)]
            # print 'k= ',k,'\nbuses =',buses
            for bus in buses:
                # print '\t bus',bus
                #print 'loadS = ',bus.loadShttp://www.cise.ufl.edu/research/sparse/umfpack/
                #print 'loadI = ',bus.loadI        
                # update load current at each bus at level k.
                bus.calc_loadI()
                #print 'loadI = ',bus.loadI        
                # update shunt current
                #print 'shuntI = ',bus.shuntI        
                bus.calc_shuntI()
                # update current to branches connected to this bus from higher level.
                bus.calc_childI()
                # update branch curret at each branch ending at level k.
                # print 'totalI = ',bus.totalI        
                totalI=bus.calc_totalI()
                # print 'totalI = ',bus.totalI                        
                for branch in bus.connected_from_branch:
                    #branch.I_to=totalI
                    branch.updateBackward()
            # go to the k-1 level
            # print 
            k=k-1
    def updateVoltage(self):
        """ Performing forward sweep level by level. starting from level 1
        """
        # calculate max level
        max_level=max(self.case.bus_level.values())
        # start from level 1
        k=1
        while k<=max_level:
            # pick a bus
            for bus in self.getBusesAtLevel(k):
                # pick a branch having the bus as its from_bus
                for branch in bus.connected_from_branch:
                    # update voltage at from_bus
                    branch.updateEfrom()
                    # update voltages at to_bus
                    branch.updateEto() # line and transformer may different implementation here
                    # set the bus voltage to branch.E_to
                    branch.to_bus.E=branch.E_to
            # after all buses are updated, go to next level
            k=k+1
    def oneIteration(self):
        self.updateCurrent()
        self.updateVoltage()
        self.presentE=np.array([b.E for b in self.case.buses])
        #print np.max(np.abs(self.oldE-self.presentE))
   
    def solve(self):
        repeat= True
        k=0
        while repeat and k<self.max_iter:
            self.oldE=np.array([b.E for b in self.case.buses])
            self.oneIteration()
            converged=self.checkConvergence()
            repeat=not converged
            k=k+1
        if converged:
            print 'Conververged after ',k,'iterations.'
        if not converged:
            print 'Sorry! I gave up after ',self.max_iter,' iterations.'
    def checkConvergence(self):
        mis=self.oldE-self.presentE
        F = np.r_[mis.real, mis.imag]
        normF = np.linalg.norm(F, np.Inf)
        if normF < self.tolerance:
            converged = True
        else:
            converged = False
            if self.verbose:
                logger.info("Difference: %.3f" % (normF - self.tolerance))
        return converged
        
    def printE(self):
        for bus in self.case.buses:
            print bus, bus.E
    
### test 
### print branch1.I_line
##bs=bfs_pf(kasus)
##
##
### print branch1.I_line
###bs.updateVoltage()
### print branch1.I_line
