import pylon as pn
import networkx as nx
from distsys import *
from pylon.util import _Named
class Traverser(_Named):
    """ The backward and forward sweep is not implemented here. Instead it is
        implemented in this clase descendant named bfs_pf
        25 Feb 2011 19.32
        There are two ways to implement the sweeps:
        1. per-level basis : not tried
        2. per-lateral basis : tried but not finished.
        which one to implement?
        25 Feb 2011 09.53 AM
    """
    def __init__(self,case):
        self.case=case
        self.firstComps=None 
        self.nextComp=None
        self.lastComp=None
        self.currentComp=None
    def determineFeeder(self):
        """ tested and working for single end point. check it for multiendpoints.
        """
        feeders=[]
        self.firstComps=[bus for bus,path in nx.shortest_path(self.case.G, source=self.case.rootbus[0]).iteritems() if bus.typeG==endNode]
        hitung=0
        for fc in self.firstComps:
            f=self.oneFeeder(fc)
            feeders.append(f)
        return feeders
    
    def oneFeeder(self,fc):
        feeder=[]
        while (fc.typeG is not (forkNode or rootNode)) and fc!=[]:
            line=(fc.connected_from_branch)
            for line1 in line:
                feeder.append(line1)
            fc= self.case.G.predecessors(fc)
            if fc==[]:
                break
            for fcc in fc:
                fc=fcc
        return feeder
    def findFromBus(self,branch):
        return branch.from_bus
    def getBusesAtLevel(self,level):
        return [bus for bus in self.case.buses if self.case.bus_level[bus]==level and level>0]
    def getBranchesEndAtLevel(self,level):
        return [bus.connected_from_branch[0] for bus in self.case.buses if self.case.bus_level[bus]==level and level>0] 
    

#### ............................... Test ................................
##T=Traverser(kasus)
##FeederList=T.determineFeeder()
