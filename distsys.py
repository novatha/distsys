import pylon as pn
import networkx as nx
import matplotlib.pyplot as pl
import numpy as np
import scipy as sc
from scipy.sparse import hstack, vstack, csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, splu
"""
b1=pn.Bus()
b2=pn.Bus()
b3=pn.Bus()
G=nx.DiGraph()

G.add_node(0)
G.add_node(1)
G.add_node(2)
G.add_edges_from([(0,1),(1,2)])
br1=pn.Branch(0,1)
br1=pn.Branch(1,2)

nx.draw(G)
pl.show()
"""
# General constant
pi=np.pi
a=np.cos(2*pi/3)+1j*np.sin(2*pi/3)
A=np.array([1,a,a*a])
rootNode="root node"
forkNode="fork node"
endNode="end node"
normalNode="normal node"
satu=np.ones(3)*A 
class Bus3(pn.Bus):
    """Three phase bus model
       I think that every bus should have self.kvbase attribut.
    """
    def __init__(self,voltage=satu,busa=None,busb=None,busc=None,base_kv=1.0):
        pn.Bus.__init__(self)
        self.Busa=busa
        self.Busb=busb
        self.Busc=busc
        self.base_kv=base_kv
        self.typeG=normalNode
        self.loadS=np.zeros(3,dtype=complex)   # 
        self.sourceS=np.zeros(3) # please assume generation as negative load
        self.shuntS=np.zeros(3)
        self.E=satu*A# balanced unity voltages
        self.connected_from_branch=[] # self is connected to root from branch
        self.connected_to_branch=[] # self is connected 
        self.loadI=np.zeros(3)
        self.sourceI=np.zeros(3)
        self.childI=np.zeros(3)
        self.shuntI=np.zeros(3)
        self.totalI=np.zeros(3)
    def __repr__(self):
        """ A function that modifies "print this_bus" result
        """
        return self.name
    
    def set_load(self,S):
        a,b,c=np.asarray(S)
        self.loadS=np.asarray(S)
    def get_load(self):
        
        return self.loadS
    
    def calc_loadI(self):
        """A function that calculate current of load directly connected 
           to this bus 
        """
        self.loadI=np.conj(self.loadS/self.E)
        return self.loadI
    def calc_sourceI(self):
        """ A function that calculate current from generator connected to this
            bus.
        """
        self.sourceI=-np.conj(self.sourceS/self.E) # assumed a negative load
        return self.sourceI
    def calc_shuntI(self):
        self.shuntI=np.conj(self.shuntS/self.E)
        return self.shuntI
    def calc_childI(self):
        """ A function that calculate currents flowing from this bus to all
            other child buses. 
        """
        self.childI=np.sum([b.I_from for b in self.connected_to_branch ],axis=0)
        return self.childI
    def calc_totalI(self):
        self.totalI=self.calc_loadI()+self.calc_sourceI()+self.calc_shuntI()+self.calc_childI()
        return self.totalI
    def addConnectedToBranch():
        pass
class Line3(pn.Branch):
    """ Three phase line model
        I think every line should have self.zbase atrribut
    """
    def __init__(self, from_bus, to_bus, Z=np.ones([3, 3]), Y=np.zeros([3, 3]),base_kv=1.0):
        pn.Branch.__init__(self,from_bus,to_bus)
        self.Z=Z
        self.Y=Y
        self.base_kv=base_kv
        self.E_from = from_bus.E#np.zeros(3)
        self.E_to=to_bus.E#np.zeros(3)
        self.I_from = np.zeros(3)
        self.I_to=np.zeros(3)
        #self.Iline=np.zeros(3)
        self.I_line=self.I_to+np.dot(self.Y,self.E_to)
        self.line_drop=np.zeros(3)
        #self.A=np.eye(3, 3) + np.dot(Z, Y/2)
        #self.B=Z
        #self.C=Y+np.dot(np.dot(Y, Z,), Y)
        #self.D=np.eye(3, 3) + np.dot(Y/2, Z)
    def getEfrom(self):
        result=self.E_to + np.dot(self.Z,self.Iline)
        #result=np.dot(self.A, self.E_to) + np.dot(self.B, self.I_to)
        return result
    def updateEfrom(self):
        self.E_from=self.from_bus.E
        return self.E_from
    def updateIfrom(self):
        self.I_from=self.I_line +np.dot(self.Y/2,self.E_from)
    def updateIline(self):
        self.I_line=self.I_to + .5*np.dot(self.Y,self.E_to)
    def updateIto(self):
        self.I_to=self.to_bus.totalI
        return self.I_to  
    def updateEto(self):
        self.E_to=self.E_from - np.dot(self.Z,self.I_line)
        self.to_bus.E=self.E_to
        return self.E_to
    def updateBackward(self):
        self.updateIto()
        self.updateIline()
        self.updateIfrom()
    def updateForward(self):
        self.updatedEfrom()
        self.updateEto()
    def updateLineDrop(self):
        self.line_drop=np.dot(self.Z,self.I_line)
        return self.line_drop
    def __repr__(self):
        strs=self.name+'\n'
        return strs
#        str2='E_from = ' +str(self.E_from)+ ' I_from = '+str(self.I_from)+'\n'
#        str3='I_to   = ' +str(self.E_to) + ' I_to   = ' + str(self.I_to)+'\n'    
    def backwardSweep(self):
        pass

#        strs=str1+str2+str3
        return strs
    

class Trafo3(Line3):
    """ A model for transformer based on:
    
        Xiao, P.; Yu, D. & Yan, W. A unified three-phase transformer model 
        for distribution load flow calculations Power Systems, 
        IEEE Transactions on, 2006, 21, 153 - 159
        
    """
    ## please evaluate pu conversion in tranformator.
    ## now updatebakcward
    
    def __init__(self,from_bus=Bus3(),to_bus=Bus3(),yt=0,base_mva= 1.0,base_kv=12.0,case_base_mva=1.0,case_base_kv=12.0,E_high= 12.47, E_low=4.16, type_connection='YgYg',type='Step-Down' ):
        Line3.__init__(self,from_bus,to_bus)
        # constants
        self.YI=np.eye(3)*yt
        self.YII=np.array([[2.,-1.,-1.],[-1.,2.,-1.],[-1.,-1.,2.]])*yt/3
        self.YIII=np.array([[-1.,1.,0.],[0.,-1.,1.],[1.,0.,-1.]])*yt/np.sqrt(3)
        self.type=type
        self.yt=yt
        self.type_connection=type_connection
        self.base_mva=base_mva
        self.base_kv=base_kv
        self.E_high=E_high
        self.E_low=E_low
        self.to_bus.base_kv=self.from_bus.base_kv*self.E_low/self.E_high
        if self.type_connection=='YgYg': # tested. there are more to model
            self.Ypp=self.YI
            self.Yps=-self.YI
            self.Ysp=-self.YI
            self.Yss=self.YI
        if self.type_connection=='YgY':
            self.Ypp=self.YII
            self.Yps=-self.YII
            self.Ysp=-self.YII
            self.Yss=self.YII
        if self.type_connection=='YgD':
            self.Ypp=self.YI
            self.Yps=self.YIII
            self.Ysp=self.YIII.transpose()
            self.Yss=self.YII
            
        self.YT=np.vstack([np.hstack([self.Ypp,self.Yps]),np.hstack([self.Ysp,self.Yss])])
    def updateIto(self):
        self.I_to=-self.to_bus.totalI # the minus sign is important because by definition it is the injected current to the secondary side.
        return self.I_to
    def updateIline(self):
        self.I_line =self.I_from
        return self.I_line
    def updateIfrom(self):
        self.updateEfromBackward()
        self.I_from= np.dot(self.Ypp,self.from_bus.E)+np.dot(self.Yps,self.E_to)
    def updateBackward(self):
        self.updateIto()
        # not relevan for a transformer self.updateIline()
        # updating from_bus voltage
        if self.type_connection=='YgYg':
            #print 'Transformer type is',self.type, ' and ',self.type_connection,' connected.'
            #self.updateEfromBackward()
            # updating primary current
            self.updateIfrom()
            self.updateIline()
            self.from_bus.totalI=self.I_from
        else:
            # calculate zero componen of voltages
            self.E_from0=np.sum(self.E_from)/3
            self.E_to0=np.sum(self.E_t0)/3
            # modify admitanve matrices
            self.E_original=self.E
            self.E_from[2]=np.ones(3)
            self.E_to_original=self.E_to
            self.E_to[2]=np.ones(3)
            self.Ysp[2]=np.ones(3)
            self.Yss[2]=np.zeros(3)
            
            pass
       
    def updateEfromBackward(self):
        b=csr_matrix(self.I_to - np.dot(self.Yss, self.E_to))
        a=csr_matrix(self.Ysp)
        self.from_bus.E=spsolve(a,b)
        return self.from_bus.E# voltage is not update in backward sweep
    def updateEto(self):
        b=csr_matrix(self.I_to-np.dot(self.Ysp,self.from_bus.E))
        a=csr_matrix(self.Yss)
        self.E_to=spsolve(a,b)
        self.to_bus.E=self.E_to
        return self.E_to
class RootBus3(Bus3):
    def __init__(self,busa=None,busb=None,busc=None,base_kv=1):
        Bus3.__init__(self,busa=None,busb=None,busc=None,base_kv=base_kv)

class EndBus3(Bus3):
    def __init__(self,busa=None,busb=None,busc=None):
        Bus3.__init__(self,busa=None,busb=None,busc=None)
        
class ForkBus3(Bus3):
    def __init__(self,busa=None,busb=None,busc=None):
        Bus3.__init__(self,busa=None,busb=None,busc=None)
class Generator3(pn.Generator):
    """Three phase Generator Model
    """
    def __init__(self, bus3):
        pn.Generator.__init__(self,bus3)
        self.E=np.ones(3)*A
        self.gena=pn.Generator(bus3.Busa) 
        self.genb=pn.Generator(bus3.Busb)
        self.genc=pn.Generator(bus3.Busc)
    
class Case3(pn.Case):
    """This class models a power system studied. It uses networkx packacge to 
       build its topology. I am thinking about using it to implement BFS power flow.
    """
    def __init__(self,name=None, base_mva=1.0, base_kv=20.0, buses=None, branches=None,  generators=None):
        
    #pn.Case.__init__(self,buses=buses,branches=branches)
        self.buses=buses
        self.base_mva=base_mva
        self.base_kv=base_kv
        self.base_z=base_kv**2/base_mva
        self.branches=branches
        self.generators=generators    
        self.G=nx.DiGraph()
        self._buildGraph()
        self.rootbus=[b for b in self.buses if self.G.pred[b] == {}]
        self.forkbuses=[b for b in self.buses if self.G.degree()[b]>2]
        self.normalbuses=[b for b in self.buses if self.G.degree()[b]==2]    
        self.endbuses=[b for b in self.buses if self.G.succ[b] == {}]
        self.updateTypeG()
        self.bus_level=nx.shortest_path_length(self.G,source=self.rootbus[0])
        self.updateBusConnBranch()
        self.buses_loadJ=self.updateBusLoadCurrent()    
        #nx.draw(self.G)
        #pl.show()
    def _buildGraph(self):
        """ add buses as nodes  and branches as edges. updates connected
        """
        self.G.add_nodes_from([bus for bus in self.buses])
        self.G.add_edges_from([(b.from_bus,b.to_bus) for b in self.branches])

    def updateBusConnBranch(self):
        for a in self.branches:
            a.from_bus.connected_to_branch.append(a)
            a.to_bus.connected_from_branch.append(a)
    # as edges are added, there should be a procedure to update bus.connected_to and bus.connected_from
    def updateTypeG(self):
        for b in self.buses:
            if self.rootbus.__contains__(b):
                b.typeG=rootNode
            if self.forkbuses.__contains__(b):
                b.typeG=forkNode
            if self.endbuses.__contains__(b):
                b.typeG=endNode
    def updateBusLoadCurrent(self):
       return [b.calc_loadI() for b in self.buses]
            
class RadNet(nx.DiGraph):
    def __init__(self):
            self.type = "Radial Distribution Network"
            nx.DiGraph.__init__(self)

def P_to_S(P,pf):
    theta=np.arccos(pf)
    Samp=P/pf
    S=Samp*np.cos(theta)+1j*Samp*np.sin(theta)
    return S

def polarToRect(magnitude,degree):
    result=magnitude*(np.cos(degree)+1j*np.sin(degree))
    return result

#-----------------------------------------------------------------------------#
### Test Line3
##zd=[[0.4013+1.4133*1j, 0.0953+0.8515*1j, 0.0953+0.7266*1j],
##    [0.0953+0.8515*1j, 0.4013+1.4133*1j, 0.0953 +0.7802*1j],
##    [0.0953+0.7266*1j, 0.0953+0.7802*1j, 0.4013+1.4133*1j]
##    ]
###zd=np.eye(3,3)
##zd=np.array(zd)
##Z1=zd*2000/5280
##Z2=zd*2500/5280
##
### Test Case3
##busa=pn.Bus()
##bus1=RootBus3(busa,busa,busa)
##bus2=Bus3(busa,busa,busa)
##bus3=Bus3(busa,busa,busa)
##bus4=Bus3(busa,busa,busa)
##Pa=1.8
##pfa=0.9
##Pb=1.8
##pfb=0.9
##Pc=1.8
##pfc=0.9
##Sa=P_to_S(Pa,pfa)
##Sb=P_to_S(Pb,pfb)
##Sc=P_to_S(Pc,pfc)
##S3=np.array([Sa,Sb,Sc])/3
##bus2.loadS=S3
##bus3.loadS=S3
##bus4.loadS=S3*.5
###bus4=Bus3(busa,busa,busa)
##bus123=[bus1,bus2,bus3,bus4]
##gen1=Generator3(bus1)
##branch1=Line3(bus1,bus2,Z1)
##branch2=Line3(bus2,bus3,Z2)
##branch3=Line3(bus2,bus4,Z2)
##branch12=[branch1,branch2,branch3]
##kasus=Case3(buses=bus123, branches=branch12,generators=[gen1])

