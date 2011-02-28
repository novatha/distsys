from bfs_pf2 import *
MVAbase=1
kVbase=7.2
Zbase=kVbase**2/MVAbase#in the present of transformer, zbase as well as kVbase on both side of transformer are differents.
# Test Line3 
# this is impedance for three wire configuration
zd=[[0.4013+1.4133*1j, 0.0953+0.8515*1j, 0.0953+0.7266*1j],
    [0.0953+0.8515*1j, 0.4013+1.4133*1j, 0.0953 +0.7802*1j],
    [0.0953+0.7266*1j, 0.0953+0.7802*1j, 0.4013+1.4133*1j]
    ]
# this is impedance for 4-wire configuration
zd=[[0.4576+1.078*1j, 0.1559+0.5017*1j, 0.1535+0.3849*1j],
    [0.1559+0.5017*1j, 0.4666+1.0482*1j, 0.158+0.4236*1j],
    [0.1535+0.3849*1j, 0.158+0.4236*1j, 0.4615+1.0651*1j]
    ]

#zd=np.eye(3,3)
zd=np.array(zd)
Z1=zd*2000/(5280*Zbase)


# Test Case3
busa=pn.Bus()
bus1=RootBus3(busa,busa,busa,base_kv=kVbase)
bus2=Bus3(busa,busa,busa,base_kv=kVbase)
bus3=Bus3(busa,busa,busa,base_kv=kVbase)
bus4=Bus3(busa,busa,busa,base_kv=kVbase)
Pa=1.8
pfa=0.9
Pb=1.8
pfb=0.9
Pc=1.8
pfc=0.9
Sa=P_to_S(Pa,pfa)
Sb=P_to_S(Pb,pfb)
Sc=P_to_S(Pc,pfc)
S3=np.array([Sa,Sb,Sc])
bus2.loadS=0
bus3.loadS=0
bus4.loadS=S3
#bus4=Bus3(busa,busa,busa)
bus123=[bus1,bus2,bus3,bus4]
gen1=Generator3(bus1)
# line data
branch1=Line3(bus1,bus2,Z1)
branch2=Line3(bus2,bus3,Z1)

#transformer data
t_connection='Step-Down'
t_MVA=6.0
t_kVLL_high=12.47
t_zbase= t_kVLL_high**2/t_MVA
t_kVLL_low=4.16
t_R=0.01 # pu
t_X=0.06# pu
t_a=(MVAbase/t_MVA)*(t_kVLL_high/kVbase)**2 # adjustment factor
t_z=(t_R + 1j* t_X)*t_a
t_y=1/t_z
trafo1=Trafo3(bus2,bus3,yt=t_y,base_mva=t_MVA)
# line data
#new kva_base
kVbase3=t_kVLL_low*kVbase/t_kVLL_high
zbase3=kVbase3**2/MVAbase
Z2=zd*2500/(5280*zbase3)
bus4.base_kv=kVbase3
branch3=Line3(bus3,bus4,Z2,base_kv=kVbase3)
# branchess
branch12=[branch1,trafo1,branch3]
kasus=Case3(base_mva=MVAbase, base_kv=kVbase, buses=bus123, branches=branch12,generators=[gen1])



## ............................... Test ................................
T=Traverser(kasus)
FeederList=T.determineFeeder()
# test 
# print branch1.I_line
bs=bfs_pf(kasus)
bs.solve()
#print np.abs(bs.presentE)

# print branch1.I_line
#bs.updateVoltage()
# print branch1.I_line
E=[b.E*b.base_kv for b in kasus.buses]
for bb in E:
    print
    for bbb in bb:
        print "%.0f / %.1f " % (np.abs(bbb*1000),np.angle(bbb)*180/np.pi)
    print
I=[br.I_line for br in kasus.branches]

#--------------------------------IEEE 4 bus feeder----------------------------
Im=np.array([347.9,323.7,336.8])
Iang=np.array([-34.9,-154.2,85.0])*np.pi/180
I12=polarToRect(Im,Iang)
Vm=np.array([1918,2061,1981])
Vang=np.array([-9.1,-128.3,110.9])
V4=polarToRect(Vm,Vang)