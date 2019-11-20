#!/usr/bin/env python
# coding: utf-8

# In[1]:


#  G.V.S.Sri Ram, 15CE31015
#  Department of Civil Engineering, IIT, Kharagpur


# In[2]:


import math
import numpy as np
import random


# In[3]:


class Node(object):
    def __init__(self,**kwargs):
        for key,value in kwargs.items():
            vars(self)[key]=value


# In[4]:


class Element(object):
    def __init__(self,**kwargs):
        for key,value in kwargs.items():
            vars(self)[key]=value


# In[5]:


class ElementTrussXY(Element):
    def __init__(self,**kwargs):
        Element.__init__(self,**kwargs)
                
        self.mdof=4        
        self.ke=np.zeros((4,4),float)
        self.ok=np.zeros((4,4),float)
        self.ot=np.zeros((4,4),float)
        self.jdofv=np.zeros((4),int)

                         
    def computeStiffnessMatrix(self,NODE):
        
        self.c=np.zeros((self.mdof,self.mdof),float)
        
        for i in range(self.mdof):
            for j in range(self.mdof):
                self.c[i,j]=0.0

        n1=self.nodes[0]
        n2=self.nodes[1]
        self.jdofv[0]=2*n1-1
        self.jdofv[1]=2*n1
        self.jdofv[2]=2*n2-1
        self.jdofv[3]=2*n2
                         
        xx=NODE[n2].x-NODE[n1].x
        yy=NODE[n2].y-NODE[n1].y

        A=self.A
 
        E=self.E

        L=math.sqrt((xx*xx+yy*yy))
        for i in range(4):
            for j in range(4):
                self.ok[i,j]=0.0
                self.ot[i,j]=0.0
            self.ot[i,i]=1.0    

        C = xx / L
        S = yy / L
        self.ot[0,0]=C
        self.ot[0,1]=S
        self.ot[1,0]=-S
        self.ot[1,1]=C
        self.ot[2,2]=C
        self.ot[2,3]=S
        self.ot[3,2]=-S
        self.ot[3,3]=C

        AE_by_L = A * E / L

        self.ke[0,0]=      AE_by_L
        self.ke[0,2]=     -AE_by_L
    
        self.ke[2,0]=     -AE_by_L
        self.ke[2,2]=      AE_by_L

        for i in range(4):
            for j in range(4):
                self.c[i,j]=0.0
                for k in range(4):
                    self.c[i,j]=self.c[i,j]+self.ot[k,i]*self.ke[k,j]

        for i in range(4):
            for j in range(4):
                self.ok[i,j]=0.0
                for k in range(4):
                    self.ok[i,j]=self.ok[i,j]+self.c[i,k]*self.ot[k,j]

        del self.c

    def assembleStiffness(self,gk):
        
        for i in range(self.mdof):
            ii=self.jdofv[i]
            if ii > 0:
                for j in range(self.mdof):
                    jj=self.jdofv[j]
                    if jj > 0:
                       gk[ii-1,jj-1]=gk[ii-1,jj-1]+self.ok[i,j]
    
    def computeLength(self,NODE):
        n1=self.nodes[0]
        n2=self.nodes[1]
        N1=NODE[n1]
        N2=NODE[n2]
        self.l = math.sqrt((N1.x - N2.x)**2+(N1.y - N2.y)**2)
    
    
    def computeMemberForces(self,NODE):
            self.computeStiffnessMatrix(NODE)
            mfg=[0,0,0,0]
            mfl=[0,0,0,0]
            disp=[0,0,0,0]
            n1=self.nodes[0]
            n2=self.nodes[1]
            N1=NODE[n1]
            N2=NODE[n2]
            disp[0]=N1.Dx
            disp[1]=N1.Dy
            disp[2]=N2.Dx
            disp[3]=N2.Dy
           
            for i in range(4):
                mfg[i]=0.0
                for j in range(4):
                    mfg[i]=mfg[i]+self.ok[i,j]*disp[j]
                    
            for i in range(4):
                mfl[i]=0.0
                for j in range(4):
                    mfl[i]=mfl[i]+self.ot[i,j]*mfg[j]

            self.mfl=mfl
    


# In[6]:


class Structure(object):
   
    def __init__(self,**kwargs):
        self.title='Untitled'
        self.numnode=0
        self.numelem=0
        self.NODE=dict()
        self.ELEM=dict()
        self.NODE_LIST=list()
        self.ELEM_LIST=list()

    def node(self, **kwargs):
        if 'NODE' not in vars(self):
            self.NODE=dict()
        if 'NODE_LIST' not in vars(self):
            self.NODE_LIST=list()
        
        if 'nid' in kwargs:
            nid=kwargs['nid']
            self.NODE[nid]=Node(**kwargs)
            self.NODE_LIST.append(nid)
            self.numnode=self.numnode+1
                

    def element(self, **kwargs):
        if 'ELEM' not in vars(self):
            self.ELEM=dict()
        if 'ELEM_LIST' not in vars(self):
            self.ELEM_LIST=list()

        if 'eid' in kwargs:
            eid=kwargs['eid']
            if 'etype' in kwargs:
                self.etype=kwargs['etype']
            if self.etype == 'TrussXY':
               self.ELEM[eid]=ElementTrussXY(**kwargs)
            
            
            self.ELEM_LIST.append(eid)
            self.numelem=self.numelem+1
            
    def solve(self,**kwargs):
        self.ndof=2*self.numnode
        self.gk=np.zeros((self.ndof,self.ndof),float)
        self.gp=np.zeros((self.ndof,1),float)

        for eid in self.ELEM_LIST:
            self.ELEM[eid].computeStiffnessMatrix(self.NODE)
            self.ELEM[eid].assembleStiffness(self.gk)

        for nid in self.NODE_LIST:
            N=self.NODE[nid]
            i1=2*nid-1
            i2=2*nid
            
            if 'Fx' in vars(N):
                self.gp[i1-1,0]=N.Fx

            if 'Fy' in vars(N):
                self.gp[i2-1,0]=N.Fy

            if 'idx' in vars(N):
                idx=N.idx
                if idx == 1:
                    for i in range(self.ndof):
                        self.gk[i,i1-1]=0.0
                        self.gk[i1-1,i]=0.0
                    self.gk[i1-1,i1-1]=1.0

            if 'idy' in vars(N):
                idy=N.idy
                if idy == 1:
                    for i in range(self.ndof):
                        self.gk[i,i2-1]=0.0
                        self.gk[i2-1,i]=0.0
                    self.gk[i2-1,i2-1]=1.0

        self.disp=np.linalg.solve(self.gk,self.gp)

        for nid in self.NODE_LIST:
            i1=2*nid-1
            i2=2*nid
            self.NODE[nid].Dx=self.disp[i1-1,0]
            self.NODE[nid].Dy=self.disp[i2-1,0]

        for eid in self.ELEM_LIST:
            self.ELEM[eid].computeMemberForces(self.NODE)
            self.ELEM[eid].computeLength(self.NODE)
            
    def showStructure(self,**kwargs):
        for key,value in kwargs.items():
            vars(self)[key]=value
        plt.axis((-5.0,55.0,-5.0,55.0))
        ax=plt.gca()
        plt.axis('off')
        for eid in self.ELEM_LIST:
            n1=self.ELEM[eid].nodes[0]
            n2=self.ELEM[eid].nodes[1]
            N1=self.NODE[n1]
            N2=self.NODE[n2]
            p1=[N1.x,N1.y]
            p2=[N2.x,N2.y]
#            print(eid,n1,n2,p1,p2)
            l=mlines.Line2D([N1.x,N2.x],[N1.y,N2.y])
            ax.add_line(l)
            
        plt.show()
            


# In[7]:


def sri_truss_bridge_xy():
        a=5.0
        h=7.0
        L=6*a
        A_bot=0.1
        A_top=0.1
        A_diag=0.1
        A_vert=0.1
        Iz_bot=0.5e-05
        Iz_top=0.5e-05
        Iz_diag=0.5e-05
        Iz_vert=0.5e-05
        E=2.0e10
        rho=7850
        m_bar=rho*A_bot
        pstr=Structure(etype='TrussXY',title="Truss Bridge - span 30 m")
        
        pstr.node(nid=0,tagid='L0',x=0.0,  y=0, idx=1, idy=1 )
        pstr.node(nid=1,tagid='L1',x=a, y=0, Fy=-250000.0 )
        pstr.node(nid=2,tagid='L2',x=2*a, y=0, Fy=-250000.0 )
        pstr.node(nid=3,tagid='L3',x=3*a, y=0, Fy=-250000.0 )
        pstr.node(nid=4,tagid='L4',x=4*a, y=0, Fy=-250000.0 )
        pstr.node(nid=5,tagid='L5',x=5*a, y=0, Fy=-250000.0 )
        pstr.node(nid=6,tagid='L6',x=6*a, y=0, idy=1 )
        
        pstr.node(nid=7,tagid='U1',x=a, y=h )
        pstr.node(nid=8,tagid='U2',x=2*a, y=h )
        pstr.node(nid=9,tagid='U3',x=3*a, y=h )
        pstr.node(nid=10,tagid='U4',x=4*a, y=h )
        pstr.node(nid=11,tagid='U5',x=5*a, y=h )

        pstr.element(eid=1,  tagid='L0-L1',etype='TrussXY',nodes=(0,1),
                     A=A_bot,E=E,Iz=Iz_bot,rho=rho,m_bar=m_bar)
        pstr.element(eid=2,  tagid='L1-L2',nodes=(1,2),
                     A=A_bot,E=E,Iz=Iz_bot,rho=rho,m_bar=m_bar)
        pstr.element(eid=3,  tagid='L2-L3',nodes=(2,3),
                     A=A_bot,E=E,Iz=Iz_bot,rho=rho,m_bar=m_bar)
        pstr.element(eid=4,  tagid='L3-L4',nodes=(3,4),
                     A=A_bot,E=E,Iz=Iz_bot,rho=rho,m_bar=m_bar)
        pstr.element(eid=5,  tagid='L4-L5',nodes=(4,5),
                     A=A_bot,E=E,Iz=Iz_bot,rho=rho,m_bar=m_bar)
        pstr.element(eid=6,  tagid='L5-L6',nodes=(5,6),
                     A=A_bot,E=E,Iz=Iz_bot,rho=rho,m_bar=m_bar)

        
        pstr.element(eid=7,  tagid='U1-U2',nodes=(7,8),
                     A=A_top,E=E,Iz=Iz_top,rho=rho,m_bar=m_bar)
        pstr.element(eid=8, tagid='U2-U3',nodes=(8,9),
                     A=A_top,E=E,Iz=Iz_top,rho=rho,m_bar=m_bar)
        pstr.element(eid=9, tagid='U3-U4',nodes=(9,10),
                     A=A_top,E=E,Iz=Iz_top,rho=rho,m_bar=m_bar)
        pstr.element(eid=10, tagid='U4-U5',nodes=(10,11),
                     A=A_top,E=E,Iz=Iz_top,rho=rho,m_bar=m_bar)

        
        pstr.element(eid=11, tagid='L0-U1',nodes=(0,7),
                     A=A_diag,E=E,Iz=Iz_diag,rho=rho,m_bar=m_bar)
        pstr.element(eid=12, tagid='L2-U1',nodes=(2,7),
                     A=A_diag,E=E,Iz=Iz_diag,rho=rho,m_bar=m_bar)
        pstr.element(eid=13, tagid='L2-U3',nodes=(2,9),
                     A=A_diag,E=E,Iz=Iz_diag,rho=rho,m_bar=m_bar)
        pstr.element(eid=14, tagid='L4-U3',nodes=(4,9),
                     A=A_diag,E=E,Iz=Iz_diag,rho=rho,m_bar=m_bar)
        pstr.element(eid=15, tagid='L4-U5',nodes=(4,11),
                     A=A_diag,E=E,Iz=Iz_diag,rho=rho,m_bar=m_bar)
        pstr.element(eid=16, tagid='L6-U5',nodes=(6,11),
                     A=A_diag,E=E,Iz=Iz_diag,rho=rho,m_bar=m_bar)
 
        
        pstr.element(eid=17, tagid='L1-U1',nodes=(1,7),
                     A=A_vert,E=E,Iz=Iz_vert,rho=rho,m_bar=m_bar)
        pstr.element(eid=18, tagid='L2-U2',nodes=(2,8),
                     A=A_vert,E=E,Iz=Iz_vert,rho=rho,m_bar=m_bar)
        pstr.element(eid=19, tagid='L3-U3',nodes=(3,9),
                     A=A_vert,E=E,Iz=Iz_vert,rho=rho,m_bar=m_bar)
        pstr.element(eid=20, tagid='L4-U4',nodes=(4,10),
                     A=A_vert,E=E,Iz=Iz_vert,rho=rho,m_bar=m_bar)
        pstr.element(eid=21, tagid='L5-U5',nodes=(5,11),
                     A=A_vert,E=E,Iz=Iz_vert,rho=rho,m_bar=m_bar)

        
        return pstr

#        pstr.solve()
        
#        print('Nodes')
#        for nid in pstr.NODE_LIST:
#            N=pstr.NODE[nid]
#            print(nid,N.x,N.y,N.Dx,N.Dy)
            
#        print('Elements')

#        for eid in pstr.ELEM_LIST:
#            E=pstr.ELEM[eid]
#            print(eid,E.nodes,E.A,E.E,E.mfl[0])


# In[8]:


# Prepare model

pstr=sri_truss_bridge_xy()

pstr.solve()

print ('Nodes')

for nid in pstr.NODE_LIST:
    N=pstr.NODE[nid]
    print (vars(N))

print('Nodes')
for nid in pstr.NODE_LIST:
    N=pstr.NODE[nid]
    print(N.nid,N.x,N.y)
    
print('Elements')

for eid in pstr.ELEM_LIST:
    E=pstr.ELEM[eid]
    print(E.eid,E.nodes,E.A,E.E,E.mfl[0])
    


# In[9]:


import matplotlib.pyplot as plt
import matplotlib.lines as mlines

pstr.showStructure()


# In[10]:


max_stress=120*10**6
min_A=0.025

#step-1
hms=20
hmcr=0.8
par=0.3
n_max=30000


# In[11]:


def create_hm(min_A,j):
    x=random.uniform(min_A,1)
    E=pstr.ELEM[pstr.ELEM_LIST[j]]
    if abs(E.mfl[0]/x) < max_stress:
        return x
    else:
        create_hm(min_A,j)


# In[12]:


#step-2
hm=[]
for i in range(hms):
    arr=[]
    for j in range(21):
        x=create_hm(min_A,j)
        arr.append(x)
            
    hm.append(arr)
    
print(hm)


# In[13]:


def get_weight(hm,hms):
    w=[]
    for i in range(hms):
        sum_row=0
        for j in range(21):
            E=pstr.ELEM[pstr.ELEM_LIST[j]]
            sum_row+=hm[i][j]*E.rho*E.l
        w.append(sum_row)
    return w


# In[14]:


#step-3
def new_vector(hm,hms,hmcr,par):
    x_new=[]
    r1=random.random()
    if(r1<hmcr):
        arr=[]
        for j in range(21):
            for i in range(hms):
                arr.append(hm[i][j])
            x_new.append(random.choice(arr))
    else:
        for i in range(21):
            x_new.append(create_hm(min_A,i))
    r2=random.random()
    if(r2<par):
        for i in range(21):
            bw=0.01
            x_new[i]+=bw*random.uniform(-1,1)
    return x_new


# In[15]:


w=get_weight(hm,hms)
print(w)
worst=w.index(max(w))
weights=[min(w)]
for i in range(n_max):
    x_new=new_vector(hm,hms,hmcr,par)
    #print(x_new)
    sum_row=0
    for j in range(len(x_new)):
        E=pstr.ELEM[pstr.ELEM_LIST[j]]
        sum_row+=x_new[j]*E.rho*E.l
    #print(sum_row,"\t",max(w))
    if sum_row<max(w):
        hm[worst]=x_new
        #print(hm)
        w=get_weight(hm,hms)
        worst=w.index(max(w))
    if i==1000 || i==2500 || i==5000 || i==10000 || i==20000:
        print("hm",hm)
        print("w",w)
    weights.append(min(w))

#print(w)
#print(weights)


# In[16]:


import matplotlib.pyplot as plt

g = plt.figure()
plt.plot(weights)
plt.xlabel("No. of Iterations")
plt.ylabel("Weight (Kg)")
plt.show()

g.savefig("harmony.pdf")


# In[ ]:




