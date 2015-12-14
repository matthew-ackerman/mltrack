from sympy import *
import sys
import random
import numpy
import math

#taking it the dumb way cause I can.
def stdev (data, name):
	m=mean(data, name)
	SS=0
	for x in data.data[name]:
		SS+=(x-m)**2
	return math.sqrt(SS/float(len(data.data[name]) ))

def mean (data, name):
	SUM=0
	for x in data.data[name]:
		SUM+=x
	return float(SUM)/float(len(data.data[name]) )
def lmean (data, name):
	SUM=0
	for x in data.data[name]:
		SUM+=ln(x)
	return float(SUM)/float(len(data.data[name]) )

def lstdev (data, name):
        m=mean(data, name)
        SS=0
        for x in data.data[name]:
                SS+=(ln(x)-m)**2
        return math.sqrt(SS/float(len(data.data[name]) ))

class p:
	def __init__ (self, name, LB, UB):
		self.name=name
		self.sym=Symbol(name)
		self.val=0
		self.UB=UB
		self.LB=LB
		self.R=UB-LB
		if self.R<=0:
			print "range error: parameter "+name+" can take no values"
			quit()
class d:
	def __init__ (self, D):
		self.D=len(D)
		self.Ds=[]
		self.sym=[]
		self.names=[]
		for ps in D:
			self.Ds.append(p(ps[0], ps[1], ps[2]) )
		for ps in self.Ds:
			self.sym.append(ps.sym)
			self.names.append(ps.name)

	def read (self, filename):
		File=open(filename)
		self.data={}
		self.sets=[]
		for name in self.names:
			self.data[name]=[]
		for line in File:
			datum=[]
			line=line.strip('\n').split(',')
		#	try:
			for x in range(0, self.D):
				t=float(line[x])
				if t<=self.Ds[x].UB and t>=self.Ds[x].LB:
					datum.append(t)
					self.data[self.names[x]].append(t)
			self.sets.append(datum)
		#	except:
		#		print "error reading "+'\t'.join(line)
		File.close()
class system : 
	def __init__ (self, eqs, parms, data):
		self.eqs=eqs
		self.N=len(eqs)
		self.parms=parms
		self.data=data
		self.R=float("infinity")
		if len(parms)!=self.N:
			print "system does not have one equation for each parameter. Matt cannot maximize this."
			quit()
		#check to see if any of the equations have algabaic solutions
		#???

	def inc (self):
		self.Rval=zeros((self.N, 1) )
		self.Req=[]
		self.Jval=zeros(self.N)
		test=zeros(self.N)
		self.Jeq=[]
		for x in range(0, self.N):
			self.Jeq.append([])
			R=self.eqs[x]
			for parm in self.parms:
				R=R.subs(parm.sym, parm.val)
			self.Req.append(lambdify(self.data.sym, R) )
			for y in range(0, self.N):
				K=diff(self.eqs[x], parms[y].sym)
				for parm in self.parms:
					K=K.subs(parm.sym, parm.val)
				self.Jeq[x].append(lambdify(self.data.sym, K) )
		for X in self.data.sets:
			for x in range(0, self.N):
				for y in range(0, self.N):
					self.Jval[x,y]+=self.Jeq[x][y](*X)
				try:
					self.Rval[x]+=self.Req[x](*X)
				except:
					self.Rval[x]+=float("infinity")
		IJ=self.Jval.inv()
		self.R=0
		for x in range(0, self.N):
			test[x]=self.parms[x].val
		for x in range(0, self.N):
			for y in range(0, self.N):
				test[x]-=(IJ[x,y]*self.Rval[y])
			self.R+=abs(self.Rval[x])
		for x in range(0, self.N):
		#	if test[x]>self.parms[x].LB and test[x]<self.parms[x].UB:
			self.parms[x].val=test[x]
		#	else:
		#		print "MLE fail to converge : setting "+self.parms[x].name+" to "+str(test[x])
class eq :
	def __init__ (self, eq, parms, data):
		for x in range (0, len(parms) ):
			eq=eq.replace(parms[x].name, "parms["+str(x)+"].sym")
		eq=eq.replace(data.name, "data.sym")
		self.parms=parms
		self.data=data
		self.N=len(parms)
		exec("self.eq="+eq)

class surface:
	def __init__ (self, dim, size):
		self.dim=dim
		self.notdone=True
		self.val=[]
		self.size=size
		self.alls=self.size**self.dim
		self.flat=0
		self.at=zeros(1, dim)
		self.maximum=float("-infinity")
	def setto (self, x):
		while (len(self.val)<=self.flat):
			self.val.append(0)
		if x>self.maximum:
			self.maximum=x
			self.maxcor=self.at[:]
		self.val[self.flat]=x
	def inc (self):
		if (self.notdone):
			self.at[0]+=1
			self.flat+=1
			i=0
			while (i<self.dim):
				if (self.at[i]>=self.size):
					self.at[i]-=self.size
					self.at[i+1]+=1
				i+=1
			if self.flat==self.size**(self.dim)-1:
				self.notdone=False
	def cor(self, x):
		return self.at[x]
	def getpeeks(self):
		return [self.maxcor]
		ret=[]
		self.at=zeros(1, self.dim)
		self.flat=0
		self.notdone=True
		while (self.notdone):
			self.inc()
			good=True
			for i in range(0, self.dim):
				if self.at[i]<1 or self.at[i]>self.size-2:
					good=False
			peak=False
			if (good):
				peak=True
				for i in range(0, self.dim):
					if self.val[self.flat]<self.val[self.flat-self.size**i]:
						peak=False
					if self.val[self.flat]<self.val[self.flat+self.size**i]:
						peak=False
			if (peak):
				K=self.at[:]
				ret.append(K)
		return ret
class Leq :
	def __init__ (self, eq, parms, data, surf_size):
		for x in range (0, len(parms) ):
			eq=eq.replace(parms[x].name, "parms["+str(x)+"].sym")
		for x in range (0, len(data.names) ):
			eq=eq.replace(data.names[x], "data.sym["+str(x)+"]")
		self.parms=parms
		self.data=data
		self.xeq=[]
		self.N=len(parms)
		exec("self.eq="+eq)
		for x in range(0, self.N):
			deq=diff(self.eq, parms[x].sym)
			self.xeq.append(deq)
		self.S=system(self.xeq, parms, data)
		self.R=self.S.R
		if surf_size>0:
			surf=surface(self.N, surf_size)
			while(surf.notdone):
				for x in range(0, self.N):
					parms[x].val=float(surf.cor(x) )/surf.size*parms[x].R+parms[x].LB
				surf.setto(self.get_Ls() )
				surf.inc()
			self.peaks=surf.getpeeks()
			print len(self.peaks)
			for peak in self.peaks:
				print peak
			for x in range(0, self.N):
				parms[x].val=float(self.peaks[0][x] )/surf.size*parms[x].R+parms[x].LB

		self.L=0
		self.I=0
		self.E=0
	
	def get_L (self):
		self.L=0
		T=self.eq
		for parm in self.parms:
			T=T.subs(parm.sym, parm.val)
		self.fd=lambdify(self.data.sym, T)
		for x in self.data.sets:
			self.L+=self.fd(*x)
		return self.L
        def get_Ls (self):
                self.L=0
                T=self.eq
                for parm in self.parms:
                        T=T.subs(parm.sym, parm.val)
                self.fd=lambdify(self.data.sym, T)
                for x in self.data.sets[0:100]:
                        try:
                                self.L+=self.fd(*x)
                        except:
                                self.L+=float("-infinity")
                return self.L

	def get_f (self):
		T=self.eq
		for parm in self.parms:
			T=T.subs(parm.sym, parm.val)
		self.fd=lambdify(self.data.sym, T)

	def f (self, x):
		return self.fd(x)
	
	def inc(self):
		self.S.inc()
		self.R=self.S.R
		ret=[]
		for p in self.parms:
			ret.append(str(p.name)+":"+str(p.val))
		return '\t'.join(ret)
	

#read in the data from sys.argv[1], 1 datum per line

data=d([["X", 0, float("infinity")]])
#, ["T", 0, float("infinity") ] ])
data.read(sys.argv[1])

#declare all the parameters we will use in ML equations

A=p("A", 0.00001, 1000 )
#B=p("B", 0.00001, 1000 )
#A2=p("A2", 0.01, max(data.data) )
#B2=p("B2", min(data.data), max(data.data) )
#C=p("C", 0.01, 0.99 )
#LNorm=Leq("-ln (sqrt(2*pi)*(X*S*sqrt(T) ) )-(ln(X)-ln(S0)-(M-0.5*S**2)*T)**2/(2*S**2*T)", parms, data, 0)


parms=[A]

#our ML equation. Should apper as python wants to see it, except in quotes.

m=mean(data, "X")
v=stdev(data, "X")**2

Neg_exp=Leq("loggamma(X+A)-loggamma(A)-loggamma(X+1)+A*ln(1-B)+X*ln(B)", parms, data, 0)

B.val=ln(m**2/sqrt(v+m**2) )
A.val=sqrt(ln(1+v/m**2) )
L=LNorm

while (L.R>0.0001):
	L.get_L()
	L.inc()
print "LogNorm: "+L.inc()+"\tLnL:"+str(L.get_L() )

B.val=m
A.val=sqrt(v)
L=Norm
while (L.R>0.0001):
	L.get_L()
	L.inc()
print "Norm: "+L.inc()+"\tLnL:"+str(L.get_L() )

A.val=m
L=Poisson
while (L.R>0.0001):
	L.get_L()
	L.inc()
print "Poisson: "+L.inc()+"\tLnL:"+str(L.get_L() )

L=Log_Log
B.val=1
A.val=100

while (L.R>0.0001):
	L.get_L()
	L.inc()
print "Log_Log: "+L.inc()+"\tLnL:"+str(L.get_L() )
