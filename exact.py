# Simple python script to calculate energy specturm 
# and density of states  for a nodal loop semimetal 
# subject to magnetic field.

import numpy
from numpy import *
from numpy.linalg import *
import pickle
import argparse
import sys
from scipy import sqrt
import scipy.sparse as ss
import scipy.sparse.linalg as sl
# Some input parsing
parser = argparse.ArgumentParser()

parser.add_argument('--out'   , dest='out'   , default='DOS.dat'               ,help='Name of the output file')
parser.add_argument('--B'     , dest='B'     , default=0.1                     ,type=float,help='Magnetic field')
parser.add_argument('--Eran'  , dest='Eran'  , default='linspace(0.01,1.2,400)',help='Energy range')
parser.add_argument('--dim'   , dest='dim'   , default=600,type=int,help='discretization dimension')
parser.add_argument('--L'     , dest='L'     , default=10.0,type=float,help='discretization length')
args = parser.parse_args()


dim=args.dim
# real space discretized position operator
z=linspace(-args.L,args.L,dim) 
dz=(max(z)-min(z))/dim
# Pauli matrices
sx=array([[0,1],
         [1,0]])
sz=array([[1,0],
         [0,-1]])
sy=array([[0,-1j],
         [1j,0]])
# real space discretized momentum operator
ps=1j*(ss.diags([ones(dim-1)],[1])-ss.diags([ones(dim-1)],[-1]))/(2*dz)

ons=ones_like(z) 
# Helper function to generate eigen problem for px^2

makeEB=(lambda ee,BB: (
    -ss.kron(ss.diags([z**2*BB**(2/3)],[0]),eye(2))
    +ss.kron(ss.diags([ons],[0]),eye(2))
    -ee*ss.kron(ss.diags([ons],[0]),sx)
    -BB**(2/3)*1j*ss.kron(ps,sy)
    ).tocsc() )
# current operator is proportional to \sigma_x
curr=ss.kron(ss.diags([ons],[0]),sx)
# Calculationg the spectrum in an energy intervall
datE=[]
for eee in eval(args.Eran):
    # the eigen problem
    va,ve=eig(makeEB(eee,args.B).toarray())
    # groupvelocities
    vg=[]
    for i in range(len(va)):
        vvee=ve[:,i]
        vvaa=va[i]     
        vg.append(dot(conj(vvee),-2*sqrt(vvaa)*curr*vvee))
    vg=array(vg)
    # Finding relevant (real) solutions
    idx=(abs(imag(va))<1e-10)*(real(va)>0)
    # storing the found momentums group velocities and totl DOS at a given energy
    datE.append(dict(E=eee,va=va[idx],vg=vg[idx],dos=sum(1/abs(vg[idx]))))
# Writing output
with open(args.out, 'wb') as handle:
        pickle.dump(datE, handle, protocol=pickle.HIGHEST_PROTOCOL)
