import matplotlib
matplotlib.use('Agg')
import pybinding as pb
import numpy as np
import matplotlib.pyplot as plt
import pickle
import argparse
import sys

from matplotlib import pylab, mlab 

from pylab import *
from numpy import *

# Setting up argument parsing
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--out'   , dest='out'   , default=False             ,help='Name of the output file')
parser.add_argument('--dim'   , dest='dim'   , default=10   ,type=int    ,help='Width of the sample')
parser.add_argument('--theta' , dest='theta0', default=0    ,type=float  ,help='Angle theta of the field orientation')
parser.add_argument('--B'     , dest='B'     , default=0    ,type=float  ,help='Magnetic field strength')
parser.add_argument('--Eran'  , dest='Eran'  , default=1    ,type=float  ,help='Energy window for DOS is [-Eran,Eran]')
parser.add_argument('--broad' , dest='broad' , default=0.002,type=float  ,help='Broadening of the calculation')

args = parser.parse_args()
theta=args.theta0/180*pi # going radians

# This is a helper function responsible for Peierls substitution
# that is implementing a vector potential corresponding to a 
# magnetic field
def constant_magnetic_field(B,theta):
    @pb.hopping_energy_modifier
    def function(energy, x1, y1,z1, x2, y2,z2):
        x = 0.5 * (x1 + x2)
        y = 0.5 * (y1 + y2)
        z = 0.5 * (z1 + z2)
        
        A_y = B * (x*cos(theta)-z*sin(theta))
        
        peierls = A_y * (y1 - y2)
        return energy * np.exp(1j * 2*pi * peierls)
    return function

# This builds a model of a nodal loop semimetal on a cubic lattice
def mymodel(dim=20,theta=0,B=0.0):
    # Pauli matrices
    S1=array([[0,1],[1,0]])
    S2=array([[0,-1j],[1j,0]])
    S3=array([[1,0],[0,-1]])   
    Tx=S3       # Hopping in the x direction
    Ty=S3       # Hopping in the y direction
    Tz=S3+1j*S1 # Hopping in the z direction
    U=-5*S3     # Onsite term
    # Building the lattice model
    lat=pb.Lattice(a1=[1,0,0],a2=[0,1,0],a3=[0,0,1])
    lat.add_sublattices(('A', [0,0,0],U ))
    lat.add_hoppings(                 
                 ([1,0,0], 'A', 'A', Tx),
                 ([0,1,0], 'A', 'A', Ty),
                 ([0,0,1], 'A', 'A', Tz)
                )
    # generating the model object by combinig the lattice and the Peierls field 
    model = pb.Model(lat,pb.primitive(a1=dim, a2=dim,a3=dim),constant_magnetic_field(B=B,theta=theta))
    return model

# Calculate DOS
# We use jackson kernel here that is good for DOS
kpm   = pb.kpm(mymodel(args.dim,theta,args.B),silent=True,
               energy_range=[-13,13],kernel=pb.jackson_kernel()) 
# we approximate the bulk dos with the local DOS at the center of the lattice
# the calculation is performed at this step
dos=kpm.calc_ldos(energy=linspace(-args.Eran,args.Eran,3000),
                  position=[0,0,0],
                  broadening=args.broad)

# build up the name of the output file
if args.out==False:    
    out='-'.join(['dos' ,
                'dim',str(args.dim),
              'theta',str(args.theta0),
                'B'  ,str(args.B)])

# Dump data and meta data in output file
with open(out, 'wb') as handle:
    pickle.dump(args , handle, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(dos  , handle, protocol=pickle.HIGHEST_PROTOCOL)
