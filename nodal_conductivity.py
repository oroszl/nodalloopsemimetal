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
#parser.add_argument('--out'   , dest='out'   , default='sigma.dat'       ,help='Name of the output file')
parser.add_argument('--sigma' , dest='sigma' , default='zz'              ,help='Direction pair of the conductivity calculated')
parser.add_argument('--dim'   , dest='dim'   , default=10   ,type=int    ,help='Width of the sample')
parser.add_argument('--theta' , dest='theta0', default=0    ,type=float  ,help='Angle theta of the node orientation')
parser.add_argument('--phi'   , dest='phi0'  , default=0    ,type=float  ,help='Angle phi of the node orientation')
parser.add_argument('--B'     , dest='B'     , default=0    ,type=float  ,help='Magnetic field strength, fixed in the z direction')
parser.add_argument('--numrnd', dest='numrnd', default=1    ,type=int    ,help='Number of random vectors for the conductivity calculation')
parser.add_argument('--broad' , dest='broad' , default=0.1  ,type=float  ,help='Broadening of the conductivity calculation')


args = parser.parse_args()
theta=args.theta0/180*pi # going radians
phi=args.phi0/180*pi     # instead of degs


# This is a helper function responsible for Peierls substitution
# that is implementing a vector potential corresponding to a 
# magnetic field
def constant_magnetic_field(B):
    @pb.hopping_energy_modifier
    def function(energy, x1, y1, x2, y2):
        y = 0.5 * (y1 + y2)
        A_x = B * y
        peierls = A_x * (x1 - x2)
        return energy * np.exp(1j * 2*pi * peierls)
    return function
# This builds a model of a nodal loop semimetal on a cubic lattice
# the orientation of the loop is governed via spherical angles theta and phi
def mymodel(dim=20,phi=0,theta=0,B=0.0):
    S0=array([[1,0],[0,1]])
    S1=array([[0,1],[1,0]])
    S2=array([[0,-1j],[1j,0]])
    S3=array([[1,0],[0,-1]])
    N=zeros_like(S0)    
    Tx=-S3+1j*S1*sin(theta)*cos(phi)
    Ty=-S3+1j*S1*sin(theta)*sin(phi)
    Tz=-S3+1j*S1*cos(theta)
    U=-5*S3        
    lat=pb.Lattice(a1=[1,0,0],a2=[0,1,0],a3=[0,0,1])
    lat.add_sublattices(('A', [0,0,0],U ))
    lat.add_hoppings(                 
                 ([1,0,0], 'A', 'A', Tx),
                 ([0,1,0], 'A', 'A', Ty),
                 ([0,0,1], 'A', 'A', Tz)
                )
    model = pb.Model(lat,pb.primitive(a1=dim, a2=dim,a3=dim),constant_magnetic_field(B=B))
    return model

# Calculate DOS
kpm   = pb.kpm(mymodel(args.dim,phi,theta,args.B),silent=True,
               energy_range=[-13,13],kernel=pb.jackson_kernel()) # We use jackson kernel here that is good for DOS
dos=kpm.calc_ldos(energy=linspace(-15,15,1000),position=[0,0,0],broadening=0.05)

# Calculate conductivity
kpm   = pb.kpm(mymodel(args.dim,phi,theta,args.B),silent=True,
               energy_range=[-13,13],kernel=pb.lorentz_kernel()) # We use Lorentz kernel here that is needed for the conductivity
sigma = kpm.calc_conductivity(chemical_potential=linspace(-2, 2, 100),
                              broadening=args.broad,temperature=0,direction=args.sigma,
                              volume=args.dim**3,num_random=args.numrnd)
# make output file
out='-'.join(['sigma',args.sigma,
                'dim',str(args.dim),
              'theta',str(args.theta0),
                'phi',str(args.phi0),
                'B'  ,str(args.B)])

# Dump data and meta data in output file
with open(out, 'wb') as handle:
    pickle.dump(args , handle, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(dos  , handle, protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(sigma, handle, protocol=pickle.HIGHEST_PROTOCOL)
