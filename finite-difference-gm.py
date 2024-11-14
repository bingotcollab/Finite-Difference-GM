"""
@author: lucas
"""

import pathlib

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.special import lpmv, lambertw, factorial

import os, glob, sys
from pathlib import Path      #no headaches while working with files and paths
import argparse, json

class Metric:
    def __init__(self, metric_params: dict):
        self.alpha = metric_params['alpha']
    
    def turtle_coord(self,x):
        return x + 2*np.log(x-1)
    
    def inverse_turtle_coord(self, xtilde):
        '''uses a stirling aproximation for higher values. from x~500, the error (x*-x)/x
        is on the order 10^-7'''
        if type(xtilde) == np.ndarray:          
            lower = np.array(2*lambertw(np.exp((xtilde[xtilde<1000]-1)/2)/2)+1)
            u = np.array(xtilde[xtilde>=1000]/2 + np.log(2))
            higher = 2*(u-np.log(u)+np.log(u)/u) + 1
            return np.concatenate([lower,higher])
        else:
            return np.array(2*lambertw(np.exp((xtilde-1)/2)/2)+1)
        
    
    def ratio_fg(self,y,x):
        return ((1+self.alpha**2)*(x**2-y**2))**8/ \
        (((x**2)-(y**2) + (self.alpha**2)*(x-1)**2)**2 \
            - 4*(self.alpha**2)*(y**2)*(x**2-1))**4
        #return 1 + self.alpha**2 * y**2/(1+x**2)
    
    def g2(self,y,x):
        #correct one 
        return ((x**2 - y**2 + self.alpha**2 * (x**2 - 1))**2\
                + 4* (self.alpha**2) * (x**2) * (1 - y**2))**4\
            /((self.alpha**2 + 1)**8 * (x**2 - y**2)**8)
        #return 1 + self.alpha**2 *(1- y**2)*y**2/(1+x**2)
    
    def f2(self,y,x):
        return ((x**2 - y**2 + self.alpha**2 * (x**2 - 1))**2\
                + 4* (self.alpha**2) * (x**2) * (1 - y**2))**2/ \
                (((x**2)-(y**2) + (self.alpha**2)*(x-1)**2)**2 \
                - 4*(self.alpha**2)*(y**2)*(x**2-1))**2
    
    def A3(self,y,x):
        return ((4*(1-y**2)*self.alpha**3)*(2*(1+self.alpha**2)*x**3 + (1 - 3*self.alpha**2)*x**2\
                                    +y**2 +self.alpha**2))/ \
        ((1+self.alpha**2)*(((x**2)-(y**2) + (self.alpha**2)*(x**2-1))**2 \
        - 4*(self.alpha**2)*(y**2)*(x**2-1)))
    
    def b_integrand(self,y,x):
        return (self.g2(y, x) - 1)*self.ratio_fg(y, x)/(1-y**2)
    
    def c_im_integrand(self,y,x):
        return ((self.g2(y, x) * self.f2(y, x) * (x-1))/((x+1) * (1-y**2)))\
            *self.A3(y, x)
    
    def c_re_intengrand(self,y,x):
        return ((self.g2(y, x) * self.f2(y, x) * (x-1))/((x+1) * (1-y**2)))\
            *self.A3(y, x)*self.A3(y,x)
    
    
    

class Solver:
    def __init__(self, metric, lmax, Npoints_x=10000, Npoints_t=20000, xmin=-100, xmax=100,
                 t_step=10/10000, m=0, initial_value=None, li=0,
                 initial_velocity=None, initial_spread=2, initial_distance=0,
                 block_size=2000, observing_points=5, save_full=False, q=0):
        #initialize constant atributes
        self.metric = metric
        self.lmin = m
        self.lmax = lmax
        self.Lmodes = self.lmax - self.lmin
        self.Npoints_x = Npoints_x
        self.Npoints_t = Npoints_t
        self.xmin = xmin
        self.xmax = xmax
        self.x_step = (xmax - xmin)/Npoints_x
        self.m = m
        if t_step >= self.x_step:
            print("WARNING: time step is greater than the spatial step. \
                  Von Neumann conditions may not be satisfied")
        self.t_step = t_step
        self.block_size = block_size
        self.save_full = save_full
        self.initial_time = 0
        self.q = q

        
        #populate integration grid
        self.Xtilde, self.dXtilde = np.linspace(self.xmin, self.xmax,
                                                self.Npoints_x, retstep=True)
        self.X = self.metric.inverse_turtle_coord(self.Xtilde)
        self.T, self.dT = np.linspace(0, self.Npoints_t*self.t_step,
                                      self.Npoints_t, retstep=True)
        
        #set up the points in which to observe
        highest_x = Npoints_x*(initial_distance - xmin)//(xmax - xmin)
        self.observed_points = np.append(np.linspace(0, Npoints_x, observing_points - 1, endpoint=False), highest_x)
        self.observed_points.sort()
        self.observed_points = self.observed_points.astype(int)
        
        #initialize A, B and b matrices
        self.A = np.zeros((self.Lmodes,self.Lmodes,self.Npoints_x), dtype=float)
        self.b = np.zeros((self.Lmodes,self.Lmodes,self.Npoints_x), dtype=float)
        self.B = np.zeros((self.Lmodes,self.Lmodes,self.Npoints_x), dtype=float)
        self.C = np.zeros((self.Lmodes,self.Lmodes,self.Npoints_x), dtype=float)
        
        #populating these matrices
        for l in range(self.lmin, self.lmax):
            norm_l = self.Norm(l)              #uglier to put here but saves time
            for k in range(self.lmin, self.lmax):
                norm_k = self.Norm(k)
                for i, x in enumerate(self.X):
                    self.A[l-self.lmin,k-self.lmin,i] = norm_l*norm_k*quad(self.a_prod,-1,1,args=(x,l,k))[0]
                    if self.A[l-self.lmin,k-self.lmin,i] < 1e-12:
                        self.A[l-self.lmin,k-self.lmin,i] = 0
                    self.b[l-self.lmin,k-self.lmin,i] = norm_l*norm_k*quad(self.b_prod,-1,1,args=(x,l,k))[0]
                    if self.b[l-self.lmin,k-self.lmin,i] < 1e-12:
                        self.b[l-self.lmin,k-self.lmin,i] = 0
                    self.B[l-self.lmin,k-self.lmin,i] = ((x - 1)/(x + 1)**3) \
                        *(self.A[l-self.lmin,k-self.lmin,i]*self.Veff(l,x) + (self.m**2)*self.b[l-self.lmin,k-self.lmin,i])
#                    self.C[l-self.lmin,k-self.lmin,i] = norm_l*norm_k*quad(self.c_prod,-1,1,args=(x,l,k))[0]

        
        
        
        #initial conditions:
        self.Psi = np.zeros((self.Lmodes,self.Npoints_x,self.block_size))
        
        
        if initial_value:
            self.Psi[:,:,0] = np.array(initial_value, dtype=float)
        else:
            self.Psi[:,:,0] = np.zeros((self.Lmodes, self.Npoints_x))
            self.Psi[0,:,0] += np.ones(self.Npoints_x)*np.exp(-(self.Xtilde-initial_distance)**2/initial_spread)
            
            # experimental
            # self.Psi[li,:,0] -= 0.3*np.ones(self.Npoints_x)*np.exp(-(self.Xtilde-initial_distance)**2/(initial_spread*4))
            
            
        if initial_velocity:
            self.initial_velocity = np.array(initial_velocity, dtype=float)
        else:
            self.initial_velocity = np.zeros((self.Lmodes, self.Npoints_x))

        if initial_value:
            #NOT IMPLEMENTED
            pass
        else:
            f_dobleline = np.zeros((self.Lmodes,self.Npoints_x))
            f_dobleline[0,:] = (2/initial_spread)*((2/initial_spread)*(self.Xtilde-initial_distance) - 1) \
                *np.exp(-(self.Xtilde-initial_distance)**2/initial_spread)
            
            #experimental
            # f_dobleline[li,:] -= 0.3*(0.5/initial_spread)*((0.5/initial_spread)*(self.Xtilde-initial_distance) - 1) \
            #     *np.exp(-(self.Xtilde-initial_distance)**2/(initial_spread*4))
            
            self.Psi[:,:,1] += np.einsum('ijk,jk->ik'
                    , (np.tensordot(np.eye(self.Lmodes),np.ones(self.Npoints_x),0) \
                        - (self.t_step**2/2)*self.B),
                    self.Psi[:,:,0]) #contracts the matrix B with F for each value of x
            self.Psi[:,:,1] += self.t_step*self.initial_velocity
            self.Psi[:,:,1] += (self.t_step**2/2)*np.einsum('ijk,jk->ik', self.A,
                    f_dobleline)
            
    def times(self,x_point):
        tmax = (self.Npoints_x - x_point) if x_point >= self.Npoints_x//2 else x_point
        return np.linspace(0, tmax*self.t_step,tmax)

                    
    def Norm(self,l):
       	return np.sqrt(((2*l+1)*factorial(l-self.m))/(2*factorial(l+self.m)))
    
    def a_prod(self,y,x,l,k):
        return self.metric.ratio_fg(y,x)*lpmv(self.m,l,y)*lpmv(self.m,k,y)
        
    def b_prod(self,y,x,l,k):
        return self.metric.b_integrand(y,x)*lpmv(self.m,l,y)*lpmv(self.m,k,y)
    
    def c_prod(self,y,x,l,k):
        return ((self.q**2 * self.metric.c_re_intengrand(y, x))\
                - (2*self.q * self.m * self.metric.c_im_integrand(y,x)))\
            *lpmv(self.m,l,y)*lpmv(self.m,k,y)
    
    def Veff(self,l,x):
        return np.array(l*(l+1) + 2/(x+1))
    
    def step(self,i_t):
        #implement 4th order, implement borders
        self.Psi[:,:,i_t] = -self.Psi[:,:,i_t-2] + 2*self.Psi[:,:,i_t - 1]
        self.Psi[:,:,i_t] += -self.dT**2 * np.einsum('ijk,jk->ik', \
                                            self.B, self.Psi[:,:,i_t - 1])
        
        D2Psi = np.zeros((self.Lmodes, self.Npoints_x))
        # D2Psi[:,0] = -0 +  16*0 - 30*self.Psi[:,0,i_t-1] + 16*self.Psi[:,1,i_t-1] \
        #     - self.Psi[:,2,i_t-1]
        D2Psi[:,1] = -0 +  16*self.Psi[:,0,i_t-1] - 30*self.Psi[:,1,i_t-1] \
            + 16*self.Psi[:,2,i_t-1] - self.Psi[:,3,i_t-1]
        D2Psi[:,-1] = -self.Psi[:,-3,i_t-1] +  16*self.Psi[:,-2,i_t-1] \
            - 30*self.Psi[:,-1,i_t-1] + 16*0 - 0
        D2Psi[:,-2] = -self.Psi[:,-4,i_t-1] +  16*self.Psi[:,-3,i_t-1] \
            - 30*self.Psi[:,-2,i_t-1] + 16*self.Psi[:,-1,i_t-1] - 0
        
        
        for i_x in range(2, self.Npoints_x-2):
            D2Psi[:,i_x] = -self.Psi[:,i_x - 2,i_t-1] +  16*self.Psi[:,i_x - 1,i_t-1] \
                - 30*self.Psi[:,i_x,i_t-1] \
                + 16*self.Psi[:,i_x + 1,i_t-1] - self.Psi[:,i_x + 2,i_t-1]
        self.Psi[:,:,i_t] += (self.dT**2 / self.dXtilde**2) \
            * np.einsum('ijk,jk->ik', self.A, D2Psi)/12
    
    def solve(self, init_t=0):
        for i_t in range(init_t*self.block_size+2,self.Npoints_t):
            if i_t % self.block_size == 0:
                self._save_data(i_t)
            self.step(i_t%self.block_size)  
        self._save_data(self.Npoints_t-1)
            
    def _save_data(self, i_t):
        print("saving at step {}".format(i_t))
        
        #saving folder
        sf = Path(f'./results/alpha_{self.metric.alpha:.3e}_xmin_{self.xmin}')   #not ideal, but not worth finding a solution atm
        sf.mkdir(exist_ok=True, parents=True)
        
        #save the whole vector
        init_time = self.T[int(max(i_t-self.block_size,0))]
        end_time = self.T[i_t]
        if self.save_full:
            file = sf / f'Psi_from_t_{init_time:04.0f}_to_{end_time:04.0f}'
            np.savez_compressed(file, self.Psi)
        
        #save checkpoints
        saving_vector = np.stack(tuple(self.Psi[:,obs_pt,:] for obs_pt in self.observed_points)
                                 , axis=1)
        file = sf / f'Psi_from_t_{init_time:04.0f}_to_{end_time:04.0f}_selected_points'
        np.savez_compressed(file, saving_vector)
        file = sf / 'observation_points'
        saving_points = np.array([self.X[obs_pt] for obs_pt in self.observed_points])
        np.savez_compressed(file, saving_points)
        
        #save last two points
        file = sf / f'last_two_times_{i_t//self.block_size:04n}'
        np.savez_compressed(file, self.Psi[:,:,-2:])
        
        
    def _retrieve_last_data(self):
        file = sorted(Path(f'./results/alpha_{self.metric.alpha:.3e}').glob('last_two_times_*.npz'))[-1]
        loaded_initial_condition = np.load(file).get('arr_0')     #np savez_compressed naming convention
        self.initial_time = int(file.name.split('_')[-1].split('.')[0])
        self.Psi[:,:,:2] = loaded_initial_condition

    def load_observations(self):
        files = files = sorted(Path(f'./results/alpha_{self.metric.alpha:.3e}').glob('Psi_from_*_to_*_selected_points.npz'))
        return np.concatenate(tuple(np.load(file).get('arr_0') for file in files),axis=2)

    def energy(self,i_t):
        pass

    
def parse_arguments():
    argpars = argparse.ArgumentParser(prog='Finite-difference scheme for Gutsanaev-Manko-like metrics',
                                      description='Takes in an axisymmetric metric and calculates a scalar perturbations')
    argpars.add_argument('--params', help='parameters file path', type=pathlib.Path, default=pathlib.Path('./parameters.json'), required=False)
    argpars.add_argument('-l', '--load', action='store_true',help='load last simulated value')
    return vars(argpars.parse_args())

def load_params(params_file):
    with params_file.open(encoding='UTF-8') as pf:
        params = json.load(pf)
    return params


if __name__ == "__main__":
    args = parse_arguments()
    params = load_params(args['params'])
    # ok, so here comes some trouble. A better practice would be to pass agnostically the parameters of the metric,
    # and the metric class should handle internally. The way things are now, it's quite prone to error, one should
    # always pay attention if reloading from file if the parameters are the same

    met = Metric(params['metric'])
    s = Solver(met, params['lmax'], m=params['m'], Npoints_x=params['Npoints_x'], t_step=params['t_step'],
               xmin=params['xmin'], xmax=params['xmax'], initial_distance=params['initial_distance'],
               Npoints_t=params['Npoints_t'], li=params['li'], initial_spread=params['initial_spread'],
               block_size=params['block_size'], save_full=params['save_full'],
               observing_points=params['observing_points'],q=params['charge'])
    if not args['load']:
        s.solve()
    else:
        s._retrieve_last_data()
        s.solve(init_t=s.initial_time)
