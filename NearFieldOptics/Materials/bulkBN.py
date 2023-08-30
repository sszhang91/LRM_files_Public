#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:44:08 2020

@author: samuelmoore


To use with lightning-rod model code:
    
import bulkBN as BN_params
import NearFieldOptics.Materials.TransferMatrixMedia as TMmat

SiO2 = Mat.SiO2_300nm
BNbulk = BN_params.LRMmat_default
BNbulk11 = BN_params.LRMmat_11BN
BNbulk10 = BN_params.LRMmat_10BN

stack1TM = TMmat.LayeredMedia((BNbulk11,100e-7),(SiO2,285e-7),exit = Mat.Si)

w = np.linspace(1380,1575,150); q = np.linspace(0.01,1.2,100)*1e5
stack1_rp = stack1TM.reflection_p(w,q=q)

"""
import sys
import numpy as np
#sys.path.insert(1,'/Users/samuelmoore/Desktop/Physics/BasovLab/Experiments/LRMc')
sys.path.insert(1,'C:\Data\Github\LRM_files')
from NearFieldOptics import Materials as Mat
from NearFieldOptics import TipModels as T

#relative to c-axis
global wTO_perp
wTO_perp = 1367.0 + 0 *1j
global wLO_perp
wLO_perp = 1610.0 +0*1j
global Gamma_perp
Gamma_perp = 29.0/10 + 0*1j

global wLO_par
wLO_par = 810.0 + 0*1j
global wTO_par
wTO_par = 746.0 + 0*1j
global Gamma_par
Gamma_par = 8.0/10 + 0*1j

global einf_perp
einf_perp = 4.95
global einf_par
einf_par = 4.1

def qdBN(w,wTOz,wTOt,wLOz,wLOt,gz,gt,epsInf_z,epsInf_t,esub,l):
    psi = -1j*np.sqrt(e_lorentz(w,wTOz,wLOz,gz,epsInf_z)/e_lorentz(w,wTOt,wLOt,gt,epsInf_t))
    airterm = np.arctan(1/psi/e_lorentz(w,wTOt,wLOt,gt,epsInf_t))
    subterm = np.arctan(esub/psi/e_lorentz(w,wTOt,wLOt,gt,epsInf_t))
    return -psi*(airterm+subterm+np.pi*l)

def qdBN_default(w,esub,l):
    return qdBN(w,wTO_par,wTO_perp,wLO_par,wLO_perp,Gamma_par,Gamma_perp,einf_par,einf_perp,esub,l)
    
def e_lorentz(w,wTO,wLO,Gamma,einf):
    return einf * (1 + (wLO**2-wTO**2)/(wTO**2-w**2-1j*Gamma*w) )

def e_lorentz_isotopic10(w):
    einf = 4.9; wLO = 1620; wTO = 1386; Gamma = 0.8;
    return einf * (1 + (wLO**2-wTO**2)/(wTO**2-w**2-1j*Gamma*w) )

#defined in terms of oscillator strength.  Useful for quasi-2D limit
def e_lorentz2D(w,wTO,osc,Gamma,einf):
    return einf + osc/(wTO**2-w**2-1j*Gamma*w)

def eps_t(w):
    return e_lorentz(w,wTO_perp,wLO_perp,Gamma_perp,einf_perp)
def eps_z(w):
    return e_lorentz(w,wTO_par,wLO_par,Gamma_par,einf_par)

def LRMmat(einf_perp1,einf_par1,wLO_perp1,wLO_par1,wTO_perp1,wTO_par1,Gamma_perp1,Gamma_par1):
    oscs_perp = [einf_perp1*(wLO_perp1**2-wTO_perp1**2)]; 
    oscs_par = [einf_par1*(wLO_par1**2-wTO_par1**2)]
    gs_perp = [Gamma_perp1]; gs_par = [Gamma_par1];
    ws_perp = [wTO_perp1]; ws_par = [wTO_par1]
    eps_lpsBN=[[(osc,w_TO,g) for w_TO,osc,g in zip(ws_perp,oscs_perp,gs_perp)]]*2
    eps_lpsBN+=[[(osc,w_TO,g) for w_TO,osc,g in zip(ws_par,oscs_par,gs_par)]]
    return Mat.AnisotropicMaterial(eps_infinity=[einf_perp1,einf_perp1,einf_par1],eps_lps=eps_lpsBN)

LRMmat_default = LRMmat(einf_perp,einf_par,wLO_perp,wLO_par,wTO_perp,wTO_par,Gamma_perp,Gamma_par)
LRMmat_11BN = LRMmat(einf_perp,einf_par,1600,802,1360,738,1,2)
LRMmat_10BN = LRMmat(einf_perp,einf_par,1630,835,1386,775,1,2)