from dftpy.ions import Ions
from dftpy.field import DirectField
from dftpy.grid import DirectGrid
from dftpy.functional import LocalPseudo, Functional, TotalFunctional
from dftpy.formats import io
from dftpy.optimization import Optimization
from dftpy.mpi import sprint
from dftpy.functional.pseudo.psp import PSP
from dftpy.constants import environ
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from ase.io import read

atoms = read('../../../KS/PBE/rho_ks_Li_gbrv.xsf')
ions, rho_ks, _ = io.read_all('../../../KS/PBE/rho_ks_Li_gbrv.xsf')
grid = rho_ks.grid
PP_list = {'Li': '/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/Li/li_pbe_v1.4.uspp.F.UPF',
           'O': '/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/Li/oxide/o_pbe_v1.2.uspp.F.UPF'}
# MaxPoints=1000
# PSEUDO = LocalPseudo(grid = grid, ions=ions, PP_list=PP_list)#, MaxPoints=MaxPoints)
rho_ini =rho_ks.copy()
PSEUDO = Functional(type='PSEUDO', grid=rho_ks.grid, ions=ions, PP_list=PP_list)
core = PSEUDO.core_density ## GBRV don't have NLCC
KE = Functional(type='KEDF',name='TFvW', y=0.2)
XC = Functional(type='XC',name='PBE', core_density=core)
HARTREE = Functional(type='HARTREE')
evaluator = TotalFunctional(KE=KE, XC=XC, HARTREE=HARTREE, PSEUDO=PSEUDO)

def delta_pp(r, rcut, a):
    d = r - rcut
    b = (3*a[0]*rcut-4*a[1]*rcut**2+5*a[2]*rcut**3)/2.0
    v = b*d**2 + a[0]*d**3 + a[1]*d**4+a[2]*d**5
    v[r>rcut] = 0.0
    return v

def lpp2vloc(r_Li, r_O, v_Li, v_O, ions, grid, zval=0.0):
    engine_Li = PSP(None)
    engine_Li.r = r_Li
    engine_Li.v = v_Li
    engine_Li._zval = zval

    engine_O = PSP(None)
    engine_O.r = r_O
    engine_O.v = v_O
    engine_O._zval = zval

    pseudo = LocalPseudo(grid = grid, ions=ions, PP_list={'Li': engine_Li, 'O':engine_O})#, MaxPoints=MaxPoints)
    pseudo.local_PP()
    return pseudo._vreal


# key='Li'
grid = rho_ks.grid
rcut_Li = 2.0 # Taken from the OEPP PP cutoff radius
r_Li = np.linspace(0, rcut_Li, 100)
rcut_O = 1.25
r_O = np.linspace(0, rcut_O, 100)
a = np.zeros(6)
KE = Functional(type='KEDF',name='TFvW', y=0.2)

evaluator = TotalFunctional(KE=KE, XC=XC, HARTREE=HARTREE, PSEUDO=PSEUDO)

ext = Functional(type='EXT')
evaluator.UpdateFunctional(newFuncDict={'EXT': ext})

opt = Optimization(EnergyEvaluator=evaluator)

rho_ini = rho_ks.copy()
environ['LOGLEVEL'] = 4
def delta_rho(a):
    a_O = a[:3]  
    a_Li = a[3:]  
    
    v_O = delta_pp(r_O, rcut_O, a_O)  
    v_Li = delta_pp(r_Li, rcut_Li, a_Li)  

    v = lpp2vloc(r_Li, r_O, v_Li, v_O, ions, grid)  
    # v = lpp2vloc(r_P, v_P, ions, grid)  

    ext.v = v

    # v = delta_pp(r, rcut, a)
    # ext.v = lpp2vloc(r_Li, v, ions, grid, key)
    rho = opt.optimize_rho(guess_rho=rho_ini)
    # rho_ini[:]=rho
    diff = 0.5 * (np.abs(rho - rho_ks)).integral()
    # if i%50==0 and i>1:
    print('aa:', a, diff)
    return diff

res = minimize(delta_rho, a, method='Powell', options={'ftol': 1.0e-4})
environ['LOGLEVEL'] = 2
a_Li = res.x[3:]
key = 'Li'
r = PSEUDO.readpp.pp[key].r
zval = PSEUDO.readpp.pp[key]._zval
vl = PSEUDO.readpp.pp[key].v
v = delta_pp(r, rcut_Li, a_Li)
v += vl

engine = PSP(None)
engine.r = r
engine.v = v
engine.info['atomicnum'] = 3
engine._zval = zval
engine.write('li_pbe_tf02vw.psp8')

a_O = res.x[:3]
key_O = 'O'
r_O = PSEUDO.readpp.pp[key_O].r
zval_O = PSEUDO.readpp.pp[key_O]._zval
vl_O = PSEUDO.readpp.pp[key_O].v
v_O = delta_pp(r_O, rcut_O, a_O)
v_O += vl_O

#v=vl

engine_O = PSP(None)
engine_O.r = r_O
engine_O.v = v_O
engine_O.info['atomicnum'] = 8
engine_O._zval = zval_O
engine_O.write('o_pbe_tf02vw.psp8')