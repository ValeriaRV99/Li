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
atoms = read('/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/Li/metal/KS/LDA/rho_ks_Li_oepp.xsf')
ions, rho_ks, _ = io.read_all('/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/Li/metal/KS/LDA/rho_ks_Li_oepp.xsf')
grid = rho_ks.grid
PP_list = {'Li': '/projectsn/mp1009_1/Valeria/Batteries/Li2S_interface/OF/PP/Li_OEPP_PZ.UPF'}
# MaxPoints=1000
PSEUDO = LocalPseudo(grid = grid, ions=ions, PP_list=PP_list)#, MaxPoints=MaxPoints)
rho_ini =rho_ks.copy()
core = PSEUDO.core_density
KE = Functional(type='KEDF',name='TFvW', y=0.2)
XC = Functional(type='XC',name='LDA', core_density=core)
HARTREE = Functional(type='HARTREE')

evaluator = TotalFunctional(KE=KE, XC=XC, HARTREE=HARTREE, PSEUDO=PSEUDO)
optimization_options = {'econv' : 1e-11*ions.nat}
opt = Optimization(EnergyEvaluator=evaluator, optimization_options = optimization_options, optimization_method = 'CG')
rho = opt.optimize_rho(guess_rho=rho_ks)
delta = 0.5 * np.abs(rho_ks - rho).integral()
print('Drho_TF0.2vW', delta)

kfmax=( np.max(rho) * 3 * np.pi )**(1.0/3.0)+0.2
from dftpy.mixer.pulay import PulayMixer
from dftpy.optimization import OESCF
pulay=PulayMixer(mp=None)
opt_options = {'econv' : 1e-6}
vW = Functional(type='KEDF', name='vW')
evaluator = TotalFunctional(KE=vW, XC=XC, HARTREE=HARTREE, PSEUDO=PSEUDO)
kedf_emb = Functional(type='KEDF',name='LMGP')#, kfmax=kfmax, kfmin=1e-10, ratio=1.02)
kedf_emb.options.update({'y':0})
evaluator_emb = TotalFunctional(KEDF_EMB = kedf_emb)
opt = Optimization(EnergyEvaluator=evaluator, optimization_options = opt_options, optimization_method = 'CG')
opt = OESCF(optimization=opt, evaluator_emb=evaluator_emb, guess_rho=rho,mixer=pulay)
rho_lmgp = opt.optimize_rho(guess_rho=rho_ks, econv=1e-6)
delta = 0.5 * np.abs(rho_ks - rho_lmgp).integral()
print('Drho_LMGP', delta)

KE = Functional(type='KEDF',name='revHC')
evaluator = TotalFunctional(KE=KE, XC=XC, HARTREE=HARTREE, PSEUDO=PSEUDO)
optimization_options = {'econv' : 1e-11*ions.nat}
opt = Optimization(EnergyEvaluator=evaluator, optimization_options = optimization_options, optimization_method = 'CG')
rho_revhc = opt.optimize_rho(guess_rho=rho_ks)
delta = 0.5 * np.abs(rho_ks - rho_revhc).integral()
print('Drho_revHC', delta)

fig, axs = plt.subplots(1, 5, figsize=(10, 3))
r = np.linspace(0,ions.cell[2][2],len(rho_ks[0,0,:].ravel()))
cut = 0
axs[0].plot(rho_ks[:,cut,cut].ravel(), label='KS LDA GBRV')
axs[0].plot(rho[:,cut,cut].ravel(), label='OF LDA TF0.2vW')
axs[0].plot((rho_lmgp[:,cut,cut]).ravel(), label='OF LMGP')
axs[0].plot((rho_revhc[:,cut,cut]).ravel(), label='OF revHC')
axs[0].legend()
axs[0].set_ylabel('n(r)')
axs[0].set_xlabel('r (au)')
axs[1].matshow(rho_ks[cut,:,:])
axs[1].set_title('KS')
axs[2].matshow(rho[cut,:,:])
axs[2].set_title('OF TF0.2vW')
axs[3].matshow(rho_lmgp[cut,:,:])
axs[3].set_title('OF LMGP')
axs[4].matshow(rho_revhc[cut,:,:])
axs[4].set_title('OF revHC')
plt.tight_layout()
fig.savefig('densities.png', dpi=500)

rho.write('rho_TF02vW_Li.xsf', ions=ions)
rho_lmgp.write('rho_LMGP_Li.xsf', ions=ions)
rho_revhc.write('rho_revHC_Li.xsf', ions=ions)

def delta_pp(r, rcut, a):
    d = r - rcut
    b = (3*a[0]*rcut-4*a[1]*rcut**2+5*a[2]*rcut**3)/2.0
    v = b*d**2 + a[0]*d**3 + a[1]*d**4+a[2]*d**5
    v[r>rcut] = 0.0
    return v

def lpp2vloc(r, v, ions, grid, key, zval=0.0):
    engine = PSP(None)
    engine.r = r
    engine.v = v
    engine._zval = zval
    pseudo = LocalPseudo(grid = grid, ions=ions, PP_list={key:engine})#, MaxPoints=MaxPoints)
    pseudo.local_PP()
    return pseudo._vreal

key='Li'
grid = rho_ks.grid
rcut = 2.0 # Taken from the GBRV PP cutoff radius
r = np.linspace(0, rcut, 100)
a = np.zeros(3)
KE = Functional(type='KEDF',name='revHC')
evaluator = TotalFunctional(KE=KE, XC=XC, HARTREE=HARTREE, PSEUDO=PSEUDO)
ext = Functional(type='EXT')
evaluator.UpdateFunctional(newFuncDict={'EXT': ext})
optimization_options = {'econv' : 1e-11*ions.nat}
opt = Optimization(EnergyEvaluator=evaluator,optimization_options = optimization_options, optimization_method = 'CG')

rho_ini = rho_ks.copy()
environ['LOGLEVEL'] = 4
def delta_rho(a):
    print('Init')
    v = delta_pp(r, rcut, a)
    ext.v = lpp2vloc(r, v, ions, grid, key)
    rho = opt.optimize_rho(guess_rho=rho_ini)
    # rho_ini[:]=rho
    diff = 0.5 * (np.abs(rho - rho_ks)).integral()
    # if i%50==0 and i>1:
    print('aa:', a, diff)
    return diff
res = minimize(delta_rho, a, method='Powell', options={'ftol': 1.0e-4})
environ['LOGLEVEL'] = 2
a = res.x
key = 'Li'
r = PSEUDO.readpp.pp[key].r
zval = PSEUDO.readpp.pp[key]._zval
vl = PSEUDO.readpp.pp[key].v
v = delta_pp(r, rcut, a)
v += vl

engine = PSP(None)
engine.r = r
engine.v = v
engine.info['atomicnum'] = 3
engine._zval = zval
np.save('a', a)
engine.write('PGBRV_LDA_revHC.psp8')
