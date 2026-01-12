# Li

This repository contains scripts to generate Local pseudopotentials (LPPs) using DFTpy for Lithium in metal (Li BCC) and oxide (Li2O). 
The organization of the repository is the following:

├── metal
│   ├── KS
│   │   ├── LDA
│   │   │   ├── QEpy_RHO_Li.ipynb
│   │   │   ├── rho_ks_Li_gbrv.xsf
│   │   │   └── rho_ks_Li_oepp.xsf
│   │   └── PBE
│   │       ├── QEpy_RHO_Li.ipynb
│   │       ├── rho_ks_Li_gbrv.xsf
│   │       └── rho_ks_Li_oepp.xsf
│   └── OF
│       ├── LDA
│       │   ├── GBRV
│       │   │   ├── revHC
│       │   │   │   ├── PGBRV_LDA_revHC.psp8
│       │   │   │   └── pp.py
│       │   │   └── TF02vW
│       │   │       ├── eos.ipynb
│       │   │       ├── PGBRV_LDA_TF02vW.psp8
│       │   │       └── pp.ipynb
│       │   └── OEPP
│       │       ├── revHC
│       │       │   ├── eos.ipynb
│       │       │   ├── PGBRV_LDA_revHC.psp8
│       │       │   └── pp.py
│       │       └── TF02vW
│       │           ├── eos.ipynb
│       │           ├── POEPP_LDA_TF02vW.psp8
│       │           └── pp.ipynb
│       └── PBE
│           ├── GBRV
│           │   ├── revHC
│           │   │   └── pbe_revhc.py
│           │   └── TF02vW
│           │       ├── eos.ipynb
│           │       ├── PGBRV_PBE_TF02vW.psp8
│           │       └── pp.ipynb
│           └── OEPP
│               ├── revHC
│               │   ├── eos.ipynb
│               │   ├── POEPP_PBE_revHC.psp8
│               │   └── pp.py
│               └── TF02vW
│                   ├── eos.ipynb
│                   ├── POEPP_PBE_TF02vW.psp8
│                   └── pp.ipynb
├── oxide
│   ├── KS
│   │   ├── LDA
│   │   │   ├── QEpy_RHO_Li.ipynb
│   │   │   └── rho_ks_Li_gbrv.xsf
│   │   └── PBE
│   │       ├── QEpy_RHO_Li.ipynb
│   │       └── rho_ks_Li_gbrv.xsf
│   └── OF
│       ├── LDA
│       │   └── GBRV
│       │       └── TF02vW
│       │           ├── eos_li2o_gbrv.ipynb
│       │           ├── eos_Li.ipynb
│       │           ├── eos_Li-O.ipynb
│       │           ├── Li_Li2O_PGBRV_LDA_TF02vW.psp8
│       │           ├── Li_PGBRV_LDA_TF02vW.psp8
│       │           ├── O_Li2O_PGBRV_LDA_TF02Vw.psp8
│       │           ├── pp_Li.ipynb
│       │           └── pp_Li-O.py
│       └── PBE
│           └── GBRV
│               └── TF02vW
│                   ├── eos_li2o_gbrv.ipynb
│                   ├── eos_Li.ipynb
│                   ├── eos_Li-O.ipynb
│                   ├── Li_Li2O_PGBRV_PBE_TF02vW.psp8
│                   ├── Li_PGBRV_PBE_TF02vW.psp8
│                   ├── O_Li2O_PGBRV_PBE_TF02vW.psp8
│                   ├── pp-Li.ipynb
│                   └── pp_Li-O.py


