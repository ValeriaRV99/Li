from __future__ import annotations
import numpy as np

import os
import math
from dataclasses import dataclass
from typing import Optional, Dict, Any
import numpy as np
import xmltodict  # pip install xmltodict

def psp8_to_upf(psp8_path: str, out_path: str, element_label: Optional[str] = None,) -> None:
    """Convert PSP8 format to UPF format.
    
    Args:
        psp8_path: Path to input PSP8 file
        out_path: Path to output UPF file
        element_label: Optional element label (if None, will try to determine from file)
    """
    psp8_data = read_psp8(psp8_path)
    r = psp8_data["r"]
    vloc = psp8_data["vloc"]
    z_valence = psp8_data["z_valence"]
    z_atom = psp8_data["z_atom"]
    element = element_label or psp8_data["element"]
    rho_core = psp8_data["rho_core"]
    
    # Compute rab (derivative of r)
    rab = np.gradient(r)
    
    # Guess dx and xmin for mesh
    zmesh = z_atom
    dx, xmin = guess_dx_xmin(r, zmesh)
    
    # Determine XC functional from pspxc code
    xc_map = {1: "PZ", 2: "LDA", 3: "PBEsol", 4: "LDA"}  # Common mappings
    xc_code = psp8_data["header"].get("pspxc", 1)
    xc_func = xc_map.get(xc_code, "LDA")
    
    # Create header attributes
    header_attrs = {
        "element": element,
        "z_valence": f"{z_valence:.8f}",
        "pseudo_type": "NC",  # Norm-conserving
        "relativistic": "scalar",
        "is_ultrasoft": "false",
        "is_paw": "false",
        "core_correction": "true" if rho_core is not None else "false",
        "functional": xc_func,
        "mesh_size": str(r.size),
        "number_of_wfc": "0",
        "number_of_proj": "0",
        "has_so": "false",
    }
    
    # Write UPF v2.0.1
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w") as f:
        f.write('<UPF version="2.0.1">\n')
        
        # PP_INFO
        f.write("<PP_INFO>\n")
        f.write(f"  Converted from PSP8 format: {os.path.basename(psp8_path)}\n")
        f.write(f"  Element: {element}, Z={z_atom:.1f}, Z_valence={z_valence:.1f}\n")
        f.write(f"  XC functional: {xc_func}\n")
        f.write("</PP_INFO>\n")
        
        # PP_HEADER
        attrs = " ".join(f'{k}="{v}"' for k, v in header_attrs.items())
        f.write(f"<PP_HEADER {attrs}/>\n")
        
        # PP_MESH
        f.write(f'<PP_MESH dx="{dx:.8e}" mesh="{r.size}" xmin="{xmin:.8e}" rmax="{r[-1]:.15E}" zmesh="{zmesh:.8f}">\n')
        f.write(f'  <PP_R type="real" size="{r.size}" columns="6">\n')
        f.write("    " + format_array(r) + "\n")
        f.write("  </PP_R>\n")
        f.write(f'  <PP_RAB type="real" size="{rab.size}" columns="6">\n')
        f.write("    " + format_array(rab) + "\n")
        f.write("  </PP_RAB>\n")
        f.write("</PP_MESH>\n")
        
        # PP_NLCC (if present)
        if rho_core is not None:
            f.write(f'<PP_NLCC type="real" size="{rho_core.size}" columns="6">\n')
            f.write("  " + format_array(rho_core) + "\n")
            f.write("</PP_NLCC>\n")
        
        # PP_LOCAL
        f.write(f'<PP_LOCAL type="real" size="{vloc.size}" columns="6">\n')
        f.write("  " + format_array(vloc*2) + "\n") ## From Ha2Ry
        f.write("</PP_LOCAL>\n")
        
        f.write("</UPF>\n")
    
    print(f"Converted PSP8 → UPF: {psp8_path} → {out_path}")

def read_psp8(path: str) -> Dict[str, Any]:
    """Read PSP8 format pseudopotential file.
    
    Returns a dictionary with:
    - element: element symbol
    - z_atom: atomic number
    - z_valence: valence charge
    - r: radial grid
    - vloc: local potential
    - rho_core: core charge density (if present)
    - header: header information
    """
    with open(path, "r") as f:
        lines = f.readlines()
    
    # Parse header
    header = {}
    header["format"] = lines[0].strip()  # Should be "DFTpy"
    
    # Line 2: zatom, zion, pspd
    parts = lines[1].split()
    z_atom = float(parts[0])
    z_valence = float(parts[1])
    pspd = parts[2] if len(parts) > 2 else ""
    header["z_atom"] = z_atom
    header["z_valence"] = z_valence
    header["pspd"] = pspd
    
    # Line 3: pspcod, pspxc, lmax, lloc, mmax, r2well
    parts = lines[2].split()
    pspcod = int(parts[0])
    pspxc = int(parts[1])
    lmax = int(parts[2])
    lloc = int(parts[3])
    mmax = int(parts[4])
    r2well = int(parts[5]) if len(parts) > 5 else 0
    header["pspcod"] = pspcod
    header["pspxc"] = pspxc
    header["lmax"] = lmax
    header["lloc"] = lloc
    header["mmax"] = mmax
    header["r2well"] = r2well
    
    # Line 4: rchrg, fchrg, qchrg
    parts = lines[3].split()
    rchrg = float(parts[0])
    fchrg = float(parts[1])
    qchrg = float(parts[2])
    header["rchrg"] = rchrg
    header["fchrg"] = fchrg
    header["qchrg"] = qchrg
    
    # Line 5: nproj
    nproj = int(lines[4].split()[0])
    header["nproj"] = nproj
    
    # Line 6: extension_switch
    extension_switch = int(lines[5].split()[0])
    header["extension_switch"] = extension_switch
    
    # Line 7: version or format indicator
    version = int(lines[6].strip())
    header["version"] = version
    
    # Parse radial grid and local potential (lines 8+)
    r = []
    vloc = []
    
    for line in lines[7:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 3:
            idx = int(parts[0])
            r_val = float(parts[1])
            v_val = float(parts[2])
            r.append(r_val)
            vloc.append(v_val)
        elif len(parts) == 2:
            # Might be continuation of previous line
            r_val = float(parts[0])
            v_val = float(parts[1])
            r.append(r_val)
            vloc.append(v_val)
    
    r = np.array(r)
    vloc = np.array(vloc)
    
    # Check for NLCC (non-linear core correction) - would be after vloc if present
    rho_core = None
    # For now, assume no NLCC unless we detect it
    
    # Determine element from atomic number
    element = None
    # for sym, z in _ATOMIC_Z.items():
    #     if abs(z - z_atom) < 0.1:
    #         element = sym
    #         break
    if element is None:
        # Fallback: use atomic number to get element symbol
        try:
            from ase.data import chemical_symbols
            element = chemical_symbols[int(z_atom)]
        except:
            element = "X"
    
    return {
        "element": element,
        "z_atom": z_atom,
        "z_valence": z_valence,
        "r": r,
        "vloc": vloc,
        "rho_core": rho_core,
        "header": header
    }


def guess_dx_xmin(r: np.ndarray, zmesh: Optional[float]) -> tuple[float, float]:
    """Heuristic for UPF <PP_MESH> attrs; QE uses numeric R/RAB, so rough values are OK."""
    # Avoid log(0)
    rpos = r[r > 0]
    if rpos.size > 1:
        lograt = np.log(rpos[1:] / rpos[:-1])
        dx = float(np.median(lograt))
        xmin = float(np.log(rpos[0] * (zmesh if zmesh else 1.0)))
    else:
        dx, xmin = 0.0, 0.0
    return dx, xmin

def format_array(a: np.ndarray, cols: int = 6) -> str:
    """Compact numeric block with fixed columns."""
    lines = []
    for i in range(0, a.size, cols):
        chunk = a[i:i+cols]
        lines.append(" ".join(f"{v: .15E}" for v in chunk))
    return "\n".join(lines)

