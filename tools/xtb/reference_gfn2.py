"""Valori di riferimento GFN2-xTB calcolati con tblite (pip install tblite) per tests/gfn2.mjs.
Uso: python3 tools/xtb/reference_gfn2.py > tests/data/gfn2-reference.json
"""
import json
import numpy as np
from tblite.interface import Calculator
from tblite.library import ffi

A2B = 1 / 0.52917721090
MOLS = {
    "H2O": ([8, 1, 1], [[0, 0, 0], [0, -0.757, 0.587], [0, 0.757, 0.587]], 0, 0),
    "CH4": ([6, 1, 1, 1, 1], [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63], [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63]], 0, 0),
    "HCl_H2O": ([17, 1, 8, 1, 1], [[0, 0, 0], [1.28, 0, 0], [3.0, 0, 0], [3.4, 0.8, 0.3], [3.4, -0.8, 0.3]], 0, 0),
    "NaCl": ([11, 17], [[0, 0, 0], [2.36, 0, 0]], 0, 0),
    "H3O+": ([8, 1, 1, 1], [[0, 0, 0.1], [0.95, 0, -0.25], [-0.47, 0.82, -0.25], [-0.47, -0.82, -0.25]], 1, 0),
    "Cl-_H2O": ([17, 8, 1, 1], [[0, 0, 0], [3.1, 0, 0], [2.15, 0.1, 0.05], [3.4, 0.9, 0.1]], -1, 0),
    "O2_triplet": ([8, 8], [[0, 0, 0], [1.21, 0, 0]], 0, 2),
    "ZnCl2": ([30, 17, 17], [[0, 0, 0], [2.07, 0, 0], [-2.07, 0.1, 0]], 0, 0),
    "CH3I": ([6, 53, 1, 1, 1], [[0, 0, 0], [2.14, 0, 0], [-0.36, 1.03, 0], [-0.36, -0.51, 0.89], [-0.36, -0.51, -0.89]], 0, 0),
    "Fe(CO)": ([26, 6, 8], [[0, 0, 0], [1.8, 0.1, 0], [2.95, 0.1, 0.05]], 0, 2),
    "PbH4": ([82, 1, 1, 1, 1], [[0, 0, 0], [1.0, 1.0, 1.0], [-1.0, -1.0, 1.0], [-1.0, 1.0, -1.0], [1.0, -1.0, -1.0]], 0, 0),
}

def cluster():
    """HCl con 8 acque su una griglia perturbata (seme fisso): sistema in cui il termine ATM non è trascurabile."""
    rng = np.random.default_rng(3)
    Z, xyz = [17, 1], [[0, 0, 0], [1.28, 0, 0]]
    w = np.array([[0, 0, 0], [0.757, 0.587, 0], [-0.757, 0.587, 0]])
    for i in range(-1, 2):
        for j in range(-1, 2):
            if i == 0 and j == 0:
                continue
            c = np.array([i * 3.0, j * 3.0, 1.5 * ((i + j) % 2)]) + rng.normal(0, 0.2, 3)
            th = rng.uniform(0, 6.28)
            R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
            for zz, p in zip([8, 1, 1], w):
                Z.append(zz); xyz.append(list(c + R @ p))
    return Z, xyz


MOLS["HCl_8H2O"] = (*cluster(), 0, 0)
out = {"source": "tblite GFN2-xTB, accuratezza 1e-10, temperatura elettronica 300 K", "unit": "Hartree, bohr", "molecules": {}}
for name, (Z, xyz, charge, uhf) in MOLS.items():
    pos = np.array(xyz, dtype=float) * A2B
    calc = Calculator("GFN2-xTB", np.array(Z), pos, charge=charge, uhf=uhf)
    calc.set("verbosity", 0)
    calc.set("accuracy", 1e-4)
    res = calc.singlepoint()
    out["molecules"][name] = {
        "Z": Z, "pos": pos.round(10).tolist(), "charge": charge, "uhf": uhf,
        "energy": float(res.get("energy")),
        "gradient": np.asarray(res.get("gradient")).round(10).tolist(),
        "charges": np.asarray(res.get("charges")).round(8).tolist(),
        "dipole": np.asarray(res.get("dipole")).round(8).tolist(),
    }
# campo elettrico uniforme (unità atomiche): solo energia e dipolo, perché in tblite 0.7.0 il gradiente con il campo
# non è la derivata della sua stessa energia (verificato alle differenze finite); il gradiente si controlla a parte
out["field"] = {}
for name, F in [("HCl_H2O", [0.01, -0.005, 0.02]), ("H3O+", [0.0, 0.02, -0.01]), ("Fe(CO)", [0.015, 0.0, 0.0])]:
    Z, xyz, charge, uhf = MOLS[name]
    calc = Calculator("GFN2-xTB", np.array(Z), np.array(xyz, dtype=float) * A2B, charge=charge, uhf=uhf)
    calc.set("verbosity", 0)
    calc.set("accuracy", 1e-4)
    calc.add("electric-field", ffi.new("double[3]", F))
    res = calc.singlepoint()
    out["field"][name] = {"field": F, "energy": float(res.get("energy")), "dipole": np.asarray(res.get("dipole")).round(8).tolist()}
print(json.dumps(out, indent=1))
