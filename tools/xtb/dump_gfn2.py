"""Esporta i parametri GFN2-xTB e D4 in js/chem/xtb/gfn2Data.js.

Fonti: parametrizzazione GFN2-xTB distribuita con dxtb (identica a quella di tblite/xtb, Bannwarth, Ehlert, Grimme,
J. Chem. Theory Comput. 15, 1652, 2019), dati di riferimento DFT-D4 di tad-dftd4 (Caldeweyher et al., J. Chem. Phys.
150, 154122, 2019) e costanti di tad-mctc. Richiede: pip install dxtb==0.4.0 tad-dftd4 tad-mctc.
Uso: python3 tools/xtb/dump_gfn2.py > js/chem/xtb/gfn2Data.js
"""
import json, math, sys, tomllib, pathlib
import torch
import dxtb
from dxtb._src.basis.slater import slater_to_gauss
import tad_dftd4 as d4
from tad_dftd4.model import D4Model
from tad_dftd4.model.base import GA_DEFAULT, GC_DEFAULT, WF_DEFAULT
from tad_dftd4.reference import d4 as d4ref
from tad_dftd4.reference.d4.charge_gfn2 import refq as gfn2_refq
from tad_mctc.data.radii import ATOMIC_RADII, COV_D3
from tad_mctc.data.en import PAULING
# tblite converte i livelli con la costante di mctc-lib (CODATA 2010), non con quella di tad-mctc
EV2AU = 1 / 27.21138505
from tad_mctc.units.length import AA2AU

torch.set_default_dtype(torch.float64)
toml = tomllib.loads((pathlib.Path(dxtb.__file__).parent / "_src/param/gfn2/gfn2-xtb.toml").read_text())
SYM = list(toml["element"].keys())
assert len(SYM) == 86 and SYM[0] == "H" and SYM[-1] == "Rn"
LNUM = {"s": 0, "p": 1, "d": 2, "f": 3}

def r(x, d=12):
    return float(f"{x:.{d}g}")

elements = []
for z, sym in enumerate(SYM, start=1):
    e = toml["element"][sym]
    shells = []
    for k, label in enumerate(e["shells"]):
        n, l = int(label[0]), LNUM[label[1]]
        alpha, coeff = slater_to_gauss(torch.tensor(e["ngauss"][k]), torch.tensor(n), torch.tensor(l),
                                       torch.tensor(e["slater"][k]), norm=True)
        shells.append({
            "l": l, "n": n,
            "level": r(e["levels"][k] * EV2AU), "kcn": r(e["kcn"][k] * EV2AU),
            "zeta": e["slater"][k], "shpoly": e["shpoly"][k], "refocc": e["refocc"][k], "lgam": e["lgam"][k],
            "alpha": [r(a, 15) for a in alpha.tolist()], "coeff": [r(c, 15) for c in coeff.tolist()],
        })
    elements.append({
        "sym": sym, "shells": shells,
        "gam": e["gam"], "gam3": e["gam3"], "zeff": e["zeff"], "arep": e["arep"], "en": e["en"],
        "dkernel": e["dkernel"], "qkernel": e["qkernel"], "mprad": e["mprad"], "mpvcn": e["mpvcn"],
    })

# DFT-D4: polarizzabilità di riferimento già scalate con le cariche GFN2 (come D4Model(ref_charges="gfn2"))
numbers = torch.arange(1, 87)
model = D4Model(numbers, ref_charges="gfn2")
alpha = model._get_alpha()                       # (86, 7, 23)
refc = d4ref.refc[numbers]; refcn = d4ref.refcovcn[numbers]; refq = gfn2_refq[numbers]
d4el = []
for i in range(86):
    refs = []
    for k in range(alpha.shape[1]):
        if refc[i, k] <= 0: continue
        refs.append({"cn": r(refcn[i, k].item()), "count": int(refc[i, k]), "q": r(refq[i, k].item()),
                     "alpha": [r(a, 10) for a in alpha[i, k].tolist()]})
    d4el.append(refs)

g = toml["hamiltonian"]["xtb"]
out = {
    "source": "GFN2-xTB (Bannwarth, Ehlert, Grimme, JCTC 15, 1652, 2019) e DFT-D4 (Caldeweyher et al., JCP 150, 154122, 2019); "
              "parametri da dxtb " + dxtb.__version__ + ", tad-dftd4 " + d4.__version__,
    "EV2AU": EV2AU, "AA2AU": AA2AU,
    "hamiltonian": {"wexp": g["wexp"], "enscale": g["enscale"], "kshell": g["shell"]},
    "repulsion": toml["repulsion"]["effective"],
    "charge": toml["charge"]["effective"],
    "thirdorder": toml["thirdorder"]["shell"],
    "multipole": toml["multipole"]["damped"],
    "dispersion": toml["dispersion"]["d4"],
    "d4": {"ga": GA_DEFAULT, "gc": GC_DEFAULT, "wf": WF_DEFAULT,
           "k4": d4.defaults.D4_K4, "k5": d4.defaults.D4_K5, "k6": d4.defaults.D4_K6, "kcn": d4.defaults.D4_KCN,
           "zeff": d4.data.ZEFF()[numbers].tolist(), "gam": [r(x) for x in d4.data.GAM()[numbers].tolist()],
           "r4r2": [r(x) for x in d4.data.R4R2()[numbers].tolist()]},
    "covD3": [r(x) for x in COV_D3()[numbers].tolist()],
    "atomicRadii": [r(x) for x in ATOMIC_RADII()[numbers].tolist()],
    "pauling": [r(x) for x in PAULING()[numbers].tolist()],
    "elements": elements,
    "d4ref": d4el,
}
print("// File generato da tools/xtb/dump_gfn2.py: non modificare a mano.")
print("// " + out["source"])
print("export const GFN2 = " + json.dumps(out, separators=(",", ":")) + ";")
