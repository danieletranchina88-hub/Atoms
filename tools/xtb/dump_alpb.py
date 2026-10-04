"""Esporta i parametri del solvente implicito ALPB (acqua, GFN2-xTB) in js/chem/xtb/alpbData.js.

I parametri sono nei sorgenti Fortran di tblite (Ehlert, Stahn, Spicher, Grimme, J. Chem. Theory Comput. 17, 4250, 2021),
inclusi nel pacchetto sorgente di PyPI:
  curl -LO https://files.pythonhosted.org/packages/source/t/tblite/tblite-0.7.0.tar.gz && tar xzf tblite-0.7.0.tar.gz
Uso: python3 tools/xtb/dump_alpb.py tblite-0.7.0 > js/chem/xtb/alpbData.js
"""
import itertools, json, math, re, sys, pathlib

root = pathlib.Path(sys.argv[1]) / "subprojects/tblite/src/tblite/solvation"


def fortran_block(path, name):
    """Numeri del costruttore Fortran `type(...) :: name = tipo( ... )`, in ordine."""
    text = path.read_text()
    body = text[text.index(f":: {name} ="):]
    nxt = body.find("type(", 10)
    body = body[: nxt if nxt > 0 else len(body)]
    return [float(x) for x in re.findall(r"(-?\d+\.\d+(?:[eEdD][-+]?\d+)?)_wp", body)]


alpb = fortran_block(root / "data/alpb/param_alpb_water.fh", "gfn2_alpb_water")
cds = fortran_block(root / "data/cds/param_alpb_water.fh", "gfn2_alpb_water")
shift = fortran_block(root / "data/shift/param_alpb_water.fh", "gfn2_alpb_water")
assert len(alpb) == 3 + 94 and len(cds) == 1 + 94 + 94 and len(shift) == 3, (len(alpb), len(cds), len(shift))

data = (root / "data.f90").read_text()
blk = data[data.index("vdw_rad_d3(1:94) = aatoau * ["):]
blk = blk[: blk.index("]")]
vdw_d3 = [float(x) for x in re.findall(r"(\d+\.\d+)_wp", blk)]
assert len(vdw_d3) == 94

# costanti CODATA 2018 come in mctc-lib
h, me, c, alpha, e, na = 6.62607015e-34, 9.1093837015e-31, 299792458.0, 7.2973525693e-3, 1.602176634e-19, 6.02214076e23
hbar = h / (2 * math.pi)
bohr = hbar / (me * c * alpha)
hartree = me * c**2 * alpha**2
aatoau = 1 / (bohr * 1e10)
kcaltoau = 1 / (hartree * na * 1e-3 / 4.184)

# griglia di Lebedev a 230 punti (stessi parametri di tblite/mesh/lebedev.f90, ld0230)
def orbit(v):
    pts = set()
    for p in itertools.permutations(v):
        for s in itertools.product([1, -1], repeat=3):
            pts.add(tuple(round(a * b, 15) + 0.0 for a, b in zip(p, s)))
    return sorted(pts)

def oh4(a): b = math.sqrt(1 - 2 * a * a); return orbit((a, a, b))
def oh5(a): b = math.sqrt(1 - a * a); return orbit((a, b, 0.0))
def oh6(a, b): return orbit((a, b, math.sqrt(1 - a * a - b * b)))
s3 = math.sqrt(1 / 3)
orbits = [
    (orbit((1.0, 0.0, 0.0)), -0.5522639919727325e-1),
    (orbit((s3, s3, s3)), 0.4450274607445226e-2),
    (oh4(0.4492044687397611), 0.4496841067921404e-2),
    (oh4(0.2520419490210201), 0.5049153450478750e-2),
    (oh4(0.6981906658447242), 0.3976408018051883e-2),
    (oh4(0.6587405243460960), 0.4401400650381014e-2),
    (oh4(0.4038544050097660e-1), 0.1724544350544401e-1),
    (oh5(0.5823842309715585), 0.4231083095357343e-2),
    (oh5(0.3545877390518688), 0.5198069864064399e-2),
    (oh6(0.2272181808998187, 0.4864661535886647), 0.4695720972568883e-2),
]
grid = [(p, w) for pts, w in orbits for p in pts]
assert len(grid) == 230, len(grid)
assert abs(sum(w for _, w in grid) - 1) < 1e-12

out = {
    "source": "ALPB(acqua) per GFN2-xTB da tblite 0.7.0 (Ehlert et al., JCTC 17, 4250, 2021)",
    "epsilon": alpb[0], "bornScale": alpb[1], "bornOffset": alpb[2] * 0.1 * aatoau, "descreening": alpb[3:],
    "probe": cds[0] * aatoau, "tension": [g * 1e-5 for g in cds[1:95]], "hbond": [-kcaltoau * s * s for s in cds[95:]],
    "gshift": shift[2] * kcaltoau, "vdwD3": [r * aatoau for r in vdw_d3],
    "lebedev230": [list(p) + [w] for p, w in grid],
}
print("// File generato da tools/xtb/dump_alpb.py: non modificare a mano.")
print("// " + out["source"])
print("export const ALPB_WATER = " + json.dumps(out, separators=(",", ":")) + ";")
