"""LiF on the polar detector grid of grid.f90: SA-CASSCF against the IAM, for one fixed
orientation and for a partially aligned ensemble (run from the unpacked source folder, ~1 min).

    python examples/grid_vs_iam.py          -> lif_grid_vs_iam.png
"""
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.tri import Triangulation

import thistle as th

warnings.simplefilter("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
lif = th.load(os.path.join(HERE, "..", "tests", "data", "molpro", "lif_sa4_casscf.out"))
grid = dict(grid="polar", energy_kev=10, theta_max=60, rings=49, azimuths=72)    # grid.f90: 10, 60, 49, 73

rows = {
    "Fixed orientation (Li-F along $y$)":
        [lif.detector(orientation=(90, 90, 0), model=m, **grid) for m in ("ab initio", "iam")],
    r"Aligned along $y$, $\langle P_2\rangle = 0.5$":
        [lif.aligned_ensemble(axis=(0, 1, 0), p2=0.5, model=m, **grid) for m in ("ab initio", "iam")],
}
fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.4), constrained_layout=True)
for r, (title, (mol, iam)) in enumerate(rows.items()):
    diff = th.compare(mol, iam, "total")["%diff total"]
    mol.save(f"lif_row{r}_casscf.dat"); iam.save(f"lif_row{r}_iam.dat")
    a = 1 / th.BOHR
    tri = Triangulation(mol.qvec[:, 0] * a, mol.qvec[:, 1] * a)
    norm = LogNorm(min(mol["total"].min(), iam["total"].min()), max(mol["total"].max(), iam["total"].max()))
    for c, (pat, lab) in enumerate(((mol, "SA-CASSCF"), (iam, "IAM"))):
        im = axs[r, c].tripcolor(tri, pat["total"], shading="gouraud", cmap="viridis", norm=norm)
        axs[r, c].set_title(f"{lab}: total", loc="left")
    fig.colorbar(im, ax=axs[r, :2], label="total intensity (electron units)", shrink=0.9)
    m = np.nanmax(np.abs(diff))
    im = axs[r, 2].tripcolor(tri, diff, shading="gouraud", cmap="RdBu_r", norm=Normalize(-m, m))
    fig.colorbar(im, ax=axs[r, 2], label=r"100 (I $-$ I$_\mathrm{IAM}$) / I$_\mathrm{IAM}$ (%)", shrink=0.9)
    axs[r, 2].set_title("difference from IAM", loc="left")
    axs[r, 0].text(-0.28, 0.5, title, transform=axs[r, 0].transAxes, rotation=90, va="center", ha="center")
for ax in axs.flat:
    ax.set_aspect("equal")
    ax.set_xlabel(r"$q_x$ (1/$\AA$)")
    ax.set_ylabel(r"$q_y$ (1/$\AA$)")
fig.suptitle("LiF, Molpro SA4-CASSCF(6,4) state 1 - 10 keV X-rays - polar grid of grid.f90 (49 rings x 72, 2theta <= 60 deg)")
fig.savefig("lif_grid_vs_iam.png", dpi=170)
print("lif_grid_vs_iam.png written")
