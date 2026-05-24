#!/usr/bin/env python3
"""
make_plots.py - Post-procesado consolidado de los experimentos AMLCS-DA.

Recolecta los resultados de los tres metodos (EnKF_MC_obs, LETKF, LEnKF) desde
runs/ y genera TODOS los plots en una sola carpeta runs/plots/ organizada:

    runs/plots/
    ├── comparison/      comparacion de los 3 metodos: log(l2) vs paso, por variable/nivel
    ├── vs_level/        RMSE vs nivel de presion (perfil vertical), por variable
    └── single/<metodo>/ analysis vs background vs NODA, por metodo/variable/nivel

Uso:
    python3 make_plots.py                          # usa la config por defecto (3 metodos t21)
    python3 make_plots.py --runs /ruta/a/runs      # otra ubicacion de runs/
    python3 make_plots.py --mc <carpeta_mc> --letkf <carpeta_letkf> --lenkf <carpeta_lenkf>

Es independiente del directorio desde el que se llama: detecta su propia ubicacion.
No requiere Basemap.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Imports locales del proyecto (deben estar en amlcs/)
from grid_resolution import grid_resolution
from postpro_tools import postpro_tools

try:
    import seaborn as sns
    sns.set_style("darkgrid")
except Exception:
    pass

matplotlib.rcParams.update({"font.size": 13})

# ----------------------------------------------------------------------
# Configuracion
# ----------------------------------------------------------------------
MODEL_VARS = ["PSG0", "PSG1", "TG0", "TG1", "TRG0", "TRG1", "UG0", "UG1", "VG0", "VG1"]
VAR_CODES = {
    "TG0": "T_0",  "UG0": "u_0",  "VG0": "v_0",  "TRG0": "Hq_0", "PSG0": "PS_0",
    "TG1": "T_1",  "UG1": "u_1",  "VG1": "v_1",  "TRG1": "Hq_1", "PSG1": "PS_1",
}
PSLVL = [30, 100, 200, 300, 500, 700, 850, 925]

# Carpetas por defecto (las que generaste en runs/)
DEFAULT_MC    = "t21_80_0.05_30_EnKF_MC_obs_1_2_102_mask_2"
DEFAULT_LETKF = "t21_80_0.05_30_LETKF_1_2_102_mask_2"
DEFAULT_LENKF = "t21_80_0.05_30_LEnKF_1_2_102_mask_2"
DEFAULT_RES   = "t21"
DEFAULT_M     = 30


def _read(path):
    return pd.read_csv(path)


# ----------------------------------------------------------------------
# 1. Comparacion de los 3 metodos: log(l2) vs paso de asimilacion
# ----------------------------------------------------------------------
def plot_comparison(out_dir, ana_mc, ana_letkf, ana_lenkf, var, lvl, ppt):
    s = str(lvl)
    zero = pd.Series(ana_mc[s].values[0])
    d_mc    = np.log(pd.concat([zero, ana_mc[s]], ignore_index=True))
    d_letkf = np.log(pd.concat([zero, ana_letkf[s]], ignore_index=True))
    d_lenkf = np.log(pd.concat([zero, ana_lenkf[s]], ignore_index=True))
    fr = np.log(ppt.noda[var][lvl, :])

    plt.figure(figsize=(9, 4))
    plt.plot(d_mc,    label="EnKF-MC-obs")
    plt.plot(d_letkf, label="LETKF")
    plt.plot(d_lenkf, label="LEnKF")
    plt.plot(fr, color="k", label="NODA")
    plt.title(f"${VAR_CODES[var]}$ at {PSLVL[lvl]} mb")
    plt.ylabel(r"$\log(\ell_2)$")
    plt.xlabel("Assimilation step")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_dir / f"error_{var}_{lvl}.png", bbox_inches="tight", dpi=120)
    plt.close()


# ----------------------------------------------------------------------
# 2. RMSE vs nivel de presion (perfil vertical)
# ----------------------------------------------------------------------
def plot_vs_level(out_dir, ana_mc, ana_letkf, ana_lenkf, var):
    rmse_mc, rmse_letkf, rmse_lenkf = [], [], []
    pslvl_p = PSLVL[:]
    for i in range(8):
        idx = str(i)
        if ("TRG" in var) and i < 2:
            pslvl_p = PSLVL[2:]
            continue
        rmse_mc.append(np.sqrt(np.sum(ana_mc[idx] ** 2) / len(ana_mc[idx])))
        rmse_letkf.append(np.sqrt(np.sum(ana_letkf[idx] ** 2) / len(ana_letkf[idx])))
        rmse_lenkf.append(np.sqrt(np.sum(ana_lenkf[idx] ** 2) / len(ana_lenkf[idx])))

    plt.figure(figsize=(4, 6))
    plt.plot(rmse_mc, pslvl_p, label="EnKF-MC-obs", marker="o")
    plt.plot(rmse_letkf, pslvl_p, label="LETKF", marker="s")
    plt.plot(rmse_lenkf, pslvl_p, label="LEnKF", marker="^")
    plt.title(f"${VAR_CODES[var]}$")
    plt.ylabel("Pressure level (mb)")
    plt.xlabel("RMSE")
    plt.gca().invert_yaxis()
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_dir / f"verror_{var}.png", bbox_inches="tight", dpi=120)
    plt.close()


# ----------------------------------------------------------------------
# 3. Por metodo: analysis vs background vs NODA
# ----------------------------------------------------------------------
def plot_single(out_dir, analysis, background, var, lvl, ppt):
    s = str(lvl)
    zero = pd.Series(background[s].values[0])
    d_ana = np.log(pd.concat([zero, analysis[s]], ignore_index=True))
    d_bck = np.log(pd.concat([zero, background[s]], ignore_index=True))
    fr = np.log(ppt.noda[var][lvl, :])

    plt.figure(figsize=(9, 4))
    plt.title(f"${VAR_CODES[var]}$ at {PSLVL[lvl]} mb")
    plt.plot(d_ana, color="r", label="Analysis")
    plt.plot(d_bck, color="b", label="Background")
    plt.plot(fr, color="k", label="NODA")
    plt.ylabel(r"$\log(\ell_2)$")
    plt.xlabel("Assimilation step")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_dir / f"single_error_{var}_{lvl}.png", bbox_inches="tight", dpi=120)
    plt.close()


# ----------------------------------------------------------------------
# Orquestador
# ----------------------------------------------------------------------
def main():
    script_dir = Path(__file__).resolve().parent          # .../AMLCS/amlcs
    project_root = script_dir.parent                        # .../AMLCS

    p = argparse.ArgumentParser(description="Post-procesado consolidado AMLCS-DA")
    p.add_argument("--runs", default=str(project_root / "runs"),
                   help="Ruta a la carpeta runs/ (por defecto <proyecto>/runs)")
    p.add_argument("--mc", default=DEFAULT_MC)
    p.add_argument("--letkf", default=DEFAULT_LETKF)
    p.add_argument("--lenkf", default=DEFAULT_LENKF)
    p.add_argument("--res", default=DEFAULT_RES)
    p.add_argument("--M", type=int, default=DEFAULT_M)
    p.add_argument("--vars", default=None,
                   help="Lista de variables separadas por coma (por defecto todas)")
    args = p.parse_args()

    runs = Path(args.runs).resolve()
    mc_path    = runs / args.mc
    letkf_path = runs / args.letkf
    lenkf_path = runs / args.lenkf

    plots_root = runs / "plots"
    cmp_dir    = plots_root / "comparison"
    lvl_dir    = plots_root / "vs_level"
    cmp_dir.mkdir(parents=True, exist_ok=True)
    lvl_dir.mkdir(parents=True, exist_ok=True)

    variables = MODEL_VARS if args.vars is None else args.vars.strip().split(",")

    # NODA (free-run vs referencia) se calcula desde el experimento MC
    gs = grid_resolution(args.res)
    ppt = postpro_tools(args.res, gs, str(mc_path), args.M)
    ppt.compute_NODA()

    print(f"* runs      : {runs}")
    print(f"* EnKF_MC   : {args.mc}")
    print(f"* LETKF     : {args.letkf}")
    print(f"* LEnKF     : {args.lenkf}")
    print(f"* plots     : {plots_root}")
    print("")

    for var in variables:
        ana_mc    = _read(mc_path    / "results" / f"{var}_ana.csv")
        ana_letkf = _read(letkf_path / "results" / f"{var}_ana.csv")
        ana_lenkf = _read(lenkf_path / "results" / f"{var}_ana.csv")

        # 1. comparacion por nivel
        for lvl in range(8):
            if ("PSG" in var) and lvl > 0:
                break
            if ("TRG" in var) and lvl < 2:
                continue
            plot_comparison(cmp_dir, ana_mc, ana_letkf, ana_lenkf, var, lvl, ppt)

        # 2. perfil vertical (no aplica a PSG, que es 2D)
        if "PSG" not in var:
            plot_vs_level(lvl_dir, ana_mc, ana_letkf, ana_lenkf, var)

        # 3. single por metodo
        for label, mpath, ana in [
            ("EnKF_MC_obs", mc_path, ana_mc),
            ("LETKF", letkf_path, ana_letkf),
            ("LEnKF", lenkf_path, ana_lenkf),
        ]:
            single_dir = plots_root / "single" / label
            single_dir.mkdir(parents=True, exist_ok=True)
            bck = _read(mpath / "results" / f"{var}_bck.csv")
            # cada metodo usa su propio NODA
            ppt_m = postpro_tools(args.res, gs, str(mpath), args.M)
            ppt_m.compute_NODA()
            for lvl in range(8):
                if ("PSG" in var) and lvl > 0:
                    break
                if ("TRG" in var) and lvl < 2:
                    continue
                plot_single(single_dir, ana, bck, var, lvl, ppt_m)

        print(f"* ENDJ - {var} listo")

    print("\n* Todos los plots en:", plots_root)


if __name__ == "__main__":
    main()