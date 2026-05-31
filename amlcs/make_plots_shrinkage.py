#!/usr/bin/env python3
"""
make_plots.py - Post-procesado consolidado de los experimentos AMLCS-DA.

Recolecta los resultados de TODOS los metodos presentes en runs/ y genera
todos los plots en una sola carpeta runs/plots/ organizada:

    runs/plots/
    ├── comparison/      comparacion de todos los metodos: log(l2) vs paso, por variable/nivel
    ├── vs_level/        RMSE vs nivel de presion (perfil vertical), por variable
    └── single/<metodo>/ analysis vs background vs NODA, por metodo/variable/nivel

A diferencia de la version anterior (que comparaba solo 3 metodos fijos), este
script AUTO-DETECTA todas las carpetas de experimentos validas dentro de runs/.
Una carpeta se considera valida si contiene results/<var>_ana.csv.

Por defecto, en auto-deteccion solo se procesan los metodos cuyo nombre
empieza por "EnKF-MC" (EnKF-MC y sus variantes) y se excluye explicitamente
"EnKF-MC-obs". LETKF queda fuera automaticamente por no coincidir con el prefijo.
Esto se puede ajustar con --include-prefixes y --exclude.

Uso:
    python3 make_plots.py                          # auto-detecta EnKF-MC y variantes (sin LETKF ni EnKF-MC-obs)
    python3 make_plots.py --runs /ruta/a/runs      # otra ubicacion de runs/
    python3 make_plots.py --methods carpeta1,carpeta2,...   # lista explicita de carpetas (ignora el filtro)
    python3 make_plots.py --include-prefixes EnKF-MC       # prefijos a incluir en auto-deteccion
    python3 make_plots.py --exclude EnKF-MC-obs,LETKF      # etiquetas exactas a excluir en auto-deteccion
    python3 make_plots.py --ref-method <carpeta>   # carpeta usada para el NODA de comparison (por defecto la 1ra)

Es independiente del directorio desde el que se llama: detecta su propia ubicacion.
No requiere Basemap.
"""

import argparse
import re
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

DEFAULT_RES = "t21"
DEFAULT_M   = 30

# Filtro de metodos por defecto (solo aplica en auto-deteccion):
#   - se incluyen solo las etiquetas que empiezan por alguno de estos prefijos
#   - se excluyen estas etiquetas exactas
# Asi, por defecto se queda con EnKF-MC y sus variantes, fuera LETKF y EnKF-MC-obs.
DEFAULT_INCLUDE_PREFIXES = ["EnKF-MC-obs", "EnKF-Sh"]
DEFAULT_EXCLUDE_LABELS   = ["LEnKF", "LETKF"]

# Prefijo comun de los nombres de carpeta: t21_80_0.05_30_<METODO>_1_2_102_mask_2
# Lo usamos para extraer una etiqueta corta y legible para la leyenda.
_PREFIX_RE = re.compile(r"^t\d+_\d+_[\d.]+_\d+_(?P<method>.+?)_\d+_\d+_\d+_mask_\d+$")


def _read(path):
    return pd.read_csv(path)


def _short_label(folder_name):
    """Genera una etiqueta corta y legible a partir del nombre de la carpeta."""
    m = _PREFIX_RE.match(folder_name)
    raw = m.group("method") if m else folder_name
    # Normaliza: EnKF_MC_obs -> EnKF-MC-obs
    return raw.replace("_", "-")


def _is_valid_method_dir(path, sample_var="TG0"):
    """Una carpeta es metodo valido si tiene results/<var>_ana.csv para alguna var."""
    results = path / "results"
    if not results.is_dir():
        return False
    for var in MODEL_VARS:
        if (results / f"{var}_ana.csv").exists():
            return True
    return False


def _passes_filter(label, include_prefixes, exclude_labels):
    """True si la etiqueta debe conservarse segun include/exclude."""
    if exclude_labels and label in exclude_labels:
        return False
    if include_prefixes and not any(label.startswith(p) for p in include_prefixes):
        return False
    return True


def _discover_methods(runs, explicit=None, include_prefixes=None, exclude_labels=None):
    """Devuelve lista de (label, path) de los metodos a procesar.

    El filtro include/exclude SOLO se aplica en auto-deteccion. Si se pasa una
    lista explicita con --methods, se respeta exactamente esa lista.
    """
    if explicit:
        names = [n.strip() for n in explicit.split(",") if n.strip()]
        dirs = [runs / n for n in names]
        apply_filter = False
    else:
        dirs = sorted(d for d in runs.iterdir()
                      if d.is_dir() and d.name != "plots")
        apply_filter = True

    methods = []
    for d in dirs:
        if not _is_valid_method_dir(d):
            print(f"  (omitida, sin results/*_ana.csv) {d.name}")
            continue
        label = _short_label(d.name)
        if apply_filter and not _passes_filter(label, include_prefixes, exclude_labels):
            print(f"  (filtrada) {label:20s}  ({d.name})")
            continue
        methods.append((label, d))
    return methods


# ----------------------------------------------------------------------
# 1. Comparacion de N metodos: log(l2) vs paso de asimilacion
# ----------------------------------------------------------------------
def plot_comparison(out_dir, methods_ana, var, lvl, ppt):
    s = str(lvl)
    plt.figure(figsize=(9, 4))
    for label, ana in methods_ana:
        zero = pd.Series(ana[s].values[0])
        d = np.log(pd.concat([zero, ana[s]], ignore_index=True))
        plt.plot(d, label=label)
    fr = np.log(ppt.noda[var][lvl, :])
    plt.plot(fr, color="k", label="NODA")
    plt.title(f"${VAR_CODES[var]}$ at {PSLVL[lvl]} mb")
    plt.ylabel(r"$\log(\ell_2)$")
    plt.xlabel("Assimilation step")
    plt.legend(loc="best", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_dir / f"error_{var}_{lvl}.png", bbox_inches="tight", dpi=120)
    plt.close()


# ----------------------------------------------------------------------
# 2. RMSE vs nivel de presion (perfil vertical) para N metodos
# ----------------------------------------------------------------------
def plot_vs_level(out_dir, methods_ana, var):
    markers = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]
    plt.figure(figsize=(4, 6))
    for k, (label, ana) in enumerate(methods_ana):
        rmse = []
        pslvl_p = PSLVL[:]
        for i in range(8):
            idx = str(i)
            if ("TRG" in var) and i < 2:
                pslvl_p = PSLVL[2:]
                continue
            rmse.append(np.sqrt(np.sum(ana[idx] ** 2) / len(ana[idx])))
        plt.plot(rmse, pslvl_p, label=label, marker=markers[k % len(markers)])
    plt.title(f"${VAR_CODES[var]}$")
    plt.ylabel("Pressure level (mb)")
    plt.xlabel("RMSE")
    plt.gca().invert_yaxis()
    plt.legend(loc="best", fontsize=9)
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

    p = argparse.ArgumentParser(description="Post-procesado consolidado AMLCS-DA (N metodos)")
    p.add_argument("--runs", default=str(project_root / "runs"),
                   help="Ruta a la carpeta runs/ (por defecto <proyecto>/runs)")
    p.add_argument("--methods", default=None,
                   help="Lista explicita de carpetas separadas por coma (ignora include/exclude)")
    p.add_argument("--include-prefixes", default=",".join(DEFAULT_INCLUDE_PREFIXES),
                   help="Prefijos de etiqueta a incluir en auto-deteccion, separados por coma "
                        "(por defecto: EnKF-MC). Usa cadena vacia para incluir todo.")
    p.add_argument("--exclude", default=",".join(DEFAULT_EXCLUDE_LABELS),
                   help="Etiquetas exactas a excluir en auto-deteccion, separadas por coma "
                        "(por defecto: EnKF-MC-obs)")
    p.add_argument("--ref-method", default=None,
                   help="Carpeta usada para el NODA en comparison (por defecto la 1ra detectada)")
    p.add_argument("--res", default=DEFAULT_RES)
    p.add_argument("--M", type=int, default=DEFAULT_M)
    p.add_argument("--vars", default=None,
                   help="Lista de variables separadas por coma (por defecto todas)")
    args = p.parse_args()

    runs = Path(args.runs).resolve()

    include_prefixes = [s.strip() for s in args.include_prefixes.split(",") if s.strip()]
    exclude_labels   = [s.strip() for s in args.exclude.split(",") if s.strip()]

    print(f"* runs      : {runs}")
    if args.methods:
        print("* metodos explicitos (--methods), no se aplica filtro include/exclude")
    else:
        print(f"* incluir prefijos : {include_prefixes if include_prefixes else '(todos)'}")
        print(f"* excluir etiquetas: {exclude_labels if exclude_labels else '(ninguna)'}")
    print("* detectando metodos...")
    methods = _discover_methods(runs, args.methods, include_prefixes, exclude_labels)
    if not methods:
        raise SystemExit("No se encontro ningun metodo valido tras el filtro (results/*_ana.csv) en runs/.")

    plots_root = runs / "plots"
    cmp_dir    = plots_root / "comparison"
    lvl_dir    = plots_root / "vs_level"
    cmp_dir.mkdir(parents=True, exist_ok=True)
    lvl_dir.mkdir(parents=True, exist_ok=True)

    variables = MODEL_VARS if args.vars is None else args.vars.strip().split(",")

    # Metodo de referencia para el NODA de comparison (por defecto el 1ro).
    ref_path = None
    if args.ref_method:
        ref_path = runs / args.ref_method
    else:
        ref_path = methods[0][1]

    gs = grid_resolution(args.res)
    ppt = postpro_tools(args.res, gs, str(ref_path), args.M)
    ppt.compute_NODA()

    print(f"* metodos   : {len(methods)}")
    for label, path in methods:
        print(f"    - {label:20s}  ({path.name})")
    print(f"* NODA ref  : {ref_path.name}")
    print(f"* plots     : {plots_root}")
    print("")

    for var in variables:
        # Carga analysis de todos los metodos para esta variable
        methods_ana = []
        for label, path in methods:
            ana_file = path / "results" / f"{var}_ana.csv"
            if not ana_file.exists():
                print(f"  (aviso) falta {ana_file.name} en {label}, se omite para {var}")
                continue
            methods_ana.append((label, _read(ana_file)))

        if not methods_ana:
            print(f"* {var}: sin datos en ningun metodo, se salta")
            continue

        # 1. comparacion por nivel
        for lvl in range(8):
            if ("PSG" in var) and lvl > 0:
                break
            if ("TRG" in var) and lvl < 2:
                continue
            plot_comparison(cmp_dir, methods_ana, var, lvl, ppt)

        # 2. perfil vertical (no aplica a PSG, que es 2D)
        if "PSG" not in var:
            plot_vs_level(lvl_dir, methods_ana, var)

        # 3. single por metodo
        for label, path in methods:
            ana_file = path / "results" / f"{var}_ana.csv"
            bck_file = path / "results" / f"{var}_bck.csv"
            if not (ana_file.exists() and bck_file.exists()):
                continue
            single_dir = plots_root / "single" / label
            single_dir.mkdir(parents=True, exist_ok=True)
            ana = _read(ana_file)
            bck = _read(bck_file)
            # cada metodo usa su propio NODA
            ppt_m = postpro_tools(args.res, gs, str(path), args.M)
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