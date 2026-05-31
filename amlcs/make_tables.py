#!/usr/bin/env python3
"""
make_tables.py - Tabla detallada de RMSE para el paper (AMLCS-DA).

Companion de make_plots.py. En vez de figuras, calcula el RMSE promediado en el
tiempo por (metodo, variable, nivel) y lo exporta listo para el paper:

    runs/tables/
    ├── table_detailed.csv     RMSE por variable y nivel (filas) x metodo + NODA (columnas)
    └── table_detailed.tex     misma tabla en LaTeX (booktabs); negrita = mejor metodo por fila

Metrica (identica a plot_vs_level de make_plots.py):
    RMSE = sqrt( mean_t( l2(t)^2 ) )
descartando los primeros --spinup pasos antes de promediar.

Reutiliza la misma auto-deteccion y filtro que make_plots.py: por defecto solo
toma los metodos cuyo nombre empieza por "EnKF-MC" y excluye "EnKF-MC-obs"
(LETKF queda fuera por no coincidir con el prefijo).

Uso:
    python3 make_tables.py                          # auto-detecta EnKF-MC y variantes, spinup por defecto
    python3 make_tables.py --spinup 20              # descarta los primeros 20 pasos
    python3 make_tables.py --fmt "%.4f"             # formato numerico de las celdas
    python3 make_tables.py --methods c1,c2          # lista explicita (ignora include/exclude)
    python3 make_tables.py --ref-method <carpeta>   # carpeta usada para el NODA (por defecto la 1ra)
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

# Imports locales del proyecto (deben estar en amlcs/)
from grid_resolution import grid_resolution
from postpro_tools import postpro_tools

# ----------------------------------------------------------------------
# Configuracion (coherente con make_plots.py)
# ----------------------------------------------------------------------
MODEL_VARS = ["PSG0", "PSG1", "TG0", "TG1", "TRG0", "TRG1", "UG0", "UG1", "VG0", "VG1"]
VAR_CODES = {
    "TG0": "T_0",  "UG0": "u_0",  "VG0": "v_0",  "TRG0": "Hq_0", "PSG0": "PS_0",
    "TG1": "T_1",  "UG1": "u_1",  "VG1": "v_1",  "TRG1": "Hq_1", "PSG1": "PS_1",
}
PSLVL = [30, 100, 200, 300, 500, 700, 850, 925]

DEFAULT_RES    = "t21"
DEFAULT_M      = 30
DEFAULT_SPINUP = 10           # pasos descartados antes de promediar (ajustable con --spinup)
DEFAULT_FMT    = "%.3e"

DEFAULT_INCLUDE_PREFIXES = ["EnKF"]
DEFAULT_EXCLUDE_LABELS   = ["LETKF", "LEnKF"]

_PREFIX_RE = re.compile(r"^t\d+_\d+_[\d.]+_\d+_(?P<method>.+?)_\d+_\d+_\d+_mask_\d+$")


# ----------------------------------------------------------------------
# Deteccion / filtro de metodos (igual que make_plots.py)
# ----------------------------------------------------------------------
def _short_label(folder_name):
    m = _PREFIX_RE.match(folder_name)
    raw = m.group("method") if m else folder_name
    return raw.replace("_", "-")


def _is_valid_method_dir(path):
    results = path / "results"
    if not results.is_dir():
        return False
    for var in MODEL_VARS:
        if (results / f"{var}_ana.csv").exists():
            return True
    return False


def _passes_filter(label, include_prefixes, exclude_labels):
    if exclude_labels and label in exclude_labels:
        return False
    if include_prefixes and not any(label.startswith(p) for p in include_prefixes):
        return False
    return True


def _discover_methods(runs, explicit=None, include_prefixes=None, exclude_labels=None):
    if explicit:
        names = [n.strip() for n in explicit.split(",") if n.strip()]
        dirs = [runs / n for n in names]
        apply_filter = False
    else:
        dirs = sorted(d for d in runs.iterdir()
                      if d.is_dir() and d.name != "plots" and d.name != "tables")
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
# Metrica
# ----------------------------------------------------------------------
def _levels_for(var):
    """Niveles validos por variable (coherente con make_plots.py)."""
    if "PSG" in var:
        return [0]               # superficie (2D)
    if "TRG" in var:
        return [2, 3, 4, 5, 6, 7]
    return [0, 1, 2, 3, 4, 5, 6, 7]


def _time_rmse(values, spinup=0):
    """RMSE en el tiempo: sqrt(mean(l2^2)), descartando los primeros 'spinup' pasos."""
    v = np.asarray(values, dtype=float)
    if spinup > 0:
        v = v[spinup:]
    v = v[~np.isnan(v)]
    if v.size == 0:
        return np.nan
    return float(np.sqrt(np.mean(v ** 2)))


def _level_label(var, lvl):
    if "PSG" in var:
        return "surface"
    return f"{PSLVL[lvl]}"


# ----------------------------------------------------------------------
# Construccion de la tabla
# ----------------------------------------------------------------------
def build_table(methods, variables, ppt, spinup):
    """Devuelve (df, method_labels). df: filas = var@nivel, cols = metodos + NODA."""
    method_labels = [label for label, _ in methods]
    rows = []

    # cache de CSVs por (label, var)
    ana_cache = {}

    def get_ana(path, var):
        key = (str(path), var)
        if key not in ana_cache:
            f = path / "results" / f"{var}_ana.csv"
            ana_cache[key] = pd.read_csv(f) if f.exists() else None
        return ana_cache[key]

    for var in variables:
        for lvl in _levels_for(var):
            row = {
                "variable": f"${VAR_CODES[var]}$",
                "level": _level_label(var, lvl),
            }
            # cada metodo
            for label, path in methods:
                ana = get_ana(path, var)
                col = str(lvl)
                if ana is not None and col in ana.columns:
                    row[label] = _time_rmse(ana[col].values, spinup)
                else:
                    row[label] = np.nan
            # NODA (baseline comun: metodo de referencia)
            try:
                noda_series = ppt.noda[var][lvl, :]
                row["NODA"] = _time_rmse(noda_series, spinup)
            except Exception:
                row["NODA"] = np.nan
            rows.append(row)

    df = pd.DataFrame(rows, columns=["variable", "level"] + method_labels + ["NODA"])
    return df, method_labels


# ----------------------------------------------------------------------
# Exportadores
# ----------------------------------------------------------------------
def export_csv(df, path):
    df.to_csv(path, index=False)


def _fmt_num(x, fmt):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "--"
    return fmt % x


def export_latex(df, method_labels, path, fmt, caption, label_tag):
    cols = list(df.columns)                       # variable, level, <metodos...>, NODA
    ncol = len(cols)
    align = "ll" + "r" * (ncol - 2)

    lines = []
    lines.append(r"\begin{table}[!ht]")
    lines.append(r"\centering")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label_tag}}}")
    lines.append(rf"\begin{{tabular}}{{{align}}}")
    lines.append(r"\toprule")

    header = ["Variable", "Level (mb)"] + method_labels + ["NODA"]
    lines.append(" & ".join(header) + r" \\")
    lines.append(r"\midrule")

    prev_var = None
    for _, r in df.iterrows():
        # mejor (minimo) entre los metodos de esta fila, para resaltar en negrita
        vals = {m: r[m] for m in method_labels if not (isinstance(r[m], float) and np.isnan(r[m]))}
        best = min(vals, key=vals.get) if vals else None

        # agrupacion visual: deja la celda de variable vacia si se repite
        var_cell = r["variable"] if r["variable"] != prev_var else ""
        prev_var = r["variable"]

        cells = [var_cell, str(r["level"])]
        for m in method_labels:
            s = _fmt_num(r[m], fmt)
            if m == best and s != "--":
                s = r"\textbf{" + s + "}"
            cells.append(s)
        cells.append(_fmt_num(r["NODA"], fmt))
        lines.append(" & ".join(cells) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def print_console(df, method_labels, fmt):
    cols = ["variable", "level"] + method_labels + ["NODA"]
    widths = {c: max(len(c), 11) for c in cols}
    head = "  ".join(f"{c:>{widths[c]}}" for c in cols)
    print(head)
    print("-" * len(head))
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            if c in ("variable", "level"):
                cells.append(f"{str(r[c]):>{widths[c]}}")
            else:
                cells.append(f"{_fmt_num(r[c], fmt):>{widths[c]}}")
        print("  ".join(cells))


# ----------------------------------------------------------------------
# Orquestador
# ----------------------------------------------------------------------
def main():
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent

    p = argparse.ArgumentParser(description="Tabla detallada de RMSE para el paper (AMLCS-DA)")
    p.add_argument("--runs", default=str(project_root / "runs"))
    p.add_argument("--methods", default=None,
                   help="Lista explicita de carpetas separadas por coma (ignora include/exclude)")
    p.add_argument("--include-prefixes", default=",".join(DEFAULT_INCLUDE_PREFIXES),
                   help="Prefijos de etiqueta a incluir en auto-deteccion (por defecto: EnKF-MC)")
    p.add_argument("--exclude", default=",".join(DEFAULT_EXCLUDE_LABELS),
                   help="Etiquetas exactas a excluir en auto-deteccion (por defecto: EnKF-MC-obs)")
    p.add_argument("--ref-method", default=None,
                   help="Carpeta usada para el NODA (por defecto la 1ra detectada)")
    p.add_argument("--spinup", type=int, default=DEFAULT_SPINUP,
                   help=f"Pasos de spin-up a descartar antes de promediar (por defecto: {DEFAULT_SPINUP})")
    p.add_argument("--fmt", default=DEFAULT_FMT,
                   help=f'Formato numerico de las celdas (por defecto: "{DEFAULT_FMT}")')
    p.add_argument("--res", default=DEFAULT_RES)
    p.add_argument("--M", type=int, default=DEFAULT_M)
    p.add_argument("--vars", default=None,
                   help="Lista de variables separadas por coma (por defecto todas)")
    args = p.parse_args()

    runs = Path(args.runs).resolve()
    include_prefixes = [s.strip() for s in args.include_prefixes.split(",") if s.strip()]
    exclude_labels   = [s.strip() for s in args.exclude.split(",") if s.strip()]

    print(f"* runs      : {runs}")
    print(f"* spinup    : {args.spinup} pasos descartados")
    if args.methods:
        print("* metodos explicitos (--methods), no se aplica filtro include/exclude")
    else:
        print(f"* incluir prefijos : {include_prefixes if include_prefixes else '(todos)'}")
        print(f"* excluir etiquetas: {exclude_labels if exclude_labels else '(ninguna)'}")
    print("* detectando metodos...")

    methods = _discover_methods(runs, args.methods, include_prefixes, exclude_labels)
    if not methods:
        raise SystemExit("No se encontro ningun metodo valido tras el filtro (results/*_ana.csv) en runs/.")

    variables = MODEL_VARS if args.vars is None else args.vars.strip().split(",")

    # NODA: metodo de referencia (1ro por defecto)
    ref_path = (runs / args.ref_method) if args.ref_method else methods[0][1]
    gs = grid_resolution(args.res)
    ppt = postpro_tools(args.res, gs, str(ref_path), args.M)
    ppt.compute_NODA()

    print(f"* metodos   : {len(methods)}")
    for label, path in methods:
        print(f"    - {label:20s}  ({path.name})")
    print(f"* NODA ref  : {ref_path.name}")

    df, method_labels = build_table(methods, variables, ppt, args.spinup)

    tables_dir = runs / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    csv_path = tables_dir / "table_detailed.csv"
    tex_path = tables_dir / "table_detailed.tex"

    export_csv(df, csv_path)
    caption = (f"Time-averaged RMSE per variable and pressure level "
               f"(first {args.spinup} assimilation steps discarded). "
               f"Bold indicates the best method per row.")
    export_latex(df, method_labels, tex_path, args.fmt, caption, "tab:rmse_detailed")

    print("\n* Tabla detallada (RMSE promediado en el tiempo):\n")
    print_console(df, method_labels, args.fmt)
    print(f"\n* CSV   : {csv_path}")
    print(f"* LaTeX : {tex_path}")


if __name__ == "__main__":
    main()