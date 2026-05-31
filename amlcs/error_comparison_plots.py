import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from grid_resolution import grid_resolution
from postpro_tools import postpro_tools

import matplotlib

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams.update({"font.size": 14})
sns.set_style("darkgrid")

model_vars = ["PSG0", "PSG1", "TG0", "TG1", "TRG0", "TRG1", "UG0", "UG1", "VG0", "VG1"]

var_codes = {
    "TG0": "T_0",
    "UG0": "u_0",
    "VG0": "v_0",
    "TRG0": "Hq_0",
    "TG1": "T_1",
    "UG1": "u_1",
    "VG1": "v_1",
    "TRG1": "Hq_1",
    "PSG0": "PS_0",
    "PSG1": "PS_1",
}

pslvl = [30, 100, 200, 300, 500, 700, 850, 925]

# ---------------------------------------------------------------------------
#  Method registry: maps the *internal* method name (as used in amlcs_da.py
#  and therefore in the run-folder name) to a nice legend label. New methods
#  only need an entry here to be auto-detected and plotted; the plotting
#  routines never hard-code how many methods there are.
# ---------------------------------------------------------------------------
METHOD_LABELS = {
    "EnKF_MC_obs": "EnKF-MC-obs",
    "LETKF": "LETKF",
    "LEnKF": "LEnKF",
    "EnKF_Sh_Binv_MSE": "EnKF-Sh-Binv-MSE",
    "EnKF_Sh_Binv_Stein": "EnKF-Sh-Binv-Stein",
    "EnKF_Sh_Binv_DA": "EnKF-Sh-Binv-DA",
}

# Order in which methods are drawn / listed (any extra methods found that are
# not in this list are appended afterwards, alphabetically).
METHOD_ORDER = [
    "EnKF_MC_obs",
    "LETKF",
    "LEnKF",
    "EnKF_Sh_Binv_MSE",
    "EnKF_Sh_Binv_Stein",
    "EnKF_Sh_Binv_DA",
]

# A distinct style per method so plots stay readable for any method count.
_PALETTE = sns.color_palette("tab10", n_colors=max(10, len(METHOD_ORDER)))
METHOD_STYLE = {
    name: dict(color=_PALETTE[i % len(_PALETTE)],
               linestyle=["-", "--", "-.", ":"][i % 4])
    for i, name in enumerate(METHOD_ORDER)
}


def _label_for(method_name):
    return METHOD_LABELS.get(method_name, method_name.replace("_", "-"))


def _style_for(method_name):
    if method_name in METHOD_STYLE:
        return METHOD_STYLE[method_name]
    i = abs(hash(method_name)) % len(_PALETTE)
    return dict(color=_PALETTE[i], linestyle="-")


def _ordered(methods):
    """Return the detected method names in a stable, sensible order."""
    known = [m for m in METHOD_ORDER if m in methods]
    extra = sorted(m for m in methods if m not in METHOD_ORDER)
    return known + extra


# ---------------------------------------------------------------------------
#  Auto-detection of method run-folders
# ---------------------------------------------------------------------------
def detect_methods(exp_pth, setting, infla, mask, r=None, s=None):
    """Scan ``exp_pth`` for run folders belonging to this experiment.

    A run folder is named
        {setting}_{METHOD}_{r}_{s}_{infla}_mask_{mask}
    e.g.  t21_80_0.05_30_EnKF_MC_obs_1_2_102_mask_2

    Returns an ordered list of (method_name, folder_name) for every method
    whose results directory exists. The method name is whatever sits between
    the setting prefix and the trailing _{r}_{s}_{infla}_mask_{mask}.
    """
    suffix_re = re.compile(
        r"^" + re.escape(str(setting)) + r"_(?P<method>.+?)_"
        r"(?P<r>\d+)_(?P<s>\d+)_" + re.escape(str(infla)) +
        r"_mask_" + re.escape(str(mask)) + r"$"
    )

    found = {}
    if not exp_pth.exists():
        return []
    for child in sorted(exp_pth.iterdir()):
        if not child.is_dir():
            continue
        m = suffix_re.match(child.name)
        if not m:
            continue
        if r is not None and int(m.group("r")) != int(r):
            continue
        if s is not None and int(m.group("s")) != int(s):
            continue
        method = m.group("method")
        if (child / "results").exists():
            found[method] = child.name

    ordered = _ordered(found.keys())
    return [(mth, found[mth]) for mth in ordered]


# ---------------------------------------------------------------------------
#  Plotting routines (driven by however many methods were detected)
# ---------------------------------------------------------------------------
def error_comparison(plots_path, ana_by_method, var, lvl, ppt):
    """Time series of log-l2 error, one curve per detected method."""
    str_lvl = str(lvl)
    fr = np.log(ppt.noda[var][lvl, :])

    plt.figure(figsize=(9, 4))

    for method_name, ana in ana_by_method.items():
        zero = pd.Series(ana[str_lvl].values[0])
        data_level = np.log(pd.concat([zero, ana[str_lvl]], ignore_index=True))
        st = _style_for(method_name)
        plt.plot(data_level, label=_label_for(method_name),
                 color=st["color"], linestyle=st["linestyle"])

    plt.plot(fr, color="k", label="NODA")

    plt.title(f"$\\mathrm{{{var_codes[var]}}}\\ at\\ {pslvl[lvl]}\\,mb$")
    plt.ylabel(r"$\log(\mathcal{l}_2)$")
    plt.xlabel(r"$\mathrm{Assimilation\;Step}$")
    plt.legend(loc="best", prop={"size": 12})

    plt.savefig(plots_path / f"error_{var}_{lvl}.png", bbox_inches="tight")
    plt.close()


def ver_error(plots_path, ana_by_method, var):
    """RMSE vs pressure level, one curve per detected method."""
    plt.figure(figsize=(4, 6))

    for method_name, ana in ana_by_method.items():
        rmse = []
        pslvl_p = pslvl
        for i in range(8):
            idx = str(i)
            if ("TRG" in var) and i < 2:
                pslvl_p = pslvl[2:]
                continue
            rmse.append(np.sqrt(sum(ana[idx] ** 2) / len(ana[idx])))
        st = _style_for(method_name)
        plt.plot(rmse, pslvl_p, label=_label_for(method_name),
                 color=st["color"], linestyle=st["linestyle"])

    plt.title(f"${{{var_codes[var]}}}$")
    plt.ylabel(r"$\mathrm{Pressure\; Levels}\ (mb)$")
    plt.xlabel(r"$\log(RMSE)$")
    plt.legend(loc="best", prop={"size": 12})

    plt.savefig(plots_path / f"verror_{var}.png", bbox_inches="tight")
    plt.close()


def _load_analysis(method_path, var):
    return pd.read_csv(method_path / "results" / f"{var}_ana.csv")


def main_general_plotter(df_params):
    root_path = Path.cwd()
    exp_pth = root_path.parents[0] / "runs"
    plots_path = exp_pth / "plots"

    for _, row in df_params.iterrows():
        grid_res = row["resolution"]
        type = row["type"]
        variables = row["variable"]

        methods = _resolve_methods(row, exp_pth)
        if not methods:
            print("* ENDJ - WARNING: no method run-folders detected for this "
                  "row; skipping.")
            continue

        method_names = [m[0] for m in methods]
        print(f"* ENDJ - Detected {len(methods)} methods: {method_names}")

        ref_folder = methods[0][1]
        ref_path = exp_pth / ref_folder

        comp_tag = _comparison_tag(row, methods)
        comp_path = plots_path / "comparison" / comp_tag
        Path(comp_path).mkdir(parents=True, exist_ok=True)

        gs = grid_resolution(grid_res)
        ppt = postpro_tools(grid_res, gs, ref_path, 30)
        ppt.compute_NODA()

        if (not isinstance(variables, str)) and _isnan(variables):
            variables = model_vars
        else:
            variables = variables.strip().split(",")

        for var in variables:
            ana_by_method = {}
            for method_name, folder in methods:
                ana_by_method[method_name] = _load_analysis(
                    exp_pth / folder, var)

            lvl = None
            if type == "errors":
                levels = row["level"]
                if (not isinstance(levels, str)) and _isnan(levels):
                    levels = range(8)
                else:
                    levels = [int(v) for v in levels.strip().split(",")]

                for lvl in levels:
                    if ("PSG" in var) and lvl > 0:
                        break
                    if ("TRG" in var) and lvl < 2:
                        continue
                    error_comparison(comp_path, ana_by_method, var, lvl, ppt)
            elif type == "errors_vs_level":
                if "PSG" in var:
                    continue
                ver_error(comp_path, ana_by_method, var)

            print(f"* ENDJ - Plot {var} {lvl} Finished")


def _resolve_methods(row, exp_pth):
    """Resolve the list of (method_name, folder) for one config row.

    Two ways to specify an experiment, in order of precedence:

    1. Auto-detect: provide ``setting`` (e.g. t21_80_0.05_30), ``infla``
       (e.g. 102) and ``mask`` (e.g. 2). Every method run-folder matching
       that experiment is discovered automatically -- 3, 6 or any number.
       Optional ``r`` / ``s`` columns further restrict the match.

    2. Explicit paths: provide one or more ``*_path`` columns (the legacy
       ``mc_path``/``letkf_path``/``lenkf_path`` still work, as does any
       column ending in ``_path`` or named ``path``). The method name is
       inferred from the folder name.
    """
    if "setting" in row.index and isinstance(row.get("setting"), str):
        setting = row["setting"]
        infla = row["infla"]
        mask = row["mask"]
        r = row["r"] if "r" in row.index and not _isnan(row.get("r")) else None
        s = row["s"] if "s" in row.index and not _isnan(row.get("s")) else None
        return detect_methods(exp_pth, setting, infla, mask, r=r, s=s)

    methods = []
    for col in row.index:
        if col == "path" or col.endswith("_path"):
            val = row[col]
            if isinstance(val, str) and val.strip():
                folder = val.strip()
                method_name = _infer_method_from_folder(folder)
                methods.append((method_name, folder))
    seen = {}
    for name, folder in methods:
        seen.setdefault(name, folder)
    ordered = _ordered(seen.keys())
    return [(m, seen[m]) for m in ordered]


def _infer_method_from_folder(folder):
    """Pull the method name out of a run-folder name.

    Folder looks like {setting}_{METHOD}_{r}_{s}_{infla}_mask_{mask}.
    Match the longest known method name that appears in the folder.
    """
    for name in sorted(METHOD_ORDER, key=len, reverse=True):
        if f"_{name}_" in folder or folder.startswith(name + "_") \
                or folder.endswith("_" + name):
            return name
    m = re.match(r"^.*?_(?P<method>[A-Za-z].*?)_\d+_\d+_\d+_mask_\d+$", folder)
    if m:
        return m.group("method")
    return folder


def _comparison_tag(row, methods):
    """A short folder name identifying this comparison."""
    if "setting" in row.index and isinstance(row.get("setting"), str):
        return f"{row['setting']}_{row['infla']}_mask_{row['mask']}"
    return methods[0][1]


def _isnan(v):
    try:
        return bool(np.isnan(v))
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare the analysis error of every data-assimilation "
        "method found for an experiment.\n\n"
        "The config CSV can use EITHER of two formats:\n"
        "  (1) auto-detect : columns  setting,infla,mask,type,resolution,"
        "variable,level\n"
        "      (optionally r,s). Every method run-folder for that experiment "
        "is\n"
        "      discovered and plotted automatically -- 3, 6 or any number.\n"
        "  (2) explicit    : columns  <name>_path,...,type,resolution,"
        "variable,level\n"
        "      (the legacy mc_path,letkf_path,lenkf_path still work).\n",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "file", help="The name of the CSV file containing the configuration"
    )
    args = parser.parse_args()

    input_file = args.file
    print("* STARTJ - Reading input file {0}".format(input_file))
    df_params = pd.read_csv(input_file)
    main_general_plotter(df_params)
