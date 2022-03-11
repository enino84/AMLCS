import argparse
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


def single_error_plotter(analysis, background, var, lvl, ppt, plots_path):
    lvl_str = str(lvl)
    zero = pd.Series(background[lvl_str].values[0])
    data_analysis_by_level = np.log(zero.append(analysis[lvl_str]))
    data_backgroudn_by_level = np.log(zero.append(background[lvl_str]))
    fr = np.log(ppt.noda[var][lvl, :])

    plt.figure(figsize=(9, 4))
    plt.title(f"$\mathrm{{{var_codes[var]}}} \ at \ {{{pslvl[lvl]}}}mb}}$")

    plt.plot(data_analysis_by_level, color="r", label="Analysis")
    plt.plot(data_backgroudn_by_level, color="b", label="Background")
    plt.plot(fr, color="k", label="NODA")

    plt.legend()
    plt.ylabel(r"$log(\mathcal{l}_2)$")
    plt.xlabel(r"$\mathrm{Assimilation\;Step}$")
    plt.legend(loc="best", prop={"size": 14})
    plt.tight_layout()
    plt.autoscale()
    plt.savefig(plots_path / f"single_error_{var}_{lvl}.png", bbox_inches="tight")
    plt.close()


def main_general_plotter(df_params):
    # input_file = sys.argv[1]
    root_path = Path.cwd()
    exp_pth = root_path.parents[0] / "runs"

    for _, row in df_params.iterrows():
        method_path = exp_pth / row["exp_path"]

        variables = row["variable"]
        levels = row["level"]
        grid_res = row["resolution"]

        if np.isnan(variables):
            variables = model_vars
        else:
            variables = variables.strip().split(",")

        if np.isnan(levels):
            levels = range(8)
        else:
            levels = levels.strip().split(",")
            levels = [int(v) for v in levels]

        plots_path = method_path / "plots" / "errors"

        gs = grid_resolution(grid_res)
        ppt = postpro_tools(grid_res, gs, method_path, 30)
        ppt.compute_NODA()

        Path(plots_path).mkdir(parents=True, exist_ok=True)

        for var in variables:
            analysis = pd.read_csv(method_path / "results" / f"{var}_ana.csv")
            bckg = pd.read_csv(method_path / "results" / f"{var}_bck.csv")

            for lvl in levels:
                if ("PSG" in var) and lvl > 0:
                    break
                if ("TRG" in var) and lvl < 2:
                    continue
                single_error_plotter(analysis, bckg, var, lvl, ppt, plots_path)
            print(f"* ENDJ - Plot {var} {lvl} Finished")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Creates a heatmap for a defined set of parameters "
        + "using the specified configurations.\n"
        + "Remember the config file must be a CSV containing the headers:\n"
        + "setting,method,infla,mask,variable,level\n"
        + "If no level or variable is set, default values will be used",
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
