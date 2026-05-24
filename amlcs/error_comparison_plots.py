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


def error_comparison(plots_path, ana_mc, ana_letkf, ana_lenkf, var, lvl, ppt):
    str_lvl = str(lvl)
    zero = pd.Series(ana_mc[str_lvl].values[0])
    data_mc_level = np.log(zero.append(ana_mc[str_lvl]))
    data_letkf_level = np.log(zero.append(ana_letkf[str_lvl]))
    data_lenk_level = np.log(zero.append(ana_lenkf[str_lvl]))
    fr = np.log(ppt.noda[var][lvl, :])

    plt.figure(figsize=(9, 4))

    plt.plot(data_mc_level, label="EnKF-MC-obs")
    plt.plot(data_letkf_level, label="LETKF")
    plt.plot(data_lenk_level, label="LEnKF")

    plt.plot(fr, color="k", label="NODA")

    plt.title(f"$\mathrm{{{var_codes[var]}}}\ at\ {{{pslvl[lvl]}}}mb}}$")
    plt.plot(fr, color="k", label="NODA")

    plt.legend()
    plt.ylabel(r"$\log(\mathcal{l}_2)$")
    plt.xlabel(r"$\mathrm{Assimilation\;Step}$")
    plt.legend(loc="best", prop={"size": 14})

    plt.savefig(plots_path / f"error_{var}_{lvl}.png", bbox_inches="tight")
    plt.close()
    # plt.show()


def ver_error(plots_path, ana_mc, ana_letkf, ana_lenkf, var):
    rmse_mc = []
    rmse_letkf = []
    rmse_lenkf = []

    pslvl_p = pslvl
    for i in range(8):
        id = str(i)

        if ("TRG" in var) and i < 2:
            pslvl_p = pslvl[2:]
            continue

        rmse_mc.append(np.sqrt(sum(ana_mc[id] ** 2) / len(ana_mc[id])))
        err_letkf = np.sqrt(sum(ana_letkf[id] ** 2) / len(ana_letkf[id]))
        rmse_letkf.append(err_letkf)
        err_lenkf = np.sqrt(sum(ana_lenkf[id] ** 2) / len(ana_lenkf[id]))
        rmse_lenkf.append(err_lenkf)

    plt.figure(figsize=(4, 6))
    plt.plot(rmse_mc, pslvl_p, label="EnKF-MC-obs")
    plt.plot(rmse_letkf, pslvl_p, label="LETKF")
    plt.plot(rmse_lenkf, pslvl_p, label="LEnKF")

    plt.title(f"${{{var_codes[var]}}}$")

    plt.ylabel(r"$\mathrm{Pressure\; Levels}\ (mb)$")
    plt.xlabel(r"$\log(RMSE)$")
    plt.legend(loc="best", prop={"size": 14})

    plt.savefig(plots_path / f"verror_{var}.png", bbox_inches="tight")
    plt.close()


def main_general_plotter(df_params):
    # input_file = sys.argv[1]
    root_path = Path.cwd()
    exp_pth = root_path.parents[0] / "runs"
    plots_path = exp_pth / "plots"

    for _, row in df_params.iterrows():
        mc_path = exp_pth / row["mc_path"]
        letkf_path = exp_pth / row["letkf_path"]
        lenkf_path = exp_pth / row["lenkf_path"]
        variables = row["variable"]
        type = row["type"]
        grid_res = row["resolution"]

        tri_path = plots_path / "tri" / row["letkf_path"].split("LETKF")[1][1:]
        Path(tri_path).mkdir(parents=True, exist_ok=True)

        gs = grid_resolution(grid_res)
        ppt = postpro_tools(grid_res, gs, mc_path, 30)
        ppt.compute_NODA()

        if np.isnan(variables):
            variables = model_vars
        else:
            variables = variables.strip().split(",")

        for var in variables:
            ana_mc = pd.read_csv(mc_path / "results" / f"{var}_ana.csv")
            ana_letkf = pd.read_csv(letkf_path / "results" / f"{var}_ana.csv")
            ana_lenkf = pd.read_csv(lenkf_path / "results" / f"{var}_ana.csv")

            if type == "errors":
                levels = row["level"]
                if np.isnan(levels):
                    levels = range(8)
                else:
                    levels = levels.strip().split(",")
                    levels = [int(v) for v in levels]

                for lvl in levels:
                    if ("PSG" in var) and lvl > 0:
                        break
                    if ("TRG" in var) and lvl < 2:
                        continue
                    error_comparison(
                        tri_path, ana_mc, ana_letkf, ana_lenkf, var, lvl, ppt
                    )
            elif type == "errors_vs_level":
                if "PSG" in var:
                    continue
                ver_error(tri_path, ana_mc, ana_letkf, ana_lenkf, var)

            print(f"* ENDJ - Plot {var} {lvl} Finished")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Creates a heatmap for a defined set of parameters "
        + "using the specified configurations.\n"
        + "Remember the config file must be a CSV containing the headers:\n"
        + "mc_path,letkf_path,lenkf_path,resolution,type,variable,levels\n"
        + "If no level or variable is set, default values will be used "
        + "except for the errors by level, which uses all default levels",
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
