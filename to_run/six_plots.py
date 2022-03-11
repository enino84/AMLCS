import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
radii = [1, 3, 5, 7]
obs_comp = [1, 2, 3, 4, 5]

n_radii = len(radii)
n_obs_comp = len(obs_comp)


def plot_six(exp_pth, routes, var, lvl, method, inf, mask, plots_path):
    data = np.zeros((n_obs_comp, n_radii))
    i = 0
    j = 0

    for route in routes:
        ana_method = pd.read_csv(exp_pth / route / "results" / f"{var}_ana.csv")
        ana_lvl_mthd = ana_method[str(lvl)]
        data[i, j] = np.sqrt(sum(ana_lvl_mthd ** 2) / len(ana_lvl_mthd))
        i = i + 1
        if i == n_obs_comp:
            i = 0
            j = j + 1

    obs_comp_per = (np.round(1 / np.array(obs_comp) ** 2, 2) * 100).astype(int)
    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(
        np.log(data), xticklabels=radii, yticklabels=obs_comp_per, cmap="Blues"
    )

    plt.title(f"$\mathrm{{{var_codes[var]}}}\ at\ {{{pslvl[int(lvl)]}}}mb}}$")
    ax.set(
        xlabel="$r$", ylabel=r"$\mathrm{\%\;Observed\;Components}$",
    )

    plt.savefig(plots_path / f"heatmap_{method}_{var}_{lvl}_{inf}_{mask}.png")
    plt.close()


def main_six_plotter(df_params):
    root_path = Path.cwd()
    exp_pth = root_path.parents[0] / "runs"
    heat_path = exp_pth / "plots" / "heatmaps"

    for _, row in df_params.iterrows():
        setting = row["setting"]  # t21_80_0.05_30
        mthd = row["method"]  # EnKF_MC_obs
        inf = row["infla"]  # 102
        mask = row["mask"]  # 2
        variables = row["variable"]
        levels = row["level"]

        if np.isnan(variables):
            variables = model_vars
        else:
            variables = variables.strip().split(",")

        if np.isnan(levels):
            levels = range(8)
        else:
            levels = levels.strip().split(",")
            levels = [int(v) for v in levels]

        root_exp_path = f"{setting}_{mthd}"
        add_set_exp = f"{inf}_mask_{mask}"

        plts_pth = heat_path / f"{root_exp_path}_{add_set_exp}"

        Path(plts_pth).mkdir(parents=True, exist_ok=True)

        print(f"variables: {variables}")
        print(f"levels: {levels}")
        print(f"radii: {radii}")
        print(f"s: {obs_comp}")

        pths_ex = []
        for r in radii:
            for s in obs_comp:
                pths_ex.append(f"{root_exp_path}_{r}_{s}_{add_set_exp}")

        for var in variables:
            for lvl in levels:
                if ("PSG" in var) and lvl > 0:
                    break
                if ("TRG" in var) and lvl < 2:
                    continue

                plot_six(exp_pth, pths_ex, var, lvl, mthd, inf, mask, plts_pth)
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
    main_six_plotter(df_params)
