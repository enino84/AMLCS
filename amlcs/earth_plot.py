import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

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

T21_spacing = 5.625
T30_spacing = 3.75
T47_spacing = 2.5
T63_spacing = 1.875
T106_spacing = 1.125

T106_lat = np.arange(-89.1416, 90, T106_spacing)

T63_lat = np.array(
    [
        -88.572,
        -86.723,
        -84.862,
        -82.999,
        -81.135,
        -79.271,
        -77.406,
        -75.541,
        -73.676,
        -71.811,
        -69.946,
        -68.081,
        -66.216,
        -64.351,
        -62.486,
        -60.620,
        -58.755,
        -56.890,
        -55.025,
        -53.160,
        -51.294,
        -49.429,
        -47.564,
        -45.699,
        -43.833,
        -41.968,
        -40.103,
        -38.238,
        -36.372,
        -34.507,
        -32.642,
        -30.777,
        -28.911,
        -27.046,
        -25.181,
        -23.316,
        -21.450,
        -19.585,
        -17.720,
        -15.855,
        -13.989,
        -12.124,
        -10.259,
        -8.394,
        -6.528,
        -4.663,
        -2.798,
        -0.933,
        0.933,
        2.798,
        4.663,
        6.528,
        8.394,
        10.259,
        12.124,
        13.989,
        15.855,
        17.720,
        19.585,
        21.450,
        23.316,
        25.181,
        27.046,
        28.911,
        30.777,
        32.642,
        34.507,
        36.372,
        38.238,
        40.103,
        41.968,
        43.833,
        45.699,
        47.564,
        49.429,
        51.294,
        53.160,
        55.025,
        56.890,
        58.755,
        60.620,
        62.486,
        64.351,
        66.216,
        68.081,
        69.946,
        71.811,
        73.676,
        75.541,
        77.406,
        79.271,
        81.135,
        82.999,
        84.862,
        86.723,
        88.572,
    ]
)

T47_lat = np.array(
    [
        -88.100,
        -85.638,
        -83.161,
        -80.681,
        -78.200,
        -75.719,
        -73.237,
        -70.755,
        -68.272,
        -65.790,
        -63.308,
        -60.825,
        -58.343,
        -55.860,
        -53.377,
        -50.895,
        -48.412,
        -45.930,
        -43.447,
        -40.964,
        -38.482,
        -35.999,
        -33.516,
        -31.034,
        -28.551,
        -26.068,
        -23.586,
        -21.103,
        -18.620,
        -16.138,
        -13.655,
        -11.172,
        -8.689,
        -6.207,
        -3.724,
        -1.241,
        1.241,
        3.724,
        6.207,
        8.689,
        11.172,
        13.655,
        16.138,
        18.620,
        21.103,
        23.586,
        26.068,
        28.551,
        31.034,
        33.516,
        35.999,
        38.482,
        40.964,
        43.447,
        45.930,
        48.412,
        50.895,
        53.377,
        55.860,
        58.343,
        60.825,
        63.308,
        65.790,
        68.272,
        70.755,
        73.237,
        75.719,
        78.200,
        80.681,
        83.161,
        85.638,
        88.100,
    ]
)
T30_lat = np.array(
    [
        -87.16,
        -83.47,
        -79.78,
        -76.07,
        -72.36,
        -68.65,
        -64.94,
        -61.23,
        -57.52,
        -53.81,
        -50.10,
        -46.39,
        -42.68,
        -38.97,
        -35.26,
        -31.54,
        -27.83,
        -24.12,
        -20.41,
        -16.70,
        -12.99,
        -9.28,
        -5.57,
        -1.86,
        1.86,
        5.57,
        9.28,
        12.99,
        16.70,
        20.41,
        24.12,
        27.83,
        31.54,
        35.26,
        38.97,
        42.68,
        46.39,
        50.10,
        53.81,
        57.52,
        61.23,
        64.94,
        68.65,
        72.36,
        76.07,
        79.78,
        83.47,
        87.16,
    ]
)
T21_lat = np.array(
    [
        -85.761,
        -80.269,
        -74.745,
        -69.213,
        -63.679,
        -58.143,
        -52.607,
        -47.070,
        -41.532,
        -35.995,
        -30.458,
        -24.920,
        -19.382,
        -13.844,
        -8.307,
        -2.769,
        2.769,
        8.307,
        13.844,
        19.382,
        24.920,
        30.458,
        35.995,
        41.532,
        47.070,
        52.607,
        58.143,
        63.679,
        69.213,
        74.745,
        80.269,
        85.761,
    ]
)

T21_lon = np.arange(0, 360, T21_spacing)
T30_lon = np.arange(0, 360, T30_spacing)
T47_lon = np.arange(0, 360, T47_spacing)
T63_lon = np.arange(0, 360, T63_spacing)
T106_lon = np.arange(0, 360, T106_spacing)

LATS = {
    "t21": T21_lat,
    "t30": T30_lat,
    "t47": T47_lat,
    "t63": T63_lat,
    "t106": T106_lat,
}

LONS = {
    "t21": T21_lon,
    "t30": T30_lon,
    "t47": T47_lon,
    "t63": T63_lon,
    "t106": T106_lon,
}


def basemap_plot(dataset, resolution, title, path, fig_name):
    lon, lat = np.meshgrid(LONS[resolution], LATS[resolution])
    var_t = dataset
    plt.figure(figsize=(10, 8), edgecolor="w")
    m = Basemap(resolution="l", projection="robin", lat_0=lat[0][0], lon_0=lon[0][0])
    m.shadedrelief(scale=0.5)
    m.drawcoastlines(linewidth=1)
    m.drawcountries(linewidth=1)
    m.drawparallels(np.arange(-90, 90, 20))
    m.drawmeridians(np.arange(0, 357.5, 20))
    color = "jet"
    m.contourf(lon, lat, var_t, latlon=True, cmap=color, alpha=0.6)
    m.colorbar(location="bottom")
    cbar = m.colorbar(location="bottom")
    cbar.set_label(label=title, size=25)
    plt.savefig(path / f"earth_{fig_name}.png")
    plt.close()
    # plt.show()


def earth_plotter(exp_path, plot_path, times, plot_type, variables, levels, res):
    for time in times:
        free_run_path = exp_path / "free_run"
        ref_solution_path = exp_path / "snapshots"
        snapshots_path = exp_path / "time_snapshots"

        if plot_type == "free_run":
            comparison = Dataset(free_run_path / f"free_run_{time}.nc", "r")
        else:
            comparison = Dataset(
                ref_solution_path / f"reference_solution_{time}.nc", "r"
            )

        snapshot = Dataset(snapshots_path / f"xa{time}.nc", "r")

        n_vars = len(variables)
        for v in range(0, n_vars):
            var = variables[v]
            xm = comparison[var]
            xr = snapshot[var]
            if "PSG" in var:
                print(f"Plotting {var}")
                err_lev = abs(xm[:, :] - xr[:, :])
                basemap_plot(
                    err_lev,
                    res,
                    f"${var_codes[var]}$",
                    plot_path,
                    f"{time}_{plot_type}_{var_codes[var]}",
                )
            else:
                if "TRG" in var:
                    xm = xm[0, :, :, :]
                    xr = xr[0, :, :, :]

                for lvl in levels:
                    if "TRG" in var and lvl < 2:
                        continue
                    print(f"Plotting level {lvl} for {var}")
                    err_l = abs(xm[lvl, :, :] - xr[lvl, :, :])
                    title = f"$\mathrm{{{var_codes[var]}}}\ at\ {{{pslvl[lvl]}}}mb}}$"
                    basemap_plot(
                        err_l,
                        res,
                        title,
                        plot_path,
                        f"{time}_{plot_type}_{var_codes[var]}_{pslvl[lvl]}",
                    )


def main_earth_plotter(df_params):
    root_path = Path.cwd()
    root_path = root_path.parents[0] / "runs"

    for _, row in df_params.iterrows():
        exp = row["experiment_name"]
        res = row["resolution"]
        plot_type = row["type"]
        times = row["times"]
        variables = row["variable"]
        levels = row["level"]

        exp_path = root_path / exp
        plot_path = exp_path / "plots" / "earth_errors"
        Path(plot_path).mkdir(parents=True, exist_ok=True)

        if np.isnan(variables):
            variables = model_vars
        else:
            variables = variables.strip().split(",")

        if np.isnan(levels):
            levels = range(8)
        else:
            variables = variables.strip().split(",")
            variables = [int(v) for v in variables]

        times = times.strip().split(",")
        times = [int(t) for t in times]

        print(f"variables: {variables}")
        print(f"levels: {levels}")
        print(f"times: {times}")

        earth_plotter(exp_path, plot_path, times, plot_type, variables, levels, res)
        print("* ENDJ - Plot Finished")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Creates a earthplot using the specified configurations.\n"
        + "Remember the config file must be a CSV containing the headers:\n"
        + "experiment_name,resolution,type,times,level,variable\n"
        + "If no time, level or variable is set, default values will be used",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "file", help="The name of the CSV file containing the configuration"
    )
    args = parser.parse_args()

    input_file = args.file
    print("* STARTJ - Reading input file {0}".format(input_file))
    df_params = pd.read_csv(input_file)
    main_earth_plotter(df_params)
