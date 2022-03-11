from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from grid_resolution import grid_resolution


class postpro_tools:
    path = Path.cwd()
    var_names = [
        "UG0",
        "UG1",
        "VG0",
        "VG1",
        "TG0",
        "TG1",
        "PSG0",
        "PSG1",
        "TRG0",
        "TRG1",
    ]

    def __init__(self, resolution, grid_resolution_module, path_exp, M):
        self.res = resolution
        self.gs = grid_resolution_module
        self.folder_path = self.path / path_exp
        self.M = M

        self.var_resol = []
        self.error_steps = []
        self.noda = {}

        self.set_resol_variables()
        self.create_errors()

    def set_resol_variables(self):
        nvar = len(self.var_names)
        for v in range(0, nvar):
            if "PSG" in self.var_names[v]:
                self.var_resol.append([self.gs.get_resolution(self.res), 1])
            else:
                self.var_resol.append([self.gs.get_resolution(self.res), 8])

    def create_errors(self):
        for v in self.var_names:
            if "PSG" in v:
                self.error_steps.append(np.zeros((1, self.M)))
            else:
                self.error_steps.append(np.zeros((8, self.M)))

    def compute_NODA(self):
        free_run_path = self.folder_path / "free_run"
        ref_sol_path = self.folder_path / "snapshots"

        for k in range(self.M):
            free_run = Dataset(free_run_path / f"free_run_{k}.nc", "r")
            ref_sol = Dataset(ref_sol_path / f"reference_solution_{k}.nc", "r")

            var_resol = self.var_resol
            n_vars = len(self.var_names)
            for v in range(0, n_vars):
                var_data = self.error_steps[v]
                xm = free_run[self.var_names[v]]
                xr = ref_sol[self.var_names[v]]
                if "PSG" in self.var_names[v]:
                    err_lev = np.linalg.norm(
                        xm[:, :].reshape(-1,) - xr[:, :].reshape(-1, 1)
                    )
                    var_data[0, k] = err_lev
                else:
                    if "TRG" in self.var_names[v]:
                        xm = xm[0, :, :, :]
                        xr = xr[0, :, :, :]

                    lat, lon = var_resol[v][0]
                    lev = var_resol[v][1]
                    err_lev = np.zeros((lev,))
                    for l in range(0, lev):
                        err_l = np.linalg.norm(
                            xm[l, :, :].reshape(-1, 1) - xr[l, :, :].reshape(-1, 1)
                        )
                        err_lev[l] = err_l
                    var_data[:, k] = err_lev
                self.error_steps[v] = var_data

        for v in range(0, n_vars):
            var = self.var_names[v]
            self.noda[var] = self.error_steps[v]
