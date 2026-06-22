# Nicolò Grilli
# Wan Wan Mohammad
# University of Bristol
# 10 Aprile 2022

# parameter optimization procedure

from scipy.optimize import minimize
import numpy as np
import os
import csv
import datetime
import scipy
from moose_simulation import MOOSE_Simulation
from experiment_comp import Experiment_Comp

class Parameters_Optimizer:

    def __init__(self,simulation,experiment,init_variables,bnds_variables,type_of_solver,
                 log_file='optimization_history.csv'):
        self.simulation = simulation # this is a MOOSE_Simulation object
        self.experiment = experiment # this is a Experiment_Comp object
        self.init_variables = init_variables
        self.bnds_variables = bnds_variables
        self.type_of_solver = type_of_solver
        # path of the CSV that records every objective evaluation
        self.log_file = log_file
        # incrementing counter, one per objective evaluation
        self.eval_count = 0
        # extra (sub-metric) column names, frozen when the header is written
        self.extra_columns = None

    # find optimal parameters
    def find_optimal_parameters(self):
        # convert bnds_variables (which may be a dtype=object numpy array of
        # (lo, hi) pairs possibly containing None) into a clean list of
        # plain-float (lo, hi) tuples. None means unbounded -> +/- inf.
        bounds = self._clean_bounds(self.bnds_variables)
        # clamp the initial guess inside the bounds so the optimizer does not
        # start from an out-of-bounds point (some solvers reject that)
        x0 = self._clamp_to_bounds(self.init_variables, bounds)
        # Nelder-Mead and Powell only support bounds= from scipy 1.7.0 (2021).
        # The bounded methods L-BFGS-B / TNC / SLSQP / trust-constr always do.
        if self._bounds_supported(self.type_of_solver):
            risultato = minimize(self.calc_residual, x0,
                                 method=self.type_of_solver,
                                 bounds=bounds,
                                 options={'disp': True})
        else:
            print('WARNING: scipy ' + scipy.__version__ + ' does not support '
                  'bounds for solver "' + str(self.type_of_solver) + '". '
                  'Running UNBOUNDED. Upgrade scipy to >=1.7.0 or set '
                  'type_of_solver to a bounded method such as "L-BFGS-B".')
            risultato = minimize(self.calc_residual, x0,
                                 method=self.type_of_solver,
                                 options={'disp': True})
        return risultato

    # convert bounds to a clean list of (lo, hi) plain-float tuples,
    # mapping None to -inf / +inf so scipy can clip against real numbers
    @staticmethod
    def _clean_bounds(bnds):
        clean = []
        for pair in bnds:
            lo, hi = pair[0], pair[1]
            lo = -np.inf if lo is None else float(lo)
            hi =  np.inf if hi is None else float(hi)
            clean.append((lo, hi))
        return clean

    # clamp the initial guess into the bounds (element-wise)
    @staticmethod
    def _clamp_to_bounds(x0, bounds):
        lo = np.array([b[0] for b in bounds], dtype=float)
        hi = np.array([b[1] for b in bounds], dtype=float)
        return np.clip(np.asarray(x0, dtype=float), lo, hi)

    # True if the running scipy supports bounds= for the given solver.
    # Nelder-Mead / Powell need scipy >= 1.7.0; bounded methods always do.
    # Version is parsed without the optional 'packaging' dependency.
    @staticmethod
    def _bounds_supported(solver):
        gradient_free = {'nelder-mead', 'powell'}
        if str(solver).lower() in gradient_free:
            try:
                parts = scipy.__version__.split('.')
                ver = (int(parts[0]), int(parts[1]))
            except Exception:
                return True
            return ver >= (1, 7)
        return True

    # human-readable names for the parameter columns; prefer the names already
    # stored on the MOOSE_Simulation object, fall back to generic names
    def _param_names(self, n):
        names = getattr(self.simulation, 'opt_variables', None)
        if names is not None and len(names) == n:
            return list(names)
        return ['param_%d' % i for i in range(n)]

    # append one row (parameters + residual + optional sub-metrics) to the log.
    # opens/flushes/closes per row so a killed run keeps its history;
    # never raises, so a logging problem cannot abort the optimization.
    def _log_evaluation(self, x, residual, extra):
        try:
            x = np.atleast_1d(np.asarray(x, dtype=float))
            param_names = self._param_names(len(x))

            # freeze the set of extra columns on the first evaluation,
            # so the CSV schema (header width and order) stays stable
            if self.extra_columns is None:
                self.extra_columns = sorted(extra.keys()) if isinstance(extra, dict) else []

            write_header = (not os.path.exists(self.log_file)) or (os.path.getsize(self.log_file) == 0)

            with open(self.log_file, 'a', newline='') as f:
                writer = csv.writer(f)
                if write_header:
                    header = ['eval'] + param_names + ['residual', 'timestamp'] + self.extra_columns
                    writer.writerow(header)

                row = [self.eval_count]
                row += [repr(float(v)) for v in x]
                row += [repr(float(residual)), datetime.datetime.now().isoformat()]
                # emit extra metrics strictly in the frozen order; missing /
                # None metrics are written as an empty cell so columns never shift
                for col in self.extra_columns:
                    val = extra.get(col) if isinstance(extra, dict) else None
                    row.append('' if val is None else repr(float(val)))

                writer.writerow(row)
                f.flush()
                os.fsync(f.fileno())
        except Exception as exc:  # logging must never break the optimization
            print('Warning: could not write optimization log: %s' % exc)

    # calculate residual
    # this function is used by scipy minimize
    # x is the array of parameters to optimize
    # that will change at each iteration
    def calc_residual(self,x):
        self.eval_count += 1
        self.simulation.create_input_file(x)
        self.simulation.launch_simulation()
        temp_residual = self.experiment.calc_residual()
        self.experiment.plot_exp_sim()
        # optional sub-errors stashed by Cyclic.calc_residual (None for the
        # base monotonic Experiment_Comp, which never sets last_metrics)
        extra = getattr(self.experiment, 'last_metrics', None)
        self._log_evaluation(x, temp_residual, extra)
        return temp_residual
