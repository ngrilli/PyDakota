# Wan Wan
# Nicolò Grilli
# University of Bristol
# 14 March 2026

# class for cyclic data and calibration using backstress and effective stress
#instead of using cycle_period, strain_rate_elastic_strain_range for calibrating backstress and effective stress
#we use strain_eng from sim csv directly and we detect cycle from positive strain peaks (peak-to-peak segmentation)
#then we detect the elastic part by scanning the reverse branch for the most linear window
#then we find sigma0 using plastic-offset criterion  
#sim output must have minimum: time	avg_ezz	avg_szz	dt_size	gamma_avg	strain_eng	u_front_avg??
#experiment...



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from experiment_comp import Experiment_Comp

class Cyclic(Experiment_Comp):

    def __init__(self,
                 output_file_name, #raw simulation csv
                 experiment_file, #processed experiment
                 name_exp,   #["Cycle number", "Standard force"]
                 name_sim,   #["time", "sim_stress_col??"]
                 prefactors_exp_xy,  #cycle_scale, stress_scale
                 prefactors_sim_xy,  #time_scale, stress_scale
                 figure_file_name,  #base figure filename
                #backstress:simulation csv with gamma
                 sim_time_col="time",  #raw sim time column
                 sim_stress_col="avg_szz",  #raw sim stress column (MPa)
                 sim_strain_col="strain_eng",  #raw sim engineering strain column
                 sim_gamma_col="gamma_avg",  #raw sim cumulative gamma column 
                 #backstress: cycle segmentation controls (peak-to-peak)
                 peak_frac=0.98,  #fraction of strain amplitude used as peak threshold (near-max positive peaks) 
                 peak_min_sep=50,  #min index separation between peaks (prevents picking multiple peaks per cycle) 
                 min_points_per_cycle=50,  #minimum points in a cycle segment (ensures stable fitting)
                 M_TAYLOR=3.06,  #Taylor factor for resolving stresses: tau = sigma_E/M, chi = X/M 
                 #backstress: reverse yield extraction (sigma0)
                 ep_offset=0.0005,  #plastic offset threshold used to define reverse yield 0.05% [from MS.Pham thesis]
                 require_after_eps0=True,  #only search sigma0 after reverse strain crosses <= 0
                 require_after_sig0=False,  #only search sigma0 after stress crosses <= 0
                 elastic_min_pts=10,  #minimum number of points for candidate elastic window line fit
                 elastic_max_pts=40,  #maximum number of points for candidate elastic window line fit
                 elastic_max_start_shift=8,  #how far after peak we allow elastic window to start
                 min_E_GPa=10.0,  #lower physical bound for fitted E (GPa)
                 max_E_GPa=500.0,  #upper physical bound for fitted E (GPa)
                 #setting for backstress/max stress objective weights
                 w_tau=1.0,  #weight for tau_eff_resolved(gamma) residual term
                 w_chi=1.0,  #weight for chi_backstress_resolved(gamma) residual term 
                 w_sigma_max=0.0,  #weight for sigma_max (set >0 to calibrate max stress too)
                 w_sigma0=0.0,  #weight for sigma0_reverse_yield
                 w_sigmaE=0.0,  #weight for sigma_E_effective
                 w_X=0.0,  # weight for X_backstress
                 experiment_file_raw=None,  # optional raw experiment file for sigma_max vs cycle (LCF4_Aftercond.csv)
                 save_sim_cycle_csv=True): #write sim*cycle.csv for debugging 
        
        super().__init__(output_file_name,
                         experiment_file,
                         name_exp,
                         name_sim,
                         prefactors_exp_xy,
                         prefactors_sim_xy,
                         figure_file_name)
        
        self.sim_time_col = sim_time_col 
        self.sim_stress_col = sim_stress_col
        self.sim_strain_col = sim_strain_col 
        self.sim_gamma_col = sim_gamma_col 
        self.experiment_file_raw = experiment_file_raw  #store raw exp filename (can be None)
        
        #backstress
        self.peak_frac = float(peak_frac) 
        self.peak_min_sep = int(peak_min_sep) 
        self.min_points_per_cycle = int(min_points_per_cycle)  

        self.ep_offset = float(ep_offset)
        self.require_after_eps0 = bool(require_after_eps0)  
        self.require_after_sig0 = bool(require_after_sig0) 

        self.elastic_min_pts = int(elastic_min_pts) 
        self.elastic_max_pts = int(elastic_max_pts)  
        self.elastic_max_start_shift = int(elastic_max_start_shift)  
        self.min_E_GPa = float(min_E_GPa)  
        self.max_E_GPa = float(max_E_GPa)  

        self.M_TAYLOR = float(M_TAYLOR)  

        self.w_tau = float(w_tau)  
        self.w_chi = float(w_chi)  
        self.w_sigma_max = float(w_sigma_max)  
        self.w_sigma0 = float(w_sigma0)  
        self.w_sigmaE = float(w_sigmaE)  
        self.w_X = float(w_X)  

        self.save_sim_cycle_csv = bool(save_sim_cycle_csv)  # store whether to output the processed sim cycle CSV
        #plot filenames 
        self.figure_file_name_sigma_max = figure_file_name.replace(".png", "_sigma_max.png") 
        self.figure_file_name_chi_tau = figure_file_name.replace(".png", "_chi_tau.png") 

        # Option B logging: only write scales when best total improves
        self.best_total = None

    # fit straight line y = m*x + b
    def linear_fit(self, x, y):
        A = np.vstack([x, np.ones_like(x)]).T
        m, b = np.linalg.lstsq(A, y, rcond=None)[0]
        return float(m), float(b)
    
    #choose which part of the reverse branch is elastic window
    def r2(self, x, y, m, b):
        yhat = m * x + b
        ss_res = np.sum((y - yhat) ** 2) #residual sum of squares
        ss_tot = np.sum((y - np.mean(y)) ** 2)  #total sum of squares
        if ss_tot <= 0:
            return 0.0
        return float(1.0 - ss_res / ss_tot)
    
    #to make sure gamma is always increase (cumulative gamma)
    def make_gamma_cum(self, g):
        g = np.asarray(g, float)
        m = np.isfinite(g)
        if not np.any(m):
            return g
        first = np.where(m)[0][0]
        g = g - g[first]               # shift so gamma starts at 0
        g = np.maximum.accumulate(g)   # enforce monotone increase
        return g
    
    #cycle segmentation
    
    #defines where cycles start/end 
    def find_positive_strain_peaks(self, eps): #identify cycle boundary peaks in strain 
        eps = np.asarray(eps, float)
        N = len(eps)
        if N < 3:
            return []

        eps_amp = float(np.nanmax(np.abs(eps)))  #estimate strain amplitude
        if not np.isfinite(eps_amp) or eps_amp <= 0:
            return []
        #noisy data exp 
        thr = self.peak_frac * eps_amp   #threshold to determine the max strain 
        cand = []
        for i in range(1, N - 1):     #check local maxima
            if eps[i] >= thr and eps[i] >= eps[i - 1] and eps[i] >= eps[i + 1]:
                cand.append(i)
        if not cand:
            return []
        peaks = [cand[0]]
        for idx in cand[1:]:
            if idx - peaks[-1] >= self.peak_min_sep: #far enough then consider next cycle 
                peaks.append(idx)
            else: #if two peaks are too close, keep the higher one
                if eps[idx] > eps[peaks[-1]]:
                    peaks[-1] = idx
        return peaks
    
    #defines one cycle
    def build_segments(self, peaks, N):
        segs = []
        for a, b in zip(peaks[:-1], peaks[1:]):
            b_incl = min(b + 1, N)
            if (b_incl - a) >= self.min_points_per_cycle: #if too short not a full cycle 
                segs.append((a, b_incl))
        return segs
    
    #define the max sress inside each cycle 
    def tensile_peak_index(self, eps_seg):
        emax = float(np.nanmax(eps_seg))
        if not np.isfinite(emax):
            return 0
        thr = self.peak_frac * emax

        cand = []
        for i in range(1, len(eps_seg) - 1):
            if eps_seg[i] >= thr and eps_seg[i] >= eps_seg[i - 1] and eps_seg[i] >= eps_seg[i + 1]:
                cand.append(i)
        if cand:
            return int(cand[0])

        idx = np.where(eps_seg >= thr)[0] #some data has multiple points at max strain so fallback
        if idx.size:
            return int(idx[0])

        return int(np.argmax(eps_seg))
    
     #--Elastic detection on reverse branch and reverse yield for sigma0--
     ## Elastic detection method (reverse/unloading branch)
    def best_reverse_elastic_fit(self, eps_seg, sig_seg, i_peak):
        N = len(sig_seg)  # segment length
        if i_peak >= N - (self.elastic_min_pts + 2):  # if too few points after peak
            return np.nan, np.nan, None  
        best = None  # store best window
        i0_base = i_peak + 1  
        for shift in range(0, self.elastic_max_start_shift + 1):  
            i0 = i0_base + shift  
            if i0 >= N - self.elastic_min_pts:  # if not enough points for smallest window
                break  
            for w in range(self.elastic_min_pts, self.elastic_max_pts + 1):  # scan window sizes
                j1 = i0 + w  
                if j1 > N:  # if window exceeds segment
                    break  
                x = eps_seg[i0:j1]  # elastic strain
                y = sig_seg[i0:j1]  #elastic stress
                m = np.isfinite(x) & np.isfinite(y) 
                if np.sum(m) < self.elastic_min_pts:  # ensure enough points
                    continue

                x = x[m]  
                y = y[m] 

                E, b = self.linear_fit(x, y)  #fit line sigma 
                if not np.isfinite(E): 
                    continue  #skip
                E_GPa = E / 1000.0  #convert MPa
                if not (self.min_E_GPa <= E_GPa <= self.max_E_GPa):  #incase doesnt make sense 
                    continue  # skip

                R2 = self.r2(x, y, E, b) 
                if best is None or R2 > best["R2"]:  
                    best = {"E": E, "b": b, "R2": R2, "i0": i0, "w": w} 

        if best is None:  #debugging
            return np.nan, np.nan, None  

        info = {"i0_rev": int(best["i0"]), "w_elastic": int(best["w"]), "R2": float(best["R2"])} 
        return float(best["E"]), float(best["b"]), info  #return best modulus, intercept
    
    #detect sigma0 reverse yield
    def estimate_sigma0(self, eps_seg, sig_seg, i_peak):
        E, b, info = self.best_reverse_elastic_fit(eps_seg, sig_seg, i_peak)
        if not np.isfinite(E): #debugging for some cycles too narrow 
            return np.nan, np.nan, np.nan

        eps_p = eps_seg - (sig_seg - b) / E #estimate plastic strain relative to the fitted elastic line
        eps_p_peak = float(eps_p[i_peak])

        i0 = info["i0_rev"]
        eps_rev = eps_seg[i0:]
        sig_rev = sig_seg[i0:]
        eps_p_rev = eps_p[i0:]

        start_k = 0

        #constraint only start searching once strain has crossed zero on the reverse path
        if self.require_after_eps0:
            k0 = np.where(eps_rev <= 0.0)[0]
            start_k = int(k0[0]) if k0.size else 0

        # make sure only after stress crosses zero
        if self.require_after_sig0:
            ks = np.where(sig_rev <= 0.0)[0]
            if ks.size:
                start_k = max(start_k, int(ks[0]))

        delta = np.abs(eps_p_rev - eps_p_peak) #plastic offset relative to tensile peak
        hit = np.where(delta[start_k:] >= self.ep_offset)[0]
        if hit.size == 0: #debugging
            return np.nan, float(E), float(b)

        j = int(hit[0] + start_k) #recover full index
        return float(sig_rev[j]), float(E), float(b)

    #--simulation backstress/effective stress process to csv-- 
    
    def process_simulation_to_cycle_df(self, sim_df):
        t = self.prefactors_sim_xy[0] * sim_df[self.sim_time_col].to_numpy()  
        sig = self.prefactors_sim_xy[1] * sim_df[self.sim_stress_col].to_numpy() 
        eps = sim_df[self.sim_strain_col].to_numpy() 
        g = sim_df[self.sim_gamma_col].to_numpy()  # read simulation cumulative gamma
        #sort by time 
        order = np.argsort(t) 
        t = t[order]  
        sig = sig[order]  
        eps = eps[order] 
        g = g[order] 

        g = self.make_gamma_cum(g) #gamma_cum
        #detect cycle peaks in strain
        peaks = self.find_positive_strain_peaks(eps) 
        segs = self.build_segments(peaks, len(eps))  

        rows = [] 
        # loop over all cycle segments
        for c, (a, b) in enumerate(segs, start=1):  
            eps_seg = eps[a:b] #segment strain
            sig_seg = sig[a:b]  #segment stress
            g_seg = g[a:b]  #segment gamma

            i_peak = self.tensile_peak_index(eps_seg)  
            sigma_max = float(sig_seg[i_peak])  # peak stress at tensile peak

            sigma0, E_fit, b_fit = self.estimate_sigma0(eps_seg, sig_seg, i_peak) #estimate reverse yield stress sigma0

            if np.isfinite(sigma0):  
                sigma_E = 0.5 * (sigma_max + sigma0)  #effective (isotropic-like) part of stress
                X = 0.5 * (sigma_max - sigma0)  #backstress (kinematic-like) part of stress
                tau = sigma_E / self.M_TAYLOR  #resolved effective stress
                chi = X / self.M_TAYLOR  #resolved backstress
            else:  #debugging
                sigma_E = np.nan  
                X = np.nan  
                tau = np.nan 
                chi = np.nan 

            gamma_start = float(g_seg[0])  #gamma at cycle start
            gamma_end = float(g_seg[-1])  #gamma at cycle end
            gamma_cum = gamma_end 
            delta_gamma = gamma_end - gamma_start  #per-cycle gamma increment

            rows.append({  
                "cycle": int(c),  #cycle number in processed CSV
                "sigma_max": float(sigma_max),  #sigma_max in this cycle
                "sigma0_reverse_yield": float(sigma0),  #reverse yield stress in this cycle
                "sigma_E_effective": float(sigma_E),  #effective stress component
                "X_backstress": float(X),  #backstress component
                "tau_eff_resolved": float(tau),  #resolved effective stress
                "chi_backstress_resolved": float(chi),  #resolved backstress
                "E_fit_MPa": float(E_fit),  #fitted E from reverse elastic part 
                "b_fit_MPa": float(b_fit),  #fitted intercept from reverse elastic part 
                "delta_gamma": float(delta_gamma),  #gamma increment this cycle
                "gamma_cum": float(gamma_cum),  #gamma used as comparison coordinate
                "gamma_start": float(gamma_start), 
                "gamma_end": float(gamma_end),  
            })

        sim_cycle_df = pd.DataFrame(rows)  
        #debugging we save the per-cycle sim file
        if self.save_sim_cycle_csv: 
            out_csv = self.output_file_name.replace(".csv", "_cycle_metrics.csv")  
            sim_cycle_df.to_csv(out_csv, index=False)  

        return sim_cycle_df 
    
    #Get experiment sigma_max vs cycle from processed experiment file
    def get_exp_sigma_cycle(self, exp_df):
        if self.experiment_file_raw is not None:
            exp_raw_df = pd.read_csv(self.experiment_file_raw)  # read raw experimental csv

            if ("Cycle number" in exp_raw_df.columns) and ("Standard force" in exp_raw_df.columns): #check
                cyc = np.rint(exp_raw_df["Cycle number"].to_numpy()).astype(int)
                sig = exp_raw_df["Standard force"].to_numpy(dtype=float)
                cycles_u = np.unique(cyc)  #unique cycle numbers in experiment
                x_exp = []
                y_exp = []
                # loop over each cycle number
                for c in cycles_u:
                    mm = (cyc == c) #mask all raw points belonging to this cycle
                    if np.any(mm):
                        x_exp.append(int(c))
                        y_exp.append(float(np.max(sig[mm])))

            # convert to arrays
            x_exp = np.asarray(x_exp, dtype=float)
            y_exp = np.asarray(y_exp, dtype=float)

            # sort by cycle number
            order = np.argsort(x_exp)
            return x_exp[order], y_exp[order]

        if "sigma_max" not in exp_df.columns:
            return None, None
        
        y_exp = exp_df["sigma_max"].to_numpy(dtype=float)
        if "cycle" in exp_df.columns:
            x_exp = exp_df["cycle"].to_numpy(dtype=float)
        elif "Cycle" in exp_df.columns:
            x_exp = exp_df["Cycle"].to_numpy(dtype=float)
        else:
            x_exp = np.arange(1, len(y_exp) + 1, dtype=float)

        order = np.argsort(x_exp)
        return x_exp[order], y_exp[order]


    #Get simulation sigma_max vs cycle
    def get_sim_sigma_cycle(self, sim_df):
        if ("cycle" not in sim_df.columns) or ("sigma_max" not in sim_df.columns): #debugging
            return None, None
        
        x_sim = sim_df["cycle"].to_numpy(dtype=float) # extract cycle
        y_sim = sim_df["sigma_max"].to_numpy(dtype=float) #extract sigma_max
        order = np.argsort(x_sim)
        return x_sim[order], y_sim[order]
            
    #--backstress/effectivestress diagmostic quality check to check sigma0, youngs modulus-- 
    def save_plots(self, exp_df, sim_df):
        x_exp, y_exp = self.get_exp_sigma_cycle(exp_df)
        x_sim, y_sim = self.get_sim_sigma_cycle(sim_df)

        # if both exist, make sigma_max plot
        if (x_exp is not None) and (x_sim is not None):
            max_sim_cycle = np.max(x_sim)

            m_exp = (x_exp <= max_sim_cycle)
            x_exp_plot = x_exp[m_exp]
            y_exp_plot = y_exp[m_exp]

            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(x_sim, y_sim, linewidth=2.5, label="simulation")
            ax.plot(x_exp_plot, y_exp_plot, marker="o", linewidth=2.0, label="experiment")
            ax.set_xlabel("Cycle number")
            ax.set_ylabel("Max stress (MPa)")
            ax.grid(True, alpha=0.3)
            ax.legend()
            fig.tight_layout()
            fig.savefig(self.figure_file_name_sigma_max, dpi=200)
            plt.close(fig)

        #check if all columns exist for combined chi/tau plot 
        have_cols = (
            ("gamma_cum" in exp_df.columns) and ("gamma_cum" in sim_df.columns) and
            ("chi_backstress_resolved" in exp_df.columns) and ("chi_backstress_resolved" in sim_df.columns) and
            ("tau_eff_resolved" in exp_df.columns) and ("tau_eff_resolved" in sim_df.columns)
        )

        #then plot
        if have_cols:
            g_exp = exp_df["gamma_cum"].to_numpy(dtype=float)
            chi_exp = exp_df["chi_backstress_resolved"].to_numpy(dtype=float)
            tau_exp = exp_df["tau_eff_resolved"].to_numpy(dtype=float)
            g_sim = sim_df["gamma_cum"].to_numpy(dtype=float)
            chi_sim = sim_df["chi_backstress_resolved"].to_numpy(dtype=float)
            tau_sim = sim_df["tau_eff_resolved"].to_numpy(dtype=float)
            oe = np.argsort(g_exp) #sort experiment gamma
            os = np.argsort(g_sim) # sort simulation gamma

            #reorder experiment arrays
            g_exp = g_exp[oe]
            chi_exp = chi_exp[oe]
            tau_exp = tau_exp[oe]

            #reorder simulation arrays
            g_sim = g_sim[os]
            chi_sim = chi_sim[os]
            tau_sim = tau_sim[os]

            # create figure
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.plot(g_sim, chi_sim, marker="^", linewidth=0, markersize=6, label="SIM: χ = X/M")
            ax.plot(g_sim, tau_sim, marker="v", linewidth=0, markersize=6, label="SIM: τ = σE/M")
            ax.plot(g_exp, chi_exp, marker="+", linewidth=0, markersize=10, label="EXP: χ = X/M")
            ax.plot(g_exp, tau_exp, marker="+", linewidth=0, markersize=10, label="EXP: τ = σE/M")
            ax.set_xlabel("Cumulative slip γ")
            ax.set_ylabel("Resolved stress (MPa)")
            ax.grid(True, alpha=0.3)
            ax.legend()
            fig.tight_layout()
            fig.savefig(self.figure_file_name_chi_tau, dpi=200)
            plt.close(fig)
    
    #calculate the difference between experiment and simulation: override
    def calc_residual(self):
        #read processed experimental dataframe
        exp_df = pd.read_csv(self.experiment_file)
        #read raw simulation dataframe
        sim_raw_df = pd.read_csv(self.output_file_name)
        sim_df = self.process_simulation_to_cycle_df(sim_raw_df) #process simulation into one row per cycle

        # save processed simulation dataframe with fixed filename
        sim_df.to_csv("resolved_stresses_gamma_simulation.csv", index=False)

        # write diagnostic plots
        self.save_plots(exp_df, sim_df)

        # require gamma columns for chi/tau comparison
        if ("gamma_cum" not in exp_df.columns) or ("gamma_cum" not in sim_df.columns):
            return 1.0e30

        # read experiment gamma 7 simulation
        g_exp = exp_df["gamma_cum"].to_numpy()
        g_sim = sim_df["gamma_cum"].to_numpy()

        # require enough points
        if g_exp.size < 5 or g_sim.size < 5:
            return 1.0e30

        oe = np.argsort(g_exp)
        os = np.argsort(g_sim)

        # reorder gamma arrays
        g_exp = g_exp[oe]
        g_sim = g_sim[os]

        # read experiment chi if present
        chi_exp = exp_df["chi_backstress_resolved"].to_numpy()[oe] if "chi_backstress_resolved" in exp_df.columns else None

        # read experiment tau if present
        tau_exp = exp_df["tau_eff_resolved"].to_numpy()[oe] if "tau_eff_resolved" in exp_df.columns else None

        # read simulation chi if present
        chi_sim = sim_df["chi_backstress_resolved"].to_numpy()[os] if "chi_backstress_resolved" in sim_df.columns else None

        # read simulation tau if present
        tau_sim = sim_df["tau_eff_resolved"].to_numpy()[os] if "tau_eff_resolved" in sim_df.columns else None

        # read experiment sigma0 if present
        sigma0_exp = exp_df["sigma0_reverse_yield"].to_numpy()[oe] if "sigma0_reverse_yield" in exp_df.columns else None

        # read simulation sigma0 if present
        sigma0_sim = sim_df["sigma0_reverse_yield"].to_numpy()[os] if "sigma0_reverse_yield" in sim_df.columns else None

        # read experiment sigma_E if present
        sigmaE_exp = exp_df["sigma_E_effective"].to_numpy()[oe] if "sigma_E_effective" in exp_df.columns else None

        # read simulation sigma_E if present
        sigmaE_sim = sim_df["sigma_E_effective"].to_numpy()[os] if "sigma_E_effective" in sim_df.columns else None

        # read experiment X if present
        X_exp = exp_df["X_backstress"].to_numpy()[oe] if "X_backstress" in exp_df.columns else None

        # read simulation X if present
        X_sim = sim_df["X_backstress"].to_numpy()[os] if "X_backstress" in sim_df.columns else None

        # compute overlap gamma lower bound
        gmin = max(float(np.nanmin(g_exp)), float(np.nanmin(g_sim)))

        # compute overlap gamma upper bound
        gmax = min(float(np.nanmax(g_exp)), float(np.nanmax(g_sim)))

        # mask experiment points inside overlap window
        m = (g_exp >= gmin) & (g_exp <= gmax) & np.isfinite(g_exp)

        # require enough overlap points
        if np.count_nonzero(m) < 5:
            return 1.0e30

        # gamma grid used for interpolation/comparison
        g_use = g_exp[m]

        # interpolate simulation chi onto experiment gamma grid
        chi_sim_on = np.interp(g_use, g_sim, chi_sim) if chi_sim is not None else None

        # interpolate simulation tau onto experiment gamma grid
        tau_sim_on = np.interp(g_use, g_sim, tau_sim) if tau_sim is not None else None

        # interpolate simulation sigma0 onto experiment gamma grid
        sigma0_sim_on = np.interp(g_use, g_sim, sigma0_sim) if sigma0_sim is not None else None

        # interpolate simulation sigma_E onto experiment gamma grid
        sigmaE_sim_on = np.interp(g_use, g_sim, sigmaE_sim) if sigmaE_sim is not None else None

        # interpolate simulation X onto experiment gamma grid
        X_sim_on = np.interp(g_use, g_sim, X_sim) if X_sim is not None else None

        # get experiment sigma_max vs cycle
        x_exp_sig, y_exp_sig = self.get_exp_sigma_cycle(exp_df)

        # get simulation sigma_max vs cycle
        x_sim_sig, y_sim_sig = self.get_sim_sigma_cycle(sim_df)

        # initialize sigma_max residual vector
        dsig = None

        # if both sigma_max curves exist, compare on common cycle numbers
        if (x_exp_sig is not None) and (x_sim_sig is not None):
            exp_map = {int(round(x)): float(y) for x, y in zip(x_exp_sig, y_exp_sig)} # build experiment map cycle->sigma_max
            sim_map = {int(round(x)): float(y) for x, y in zip(x_sim_sig, y_sim_sig)}
            common = sorted(set(exp_map.keys()).intersection(set(sim_map.keys())))

            #require at least a few cycles
            if len(common) >= 3:
                exp_vec = np.array([exp_map[c] for c in common], dtype=float)
                # build simulation sigma_max vector on common cycles
                sim_vec = np.array([sim_map[c] for c in common], dtype=float)
                # compute sigma_max residual vector
                dsig = exp_vec - sim_vec

        #determine whether this is sigma-only mode
        calib_sigma_only = (self.w_sigma_max > 0.0) and (self.w_chi == 0.0) and (self.w_tau == 0.0)

        #choose what the base Experiment_Comp plot should show
        if calib_sigma_only and (x_exp_sig is not None) and (x_sim_sig is not None):
            self.exp_x_data = x_exp_sig
            self.exp_y_data = y_exp_sig
            self.sim_x_data = x_sim_sig
            self.sim_y_data = y_sim_sig
            self.name_exp = ["Cycle number", "sigma_max"]
        else:
            #otherwise show chi vs gamma if available
            if (chi_exp is not None) and (chi_sim_on is not None):
                self.exp_x_data = g_use
                self.exp_y_data = chi_exp[m]
                self.sim_x_data = g_use
                self.sim_y_data = chi_sim_on
                self.name_exp = ["gamma_cum", "chi_backstress_resolved"]

        # initialize list of residual blocks
        r_parts = []

        #append residual if active
        if self.w_chi > 0.0 and (chi_exp is not None) and (chi_sim_on is not None):
            r_parts.append(np.sqrt(self.w_chi) * (chi_exp[m] - chi_sim_on))

        if self.w_tau > 0.0 and (tau_exp is not None) and (tau_sim_on is not None):
            r_parts.append(np.sqrt(self.w_tau) * (tau_exp[m] - tau_sim_on))

        if self.w_sigma_max > 0.0 and (dsig is not None):
            r_parts.append(np.sqrt(self.w_sigma_max) * dsig)

        if self.w_sigma0 > 0.0 and (sigma0_exp is not None) and (sigma0_sim_on is not None):
            r_parts.append(np.sqrt(self.w_sigma0) * (sigma0_exp[m] - sigma0_sim_on))

        if self.w_sigmaE > 0.0 and (sigmaE_exp is not None) and (sigmaE_sim_on is not None):
            r_parts.append(np.sqrt(self.w_sigmaE) * (sigmaE_exp[m] - sigmaE_sim_on))

        if self.w_X > 0.0 and (X_exp is not None) and (X_sim_on is not None):
            r_parts.append(np.sqrt(self.w_X) * (X_exp[m] - X_sim_on))

        if not r_parts:
            return 1.0e30

        r = np.concatenate(r_parts)
        total = float(np.linalg.norm(r))

        #compute RMS
        rms_sigma = float(np.sqrt(np.mean(dsig ** 2))) if dsig is not None else None
        rms_chi = float(np.sqrt(np.mean((chi_exp[m] - chi_sim_on) ** 2))) if (chi_exp is not None and chi_sim_on is not None) else None
        rms_tau = float(np.sqrt(np.mean((tau_exp[m] - tau_sim_on) ** 2))) if (tau_exp is not None and tau_sim_on is not None) else None

        #return scalar objective to optimizer
        return total