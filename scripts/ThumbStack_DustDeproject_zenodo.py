import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import corner


import scipy
import scipy.constants
from scipy.optimize import minimize

import emcee

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

# Fundamental Constants
Tcmb = 2.726   # K
h = 6.63e-34   # SI
kB = 1.38e-23  # SI
c = scipy.constants.c

# CIB Deprojection constants (fiducial)
beta = 1.6
T_CIB = 10.7

# CIB and tSZ frequency dependence functions (in temperature units)
def f(nu):
    """frequency dependence for tSZ temperature
    """
    x = h*nu/(kB*Tcmb)
    return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.

def dB_dT(nu, T):
    x = h * nu / (kB * T)
    numerator = 2 * h * nu**3 / c**2 * x * np.exp(x)
    denominator = T * (np.exp(x) - 1)**2
    return numerator / denominator

def f_CIB(nu, beta, Tcib):
    term1 = (nu)**(3 + beta)/(np.exp(h*nu/(kB*Tcib)) - 1)    
    return term1 * dB_dT(nu, Tcmb)

    # Load data
pathProfiles = '/pscratch/sd/r/rhliu/projects/ThumbStack/profiles/'
pathCov = '/pscratch/sd/r/rhliu/projects/ThumbStack/covariances/ACT_DR5/'
pathFig = '../figures/act_dr5_cib_removal/'

filterType = 'ringring2'

df = pd.read_csv(pathProfiles + 'DR5_profiles_T' + '_' + filterType + '.csv', index_col=0)

# Define column keys for the data
DESI_keys = ['DESI_pz1', 'DESI_pz2', 'DESI_pz3', 'DESI_pz4']
# ACT_keys = ['act_dr5_f090_T_2', 'act_dr5_f150_T_2', 'act_dr5_f220_T_2']
ACT_keys = ['act_dr5_f90_ringring2', 'act_dr5_f150_ringring2', 'act_dr5_f220_ringring2']

CMB_namepublic = ['ACT DR5 (90GHz $\Delta T$ map)', 
                  'ACT DR5 (150GHz $\Delta T$ map)',
                  'ACT DR5 (220GHz $\Delta T$ map)']

freq = np.array([90e9, 150e9, 220e9]) # Hz


# Data looks right, so now we can do fitting, first define the functions to fit to:

def CIB_scale(beta, nu, T_CIB=T_CIB):
    return (f_CIB(nu, beta, T_CIB) / f_CIB(220e9, beta, T_CIB))  

def T_to_y(nu):
    factor = (180.*60./np.pi)**2 # unit conversion from sr to arcmin^2
    return 1.0 / (Tcmb * f(nu) * 1.e6)

# Chi-squared calculator
def chi2(beta, y1, y2, y3, C1, C2, C3, nu1, nu2):
    H1 = T_to_y(nu1)
    H2 = T_to_y(nu2)
    G1 = CIB_scale(beta, nu1)
    G2 = CIB_scale(beta, nu2)

    # Residual
    r = H1 * (y1 - G1 * y3) - H2 * (y2 - G2 * y3)

    # Covariance: scaled by squared H terms, which ensures positive definiteness
    Cr = H1**2 * C1 + H2**2 * C2 + (H1 * G1 - H2 * G2)**2 * C3
    Cr_inv = np.linalg.inv(Cr)

    return r.T @ Cr_inv @ r

# Emcee sampling priors and posteriors:

def log_prior(beta):
    if 0.01 < beta < 10.0:
        return 0.0  # log(1)
    return -np.inf
    
def log_likelihood(beta, y1, y2, y3, C1, C2, C3, nu1, nu2):
    return -0.5 * chi2(beta, y1, y2, y3, C1, C2, C3, nu1, nu2)

def log_posterior(beta, y1, y2, y3, C1, C2, C3, nu1, nu2):
    lp = log_prior(beta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(beta, y1, y2, y3, C1, C2, C3, nu1, nu2)

# Now do emcee chain for each of the four bins:
beta_samples = []
beta_medians = []
beta_opt = []
ndim = 1       # one parameter: beta
nwalkers = 32  # number of MCMC walkers
nsteps = 5000  # total steps to run
p0 = np.random.uniform(0.3, 3.9, size=(nwalkers, ndim))

RAp = np.array(df['RApArcmin'])
mask = RAp >=3.0

for i, DESI_key in enumerate(DESI_keys):
    y1 = df[DESI_key + '_' + ACT_keys[0]].to_numpy()[mask]
    y2 = df[DESI_key + '_' + ACT_keys[1]].to_numpy()[mask]
    y3 = df[DESI_key + '_' + ACT_keys[2]].to_numpy()[mask]
    
    C1 = np.load(pathCov + DESI_key + '_' + ACT_keys[0] + '.npy')[2:, 2:]
    C2 = np.load(pathCov + DESI_key + '_' + ACT_keys[1] + '.npy')[2:, 2:]
    C3 = np.load(pathCov + DESI_key + '_' + ACT_keys[2] + '.npy')[2:, 2:]
    
    nu1 = freq[0]
    nu2 = freq[1]

    # Define sampler
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_posterior,
        args=(y1, y2, y3, C1, C2, C3, nu1, nu2)
    )
    
    # Run MCMC
    sampler.run_mcmc(p0, nsteps, progress=True)
    
    # Flatten and discard burn-in
    samples = sampler.get_chain(discard=1000, thin=10, flat=True)

    # Estimate median and confidence interval
    beta_median = np.median(samples)
    beta_lower, beta_upper = np.percentile(samples, [16, 84])
    print(DESI_key)
    print(f"β = {beta_median:.3f} (+{beta_upper - beta_median:.3f}/-{beta_median - beta_lower:.3f})")
    beta_samples.append(samples)
    beta_medians.append(beta_median)

    # --- Minimize chi^2 ---
    res = minimize(
        lambda beta: chi2(beta[0], y1, y2, y3, C1, C2, C3, nu1, nu2),
        x0=[1.0],  # initial guess for beta
        method='L-BFGS-B', 
        bounds=[(0.01, 10.0)])
    
    # --- Output results ---
    if res.success:
        beta_best = res.x[0]
        beta_opt.append(beta_best)
        print(f"Best-fit β = {beta_best:.4f}")
        print(f"Minimum χ² = {res.fun:.4f}")
    else:
        print("Optimization failed:", res.message)

def compute_deproj_profile(beta, y, y3, Cy, C3, nu):
    H = T_to_y(nu)
    G = CIB_scale(beta, nu)

    factor = (180.*60./np.pi)**2 # unit conversion from sr to arcmin^2

    # Compute f
    f = H * (y - G * y3)
    
    # Propagate covariance
    Cf = H**2 * (Cy + G**2 * C3)
    
    # 1σ uncertainty per element
    sigma_f = np.sqrt(np.diag(Cf))
    
    return f, sigma_f, Cf

beta = beta_median
RAp = np.array(df['RApArcmin'])
factor = (180.*60./np.pi)**2 # unit conversion from sr to arcmin^2
df2 = pd.DataFrame()
df2['radii'] = RAp
cov_dict = {}

for i, DESI_key in enumerate(DESI_keys):
    for j, ACT_key in enumerate(ACT_keys):
        if j==2:
            continue
        name = DESI_key + '_' + ACT_key
        name2 = DESI_key + '_' + ACT_keys[2]

        nu = freq[j]
        # beta = beta_medians[i]
        beta = beta_opt[i]
        Covar = np.load(pathCov + name + '.npy')
        C_220 = np.load(pathCov + name2 + '.npy')

        profile, profile_err, profile_cov = compute_deproj_profile(beta, df[name].to_numpy(), df[name2].to_numpy(), Covar, C_220, nu)
        df2[name] = profile
        df2[name + '_err'] = profile_err
        cov_dict[name] = profile_cov / factor**2

df_save_title = '../zenodo/fig15.csv'
df2.to_csv(df_save_title)
np.savez('../zenodo/fig15_cov.npz', **cov_dict)
