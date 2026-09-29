"""Package the ringring stacked profiles and covariances for sharing.

Reads the ThumbStack outputs written by ThumbStack_AllProfiles_ringring.py for the
published ACT DR6.02 NILC fiducial Compton-y map and writes, into ../data/:

  ThumbStack_ringring_ACTDR6.02_fiducial_profiles.csv
      RApArcmin plus pz<N>_act_dr602_fiducial and pz<N>_act_dr602_fiducial_err
  ThumbStack_ringring_ACTDR6.02_fiducial_cov.npz
      one 9 x 9 bootstrap covariance per redshift bin, same keys, plus 'provenance'

Everything is in y * arcmin^2 (covariances in (y * arcmin^2)^2), so the errors and
the covariance are mutually consistent: err = sqrt(diag(cov)).

Note: the earlier zenodo release used y * arcmin^2 for the profiles but left the
covariances in sr^2. This file does not follow that; the units here are uniform.

Run from scripts/:  python make_ringring_data_release.py
"""

import os
import subprocess
import numpy as np
import pandas as pd

PATH_OUT = "/pscratch/sd/r/rhliu/projects/ThumbStack/output/thumbstack/"
MAP_NAME = "act_dr602_fiducial"
MAP_FILE = "act-planck_dr6.02_nilc_ComptonY.fits"
MASK_FILE = "wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits"
FILTER_TYPE = "ringring"
EST = "tsz_uniformweight"
Z_BINS = ["pz1", "pz2", "pz3", "pz4"]
DATA_DIR = "../data/"
STEM = "ThumbStack_ringring_ACTDR6.02_fiducial"

FACTOR = (180. * 60. / np.pi)**2   # sr -> arcmin^2


def stackDir(zbin: str) -> str:
    """Return the ThumbStack output directory for one redshift bin."""
    return PATH_OUT + "DESI_" + zbin + "_" + MAP_NAME + "_" + FILTER_TYPE + "/"


def gitCommit() -> str:
    """Return the current git commit, flagging this script if it is not committed there."""
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL).decode().strip()
        dirty = subprocess.check_output(
            ["git", "status", "--porcelain", "--", "make_ringring_data_release.py"],
            stderr=subprocess.DEVNULL).decode().strip()
        return commit + (" (this script uncommitted at the time of writing)" if dirty else "")
    except (subprocess.CalledProcessError, OSError):
        return "unknown"


columns = {}
covariances = {}
nObj = {}
sanityChecks = []
rApArcmin = None

for zbin in Z_BINS:
    d = stackDir(zbin)
    key = zbin + "_" + MAP_NAME

    # measured profile: R [arcmin], stacked y [sr], analytic sigma [sr]
    measured = np.genfromtxt(d + FILTER_TYPE + "_" + EST + "_measured.txt")
    cov = np.genfromtxt(d + "cov_" + FILTER_TYPE + "_" + EST + "_bootstrap.txt")

    # guard against a pipeline format change silently broadcasting
    if measured.ndim != 2 or measured.shape[1] != 3:
        raise ValueError("unexpected measured.txt shape %s for %s" % (measured.shape, key))
    if cov.shape != (measured.shape[0], measured.shape[0]):
        raise ValueError("unexpected covariance shape %s for %s" % (cov.shape, key))

    if rApArcmin is None:
        rApArcmin = measured[:, 0]
    elif not np.allclose(rApArcmin, measured[:, 0]):
        raise ValueError("aperture radii differ between redshift bins")

    profile = measured[:, 1] * FACTOR                 # y * arcmin^2
    covArcmin = cov * FACTOR**2                       # (y * arcmin^2)^2
    sigma = np.sqrt(np.diag(covArcmin))

    if not (np.isfinite(profile).all() and np.isfinite(covArcmin).all()):
        raise ValueError("non-finite values for " + key)
    if not np.allclose(covArcmin, covArcmin.T):
        raise ValueError("covariance is not symmetric for " + key)
    if not np.all(np.linalg.eigvalsh(covArcmin) > 0.):
        raise ValueError("covariance is not positive definite for " + key)

    # analytic errors are an independent estimate; they should agree closely
    analytic = measured[:, 2] * FACTOR
    ratio = sigma / analytic
    check = "%s: bootstrap/analytic sigma ratio %.3f - %.3f" % (key, ratio.min(), ratio.max())
    print(check)
    sanityChecks.append(check)

    columns[key] = profile
    columns[key + "_err"] = sigma
    covariances[key] = covArcmin
    nObj[zbin] = int(pd.read_csv(d + "overlap_flag.txt", header=None).to_numpy().sum())

# profiles -> CSV (column order: R, then each bin's value and error)
frame = pd.DataFrame({"RApArcmin": rApArcmin})
for zbin in Z_BINS:
    key = zbin + "_" + MAP_NAME
    frame[key] = columns[key]
    frame[key + "_err"] = columns[key + "_err"]

if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)
pathCsv = DATA_DIR + STEM + "_profiles.csv"
frame.to_csv(pathCsv)
print("wrote " + pathCsv)

provenance = "\n".join([
    "ThumbStack stacked tSZ profiles, ring-ring aperture photometry filter.",
    "Compton-y map : " + MAP_FILE + " (published ACT DR6.02 NILC, LAMBDA)",
    "Mask          : " + MASK_FILE,
    "Catalogues    : DESI LRG photometric redshift bins pz1-pz4",
    "Filter        : " + FILTER_TYPE + ", inner hole 1.0 arcmin, outer ring to R*sqrt(2)",
    "Radii         : " + ", ".join("%.1f" % r for r in rApArcmin) + " arcmin",
    "Estimator     : " + EST + " (no hit map, so uniform weighting)",
    "Units         : profiles and errors in y*arcmin^2; covariances in (y*arcmin^2)^2",
    "Errors        : err = sqrt(diag(cov)), from 10000 bootstrap resamples",
    "The apertures overlap (outer ring reaches R*sqrt(2)), so neighbouring radial bins",
    "are correlated: use the covariance for any fitting, not the diagonal errors alone.",
    "Objects overlapping the map, before the point-source and outlier cuts: " +
    ", ".join("%s %d" % (z, nObj[z]) for z in Z_BINS),
    "Checks        : " + "; ".join(sanityChecks),
    "Produced by scripts/make_ringring_data_release.py at commit " + gitCommit(),
])
pathNpz = DATA_DIR + STEM + "_cov.npz"
np.savez(pathNpz, provenance=np.array(provenance), **covariances)
print("wrote " + pathNpz)
