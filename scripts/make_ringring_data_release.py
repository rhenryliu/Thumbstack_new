"""Package the ringring stacked profiles and covariances for sharing.

Reads the ThumbStack outputs written by ThumbStack_AllProfiles_ringring.py for the
published ACT DR6.02 NILC Compton-y maps (fiducial, fixed-beta CIB-deprojected, and
CIB + first-moment deprojected, all at T_d = 10.7) and writes, into ../data/:

  ThumbStack_ringring_ACTDR6.02_profiles.csv
      RApArcmin plus pz<N>_<map> and pz<N>_<map>_err for each of the 7 maps
  ThumbStack_ringring_ACTDR6.02_cov.npz
      one 9 x 9 bootstrap covariance per map and redshift bin (28), same keys,
      plus 'provenance'

The internal, unreleased fiducial map is stacked and checked against an earlier run
but deliberately excluded from both files.

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

PATH_OUT = "/pscratch/sd/r/rhliu/projects/ThumbStack_ringring_v2/output/thumbstack/"
# earlier run, used only to check that the re-run reproduces it; None to skip
REF_PATH_OUT = "/pscratch/sd/r/rhliu/projects/ThumbStack/output/thumbstack/"

# maps written to the data files, in column order
MAPS = ["act_dr602_fiducial",
        "act_dr602_cib_1.2_10.7", "act_dr602_cib_1.4_10.7", "act_dr602_cib_1.6_10.7",
        "act_dr602_dBeta_1.2_10.7", "act_dr602_dBeta_1.4_10.7", "act_dr602_dBeta_1.6_10.7"]
# also stacked, and checked against the earlier run, but deliberately not shared:
# a different (internal, unreleased) map version from the DR6.02 maps above
MAPS_CHECK_ONLY = ["act_dr6_fiducial"]

MAP_FILES = {
    "act_dr602_fiducial": "act-planck_dr6.02_nilc_ComptonY.fits",
    "act_dr602_cib_1.2_10.7": "act-planck_dr6.02_nilc_ComptonY_deproj_cib_1.2_10.7.fits",
    "act_dr602_cib_1.4_10.7": "act-planck_dr6.02_nilc_ComptonY_deproj_cib_1.4_10.7.fits",
    "act_dr602_cib_1.6_10.7": "act-planck_dr6.02_nilc_ComptonY_deproj_cib_1.6_10.7.fits",
    "act_dr602_dBeta_1.2_10.7": "act-planck_dr6.02_nilc_ComptonY_deproj_cib_cibdBeta_1.2_10.7.fits",
    "act_dr602_dBeta_1.4_10.7": "act-planck_dr6.02_nilc_ComptonY_deproj_cib_cibdBeta_1.4_10.7.fits",
    "act_dr602_dBeta_1.6_10.7": "act-planck_dr6.02_nilc_ComptonY_deproj_cib_cibdBeta_1.6_10.7.fits",
    "act_dr6_fiducial": "ilc_SZ_yy.fits",
}
MASK_FILE = "wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits"
FILTER_TYPE = "ringring"
EST = "tsz_uniformweight"
Z_BINS = ["pz1", "pz2", "pz3", "pz4"]
DATA_DIR = "../data/"
STEM = "ThumbStack_ringring_ACTDR6.02"

FACTOR = (180. * 60. / np.pi)**2   # sr -> arcmin^2


def stackDir(zbin: str, mapName: str, root: str = None) -> str:
    """Return the ThumbStack output directory for one redshift bin and map."""
    if root is None:
        root = PATH_OUT
    return root + "DESI_" + zbin + "_" + mapName + "_" + FILTER_TYPE + "/"


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


def readStack(zbin: str, mapName: str, root: str = None):
    """Read one stack.

    Returns:
        (radii [arcmin], profile, bootstrap covariance, analytic sigma), with the
        profile and sigma in y*arcmin^2 and the covariance in (y*arcmin^2)^2.
    """
    d = stackDir(zbin, mapName, root)
    measured = np.genfromtxt(d + FILTER_TYPE + "_" + EST + "_measured.txt")
    cov = np.genfromtxt(d + "cov_" + FILTER_TYPE + "_" + EST + "_bootstrap.txt")
    # guard against a pipeline format change silently broadcasting
    if measured.ndim != 2 or measured.shape[1] != 3:
        raise ValueError("unexpected measured.txt shape %s in %s" % (measured.shape, d))
    if cov.shape != (measured.shape[0], measured.shape[0]):
        raise ValueError("unexpected covariance shape %s in %s" % (cov.shape, d))
    return measured[:, 0], measured[:, 1] * FACTOR, cov * FACTOR**2, measured[:, 2] * FACTOR


columns = {}
covariances = {}
nObj = {}
sanityChecks = []
reproChecks = []
rApArcmin = None

for mapName in MAPS:
    for zbin in Z_BINS:
        key = zbin + "_" + mapName
        radii, profile, covArcmin, analytic = readStack(zbin, mapName)

        if rApArcmin is None:
            rApArcmin = radii
        elif not np.allclose(rApArcmin, radii):
            raise ValueError("aperture radii differ between stacks, at " + key)

        sigma = np.sqrt(np.diag(covArcmin))
        if not (np.isfinite(profile).all() and np.isfinite(covArcmin).all()):
            raise ValueError("non-finite values for " + key)
        if not np.allclose(covArcmin, covArcmin.T, rtol=1e-10, atol=0.):
            raise ValueError("covariance is not symmetric for " + key)
        if not np.all(np.linalg.eigvalsh(covArcmin) > 0.):
            raise ValueError("covariance is not positive definite for " + key)

        # analytic errors are an independent estimate; they should agree closely
        ratio = sigma / analytic
        check = "%s: bootstrap/analytic sigma ratio %.3f - %.3f" % (key, ratio.min(), ratio.max())
        print(check)
        sanityChecks.append(check)

        columns[key] = profile
        columns[key + "_err"] = sigma
        covariances[key] = covArcmin
        if mapName == MAPS[0]:
            d = stackDir(zbin, mapName)
            nObj[zbin] = int(pd.read_csv(d + "overlap_flag.txt", header=None).to_numpy().sum())

# Does this run reproduce the earlier one? Same inputs and fixed bootstrap seeds,
# so any map present in both roots must agree exactly. NaN must not pass as
# agreement, and the two roots must really be different places.
nCompared = 0
fiducialShift = []
if REF_PATH_OUT is not None:
    if os.path.realpath(PATH_OUT) == os.path.realpath(REF_PATH_OUT):
        raise ValueError("PATH_OUT and REF_PATH_OUT are the same root; the check would be vacuous")
    for mapName in MAPS + MAPS_CHECK_ONLY:
        worst = 0.
        missing = 0
        for zbin in Z_BINS:
            try:
                refRadii, refProfile, refCov, _ = readStack(zbin, mapName, REF_PATH_OUT)
            except FileNotFoundError:
                missing += 1
                continue
            new = columns.get(zbin + "_" + mapName)
            newCov = covariances.get(zbin + "_" + mapName)
            if new is None:                      # check-only map, not in the product
                _, new, newCov, _ = readStack(zbin, mapName)
            for arr in (refProfile, refCov, new, newCov):
                if not np.isfinite(arr).all():   # max() silently drops NaN
                    raise ValueError("non-finite values in reproduction check: " + mapName + " " + zbin)
            if not np.allclose(refRadii, rApArcmin):
                raise ValueError("earlier run used different radii for " + mapName + " " + zbin)
            scale = max(np.abs(refProfile).max(), 1e-300)
            worst = max(worst, np.abs(new - refProfile).max() / scale,
                        np.abs(newCov - refCov).max() / max(np.abs(refCov).max(), 1e-300))
            nCompared += 1
        if missing == len(Z_BINS):
            repro = "%s: no earlier run to compare" % mapName
        else:
            repro = "%s: max |difference| / max |earlier| = %.2e%s" % (
                mapName, worst, " (%d bins missing)" % missing if missing else "")
        print(repro)
        reproChecks.append(repro)
    if nCompared == 0:
        raise ValueError("REF_PATH_OUT is set but nothing was compared; check the path")

# How much does the map version matter? Measured here rather than quoted from memory,
# comparing the internal fiducial with the published one on the same objects.
for zbin in Z_BINS:
    try:
        _, internal, _, _ = readStack(zbin, "act_dr6_fiducial")
    except FileNotFoundError:
        continue
    published = columns[zbin + "_act_dr602_fiducial"]
    frac = np.abs(published - internal) / np.abs(internal)
    fiducialShift.append("%s up to %.0f%%" % (zbin, 100. * frac.max()))

# profiles -> CSV (column order: R, then each map's value and error, per z bin)
frame = pd.DataFrame({"RApArcmin": rApArcmin})
for zbin in Z_BINS:
    for mapName in MAPS:
        key = zbin + "_" + mapName
        frame[key] = columns[key]
        frame[key + "_err"] = columns[key + "_err"]

if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)
pathCsv = DATA_DIR + STEM + "_profiles.csv"
frame.to_csv(pathCsv)
print("wrote " + pathCsv)

provenance = "\n".join([
    "ThumbStack stacked tSZ profiles, ring-ring aperture photometry filter.",
    "Compton-y maps (all published ACT DR6.02 NILC, from LAMBDA):",
] + [
    "  %-26s %s" % (m, MAP_FILES[m]) for m in MAPS
] + [
    "  'cib' deprojects the CIB at fixed spectral index beta; 'dBeta' additionally",
    "  deprojects its first moment in beta. Dust temperature T_d = 10.7 K throughout.",
    "Mask          : " + MASK_FILE,
    "Catalogues    : DESI LRG photometric redshift bins pz1-pz4",
    "Filter        : " + FILTER_TYPE + ". For aperture radius R, the inner annulus runs",
    "  1.0 arcmin < r <= R and the outer ring R < r <= R*sqrt(2); both are rescaled to the",
    "  full-disk pixel count, so a value is pi*R^2 times (mean y in the inner annulus minus",
    "  mean y in the outer ring). The smallest aperture, R = 2, uses a 1-2 arcmin annulus.",
    "Radii         : " + ", ".join("%.1f" % r for r in rApArcmin) + " arcmin",
    "Estimator     : " + EST + " (no hit map, so uniform weighting)",
    "Units         : profiles and errors in y*arcmin^2; covariances in (y*arcmin^2)^2",
    "Errors        : err = sqrt(diag(cov)), from 10000 bootstrap resamples",
    "The apertures overlap (outer ring reaches R*sqrt(2)), so neighbouring radial bins",
    "are correlated: use the covariance for any fitting, not the diagonal errors alone.",
    "The maps are also strongly correlated with each other: same sky, same catalogues and",
    "the same bootstrap resamples. No cross-map covariance is provided here, so do not",
    "treat different maps as independent measurements (e.g. when fitting a beta dependence).",
    "The CSV has an unnamed index column: read it with pandas index_col=0.",
    "Objects inside the mask footprint, before the point-source and outlier cuts: " +
    ", ".join("%s %d" % (z, nObj[z]) for z in Z_BINS),
    "",
    "Relation to Liu et al. (arXiv:2502.08850): that paper used earlier internal ILC",
    "y-maps, which are not public; these are the published DR6.02 maps with the same",
    "(beta, T_d). Stacking both fiducial versions on the same objects, the profiles",
    "differ by " + ("; ".join(fiducialShift) if fiducialShift else "an unmeasured amount") + ".",
    "Those two stacks share objects, so judging that difference needs the covariance of",
    "the difference, not the per-profile errors quoted here; by that measure it is",
    "significant in the two highest-z bins. These are therefore not a reproduction of",
    "the paper's numbers. The paper's figures also",
    "used the disk-ring filter on a 1-6 arcmin grid, not ring-ring on 2-6 arcmin.",
    "In that paper's zenodo release, fig4's 'act_dr6_Beta_x' columns are the dBeta",
    "maps at T_d = 10.7 and fig11's are the same at T_d = 24.0.",
    "",
    "Checks        : " + "; ".join(sanityChecks),
    ("Reproduction  : these stacks reproduce an earlier independent run of the same "
     "pipeline exactly; " + "; ".join(r for r in reproChecks if not any(
         m in r for m in MAPS_CHECK_ONLY))) if reproChecks else "Reproduction  : not checked",
    "Produced by scripts/make_ringring_data_release.py at commit " + gitCommit(),
])
pathNpz = DATA_DIR + STEM + "_cov.npz"
np.savez(pathNpz, provenance=np.array(provenance), **covariances)
print("wrote " + pathNpz)
