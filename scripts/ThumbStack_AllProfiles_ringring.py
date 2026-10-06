import sys
sys.path.append('../src/')


from importlib import reload
import universe
reload(universe)
from universe import *

import mass_conversion
reload(mass_conversion)
from mass_conversion import *

import catalog
reload(catalog)
from catalog import *

import thumbstack
reload(thumbstack)
from thumbstack import *

import cmb
reload(cmb)
from cmb import *
# from headers import *
from cmbMap import *
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

import json
import pprint

##################################################################################
# Parameters

filterType = 'ringring'
T_CIB = '10.7'
# T_CIB = '24.0'

# Aperture radii [arcmin]. For ringring, rApMinArcmin must be > rApInnerRad,
# otherwise the inner ring of the first aperture is empty and the filter is NaN.
rApMinArcmin = 2.
rApMaxArcmin = 6.
nRAp = 9
rApInnerRad = 1.

# Optional JSON overrides (for smoke tests), e.g.
# python ThumbStack_AllProfiles_ringring.py '{"nobj": 20000, "suffix": "_smoketest",
#     "maps": ["act_dr6_fiducial"], "catalogs": ["DESI_pz1"]}'
param_Dict = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
print('Parameters:')
pprint.pprint(param_Dict)
nObj = param_Dict.get('nobj', None)            # keep only the first nObj objects per catalog
SUFFIX = param_Dict.get('suffix', '')          # appended to the ThumbStack output names
mapSubset = param_Dict.get('maps', None)       # list of CMB map names to keep
catalogSubset = param_Dict.get('catalogs', None)  # list of catalog names to keep

# Output root; ThumbStack writes to <pathOut>output/thumbstack/<name>/. A new root
# keeps a run from overwriting an earlier one.
pathOut = param_Dict.get('path_out', "/pscratch/sd/r/rhliu/projects/ThumbStack_ringring_v2/")
if not pathOut.endswith('/'):
    pathOut += '/'  # ThumbStack builds its path as pathOut + "output/thumbstack/"
if pathOut == '/' or os.path.realpath(pathOut) == os.path.realpath(
        "/pscratch/sd/r/rhliu/projects/ThumbStack/"):
    sys.exit("path_out would overwrite an earlier run's root (or is empty): " + pathOut)

# ACT DR6 maps and masks (scratch; backed up on CFS m4031)
MAP_DIR = "/pscratch/sd/r/rhliu/projects/ThumbStack/ACT_DR6/"
pathMask = MAP_DIR + 'wide_mask_GAL070_apod_1.50_deg_wExtended.fits'
pathMask2 = MAP_DIR + 'wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'

# Old ThumbStack repo holds the DESI_pz catalogs (not present in ThumbStack_new)
pathCatalogs = "/global/homes/r/rhliu/projects/ThumbStack/output/catalog/"
pathCatalogs = param_Dict.get('catalog_dir', pathCatalogs)  # e.g. a subsampled copy for smoke tests
if not pathCatalogs.endswith('/'):
    pathCatalogs += '/'  # Catalog builds its path as pathOut + name

plot_Path = ("./figures/ThumbStack_AllPlots_" + filterType + "_dbeta_"
             + T_CIB.replace('.', '') + "_v2" + SUFFIX + ".pdf")
# distinct path per job, so concurrent jobs don't race on one file
plot_Path = param_Dict.get('plot_path', plot_Path)

print('Resolved paths:')
print('  output root  : ' + pathOut)
print('  catalogues   : ' + pathCatalogs)
print('  summary plot : ' + plot_Path)

##################################################################################

nProc = 64  # 1 haswell node on cori
##################################################################################
# CMB maps: (path, mask, convert_y, name, public name, nu)
# 'act_dr602_*' are the published ACT DR6.02 NILC maps from LAMBDA
# (see act_dr6.02_nilc_wget.sh). 'act_dr6_fiducial' is the older internal
# ILC y-map; it is a different map version from the DR6.02 fiducial
# (4' block correlation 0.96, and a half-pixel Dec offset in the CAR grid).
# 'cib' deprojects the CIB at a fixed spectral index beta; 'dBeta' additionally
# deprojects its first moment in beta. Both at dust temperature T_CIB.

CMB_params = [(MAP_DIR + "ilc_SZ_yy.fits",
               pathMask2, False, 'act_dr6_fiducial', 'ACT DR6 (fiducial, internal)', 93.e9),
              (MAP_DIR + "act-planck_dr6.02_nilc_ComptonY.fits",
               pathMask2, False, 'act_dr602_fiducial', 'ACT DR6.02 (fiducial)', 93.e9),
              (MAP_DIR + "act-planck_dr6.02_nilc_ComptonY_deproj_cib_1.2_"+T_CIB+".fits",
               pathMask2, False, 'act_dr602_cib_1.2_'+T_CIB, r'ACT DR6.02 (CIB $\beta$ 1.2)', 93.e9),
              (MAP_DIR + "act-planck_dr6.02_nilc_ComptonY_deproj_cib_1.4_"+T_CIB+".fits",
               pathMask2, False, 'act_dr602_cib_1.4_'+T_CIB, r'ACT DR6.02 (CIB $\beta$ 1.4)', 93.e9),
              (MAP_DIR + "act-planck_dr6.02_nilc_ComptonY_deproj_cib_1.6_"+T_CIB+".fits",
               pathMask2, False, 'act_dr602_cib_1.6_'+T_CIB, r'ACT DR6.02 (CIB $\beta$ 1.6)', 93.e9),
              (MAP_DIR + "act-planck_dr6.02_nilc_ComptonY_deproj_cib_cibdBeta_1.2_"+T_CIB+".fits",
               pathMask2, False, 'act_dr602_dBeta_1.2_'+T_CIB, r'ACT DR6.02 (d$\beta$ 1.2)', 93.e9),
              (MAP_DIR + "act-planck_dr6.02_nilc_ComptonY_deproj_cib_cibdBeta_1.4_"+T_CIB+".fits",
               pathMask2, False, 'act_dr602_dBeta_1.4_'+T_CIB, r'ACT DR6.02 (d$\beta$ 1.4)', 93.e9),
              (MAP_DIR + "act-planck_dr6.02_nilc_ComptonY_deproj_cib_cibdBeta_1.6_"+T_CIB+".fits",
               pathMask2, False, 'act_dr602_dBeta_1.6_'+T_CIB, r'ACT DR6.02 (d$\beta$ 1.6)', 93.e9)
             ]
if mapSubset is not None:
    CMB_params = [element for element in CMB_params if element[3] in mapSubset]
    if len(CMB_params) != len(mapSubset):
        sys.exit("Unknown map name in 'maps': " + str(mapSubset))

##################################################################################
# Pre-flight check: fail fast before the slow catalog load

# full-sky 0.5' CAR, 10320 x 43200 pixels: float32 maps, float64 masks (+ FITS header)
expectedSize = {'map': 1783298880, 'mask': 3566594880}  # bytes
missing = []
for element in CMB_params:
    for path, kind in [(element[0], 'map'), (element[1], 'mask')]:
        if not os.path.exists(path):
            missing.append(path + " (missing)")
        elif os.path.getsize(path) != expectedSize[kind]:
            missing.append(path + " (size " + str(os.path.getsize(path)) + ", expected " + str(expectedSize[kind]) + ")")
if len(missing) > 0:
    sys.exit("Pre-flight check failed:\n" + "\n".join(sorted(set(missing))))

##################################################################################
##################################################################################

# cosmological parameters
u = UnivMariana()

# M*-Mh relation
massConversion = MassConversionKravtsov14()
# massConversion.plot()

##################################################################################
# Galaxy Catalogs (from DESI)

print("Read galaxy catalogs")
tStart = time()

catalogNames = ["DESI_pz1", "DESI_pz2", "DESI_pz3", "DESI_pz4"]
if catalogSubset is not None:
    catalogNames = [name for name in catalogNames if name in catalogSubset]
    if len(catalogNames) != len(catalogSubset):
        sys.exit("Unknown catalog name in 'catalogs': " + str(catalogSubset))
catalogs = {}
for i, name in enumerate(catalogNames):
    catalogs[name] = Catalog(u, massConversion, name=name, nameLong="DESI pz bin " + name[-1],
                             save=False, nObj=nObj, pathOut=pathCatalogs)

tStop = time()
print("took "+str(round((tStop-tStart)/60., 2))+" min")

###################################################################################
# Read CMB maps

CMB_pathlist = []
CMB_masklist = []
CMB_convert = []
CMB_name = []
CMB_namepublic = []
CMB_nu = []
for element in CMB_params:
    CMB_pathlist.append(element[0])
    CMB_masklist.append(element[1])
    CMB_convert.append(element[2])
    CMB_name.append(element[3])
    CMB_namepublic.append(element[4])
    CMB_nu.append(element[5])


filterTypes = [filterType] * len(CMB_params)

cmbMap_list = []

for i, path in enumerate(CMB_pathlist):
    cmap = cmbMap(path,
                  pathMask=CMB_masklist[i],
                  pathHit=None,
                  nu=CMB_nu[i], unitLatex=r'y', convert_y=CMB_convert[i],
                  name=CMB_name[i])
    cmbMap_list.append(cmap)

catalogKeys = catalogs.keys()

###################################################################################
# Stacking
save = True
# save = False
ts_list = [[] for _ in range(len(CMB_masklist))]


for key in list(catalogKeys):
    catalog = catalogs[key]

    for i, cmap in enumerate(cmbMap_list):

        tStart = time()
        ts = ThumbStack(u, catalog,
                        cmap.map(),
                        cmap.mask(),
                        cmap.hit(),
                        catalog.name + '_' + cmap.name + '_' + filterTypes[i] + SUFFIX,
                        nameLong=None,
                        save=save,
                        nProc=nProc,
                        filterTypes=filterTypes[i],
                        doMBins=False,
                        doBootstrap=True,
                        # doStackedMap=True,
                        doVShuffle=False,
                        cmbNu=cmap.nu,
                        cmbUnitLatex=cmap.unitLatex,
                        pathOut=pathOut,
                        rApMinArcmin=rApMinArcmin,
                        rApMaxArcmin=rApMaxArcmin,
                        nRAp=nRAp,
                        rApInnerRad=rApInnerRad)
        ts_list[i].append(ts)
        tStop = time()
        print("stack " + ts.name + " took " + str(round((tStop-tStart)/60., 2)) + " min")

###################################################################################

# Next for plotting:


# Parameters
factor = (180.*60./np.pi)**2
est = 'tsz_uniformweight'

# Plotting


###############################
fig, subplots = plt.subplots(2, 2, figsize=(8,8), sharex='col', sharey='row')
subplots = subplots.ravel()

for i, key in enumerate(list(catalogKeys)):
    ax = subplots[i]
    for j in range(len(ts_list)):
        tsj = ts_list[j][i]
        filterType = filterTypes[j]

        ax.errorbar(tsj.RApArcmin, factor * tsj.stackedProfile[filterType+"_"+est],
                    factor * tsj.sStackedProfile[filterType+"_"+est],
                    label=CMB_namepublic[j], lw=2, ls='-.', capsize=6)

    ax.set_title(key)
    ax.grid()
    if i>=2:
        ax.set_xlabel(r'$R$ [arcmin]')
    if i==0 or i==2:
        ax.set_ylabel(r'Compton Y-parameter $[\mathrm{arcmin}^2]$')

ax.legend(fontsize=10, labelspacing=0.1)
plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()
os.makedirs(os.path.dirname(plot_Path) or '.', exist_ok=True)
fig.savefig(plot_Path, dpi=100) # bbox_inches='tight')

print('Done!!!')
