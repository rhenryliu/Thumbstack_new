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

from computeProfiles import computeProfiles
from astropy.table import Table

# Get JSON Parameters:
import json
import pprint

param_Dict = json.loads(sys.argv[1])
locals().update(param_Dict)
print('Parameters:')
pprint.pprint(param_Dict)
'''
Params:
catalog_name = catalog name for saving
save = true or false
'''

# Make Catalog


# catalog_path = '/global/cfs/cdirs/desicollab/science/c3/DESI-Lensing/desi_catalogues/v1.5/BGS_BRIGHT_clustering.dat.fits'
catalog_path = '/global/cfs/cdirs/desicollab/science/c3/DESI-Lensing/desi_catalogues/v1.5/' + catalog_name + '.fits'

cat = Table.read(catalog_path)


massConversion = MassConversionKravtsov14()
# Converting a flat virial mass to Stellar Mass, this is what we're using for every object for now.
MStellar = massConversion.fmVirTomStar(2e13)
u = UnivMariana()

try:
    df = cat.to_pandas()
except:
    cat.remove_columns(['BITWEIGHTS'])
    df = cat.to_pandas()

df2 = df.loc[:, ('TARGETID','RA', 'DEC', 'Z')]

# colnames = ['coordX', 'coordY', 'coordZ', 'dX', 'dY', 'dZ', 'dXKaiser', 'dYKaiser', 'dZKaiser',
#             'vX', 'vY', 'vZ', 'vR', 'vTheta', 'vPhi']

# for col in colnames:
#     df2[col]=0

df2['Mstellar'] = MStellar
df2['Mvir'] = massConversion.fmStarTomVir(MStellar)
# df2['integratedKSZ'] = 0
# addHaloMass(df2, u, massConversion)
# addIntegratedY(df2, u)
# df2.to_csv(r'/pscratch/sd/r/rhliu/projects/ThumbStack/catalogs/tempcatalog.txt', header=None, index=None, sep=' ', mode='w')

CMB_path = "/pscratch/sd/r/rhliu/projects/ThumbStack/ACT_DR6/" + "ilc_SZ_yy.fits"
CMB_mask = '/pscratch/sd/r/rhliu/projects/ThumbStack/ACT_DR6/wide_mask_GAL070_apod_1.50_deg_wExtended.fits'
CMB_mask = '/pscratch/sd/r/rhliu/projects/ThumbStack/ACT_DR6/wide_mask_GAL070_apod_1.50_deg_wExtended_srcfree_Will.fits'
CMB_convert = False
CMB_name = 'act_dr6_fiducial'
CMB_namepublic = 'ACT DR6 (fiducial)'
CMB_nu = 93.e9


cmap = cmbMap(CMB_path,
              pathMask=CMB_mask,
              pathHit=None,
              nu=CMB_nu, unitLatex=r'y', 
              convert_y=CMB_convert,
              name=CMB_name)

nProc = 64
# save = True
# save = False
filterType = 'diskring'

len_df = len(df)
# catalog_name = 'BGS_BRIGHT_clustering.dat'
df3 = df2

# catalog_name = 'BGS_BRIGHT_clustering_part1'
# df3 = df2.iloc[:int(len_df/2)]

# catalog_name = 'BGS_BRIGHT_clustering_part2'
# df3 = df2.iloc[int(len_df/2):]

# print(df2.shape)
print(df3.shape)

stack_catalogue = make_Catalog(u, massConversion, df3, catalog_name)

ts = ThumbStack(u, stack_catalogue, 
                cmap.map(), 
                cmap.mask(), 
                cmap.hit(), 
                catalog_name + '_' + cmap.name,
                nameLong=None, 
                save=save, 
                nProc=nProc,
                filterTypes=filterType,
                doMBins=False, 
                doBootstrap=True,
                # doStackedMap=True,
                doVShuffle=False, 
                cmbNu=cmap.nu, 
                cmbUnitLatex=cmap.unitLatex,
                rApMinArcmin=1.,
                rApMaxArcmin=30.,
                nRAp=256,)


# Now to save the catalogues (as a pd dataframe, combine later for fits table file)
allProfiles = computeProfiles(ts, 'diskring')
mask = ts.catalogMask(overlap=True, psMask=True, filterType='diskring', mVir=None)

df3['Mask'] = mask.astype(int)
NObj = df3.shape[0]

for i, R in enumerate(ts.RApArcmin):
    print(str(R))
    stacks = np.zeros(NObj)
    stacks[mask] = allProfiles[:, i]

    df3[str(R)] = stacks

df_save = df3.drop(columns=['Mstellar', 'Mvir'])
# df_save.to_csv(r'/pscratch/sd/r/rhliu/projects/ThumbStack/catalogs/' + catalog_name + '.csv', index=None, mode='w')
catalog_fits = Table.from_pandas(df_save)
catalog_fits.write('/pscratch/sd/r/rhliu/projects/ThumbStack/catalogs/for_Sven/' + catalog_name + '_new_tSZ.fits', format='fits', overwrite=True)

print('Done!!!')