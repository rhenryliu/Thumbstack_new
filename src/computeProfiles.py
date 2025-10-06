from headers import *

##################################################################################
##################################################################################


def computeProfiles(ts, filterType, est='tsz_uniformweight', iBootstrap=None, iVShuffle=None, tTh='',  mVir=None, z=[0., 100.], mask=None):
    """Returns the estimated profile and its uncertainty for each aperture.
    est: string to select the estimator
    iBootstrap: index for bootstrap resampling
    iVShuffle: index for shuffling velocities
    tTh: to replace measured temperatures by a theory expectation
    ts: option to specify another thumbstack object
    """

    # tStart = time()

    print(("- Compute stacked profile: "+filterType+", "+est+", "+tTh))

    # compute stacked profile from another thumbstack object
    if ts is None:
        ts = self
    if mVir is None:
        mVir = [ts.mMin, ts.mMax]

    # select objects that overlap, and reject point sources
    if mask is None:
        mask = ts.catalogMask(overlap=True, psMask=True,
                              filterType=filterType, mVir=mVir, z=z)

#      tMean = ts.meanT[filterType].copy()

    # temperatures [muK * sr]
    if tTh == '':
        t = ts.filtMap[filterType].copy()  # [muK * sr]
    elif tTh == 'tsz':
        # expected tSZ signal
        # AP profile shape, between 0 and 1
        sigma_cluster = 3.  # 1.5  # arcmin
        shape = ts.ftheoryGaussianProfile(
            sigma_cluster)  # between 0 and 1 [dimless]
        # multiply by integrated y to get y profile [sr]
        t = np.column_stack(
            [ts.Catalog.integratedY[:] * shape[iAp] for iAp in range(ts.nRAp)])
        # convert from y profile to dT profile if needed
        if ts.cmbUnitLatex == r'$\mu$K':
            nu = ts.cmbNu   # Hz
            Tcmb = 2.726   # K
            h = 6.63e-34   # SI
            kB = 1.38e-23  # SI

            def f(nu):
                """frequency dependence for tSZ temperature
                """
                x = h*nu/(kB*Tcmb)
                return x*(np.exp(x)+1.)/(np.exp(x)-1.) - 4.
            # t *= 2. * f(nu) * Tcmb * 1.e6  # [muK * sr]
            t *= f(nu) * Tcmb * 1.e6  # [muK * sr]
    t = t[mask, :]
#     tMean = tMean[mask,:]
    # -v/c [dimless]
    v = -ts.Catalog.vR[mask] / 3.e5
    v -= np.mean(v)

    # true filter variance for each object and aperture,
    # valid whether or not a hit count map is available
    s2Full = ts.filtVarTrue[filterType][mask, :]
    # Variance from hit count (if available)
    s2Hit = ts.filtHitNoiseStdDev[filterType][mask, :]**2
    # print "Shape of s2Hit = ", s2Hit.shape
    # halo masses
    m = ts.Catalog.Mvir[mask]

    if iBootstrap is not None:
        # make sure each resample is independent,
        # and make the resampling reproducible
        np.random.seed(iBootstrap)
        # list of overlapping objects
        nObj = np.sum(mask)
        # print "sample "iBootstrap, ";", nObj, "objects overlap with", ts.name
        I = np.arange(nObj)
        # choose with replacement from this list
        J = np.random.choice(I, size=nObj, replace=True)
        #
        t = t[J, :]
        # tMean = tMean[J,:]
        v = v[J]
        s2Hit = s2Hit[J, :]
        s2Full = s2Full[J, :]
        m = m[J]

    if iVShuffle is not None:
        # make sure each shuffling is independent,
        # and make the shuffling reproducible
        np.random.seed(iVShuffle)
        # list of overlapping objects
        nObj = np.sum(mask)
        I = np.arange(nObj)
        # shuffle the velocities
        J = np.random.permutation(I)
        #
        v = v[J]

    # tSZ: uniform weighting
    if est == 'tsz_uniformweight':
        weights = np.ones_like(s2Hit)
        norm = 1./np.sum(weights, axis=0)
    # tSZ: detector-noise weighted (hit count)
    elif est == 'tsz_hitweight':
        weights = 1./s2Hit
        norm = 1./np.sum(weights, axis=0)
    # tSZ: full noise weighted (detector noise + CMB)
    elif est == 'tsz_varweight':
        weights = 1./s2Full
        norm = 1./np.sum(weights, axis=0)


    # tStop = time()
    # print "stacked profile took", tStop-tStart, "sec"

    # return the stacked profiles
    # stack = norm * np.sum(t * weights, axis=0)
    # sStack = norm * np.sqrt(np.sum(s2Full * weights**2, axis=0))
    return t