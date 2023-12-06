## Chapter 4.1.2
```Python
with pm.Model() as model:
  sigma = 3.3
  #Prior Distribution
  pm.Uniform('m', lower = 0, upper = 33)
  #Data Model
  pm.Normal('obsm', mu = model['m'], sigma = sigma, observed = obsm)
```

## Chapter 4.2.1
```Python
with pm.Model() as model:
    pm.Uniform("s", lower = 0, upper = 1.0e7)
    pm.Poisson("obsn", mu = model['s'], observed = obsn)
```

## Chapter 4.2.2
```Python
with pm.Model() as model:
    #uniform prior on s
    pm.Uniform("s", lower = 0, upper = 1.0e7)
    pm.Poisson("obsn", mu = model['s'], observed = obsn)
```

## Chapter 4.2.3
```Python
with pm.Model() as model:
    pm.Beta("f", alpha = 1, beta = 1)
    pm.Binomial("obsn", n = n, p = model['f'], observed = obsn)
```

## Chapter 5.1.2
```Python
with pm.Model() as model:
    pm.LogNormal("s", mu = 4.605, sigma = 10)
    for i in range(0, len(obsn)):
        pm.Poisson("obsn[" + str(i) + "]", mu = model['s'], observed = obsn)
```

## Chapter 5.1.4
```Python
with pm.Model() as model:
    #Gauss prior with large sigma
    pm.TruncatedNormal("s", mu = 0 , sigma = 10, lower = 0)
    #Uniform
    pm.Uniform("s2", lower = 0, upper = 1.0e7)
    #Gauss likelihood
    pm.Poisson("obsn", mu = model['s'], observed = obsn)
    #Uniform likelihood
    pm.Poisson("obsn2", mu = model['s2'], observed = obsn)
```

## Chapter 5.4.1
```Python
#we use zero trick, because alpha < -1 makes gamma function undefined
zeros = 0
C = 10

with pm.Model() as model:
  #prior, using the zero trick
  pm.Uniform("lgM", lower = 12, upper = 18)
  phi = 2.3025*0.4*(alpha + 1)*(model["lgM"] - lgMstar) + 10**(0.4*(model["lgM"] - lgMstar))
  pm.Poisson("zeros", mu = phi + C, observed = zeros)
  #likelihood
  pm.Normal("obslgM", mu = model["lgM"], sigma = 0.3, observed = obslgM)
```

## Chapter 5.4.2
```Python
with pm.Model() as model:
    pm.Uniform("zphot", lower = 0, upper = 2.0)
    pm.Normal("zspec", mu = model['zphot'], sigma = (0.03*(1 + model['zphot'])**3))
```

## Chapter 6.1.1
```Python
with pm.Model() as model:
  pm.Normal("vcent", mu = 10000, sigma = 1000)
  pm.Uniform("sigma_clus", lower = 0, upper = 5000)

  for i in range(0, len(obsv)):
      pm.Normal("v[" + str(i) + "]", mu = model['vcent'], sigma = model['sigma_clus'])
      pm.Normal("obsv[" + str(i) + "]", mu = model['v[' + str(i) + ']'], sigma = sigma_v[i], observed = obsv[i])
```

## Chapter 6.1.2
```Python
with pm.Model() as model:
  pm.Uniform("s", lower = 0, upper = 1.0e7)
  pm.Uniform("bkg", lower = 0, upper = 1.0e7)
  pm.Poisson("obstot", mu = (model['s'] + model['bkg']/C), observed = obstot)
  pm.Poisson("obsbkg", mu = model['bkg'], observed = obsbkg)
```

## Chapter 6.1.3
```Python
with pm.Model() as model:
    pm.Uniform("nbkg", lower = 1, upper = 1e7)
    pm.Uniform("nclus", lower = 1, upper = 1e7)
    pm.Beta("fbkg", alpha = 1, beta = 1)
    pm.Beta("fclus", alpha = 1, beta = 1)
    f = (model['fbkg']*model['nbkg']/C + model['fclus']*model['nclus'])/(model['nbkg']/C + model['nclus'])

    obsnbkg = pm.Poisson("obsnbkg", mu = model['nbkg'], observed = obsnbkg)
    obsntot = pm.Poisson("obsntot", mu = (model['nbkg']/C + model['nclus']), observed = obsntot)
    obsbluebkg = pm.Binomial("obsbluebkg", p = model['fbkg'], n = model['obsnbkg'], observed = obsbluebkg)
    obsnbluetot = pm.Binomial("obsnbluetot", p = f, n = model['obsntot'], observed = obsbluetot)
```

## Chapter 6.1.4
```Python
with pm.Model() as model:
    pm.Uniform("S", lower = 0, upper = 1.0e7)
    pm.Uniform("HR", lower = -1, upper = 1)
    H = 2*model['S']/(1 - model['HR']) - model['S']
    pm.Uniform("bkgH", lower = 0, upper = 1.0e7)
    pm.Uniform("bkgS", lower = 0, upper = 1.0e7)

    pm.Poisson("obsbkgH", mu = model['bkgH'], observed = obsbkgH)
    pm.Poisson("obsbkgS", mu = model['bkgS'], observed = obsbkgS)
    pm.Poisson("obstotH", mu = (H + model['bkgH']/CH) observed = obstotH)
    pm.Poisson("obstotS", mu = (model['S'] + model['bkgS']/CS), observed = obstotS)
```

## Chapter 6.1.5
```Python
with pm.Model() as model:
    alpha = [1, 1, 1]

    tot = pm.Uniform("tot", lower = 0, upper = 1.0e7)

    p = pm.Dirichlet("p", alpha)
    S = p[0]*tot
    M = p[1]*tot
    H = p[2]*tot

    pm.Deterministic('HR1', (M - S)/(H + M + S))
    pm.Deterministic('HR2', (H - M)/(H + M + S))

    pm.Uniform("pkgS", upper = 0, lower = 1.0e7)
    pm.Uniform("bkgM", upper = 0, lower = 1.0e7)
    pm.Uniform("bkgH", upper = 0, lower = 1.0e7)

    pm.Poisson("obsbkgS", mu = model['bkgS'], observed = obsbkgS)
    pm.Poisson("obsbkgM", mu = model['bkgM'], observed = obsbkgM)
    pm.Poisson("obsbkgH", mu = model['bkgH'], observed = obsbkgH)
    pm.Poisson("obstotS", mu = (S + model['bkgS']/CS), observed = obstotS)
    pm.Poisson("obstotM", mu = (M + model['bkgM']/CM), observed = obstotM)
    pm.Poisson("obstotH", mu = (H + model['bkgH']/CH), observed = obstotH)
```

## Chapter 6.3.1
```Python
with pm.Model() as model:

  for i in range(0, len(obstot)):
    pm.Uniform("nbkg[" + str(i) + "]", lower = 1, upper = 1e7)
    pm.Uniform("nclus[" + str(i) + "]", lower = 0, upper = 1e7)
    pm.LogNormal("nbkgind[" + str(i) + "]", mu = np.log(model['nbkg[' + str(i) + ']']), sigma = 0.2)
    pm.Poisson("obsbkg[" + str(i) + "]", mu = model['nbkg[' + str(i) + ']'], observed = obsbkg[i])
    pm.Poisson("obstot[" + str(i) + "]", mu = (model['nclus[' + str(i) + ']'] + model['nbkgind[' + str(i) + ']']/nbox[i]), observed = obstot[i])

    #for plotting pourposes only
    #pm.Deterministic("lgLx[" + str(i) + "]", np.log(model['nclus[' + str(i) + ']'])/2.30258 + C[i])
```

## Chapter 6.3.2
```Python
with pm.Model() as model:
  pm.Normal("H0", mu = 70.1, sigma = 1.3)
  pm.Uniform("omegam", lower = 0, upper = 1)
  pm.Normal("omegabh2", mu = 0.0462, sigma = 0.0012)
  omegab = model['omegabh2']*(model['H0']/70)**-2
  fb = omegab/model['omegam']
  pm.Normal("corr", mu = 0.874, sigma = 0.023)
  fbcorr = fb*model['corr']
  pm.Uniform("intr_scatt", lower = 0.0, upper = 0.1)
  pm.Normal("fcoldonfgas", mu = 0.18, sigma = 0.05)

  for i in range(0, len(obsfb)):
    pm.Normal("fbclus[" + str(i) + "]", mu = fbcorr, sigma = model['intr_scatt'])
    pm.Normal("obsfgas[" + str(i) + "]", mu = (model["fbclus[" + str(i) + "]"]*(1 - model['fcoldonfgas'])), sigma = errfb[i], observed = obsfb[i])
```

## Chapter 6.3.3
```Python
with pm.Model() as model:
    for i in range(0, len(obstot1)):
        pm.Uniform("nbkg[" + str(i) + "]", lower = 0, upper = 1e7)
        pm.Uniform("nclus2[" + str(i) + "]", lower = 0, upper = 1e7) # taken indep for the time being
        pm.Uniform("nclus1[" + str(i) + "]", lower = model['nclus2[' + str(i) + ']']/24, upper = 1e7)
        pm.LogNormal("nbkgind[" + str(i) + "]", mu = np.log(model['nbkg[' + str(i) + ']']), sigma = 1/np.sqrt(1/0.2/0.2))
        pm.Poisson("obsbkg[" + str(i) + "]", mu = model['nbkg[' + str(i) + ']'], observed = obsbkg[i])
        # aperture 1 (small)
        pm.Poisson("obstot1[" + str(i) + "]", mu = (model['nclus1[' + str(i) + ']'] + model['nbkgind[' + str(i) + ']']/nbox[i]/25), observed = obstot1[i])
        # aperture 2 (large minus small)
        pm.Poisson("obstot2[" + str(i) + "]", mu = (model['nclus2[' + str(i) + ']'] + model['nbkgind[' + str(i) + ']']/nbox[i]*24/25), observed = obstot2[i])
```

## Chapter 6.3.4.2
```Python
with pm.Model() as model:
    pm.TruncatedNormal("bkg_lihe", mu = 0.03, sigma = 0.02, lower = 0)
    pm.Exponential("bkg_fastn", lam = 200)
    pm.Exponential("bkg_fastn2", lam = 58)
    pm.TruncatedNormal("bkg_untagmu", mu = 0.011, sigma = 0.001, lower = 0)
    pm.TruncatedNormal("bkg_acccoin", mu = 0.080, sigma = 0.001, lower = 0)
    pm.Exponential("bkg_timcor", lam = 88)
    pm.Exponential("bkg_fot", lam = 765)
    pm.TruncatedNormal("bkg_spontfi", mu = 0.003, sigma = 0.0003, lower = 0)
    pm.TruncatedNormal("bkg_scin", mu = 0.021, sigma = 0.002, lower = 0)
    pm.Exponential("bkg_buff", lam = 37)
    pm.Deterministic("bkg", (model['bkg_lihe'] + model['bkg_fastn'] + model['bkg_fastn2'] + model['bkg_untagmu'] + model['bkg_acccoin'] + model['bkg_timcor'] + model['bkg_fot'] + model['bkg_spontfi'] + model['bkg_scin'] + model['bkg_buff']))

    pm.Normal("bkg_react1", mu = 5.0, sigma = 0.3)
    pm.Normal("bkg_react2", mu = 9.4, sigma = 0.6)
    bkg_tot = 2.526*bkg + bkg_react1 + bkg_react2

    pm.Uniform("s", lower = 0, upper = 100)
    pm.Poisson("obstot", mu = (model['s'] + model['bkg_tot']), observed = obstot)
```

## Chapter 6.3.4.3
```Python
with pm.Model() as model:
    pm.Uniform("s", lower = 0, upper = 100)
    pm.Normal("bkg", mu = 0.31, sigma = 0.05)
    pm.Normal("bkg_react", mu = 5.0, sigma = 0.3)
    bkg_tot = model['bkg'] + model['bkg_react']

    pm.Poisson("obstot", mu = (model['s'] + model['bkg_tot']), observed = obstot)
```

## Chapter 6.3.5.1
```Python
with pm.Model() as model:
    pm.Uniform("s", lower = 0, upper = 10)
    pm.Uniform("bkg", lower = 0, upper = 10)

    pm.Poisson("obstot", mu = (model['s'] + model['bkg']), observed = obstot)
```

## Chapter 6.3.5.2
```Python
with pm.Model() as model:
    pm.Uniform("s", lower = 0, upper = 20)
    pm.Exponential("bkg", lam = 1.45)
    pm.Poisson("obstot", mu = (model['s'] + model['bkg']), observed = obstot)
```

## Chapter 7.2
```Python
with pm.Model() as model:
    sigma_intr = 3
    pm.Uniform("xcent", lower = -10, upper = 10)

    for i in range(0, len(x)):
        pm.TruncatedNormal("x[" + str(i) + "]", mu = model['xcent'], sigma = sigma_intr, lower = -1, observed = x[i])

with pm.Model() as model:
    pm.Uniform("sigma_intr", lower = 0, upper = 10)
    xcent = 0

    for i in range(0, len(data['x'])):
        pm.TruncatedNormal("x[" + str(i) + "]", mu = xcent, sigma = model['sigma_intr'], lower = -1, observed = x[i])
```

## Chapter 7.4
```Python
with pm.Model() as model:
    sigma_intr = 3
    pm.Uniform("xcent", lower = -10, upper = 10)

    for i in range(0, len(obsx)):
        pm.TruncatedNormal("x[" + str(i) + "]", mu = model['xcent'], sigma = sigma_intr, lower = 3)

        pm.Normal("obsx[" + str(i) + "]", mu = model['x[' + str(i) + ']'], sigma = err[i], observed = obsx[i])
```

## Chapter 8.1.3.1
```Python
with pm.Model():
    covar = [[1, 0.5], [0.5, 1]]
    mu = [0, 0]

    val = pm.MvNormal("val", mu = mu, cov = covar)
```

## Chapter 8.2
```Python
with pm.Model() as model:
  eff = [None]*len(nrec)

  pm.Uniform("A", lower = 0, upper = 1)
  pm.Uniform("B", lower = 0, upper = 1)
  pm.Uniform("mu", lower = 0, upper = 100)
  pm.Uniform("sigma", lower = 0, upper = 100)

  for i in range (0, len(nrec)):
    eff[i] = model['A'] + (model['B'] - model['A'])*0.5*(1 + pm.math.erf(((E[i] - model['mu'])/model['sigma']) / pm.math.sqrt(2)))

    pm.Binomial("nrec[" + str(i) + "]", p = eff[i], n = ninj[i], observed = nrec[i])
    pm.Binomial("nrec_rep[" + str(i) + "]", p = eff[i], n = ninj[i], observed = nrec[i])
```

## Chapter 8.3
```Python
with pm.Model() as model:
    z = [None]*len(data['obsy'])

    pm.Uniform("intrscat", lower = 0, upper = 3)
    pm.Normal("alpha", mu = 0.0, sigma = 100)
    pm.StudentT("beta", mu = 0, nu = 1, sigma = 1)

    for i in range (0, len(k)):
        #modeling ordinate vs x
        z[i] = model['alpha'] + 0.1 + model['beta']*(k[i] - 0.03)
        #modeling ordinate
        pm.Normal("y[" + str(i) + "]", mu = z[i], sigma = model['intrscat'])
        pm.Normal("obsy[" + str(i) + "]", mu = model['y[' + str(i) + ']'], sigma = errobsy[i], observed = obsy[i])
```

## Chapter 8.4
```Python
with pm.Model() as model:
  #reduced range to improve convergence
  pm.StudentT("a", mu = 0, nu = 1, sigma = 1)
  pm.Normal("b", mu = 8.0, sigma = 10)
  #simplified prior
  pm.Uniform("intrscat", lower = 0, upper = 10)

  for i in range (0, len(obsx)):
    pm.Uniform("x[" + str(i) + "]", lower = 1, upper = 4)

    pm.Normal("y[" + str(i) + "]", mu = (model['b'] + model['a']*(model['x[' + str(i) + ']'] - 2.3)), sigma = model['intrscat'])
    pm.Normal("obsx[" + str(i) + "]", mu = model['x[' + str(i) + ']'], sigma = errx[i], observed = obsx[i])
    pm.Normal("obsy[" + str(i) + "]", mu = model['y[' + str(i) + ']'], sigma = erry[i], observed = obsy[i])
```

## Chapter 8.5
```Python
with pm.Model() as model:
  f = [None]*len(obsntot)
  fclus = [None]*len(obsntot)

  pm.Normal("alpha", mu = 0, sigma = 10)
  pm.Normal("beta", mu = 0, sigma = 10)
  pm.Normal("gamma", mu = 0, sigma = 10)
  pm.Normal("zeta", mu = 0, sigma = 10)
  pm.Normal("lgfclus0", mu = 0, sigma = 10)

  for i in range (0, len(obsntot)):
    pm.Uniform("nbkg[" + str(i) + "]", lower = 1, upper = 10000)
    pm.Beta("fbkg[" + str(i) + "]", alpha = 1, beta = 1)
    pm.Uniform("nclus[" + str(i) + "]", lower = 1, upper = 1000)

    fclus[i] = pm.invlogit(model['lgfclus0'] + model['alpha']*np.log(r200[i]/0.25) + model['beta']*(lgM[i] - 11) + model['gamma']*(z[i] - 0.3) + model['zeta']*(lgM[i] - 11)*(z[i] - 0.3))

    f[i] = (model['fbkg[' + str(i) + ']']*model['nbkg[' + str(i) + ']']/C[i] + fclus[i]*model['nclus[' + str(i) + ']'])/(model['nbkg[' + str(i) + ']']/C[i] + model['nclus[' + str(i) + ']'])

    pm.Poisson("obsnbkg[" + str(i) + "]", mu = model['nbkg[' + str(i) + ']'], observed = obsnbkg[i])
    pm.Poisson("obsntot[" + str(i) + "]", mu = (model['nbkg[' + str(i) + ']']/C[i] + model['nclus[' + str(i) + ']']), observed = obsntot[i])

    pm.Binomial("obsnbluebkg[" + str(i) + "]", p = model['fbkg[' + str(i) + ']'], n = model['obsnbkg[' + str(i) + ']'], observed = obsnbluebkg[i])
    pm.Binomial("obsnbluetot[" + str(i) + "]", p = f[i], n = model['obsntot[' + str(i) + ']'], observed = obsnbluetot[i])
```

## Chapter 8.6
```Python
with pm.Model() as model:
    pm.Uniform("intrscat", lower = 0, upper = 10)
    pm.Normal("alpha", mu = 0.0, sigma = 100)
    pm.StudentT("beta", mu = 0, nu = 1, sigma = 1)
    zptB = 24

    sB = [None]*len(obstotB)
    sS = [None]*len(obstotB)

    for i in range (0, len(obstotB)):
        pm.Uniform("bkgB[" + str(i) + "]", lower = 0, upper = 1e7)
        pm.Uniform("bkgS[" + str(i) + "]", lower = 0, upper = 1e7)
        pm.Uniform("magB[" + str(i) + "]", lower = 18, upper = 25)
        pm.Normal("lgfluxS[" + str(i) + "]", mu = ((model['magB[' + str(i) + ']'] - 22)*model['beta'] + model['alpha']),  sigma = model['intrscat'])
        sB[i] = pow(10, (zptB - model['magB[' + str(i) + ']'])/2.5)
        sS[i] = pow(10, model['lgfluxS[' + str(i) + ']'] - zptS[i])

        pm.Poisson("obstotB[" + str(i) + "]", mu = (sB[i] + model['bkgB[' + str(i) + ']']/CB[i]), observed = obstotB[i])
        pm.Poisson("obsbkgB[" + str(i) + "]", mu = (model['bkgB[' + str(i) + ']']), observed = obsbkgB[i])
        pm.Poisson("obstotS[" + str(i) + "]", mu = (sS[i] + model['bkgS[' + str(i) + ']']/CS[i]), observed = obstotS[i])
        pm.Poisson("obsbkgS[" + str(i) + "]", mu = (model['bkgS[' + str(i) + ']']), observed = obsbkgS[i])
```

## Chapter 8.7
```Python
n = [None]*len(obsn)
lgM = [None]*len(obsn)
phi = [None]*len(obsn)

zeros = obsn - obsn
C = 10

with pm.Model() as model:
  pm.StudentT("alpha", mu = 0, nu = 1, sigma = 1)
  pm.Normal("beta", mu = 14.4, sigma = 3)

  for i in range (0, len(obsn)):
    #zero trick for the schechter function (alpha = -2)
    pm.Uniform("lgn[" + str(i) + "]", lower = -1, upper = 3)
    n[i] = pow(10, model['lgn[' + str(i) + ']'])
    #-log likelihood
    phi[i] = n[i]/10**2 - (-2 + 1 - 1)*model['lgn[' + str(i) + ']']
    pm.Poisson("zeros[" + str(i) + "]", mu = (phi[i] + C), observed = zeros[i])

    lgM[i] = model['alpha']*(model['lgn[' + str(i) + ']'] - 1.5) + model['beta']
    pm.Poisson("obsn[" + str(i) + "]", mu = n[i], observed = obsn[i])
    pm.Normal("obslgMm[" + str(i) + "]", mu = lgM[i], sigma = err[i], observed = obslgMm[i])
```

## Chapter 8.9
```Python
with pm.Model() as model:
    lgnm = [None]*len(obstot)
    #simplified prior
    pm.Uniform("intrscat", lower = 0, upper = 5)
    pm.Normal("alpha", mu = 0.0, sigma = 100)
    pm.StudentT("beta", mu = 0, nu = 1, sigma = 1)

    for i in range (0, len(obstot)):
        #modeling data selection
        pm.Normal("lgM[" + str(i) + "]", mu = 14.5, sigma = 0.33)
        #decomment for ignoring data selection
        #pm.Uniform("lgM[" + str(i) + "]", lower = 13, upper = 16)
        #modeling mass -n relation
        lgnm[i] = model['alpha'] + 1.5 + model['beta']*(model['lgM[' + str(i) + ']'] - 14.5)
        pm.Normal("lgn[" + str(i) + "]", mu = lgnm[i], sigma = model['intrscat'])
        pm.Uniform("nbkg[" + str(i) + "]", lower = 0, upper = 3000)

        pm.Poisson("obsbkg[" + str(i) + "]", mu = model['nbkg[' + str(i) + ']'], observed = obsbkg[i])
        pm.Poisson("obstot[" + str(i) + "]", mu = (model['nbkg[' + str(i) + ']']/C[i] + 10**model['lgn[' + str(i) + ']']), observed = obstot[i])
        pm.Normal("obslgM[" + str(i) + "]", mu = model['lgM[' + str(i) + ']'], sigma = errlgM[i], observed = obslgM[i])
```

## Chapter 8.10
```Python
nu = 6
obsvarlgM = pow(errlgM, 2)

with pm.Model() as model:
  z = [None]*len(obstot)
  #simplified prior
  pm.Uniform("intrscat", lower = 0, upper = 5)
  pm.Normal("alpha", mu = 0.0, sigma = 100)
  pm.StudentT("beta", mu = 0, nu = 1, sigma = 1)

  for i in range (0, len(obstot)):
    pm.Uniform("n[" + str(i) + "]", lower = 0, upper = 3000)
    pm.Uniform("nbkg[" + str(i) + "]", lower = 0, upper = 3000)
    z[i] = model['alpha'] + 14.5 + model['beta']*(np.log(model['n[' + str(i) + ']'])/2.30258 - 1.5)
    pm.Normal("lgM[" + str(i) + "]", mu = z[i], sigma = model['intrscat'])

    pm.Gamma("precy[" + str(i) + "]", alpha = 1.0e-5, beta = 1.0e-5)
    pm.Normal("obslgM[" + str(i) + "]", mu = model['lgM[' + str(i) + ']'], sigma = 1/pm.math.sqrt(model['precy[' + str(i) + ']']), observed = obslgM[i])
    pm.Gamma("obsvarlgM[" + str(i) + "]", alpha = (0.5*nu), beta = (0.5*nu*model['precy[' + str(i) + ']']), observed = obsvarlgM[i])

    pm.Poisson("obsbkg[" + str(i) + "]", mu = model['nbkg[' + str(i) + ']'], observed = obsbkg[i])
    pm.Poisson("obstot[" + str(i) + "]", mu = (model['nbkg[' + str(i) + ']']/C[i] + model['n[' + str(i) + ']']), observed = obstot[i])
```

## Chapter 8.11
```Python
with pm.Model() as model:
  z1 = [None]*len(x1)
  z2 = [None]*len(x2)
  z3 = [None]*len(x3)

  pm.Uniform("intrscat", lower = 0, upper = 3, shape = 3)
  pm.Normal("alpha", mu = 0.0, sigma = 100, shape = 3)
  pm.StudentT("beta", mu = 0, nu = 1, sigma = 1)

  for i in range(0, len(x1)):
    #modeling y - x
    z1[i] = model['alpha'][0] + 0.1 + model['beta']*(x1[i] - 0.03)
    #modeling y
    pm.Normal("y1[" + str(i) + "]", mu = z1[i], sigma = model['intrscat'][0])
    pm.Normal("obsy1[" + str(i) + "]", mu = model['y1[' + str(i) + ']'], sigma = err.obsy1[i], observed = obsy1[i])

  for i in range(0, len(x2)):
    #modeling y - x
    z2[i] = model['alpha'][1] + 0.1 + model['beta']*(x2[i] - 0.03)
    #modeling y
    pm.Normal("y2[" + str(i) + "]", mu = z2[i], sigma = model['intrscat'][1])
    pm.Normal("obsy2[" + str(i) + "]", mu = model['y2[' + str(i) + ']'], sigma = err.obsy2[i], observed = obsy2[i])

  for i in range(0, len(x3)):
    #modeling y - x
    z3[i] = model['alpha'][2] + 0.1 + model['beta']*(x3[i] - 0.03)
    #modeling y
    pm.Normal("y3[" + str(i) + "]", mu = z3[i], sigma = model['intrscat'][2])
    pm.Normal("obsy3[" + str(i) + "]", mu = model['y3[' + str(i) + ']'], sigma = err.obsy3[i], observed = obsy3[i])
```

## Chapter 8.12.2.6
```Python
with pm.Model() as model:
    ft = [None]*len(modeT)
    Abcor = [None]*len(modeT)
    Abcor_rep = [None]*len(modeT)

    pm.Uniform("Abz02", lower = 0, upper = 1)
    pm.StudentT("alpha", mu = 0, nu = 1, sigma = 1)
    pm.Uniform("tau", lower = 1, upper = 100)
    pm.Uniform("intrscat", lower = 0, upper = 1)
    pm.Normal("factor", mu = (0.77 - 1), sigma = 0.065, transform = Interval(-1, None))

    for i in range (0, len(modeT)):
        pm.Uniform("T[" + str(i) + "]", lower = 1, upper = 20)
        ft[i] = model['Abz02']*(1 - np.exp(-t[i]/model['tau']))/(1 - np.exp(-11/model['tau']))
        pm.LogNormal("Ab[" + str(i) + "]", mu = np.log(ft[i]), sigma = model['intrscat'])

        Abcor[i] = model['Ab[' + str(i) + ']']*pow((model['T[' + str(i) + ']']/5), model['alpha'])*(1 + model['factor']*tid[i])
        pm.Normal("modeAb[" + str(i) + "]", mu = Abcor[i], sigma = sigmaAb[i], observed = modeAb[i])
        pm.LogNormal("modeT[" + str(i) + "]", mu = np.log(model['T[' + str(i) + ']']), sigma = sigmaT[i], observed = modeT[i])

        # for p-value computation
        pm.LogNormal("Ab_rep[" + str(i) + "]", mu = np.log(ft[i]), sigma = model['intrscat'])
        Abcor_rep[i] = model['Ab_rep[' + str(i) + ']']*pow((model['T[' + str(i) + ']']/5), model['alpha'])*(1 + model['factor']*tid[i])
        pm.Normal("modeAb_rep[" + str(i) + "]", mu = Abcor_rep[i], sigma = sigmaAb[i])
```

## Chapter 8.12.3
```Python
with pm.Model() as model:
    Abcor = [None]*len(modeT)

    pm.Uniform("meanAb", lower = 0.0, upper = 1)
    pm.Normal("alpha", mu = -0.12, sigma = 0.09)
    pm.Normal("intrscat", mu = 0.18, sigma = 0.03)
    pm.Normal("factor", mu = -0.22, sigma = 0.045, transform = Interval(-1, None))

    for i in range(0, len(modeT)):
        pm.Uniform("T[" + str(i) + "]", lower = 1, upper = 20)
        pm.LogNormal("Ab[" + str(i) + "]", mu = np.log(model['meanAb']), sigma = model['intrscat'])
        Abcor[i] = model['Ab[' + str(i) + ']']*pow((model['T[' + str(i) + ']']/5), model['alpha'])*(1 + model['factor']*tid[i])
        pm.TruncatedNormal("modeAb[" + str(i) + "]", mu = Abcor[i], sigma = sigmaAb[i], lower = 0, observed = modeAb[i])
        pm.LogNormal("modeT[" + str(i) + "]", mu = np.log(model['T[' + str(i) + ']']), sigma = sigmaT[i], observed = modeT[i])
```

## Chapter 8.12.4
```Python
with pm.Model() as model:
  z = [None]*len(obstot)
  nclus = [None]*len(obstot)

  pm.Uniform("intrscat", lower = 0, upper = 10)
  pm.Normal("alpha", mu = 0.0, sigma = 100)
  pm.StudentT("beta", mu = 0, nu = 1, sigma = 1)

  for i in range(0, len(obstot)):
    # n200 terms
    pm.Uniform("n[" + str(i) + "]", lower = 1, upper = 3000)
    # modelling Lx -n200 relation
    z[i] = model['alpha'] + 44 + model['beta']*(np.log(model['n[' + str(i) + ']'])/2.30258 - 1.8)
    pm.Normal("lgLx[" + str(i) + "]", mu = z[i], sigma = model['intrscat'])
    # convert nclus in Lx
    nclus[i] = np.exp(2.30258*(model['lgLx[' + str(i) + ']'] - C[i]))
    # modelling x-ray protons
    pm.Uniform("nbkg[" + str(i) + "]", lower = 0, upper = 10000)
    pm.Poisson("obsbkg[" + str(i) + "]", mu = model['nbkg[' + str(i) + ']'], observed = obsbkg[i])
    pm.LogNormal("nbkgind[" + str(i) + "]", mu = np.log(model['nbkg[' + str(i) + ']']), sigma = 1/np.sqrt((1/0.2/0.2)))
    pm.Poisson("obstot[" + str(i) + "]", mu = (nclus[i] + model['nbkgind[' + str(i) + ']']/nbox[i]), observed = obstot[i])
    # modelling galaxy counts

    pm.Uniform("ngalbkg[" + str(i) + "]", lower = 0, upper = 3000)
    pm.Poisson("obsgalbkg[" + str(i) + "]", mu = model['ngalbkg[' + str(i) + ']'], observed = obsgalbkg[i])
    pm.Poisson("obsgaltot[" + str(i) + "]", mu = (model['ngalbkg[' + str(i) + ']']/Cgal[i] + model['n[' + str(i) + ']']), observed = obsgaltot[i])
```

## Chapter 9.1.1
```Python
sig = 3.3
with pm.Model() as model:
    pm.LogNormal("m", mu = 1, sigma = 1)
    pm.Normal("obsm", mu = model['m'], sigma = sig, observed = obsm)
```

## Chapter 9.1.1
```Python
with pm.Model() as model:
    logYm = [None]*len(obslogM)

    pm.Uniform("intrscat", lower = 0, upper = 2)
    pm.Normal("alpha", mu = 0.0, sigma = 100)
    pm.StudentT("beta", mu = 0, nu = 1, sigma = 1)
    #gamma = 2/3 # to hold gamma fixed
    gamma = pm.StudentT("gamma", mu = 0, nu = 1, sigma = 1)

    for i in range (0, len(obslogM)):
        pm.Uniform("logM[" + str(i) + "]", lower = 13, upper = 17)
        logYm[i] = model['alpha'] - 4.3 + model['beta']*(model['logM[' + str(i) + ']'] - 14.778) + gamma*np.log(Ez[i])/2.303
        pm.Normal("logYSZ[" + str(i) + "]", mu = logYm[i], sigma = model['intrscat'])
        pm.Normal("obslogM[" + str(i) + "]", mu = model['logM[' + str(i) + ']'], sigma = errlogM[i], observed = obslogM[i])
        pm.Normal("obslogYSZ[" + str(i) + "]", mu = model['logYSZ[' + str(i) + ']'], sigma = errlogYSZ[i], observed = obslogYSZ[i])
```
