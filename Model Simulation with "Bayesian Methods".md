# FITTING A MODEL WITH THE LIBRARY "Bayesian-Methods"

### STEP 1: import libraries:

```Python
import pymc as pm
import corner
import numpy as np
import arviz as az
from pymc.distributions.transforms import Interval
```
The last line is only needed for some of the codes

### STEP 2: setup the code for running:

For this example, we fit the model from chapter 6.1.2.
First of all, we need to import the data from the dataset:
```Python
(obstot, obsbkg, C) = np.loadtxt('/Dataset/data6.1.2.dat.R', delimiter = ",", unpack = True)
```
The full model is specified in the library.
```Python
with pm.Model() as model:
  pm.Uniform("s", lower = 0, upper = 1.0e7)
  pm.Uniform("bkg", lower = 0, upper = 1.0e7)
  pm.Poisson("obstot", mu = (model['s'] + model['bkg']/C), observed = obstot)
  pm.Poisson("obsbkg", mu = model['bkg'], observed = obsbkg)
```

### Step 3: run the simulation:

To sample the posterior distributions we use this instruction:
```Python
chain = pm.sample(return_inferencedata = True)
```

![Screenshot 2023-11-04 120046](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/e88c5b1b-d5a6-4cdb-bc30-1f92dfbb14b0)

To get a summary of the posterior you can use:
```Python
print(az.summary(chain.posterior['s']))
print(az.summary(chain.posterior['bkg']))
```

![Screenshot 2023-11-09 132621](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/3e973628-0646-4fa2-b087-e602f8e4d983)

To draw a plot of the posterior you can use:
```Python
figure = corner.corner(chain.posterior['s'][0])
figure = corner.corner(chain.posterior['bkg'][0])
```

![Screenshot 2023-11-09 132627](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/e38c3618-d413-405c-b83d-6d6eaa04038b)
