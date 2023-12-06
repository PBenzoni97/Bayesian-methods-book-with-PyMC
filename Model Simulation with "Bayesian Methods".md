# SIMULATE A MODEL WITH THE LIBRARY "BAYESIAN METHODS"
Here we can see an example of a model simulation with 'bayesianmethods'
### STEP 1: import libraries:
First of all, we need to import the libraries: 
```Python
import pymc as pm
import corner
import numpy as np
import arviz as az
from pymc.distributions.transforms import Interval
```
The last line is specifically needed for the codes from chapter 8.12.2.6 and 8.12.3
### STEP 2: setup the code for running:
We start by choosing a model: for this example, we simulate the model from chapter 6.1.2
```Python
(obstot, obsbkg, C) = np.loadtxt('/Dataset/data6.1.2.dat.R', delimiter = ",", unpack = True)

with pm.Model() as model:
  pm.Uniform("s", lower = 0, upper = 1.0e7)
  pm.Uniform("bkg", lower = 0, upper = 1.0e7)
  pm.Poisson("obstot", mu = (model['s'] + model['bkg']/C), observed = obstot)
  pm.Poisson("obsbkg", mu = model['bkg'], observed = obsbkg)

  chain = pm.sample(return_inferencedata = True)

  print(az.summary(chain.posterior['A']))
  print(az.summary(chain.posterior['B']))

  figure = corner.corner(chain.posterior['A'][0])
  figure = corner.corner(chain.posterior['B'][0])
```
First of all, we need to import the data from the dataset:
```Python
(obstot, obsbkg, C) = np.loadtxt('/Dataset/data6.1.2.dat.R', delimiter = ",", unpack = True)
```
Then we specify a model, wich will contain the priors, the likelihoods and the instructions for sampling and plotting the result:
In this case the priors are:
```Python
pm.Uniform("s", lower = 0, upper = 1.0e7)
pm.Uniform("bkg", lower = 0, upper = 1.0e7)
```
And the likelihoods are:
```Python
pm.Poisson("obstot", mu = (model['s'] + model['bkg']/C), observed = obstot)
pm.Poisson("obsbkg", mu = model['bkg'], observed = obsbkg)
```
To sample the posterior distributions we use this instruction:
```Python
chain = pm.sample(return_inferencedata = True)
```
### Step 3: run the simulation:
Now we can start the simulation.

During the execution, a green progressive bar will show up. The bar shows the progress of the sampling, the number of samples drawn, the extimated time completion time, the number of the chain and the divergences:

![Screenshot 2023-11-04 120046](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/e88c5b1b-d5a6-4cdb-bc30-1f92dfbb14b0)

The result of the execution will be a table containing the details of the obtained posterior:

![Screenshot 2023-11-09 132621](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/3e973628-0646-4fa2-b087-e602f8e4d983)

We obtain this table with these instructions:
```Python
print(az.summary(chain.posterior['s']))
print(az.summary(chain.posterior['bkg']))
```

And a plot of the curve:

![Screenshot 2023-11-09 132627](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/e38c3618-d413-405c-b83d-6d6eaa04038b)

We obtain this table with these instructions:
```Python
figure = corner.corner(chain.posterior['s'][0])
figure = corner.corner(chain.posterior['bkg'][0])
```
