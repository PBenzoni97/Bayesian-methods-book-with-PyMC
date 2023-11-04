# Bayesian-Methods
## Python translation of the codes from "Bayesian Methods for the Physical Sciences"
'bayesianmethods' contains the codes from the book "Bayesian Methods for the Physical Sciences" by Andreon and Weaver, translated from JAGS to Python using the library PyMC.

### Requirements
'bayesianmethods' requires the following libraries:
-  [pymc](https://www.pymc.io/)
-  [corner](https://corner.readthedocs.io/en/latest/)
-  [numpy](https://numpy.org/)
-  [arviz](https://python.arviz.org/)

Now we'll see a quick guide on how to use the codes from 'bayesianmethods':
### STEP 1: install libraries:
First of all, we need to install PyMC and Corner.py. In order to install them we need to run this code lines on our notebook:

![Screenshot 2023-11-04 114656](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/bf624900-6ab3-46c2-a65a-9d2d763a6cbf)

### STEP 2: import libraries:
After installing the libraries, we need to import them: 

![Screenshot 2023-11-04 114714](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/19382fc9-0665-48a3-aee1-6fb0f56fe1f7)

Note: the installation of the library 'numpy' is not required, you only need to import it.

The last line is needed specifically for the codes from chapter 8.12.2.6 and 8.12.3
### STEP 3: setup the code for running:
The codes contained in 'bayesianmethods' is not ready to be run, but it needs a few add-ons before it can be executed.

As an example, we'll look at the first code from 'bayesianmethods', the one from chapter 4.1.2.

This is how the code snippet in 'bayesianmethods' is like:

![Screenshot 2023-10-16 112015](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/60b61b60-e1d4-4192-bd57-8351139bfb68)

In order to have it ready to run we need two things: to read the observed data from the dataset and to sample and plot the distribution.

To read the observed data we need to add this line of code at the beginning of the model:

![Screenshot 2023-11-04 114100](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/79195404-4866-49ec-ab83-ff9065faa32c)

Be sure to match the name of the variable, to put the correct address of where the dataset is located on your environment and to insert the correct delimiter.

To sample and plot the posterior distributions we need to add this lines at the end of the model:

![Screenshot 2023-11-04 114756](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/0ffbe89b-26e1-4f6e-9e0e-be105aa08103)

Be sure to put the correct name of the distribution used in the model where needed.

You can also customize the sampling steps by adding the parameter "draws" and specifying the number of steps. 
(for more details look at https://www.pymc.io/projects/docs/en/stable/api/generated/pymc.sample.html)

If the model looks like this, then it's ready to run:

![Screenshot 2023-11-04 115146](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/eb69d7ab-39af-4f7c-bfc8-33c4750624e9)

### Step 4: run the simulation:
Now we can start the simulation.

During the execution, a green progressive bar will show up:

The bar shows the progress of the sampling, the number of samples drawn, the extimated time completion time, the number of the chain and the divergences

![Screenshot 2023-11-04 120046](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/e88c5b1b-d5a6-4cdb-bc30-1f92dfbb14b0)

The result of the execution will be a table containing the details of the obtained posterior:

![Screenshot 2023-11-04 120155](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/796333a9-a1a2-4025-a45b-dec3bac5c4a4)

And a plot of the curve:

![Screenshot 2023-11-04 120233](https://github.com/ilBenza97/Bayesian-Methods/assets/145661415/b92fb55b-702e-4ef4-a1fd-c54f560e3393)

