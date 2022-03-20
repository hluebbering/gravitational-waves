# Gravitational Wave Analysis in Python

> `GravitationalWaves`: a Python package designed to simulate, detect, and analyze continuous gravitational wave-forms. In addition to creating simulations of gravitational waves, the package also takes in observed data for comparison in detection and analysis.

GravitationalWaves is a Python package for simulating the gravitational waveforms of binary mergers of black holes and neutron stars, computing several properties of these mergers and waveforms, and evaluating their detectability. In addition, GravitationalWaves also takes in observed data from gravitational wave detectors to compare data and predict detection rates. 




**Definition.** Gravitational waves are invisible distortions in spacetime, caused by massive events such as collisions between two black holes or neutron stars.

The importance of studying gravitational waveforms stems from the idea of detecting and using gravity to estimate other dynamical astrophysical phenomena, giving enormous potential for discovering parts of the universe that are invisible to the eye, such as black holes and other unknowns. 






[Link to the presentation document](https://github.com/hluebbering/GravitationalWaves/docs/presentation.html).

--------------------

## Installation Instructions

### Pip Install

Pip supports installing packages from a Github repository using the URI form `git+https://github.com/user/project.git@{version}`. To pip install the GravitationalWaves package from Github, run the following in the command line:


```bash
pip install -e git+https://github.com/hluebbering/GravitationalWaves.git#egg=GravitationalWaves
```

Running this command clones and installs the latest version of this package from Github.



### Manual Package Install

Most Python packages are now designed to be compatible with pip. If you have a package that’s not compatible, then you’ll need to do a manual installation.

How to manually install a Python package:

1. Download or clone this repository into your local directory.
2. Unzip the repository if it's downloaded as a zip file.
3. Open your command windows and change the working directory to the repository containing setup.py using the `cd` command.
4. Use Python to install the package using the following command:

```bash
python setup.py install
```

### Updating the Package

To update this package, run the following command:

```bash
pip install GravitationalWaves --upgrade
```


--------------------

## Dependencies

List of dependencies:
- python >= 3.7
- pip >= 21.0.0
- units
- importlib
- numba >= 0.50
- numpy >= 1.16
- astropy >= 4.0
- scipy >= 1.5.0
- matplotlib >= 3.3.2
- seaborn >= 0.11.1
- schwimmbad >= 0.3.2
- legwork >= 0.2.4
- pip:
   - GravitationalWaves >= 0.1 


--------------------

## Usage Examples

1. [Instantiate a Source Class](https://github.com/hluebbering/GravitationalWaves/blob/main/examples/01_InstantiateSourceClass.ipynb) 
2. [Calculate Signal-to-Noise Ratio (SNR)](https://github.com/hluebbering/GravitationalWaves/blob/main/examples/02_CalculateSNR.ipynb)
3. [Plot Source Parameters Distribution](https://github.com/hluebbering/GravitationalWaves/blob/main/examples/03_PlotSourceDistribution.ipynb)
4. [Plot Sensitivity Curve](https://github.com/hluebbering/GravitationalWaves/blob/main/examples/04_Visualizations.ipynb)
5. [Simulate Gravitational Waves](https://github.com/hluebbering/GravitationalWaves/blob/main/examples/05_SimulateGravitationalWaves.ipynb)


### Demo 1. Single source SNR calculation


The most basic use case for GravitationalWaves is to calculate the signal-to-noise ratio of a single stellar-mass binary system. Using the package's source module, we first generate a source class and then calculate its SNR.


```python
import GravitationalWaves as gw
import astropy.units as u
source = gw.source.Source(m_1 = 11 * u.Msun,
                          m_2 = 11 * u.Msun,
                          ecc = 0.3,
                          f_orb = 1e-4 * u.Hz,
                          dist = 9 * u.kpc,
                          interpolate_g = False)
                          
source.get_snr()
```

For this example, GravitationalWaves checks whether the source is eccentric/circular and evolving/stationary, and chooses the fastest method to accurately calculate the SNR. 


### Demo 2. Multiple source SNR calculate

In the next example, we use GravitationalWaves to calculate the detectability of a collection of sources. We first import the package, then create a random set of possible LISA sources.


```python
import GravitationalWaves.source as source
import numpy as np
import astropy.units as u

# create a random collection of sources
nVals = 1800
mass1 = np.random.uniform(0, 12, n_values) * u.Msun
mass2 = np.random.uniform(0, 12, n_values) * u.Msun
eccInit = 1 - np.random.power(6, n_values)
lumDist = np.random.normal(9, 1.5, n_values) * u.kpc
orbFreq = 12**(-6 * np.random.power(3, n_values)) * u.Hz
```


Using these random sources, we can instantiate a Source class.

```python
sources = source.Source(m_1 = mass1, m_2 = mass2, 
                        ecc = eccInit, dist = lumDist, 
                        f_orb = orbFreq)
```

Now, we can calculate the SNR (signal-to-noise ratio) for these sources. This function splits the sources based on whether they are stationary/evolving and circular/eccentric.


```python
snr = sources.get_snr(verbose = True)
```

The SNR values are now stored in `sources.snr`, and we can mask any values that don’t meet some detectable threshold. In the following, we set the threshold to 7.


```python
detectableSources = sources.snr > 7
print("{} out of the {} sources are detectable".format(
  len(sources.snr[detectableSources]), nVals))
```


--------------------

## References

- https://joss.theoj.org/papers/10.21105/joss.02968
- https://joss.theoj.org/papers/10.21105/joss.03000


--------------------

## Project Directory

- [x] README.md file that gives an overview of the project
- [x] LICENSE file
- [x] setup.py file that initializes the project after it has been cloned
- [ ] doc folder that contains documentation (including the functional specification, design specification, and final project presentation)
- [x] python package folder (with the same name as the repository) that is structured as one or more python modules (e.g., with init.py files) and test files that begin with "test_".
- [x] examples folder that contains examples of using the packages




```markdown
gravitational-waves
 ┣ docs
 ┃ ┗ presentation.Rmd
 ┣ examples
 ┃ ┣ 01_InstantiateSourceClass.ipynb
 ┃ ┣ 02_CalculateSNR.ipynb
 ┃ ┣ 03_PlotSourceDistribution.ipynb
 ┃ ┣ 04_Visualizations.ipynb
 ┃ ┣ 05_SimulateGravitationalWaves.ipynb
 ┃ ┗ README.md
 ┣ GravitationalWaves
 ┃ ┣ tests
 ┃ ┃ ┣ test_psd.py
 ┃ ┃ ┣ test_source.py
 ┃ ┃ ┣ test_strain.py
 ┃ ┃ ┣ test_utils.py
 ┃ ┃ ┣ test_visualization.py
 ┃ ┃ ┣ test_wavesim.py
 ┃ ┃ ┗ __init__.py
 ┃ ┣ psd.py
 ┃ ┣ R.npy
 ┃ ┣ source.py
 ┃ ┣ strain.py
 ┃ ┣ utils.py
 ┃ ┣ visualization.py
 ┃ ┣ wavesim.py
 ┃ ┗ __init__.py
 ┣ LICENSE
 ┣ README.md
 ┗ setup.py
```


