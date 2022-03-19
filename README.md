# Gravitational Wave Analysis in Python

The main goal of the GravitationalWaves package is to simulate, detect, and analyze continuous gravitational waves in Python. In addition to creating simulations of gravitational waves, the package also takes in observed data for comparison in detection and analysis.


--------------------

## Installation Instructions

### Pip Install

Pip supports installing packages from a Github repository using the URI form `git+https://github.com/user/project.git@{version}`. To pip install the GravitationalWaves package from Github, run the following in the command line:


```
pip install -e git+https://github.com/hluebbering/GravitationalWaves.git
```

Running this command clones and installs the latest version of this package from Github.



### Manual Package Install

Most Python packages are now designed to be compatible with pip. If you have a package that’s not compatible, then you’ll need to do a manual installation.

How to manually install a Python package:

1. Download or clone this repository into your local directory.
2. Unzip the repository if it's downloaded as a zip file.
3. Open your command windows and change the working directory to the repository containing setup.py using the `cd` command.
4. Use Python to install the package using the following command:

```
python setup.py install
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


Docstrings on each module and class
Docstrings on each public method, contains parameters and return

--------------------

## Sources
- https://joss.theoj.org/papers/10.21105/joss.02968
- https://joss.theoj.org/papers/10.21105/joss.03000


--------------------

## To Do

- [x] README.md file that gives an overview of the project
- [x] LICENSE file
- [x] setup.py file that initializes the project after it has been cloned
- [ ] doc folder that contains documentation (including the functional specification, the design specification, and the final project presentation or poster)
- [ ] python package folder (with the same name as the repository) that is structured as one or more python modules (e.g., with init.py files) and test files that begin with “test_”.
- [x] examples folder that contains examples of using the packages


### Project presentation

- [ ] Background. Describe the problem or area being addressed.
- [ ] Data used. What data did you use? How was it obtained? What are its limitations?
- [ ] Use cases. How users will interact with your system in a way that addresses the problem area.
- [ ] Demo. Demonstrate your software.
- [ ] Design. Describe the components and how they interact to accomplish the use cases.
- [ ] Project Structure. Show the structure of your github repository.
- [ ] Lessons learned and future work. Focus on software engineering lessons.




gravitational-waves

* [docs/](.\gravitational-waves\docs)
  * [index.html](.\gravitational-waves\docs\index.html)
  * [present.Rmd](.\gravitational-waves\docs\present.Rmd)
* [examples/](.\gravitational-waves\examples)
  * [01_simulate_gravitational_waves.py](.\gravitational-waves\examples\01_simulate_gravitational_waves.py)
  * [02_calculate_SNR.py](.\gravitational-waves\examples\02_calculate_SNR.py)
  * [03_horizon_distance.py](.\gravitational-waves\examples\03_horizon_distance.py)
* [gravitational-waves/](.\gravitational-waves\gravitational-waves)
  * [tests/](.\gravitational-waves\gravitational-waves\tests)
    * [test_inspiral.py](.\gravitational-waves\gravitational-waves\tests\test_inspiral.py)
    * [test_source.py](.\gravitational-waves\gravitational-waves\tests\test_source.py)
    * [test_utils.py](.\gravitational-waves\gravitational-waves\tests\test_utils.py)
    * [test_visualization.py](.\gravitational-waves\gravitational-waves\tests\test_visualization.py)
    * [__init__.py](.\gravitational-waves\gravitational-waves\tests\__init__.py)
  * [inspiral.py](.\gravitational-waves\gravitational-waves\inspiral.py)
  * [source.py](.\gravitational-waves\gravitational-waves\source.py)
  * [utils.py](.\gravitational-waves\gravitational-waves\utils.py)
  * [visualization.py](.\gravitational-waves\gravitational-waves\visualization.py)
  * [__init__.py](.\gravitational-waves\gravitational-waves\__init__.py)
* [LICENSE](.\gravitational-waves\LICENSE)
* [README.md](.\gravitational-waves\README.md)
* [setup.py](.\gravitational-waves\setup.py)



