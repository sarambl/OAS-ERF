# OAS-ERF
Analysis for paper:
- Blichner, S. M., Sporre, M. K., and Berntsen, T. K.: Reduced effective radiative forcing from cloud-aerosol interactions (ERFaci) with improved treatment of early aerosol growth in an Earth System Model, Atmos. Chem. Phys. Discuss. [preprint], https://doi.org/10.5194/acp-2021-151, in review, 2021.



# Setup:
### Download:
```bash
git clone git@github.com:sarambl/OAS-ERF.git 
cd OAS-ERF/
```


### Install environment: 
```bash
# install environment etc
conda env create -f environment.yml
conda activate env_oas_erf
conda develop .

cd ../
git clone https://git.nilu.no/ebas/ebas-io.git  cd
cd ebas-io/dist/
pip install  ebas_io-3.6.1-py3-none-any.wh


```

## Download data:


## Edit settings and paths: 
Edit paths at the top of [oas_erf/constants.py](oas_erf/constants.py). 


##


## To reproduce results:
### Create Nd datasets:
Run for correct time period (edit "startyear" and "startyear" in script) and relevant cases
```bash
cd oas_erf/preprocess/
python Nd.py
python preproc_maps.py
```

### Run notebooks:
Run the notebooks in [oas_erf/notebooks](oas_erf/notebooks).