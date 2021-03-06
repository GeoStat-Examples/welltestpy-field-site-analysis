[![GS-Frame](https://img.shields.io/badge/github-GeoStat_Framework-468a88?logo=github&style=flat)](https://github.com/GeoStat-Framework)
[![Gitter](https://badges.gitter.im/GeoStat-Examples/community.svg)](https://gitter.im/GeoStat-Examples/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4139373.svg)](https://doi.org/10.5281/zenodo.4139373)

# Pumping test analysis with welltestpy

## Description

In this workflow, we analyse two pumping test campaigns on two field sites:

- ***"Horkheimer Insel"*** (Heilbronn, Germany)
- ***"Lauswiesen"*** (Tübingen, Germany)

The aim is to estimate parameters of heterogeneity from transient puming test data
and to analyse their sensitivities.

### Target parameters

- mean of log-transmissivity
- variance of log-transmissivity
- length scale of log-transmissivity
- storage

### Applied methods

The applied methods utilizing effecitive head solutions for the groundwater flow
equation under a pumping test condition are described in:

> Zech, A., Müller, S., Mai, J., Heße, F., and Attinger, S.:
> Extending Theis’ Solution: Using Transient Pumping Tests to Estimate Parameters of Aquifer Heterogeneity,
> Water Resour. Res., 52, 6156–6170, https://doi.org/10.1002/2015WR018509, 2016.

These methods were implemented in [`welltestpy`](https://github.com/GeoStat-Framework/welltestpy)
to automatically interprete pumping test data.

The underlying type-curves are implemented in [`AnaFlow`](https://github.com/GeoStat-Framework/AnaFlow).

## Data sources

The data for the "Horkheimer Insel" field site was manually taken from:

> Schad H.: Variability of hydraulic parameters in non-uniform porous media:
> experiments and stochastic modeling at different scales.
> University Tübingen; 1997. Ph.D. thesis.

The pumping data from the "Lauswiesen" field site was kindly provided by
[Dr. Carsten Leven-Pfister](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/geowissenschaften/arbeitsgruppen/angewandte-geowissenschaften/angewandte-geowissenschaften-zag/hydrogeologie/hydrogeology/carsten-leven-pfister/)
and is made available on a repository of the University of Tübingen:
[Research Data Portal FDAT](http://hdl.handle.net/10900.1/bfb0b0f7-7065-4a24-8a91-92ad8aa8fc40)

Citable as:

> Leven, C. (2020):  Pumping Tests in Wells B1-B5 at the Hydrogeological Research Site Lauswiesen,
> Eberhard Karls Universität Tübingen,
> http://hdl.handle.net/10900.1/bfb0b0f7-7065-4a24-8a91-92ad8aa8fc40

## Structure

The workflow is organized by the following structure:

- `data/`
  - contains the campaign files for both sites in the `welltestpy` format
  - contains time series for diagnostic plots
- `src/` - contains the scripts to produce the results
  - `00_wtp_plot.py` - plotting well-constellation and campaign overviews
  - `01_est_run.sh` - bash file running `02_para_estimation.py` in parallel
  - `01b_est_run.sh` - bash file running `02b_para_estimation.py` in parallel
  - `02_para_estimation.py` - estimate parameters of heterogeneity from the pumping tests
  - `02b_para_estimation.py` - estimate equivalent parameters of homogeneity the pumping tests
  - `03_postpro_results.py` - plotting the estimation results for both sites
  - `04_postpro_sensitivity.py` - plotting the sensitivity results for both sites
  - `05_est_radial_sens.sh` - bash file running `06_rad_sens.py` in parallel
  - `06_rad_sens.py` - estimate parameter sensitivites depending on the radial distance
    to the pumping well. when run in serial, results will be plotted.
  - `07_comparison_len_scale.py` - generate comparison plot for different length scales
  - `08_check_unconfined_effect.py` - generate diagnostic plots
- `results/` - all produced results


## Python environment

Main Python dependencies are stored in `requirements.txt`:

```
welltestpy==1.0.3
anaflow==1.0.1
spotpy==1.5.9
mpi4py==3.0.2
matplotlib
```

You can install them with `pip` (potentially in a virtual environment):

```bash
pip install -r requirements.txt
```


## Contact

You can contact us via <info@geostat-framework.org>.


## License

MIT © 2021
