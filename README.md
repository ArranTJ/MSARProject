EISPAC Modification for Kappa Functions
=
This project repository uses the EIS Python Analysis Code (EISPAC) cited here:

Weberg, M. J., Warren, H. P., Crump, N., & Barnes, W. (2023). EISPAC - The EIS Python Analysis Code (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.1234

Please refer to the documentation provided here when installing and using EISPAC software at https://eispac.readthedocs.io/en/latest/

The EISPAC GitHub repo can be found here: https://github.com/USNavalResearchLaboratory/eispac
_______________________________________________________________________________________________________
NOTE: that this is very much a work in progress and requires extensive cleanup and modification. 

This directory contains the modified contents of EISPAC version 0.95.0 and its features.

The key modifications are in the fitting_functions.py, and extern/__init__.py files, as well as a new fit_spectraKAPPA.py file to apply the fitting of the function to the parameters

This customised EISPAC directory (provided in EISPAC_KAPPA (eispac version 0.95.0).7z) permits the use of Kappa functions to fit to the spectra.

Produce custom template files using eispac.EISFitTemplate() to include a free parameter for kappa in addition to parameters already provided in the existing templates.

Fit the observation spectra using this template by utilising eispac.fit_spectraKAPPA().

 
