EISPAC Modification for Kappa Functions
=
This project repository uses the EIS Python Analysis Code (EISPAC) cited here:

Weberg, M. J., Warren, H. P., Crump, N., & Barnes, W. (2023). EISPAC - The EIS Python Analysis Code (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.1234

Please refer to the documentation provided here when installing and using EISPAC software at https://eispac.readthedocs.io/en/latest/

The EISPAC GitHub repo can be found here: https://github.com/USNavalResearchLaboratory/eispac
_______________________________________________________________________________________________________
NOTE: This is very much a work in progress and requires extensive cleanup and modification. 

This directory contains the modified contents of EISPAC version 0.95.0 and its features.

The key modifications are in the core/fitting_functions.py, and core/__init__.py files, as well as a new core/fit_spectraKAPPA.py file to apply the fitting of the function to the parameters

This customised EISPAC directory (provided in EISPAC_KAPPA (eispac version 0.95.0).7z) permits the use of Kappa functions to fit to the spectra.

Produce custom templates using eispac.EISFitTemplate() to include a free parameter for the kappa value in addition to parameters already provided in the existing templates.

##NOTE: 
This assumes a template containing 5 parameters (p[0],p[1],p[2],p[3],p[4]). 

Changes to the number of parameters in the template will require a change to the parameter assigned to kappa in fitting_functions.py
(i.e if the number of parameters is 5, and the kappa value is the last parameter, the kappa value is passed into fitting_functions.py as p[4]. 

If the number of parameters is larger than this, line 194: p = param\[3*n:3*n+5] must be changed to permit the additional parameter indices.)##


Once the template parameters are set up correctly, the spectra in the selected observation can be fitted with the kappa function by calling eispac.fit_spectraKAPPA().

"FittingMultiFile.py" is the primary file used to call the EISPAC fitting functions and modules. This also produces a series of plots, and saves numpy data arrays into their assigned directories.
By loading EIS data + header files into the 'filequeue' folder, FittingMultiFile.py will automatically run through and process.

Statistical information about the number of kappa fits, mean kappa values, and mean Chi^2 values for each observation is saved into .csv files under the folder 'goodkappavals'.
'GOEScurve.py' downloads and produces a GOES X-ray curve, and uses Pandas to organise and arrange the datasets of 'goodkappavals' into different arrays for each phase of the flare.
This file additionally shows subplots for counts of kappa values alongside the GOES X-ray curve from the data provided in 'goodkappavals'.
The folders Kappa1 and Kappa2 contain examples of statistical information produced during the project.
 
