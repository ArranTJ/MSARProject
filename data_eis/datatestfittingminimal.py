# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 15:43:06 2024

@author: Arran
"""

import eispac

import numpy as np
import matplotlib.pyplot as plt

from astropy.visualization import astropy_mpl_style
from astropy.io import fits

plt.style.use(astropy_mpl_style)

if __name__ == '__main__':
    # input data and template files
    data_filepath = './eis_20170906_120233.data.h5'
    # template_filepath = './fe_14_264_787.1c.template.h5'
    template_filepath = './fe_23_263_760.1c.template.h5'
    template_filepath = './fe_16_262_984.1c.template.h5'
    # read fit template
    tmplt = eispac.read_template(template_filepath)
    print(tmplt)
    
    ###Template parser to parse template data values and pass into fitting functions###
    def templateparser(tmplt):
        pvallst = []
        plimslst = []
        for i in range(len(tmplt.parinfo)):
            pvallst.append(tmplt.parinfo[i]['value'])
            plimslst.append(tmplt.parinfo[i]['limits'])
            
        return pvallst, plimslst
    # tmpltfunccall = templateparser(tmplt)
    # print(tmpltfunccall[1])
    
    # print(templateparser(tmplt)[1][1][0])
    
    pvals = templateparser(tmplt)[0]
    plims = templateparser(tmplt)[1]
    
    new_templt = eispac.EISFitTemplate(value=[pvals[0],pvals[1], pvals[2], pvals[3], 10], limits = [[0.0000, 0.0000],[plims[1][0], plims[1][1]],[plims[2][0], plims[2][1]],[plims[3][0], plims[3][1]], [0, 100]], line_ids=['Fe XVI 262.984'], wmin= 262.78399658203125, wmax=263.1000061035156)
    print(new_templt)
    # # print(pvals[0])
    # # print(plims)
    
    
    # # pvallist = []
    # # pvallst.append(tmplt.parinfo[i])
    
    # ###Defining the Kappa distribution function for fitting### 
    # def kappafunc(x, amp, cen, wid, kappa):
    #     #defining the variance-like variable indicated in Jeffreys. N. S et al (2016)
    #     variance_k = wid * (1 - 3/(2*kappa))
    #     #argument used in kappa function 
    #     arg = (1 + ((x - cen)**2)/(2 * (variance_k**2) * kappa))
        
    #     return amp * (arg) ** (-1 *kappa)
    
    # def kappafitter(x, tmplt):
    #     pvals = templateparser(tmplt)[0]
    #     plims = templateparser(tmplt)[1]
        
        
        
    #     return 
    
    # print(kappafitter(1, tmplt))
        
        
        


    # def kappafitterscipy(x, tmplt):
        
    #     np.array
        
        
    #     return
    
    

    # # Read spectral window into an EISCube
    # data_cube = eispac.read_cube(data_filepath, tmplt.central_wave)

    # #Fit the data, then save it to disk and test loading it back in
    # fit_res = eispac.core.fit_spectraKAPPA(data_cube, tmplt, ncpu='max')
    # save_filepaths = eispac.save_fit(fit_res, save_dir='cwd')
    # FITS_file = eispac.export_fits(fit_res, save_dir='cwd')
    # # load_fit = eispac.read_fit(save_filepaths[0])
    # load_fit = eispac.read_fit(save_filepaths)
    
    # print("Chi2 param?", fit_res.get_params(param_name='chi2'))   
    # # for key in fit_res.fit.keys():

    # #     print(f"{key:<15} {fit_res.fit[key]} {fit_res.fit[key]}")
    
    # chi2para = fit_res.fit['chi2']
    
    # image_file = 'eis_20170906_120233.fe_14_274_203.1c-0.wid.fits'
    
    # fits.info(image_file)
    
    # image_data = fits.getdata(image_file, ext=0)
    
    # plt.figure()
    # plt.imshow(image_data)
    # plt.colorbar()
    
    # plt.show()