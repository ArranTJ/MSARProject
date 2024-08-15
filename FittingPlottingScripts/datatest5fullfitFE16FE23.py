# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 12:26:44 2024

@author: Arran
"""

# __all__ = ['kappavalue']
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.cbook as cbook
import matplotlib.colors as colors
import astropy.units as u
from astropy.coordinates import SkyCoord, SpectralCoord
from astropy.wcs.utils import wcs_to_celestial_frame
import eispac

from astropy.visualization import astropy_mpl_style
from astropy.io import fits


#importing tools to define a kappa fit
from numpy import exp, linspace, random
from scipy.optimize import curve_fit


plt.style.use(astropy_mpl_style)


#define kappavalue to read into fittingfunctions.py
# kappavalue = 1

###Template parser to parse template data values and pass into fitting functions###
def templateparser(tmplt):
    pvallst = []
    plimslst = []
    pltdlst = []
    for i in range(len(tmplt.parinfo)):
        pvallst.append(tmplt.parinfo[i]['value'])
        pltdlst.append(tmplt.parinfo[i]['limited'])
        plimslst.append(tmplt.parinfo[i]['limits'])
        
    return pvallst, plimslst, pltdlst

###Defining the Kappa distribution function for fitting### 
# def kappafunc(x, amp, cen, wid, kappa):
#     #defining the variance-like variable indicated in Jeffreys. N. S et al (2016)
#     variance_k = wid * (1 - 3/(2*kappa))
#     #argument used in kappa function 
#     arg = (1 + ((x - cen)**2)/(2 * (variance_k**2) * kappa))
    
#     return amp * (arg) ** (-1 *kappa)

###defining kappafit for scipy curvefit###
#last arg for defining initial kappa guess when calling in the func
# def kfitter(x, y, tmplt, kappinit):
#     #parameter values
#     pval = templateparser(tmplt)[0]
#     #parameter limits
#     plims = templateparser(tmplt)[1]
#     #initial values for curve fits, inital kappa guess at end
#     init_val = [pval[0], pval[1], pval[2], kappinit]
    
#     #bounds=([0., plims[1][0] ,plims[2][0], 0.], [1e100, plims[1][1], plims[2][1], 1000])
    
#     #original curvefit vals w/ bounds
#     best_vals, covar = curve_fit(kappafunc, x, y, p0=init_val, bounds=([0., plims[1][0] ,plims[2][0], 0.], [1e100, plims[1][1], plims[2][1], 1000]))
    
#     return best_vals, covar



if __name__ == '__main__':
    # Read in the fit template and EIS observation
    # data_filepath = './data_eis/eis_20170906_114820.data.h5'
    # data_filepath = './data_eis/eis_20170906_115110.data.h5'
    # data_filepath = './data_eis/eis_20170906_115401.data.h5'
    # data_filepath = './data_eis/eis_20170906_115651.data.h5'
    data_filepath = './data_eis/eis_20170906_115942.data.h5'
    # data_filepath = './data_eis/eis_20170906_120233.data.h5'
    # data_filepath = './data_eis/eis_20170906_120524.data.h5'
    # data_filepath = './data_eis/eis_20170906_120815.data.h5'
    # data_filepath = './data_eis/eis_20170906_090025.data.h5'
    # data_filepath = './multiprocesstest/' + input()
    
    
    #testing template fitting for both lines
    
    # tmplt = eispac.EISFitTemplate(value=[ 19253.7287, 262.9824, 0.0341, 118.9812, 1601.5393,263.7539,0.0649, 342.9053], line_ids=['Fe XVI 262.984', 'Fe XXIII 263.760' ], wmin=262.78399658203125, wmax= 263.96002197265625)

    #original tmplt fit
    template_filepath1 = './fe_16_262_984.1c.template.h5'
    template_filepath2 = './fe_23_263_760.1c.template.h5'
    

    tmplt1 = eispac.read_template(template_filepath1)
    
    # #establishing new updated template for kappa parameter
    pvals1 = templateparser(tmplt1)[0]
    plims1 = templateparser(tmplt1)[1]
    pltd1 = templateparser(tmplt1)[2]
    new_tmplt1 = eispac.EISFitTemplate(value=[pvals1[0],pvals1[1], pvals1[2],pvals1[3], 1000], 
                                        limits = [[0.0000, 0.0000],[plims1[1][0], plims1[1][1]],[plims1[2][0], plims1[2][1]], [plims1[3][0], plims1[3][1]], [0, 1000]],
                                        limited = [[pltd1[0][0],pltd1[0][1]], [pltd1[1][0],pltd1[1][1]], [pltd1[2][0],pltd1[2][1]], [pltd1[3][0],pltd1[3][1]], [1,1]], 
                                        line_ids=['Fe XVI 262.984'], 
                                        wmin= 262.78399658203125, wmax=263.1000061035156)
    print(new_tmplt1)
    
    
    tmplt2 = eispac.read_template(template_filepath2)
    
    
    pvals2 = templateparser(tmplt2)[0]
    plims2 = templateparser(tmplt2)[1]
    pltd2 = templateparser(tmplt2)[2]
    new_tmplt2 = eispac.EISFitTemplate(value=[pvals2[0],pvals2[1], pvals2[2],pvals2[3], 1000], 
                                        limits = [[0.0000, 0.0000],[plims2[1][0], plims2[1][1]],[plims2[2][0], plims2[2][1]],[plims2[3][0], plims2[3][1]], [0, 1000]], 
                                        limited = [[pltd2[0][0],pltd2[0][1]], [pltd2[1][0],pltd2[1][1]], [pltd2[2][0],pltd2[2][1]], [pltd2[3][0],pltd2[3][1]], [1,1]],
                                        line_ids=['Fe XXIII 263.760'], 
                                        wmin= 263.45001220703125, wmax=263.96002197265625)
    
    print(new_tmplt2)
    

    Gaussdata_cube1 = eispac.read_cube(data_filepath, tmplt1.central_wave)
    data_cube1 = eispac.read_cube(data_filepath, new_tmplt1.central_wave)
    Gaussdata_cube2 = eispac.read_cube(data_filepath, tmplt2.central_wave)
    data_cube2 = eispac.read_cube(data_filepath, new_tmplt2.central_wave)
    
    # Select a cutout of the raster
    eis_frame = wcs_to_celestial_frame(Gaussdata_cube1.wcs)
    print(eis_frame)
    lower_left = [None, SkyCoord(Tx=500, Ty=-200, unit=u.arcsec, frame=eis_frame)]
    print(lower_left)
    upper_right = [None, SkyCoord(Tx=550, Ty=-160, unit=u.arcsec, frame=eis_frame)]

    
    # #original cutout for cropping
        
    # gauss_cutout1 = Gaussdata_cube1.crop(lower_left, upper_right)
    # gauss_cutout2 = Gaussdata_cube2.crop(lower_left, upper_right)
    
    
    # raster_cutout1 = data_cube1.crop(lower_left, upper_right)
    # raster_cutout2 = data_cube2.crop(lower_left, upper_right)
    
    
    #raster 'cutout' change to entire observation

    gauss_cutout1 = Gaussdata_cube1
    gauss_cutout2 = Gaussdata_cube2
    raster_cutout1 = data_cube1
    raster_cutout2 = data_cube2
    
    # # Fit the data and save it to disk

    # # Gaussian # # 
    fit_gauss1 = eispac.fit_spectra(gauss_cutout1, tmplt1, ncpu='max')
    save_filepathsgauss1 = eispac.save_fit(fit_gauss1, save_dir='cwd')
    
    
    fit_gauss2 = eispac.fit_spectra(gauss_cutout2, tmplt2, ncpu='max')
    save_filepathsgauss2 = eispac.save_fit(fit_gauss2, save_dir='cwd')
    
    # # Kappa # # 
    fit_res1 = eispac.fit_spectraKAPPA(raster_cutout1, new_tmplt1, ncpu='max')
    save_filepaths1 = eispac.save_fit(fit_res1, save_dir='cwd')
    
    
    fit_res2 = eispac.fit_spectraKAPPA(raster_cutout2, new_tmplt2, ncpu='max')
    save_filepaths2 = eispac.save_fit(fit_res2, save_dir='cwd')

    # # Find indices and world coordinates of max intensity
    

    # sum_data_inten1 = raster_cutout1.sum_spectra().data
    # sum_data_fullregion1 = data_cube1.sum_spectra().data  #use for full region image on subplot?
    # sum_data_inten2 = raster_cutout2.sum_spectra().data
    # sum_data_fullregion2 = data_cube2.sum_spectra().data  #use for full region image on subplot?


    sum_data_inten1 = gauss_cutout1.sum_spectra().data
    sum_data_fullregion1 = Gaussdata_cube1.sum_spectra().data  #use for full region image on subplot?
    sum_data_inten2 = gauss_cutout2.sum_spectra().data
    sum_data_fullregion2 = Gaussdata_cube2.sum_spectra().data  #use for full region image on subplot?
    
    def worldtoindex(raster, xarc, yarc):
        indexcoord = raster.wcs.world_to_array_index(SpectralCoord(2.627245243818844e-08, unit = u.m ), SkyCoord(Tx=xarc, Ty=yarc, unit=u.arcsec, frame=eis_frame))
        iy, ix = indexcoord[0], indexcoord[1]
        return iy, ix

    # print(worldtoindex(raster_cutout1, 550, -240)[0])
        
    coordx = np.linspace(500,600,20)
    coordy = np.linspace(-260,-160,20)
    


        
        
    # iy, ix = np.unravel_index(sum_data_inten1.argmax(), sum_data_inten1.shape)
    iy, ix = worldtoindex(gauss_cutout1, 549, -226)[0], worldtoindex(gauss_cutout1, 549, -226)[1]
    # iy, ix = worldtoindex(gauss_cutout1, 565, -210)[0], worldtoindex(gauss_cutout1, 565, -210)[1]
    # iy, ix = worldtoindex(gauss_cutout1, 565, -220)[0], worldtoindex(gauss_cutout1, 565, -220)[1]
    # iy, ix = worldtoindex(gauss_cutout1, 558, -235)[0], worldtoindex(gauss_cutout1, 558, -235)[1]
    # iy, ix = worldtoindex(gauss_cutout1, 530, -160)[0], worldtoindex(gauss_cutout1, 530, -160)[1]
    # iy, ix = worldtoindex(gauss_cutout1, 530, -180)[0], worldtoindex(gauss_cutout1, 530, -180)[1]
    # iy, ix = worldtoindex(gauss_cutout1, 544, -222)[0], worldtoindex(gauss_cutout1, 544, -222)[1]
    # iy, ix = worldtoindex(gauss_cutout1, 580, -200)[0], worldtoindex(gauss_cutout1, 580, -200)[1]
    
    worldcoord = gauss_cutout1.wcs.array_index_to_world(iy, ix, 0)
    # indexcoord = raster_cutout1.wcs.world_to_array_index(SpectralCoord(2.627245243818844e-08, unit = u.m ), SkyCoord(Tx=550, Ty=-226, unit=u.arcsec, frame=eis_frame))
    print(iy, ix)
    print(worldcoord)
    # print(indexcoord)
    # iy, ix = np.unravel_index(np.argsort(sum_data_inten1.ravel(), axis=0)[-30], sum_data_inten1.shape)


    ex_world_coords1 = gauss_cutout1.wcs.array_index_to_world(iy, ix, 0)[1]
    print(ex_world_coords1)
    ex_world_coords2 = gauss_cutout2.wcs.array_index_to_world(iy, ix, 0)[1]
    
    
    # y_arcsec1, x_arcsec1 = ex_world_coords1.Ty.value, ex_world_coords1.Tx.value
    # y_arcsec2, x_arcsec2 = ex_world_coords2.Ty.value, ex_world_coords2.Tx.value
    
    y_arcsec, x_arcsec = ex_world_coords1.Ty.value, ex_world_coords1.Tx.value
    
    
    # print(y_arcsec1, y_arcsec2)
    # print(x_arcsec1, x_arcsec2)
    # print(ex_world_coords1,ex_world_coords2)
    
    # Extract data profile and interpolate fit at higher spectral resolution
    
    ##Gauss Data##
    
    data_x1 = gauss_cutout1.wavelength[iy, ix, :]
    data_y1 = gauss_cutout1.data[iy, ix, :]
    data_err1 = gauss_cutout1.uncertainty.array[iy, ix, :]
    ##replace negative err. values with 0 to permit plotting for affected pixels##
    data_err1[data_err1 < 0] = 0
    # print(data_err1)

    fit_x1, fit_y1 = fit_gauss1.get_fit_profile(coords=[iy,ix], num_wavelengths=100)
    chi2fit1 = fit_gauss1.fit['chi2'][iy,ix]
    # print(chi2fit1)
    # c0_x1, c0_y1 = fit_res1.get_fit_profile(0, coords=[iy,ix], num_wavelengths=100)
    
    data_x2 = gauss_cutout2.wavelength[iy, ix, :]
    data_y2 = gauss_cutout2.data[iy, ix, :]
    data_err2 = gauss_cutout2.uncertainty.array[iy, ix, :]
    ##replace negative err. values with 0 to permit plotting for affected pixels##
    data_err2[data_err2 < 0] = 0
    fit_x2, fit_y2 = fit_gauss2.get_fit_profile(coords=[iy,ix], num_wavelengths=100)
    chi2fit2 = fit_gauss2.fit['chi2'][iy,ix]
    # print(chi2fit2)
    # c0_x2, c0_y2 = fit_res2.get_fit_profile(0, coords=[iy,ix], num_wavelengths=100)


    ##Kappa Data##


    KAPdata_x1 = raster_cutout1.wavelength[iy, ix, :]
    KAPdata_y1 = raster_cutout1.data[iy, ix, :]
    KAPdata_err1 = raster_cutout1.uncertainty.array[iy, ix, :]
    ##replace negative err. values with 0 to permit plotting for affected pixels##
    KAPdata_err1[KAPdata_err1 < 0] = 0
    # print(KAPdata_err1)

    KAPfit_x1, KAPfit_y1 = fit_res1.get_fit_profile(coords=[iy,ix], num_wavelengths=100)
    KAPchi2fit1 = fit_res1.fit['chi2'][iy,ix]
    # print(KAPchi2fit1)
    kappaparam1 = fit_res1.fit['params'][iy,ix][4]
    print(kappaparam1)
    kappaerr1 = fit_res1.fit['perror'][iy,ix]
    print(kappaerr1)
    # c0_x1, c0_y1 = fit_res1.get_fit_profile(0, coords=[iy,ix], num_wavelengths=100)
    
    KAPdata_x2 = raster_cutout2.wavelength[iy, ix, :]
    KAPdata_y2 = raster_cutout2.data[iy, ix, :]
    KAPdata_err2 = raster_cutout2.uncertainty.array[iy, ix, :]
    ##replace negative err. values with 0 to permit plotting for affected pixels##
    KAPdata_err2[KAPdata_err2 < 0] = 0
    KAPfit_x2, KAPfit_y2 = fit_res2.get_fit_profile(coords=[iy,ix], num_wavelengths=100)
    KAPchi2fit2 = fit_res2.fit['chi2'][iy,ix]
    # print(KAPchi2fit2)
    kappaparam2 = fit_res2.fit['params'][iy,ix][4]
    kappaerr2 = fit_res2.fit['perror'][iy,ix]
    print(kappaerr1)
    print(kappaparam2)
    # c0_x2, c0_y2 = fit_res2.get_fit_profile(0, coords=[iy,ix], num_wavelengths=100)
                                         
    # c1_x, c1_y = fit_res.get_fit_profile(1, coords=[iy,ix], num_wavelengths=100)
    # c2_x, c2_y = fit_res.get_fit_profile(2, coords=[iy,ix], num_wavelengths=100)
    
    # Make a multi-panel figure with the cutout and example profile

    import matplotlib.gridspec as gridspec
    
    # fig = plt.figure(figsize=[10,5])
    # plot_grid = fig.add_gridspec(nrows=2, ncols=4, wspace=0.3)
    
    #####First Figure window for intensity and cutout#####
    fig = plt.figure(figsize=[10,10])
    
    
    
    plot_grid = fig.add_gridspec(nrows=2, ncols=2, wspace=0.3)
    
    
    # data_subplt = fig.add_subplot(plot_grid[0,0])
    # data_subplt.imshow(sum_data_inten, origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt = fig.add_subplot(plot_grid[0,0])
    data_subplt.scatter(x_arcsec, y_arcsec, color='r', marker='x')
    # data_subplt.imshow(sum_data_inten, origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt.imshow(sum_data_fullregion1, norm=colors.PowerNorm(gamma=0.2) ,origin='lower', extent = Gaussdata_cube1.meta['extent_arcsec'])
    data_subplt.set_title('Region Data FeXVI\n'+' ' +gauss_cutout1.meta['mod_index']['date_obs'])
    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[0,1])
    data_subplt.scatter(x_arcsec, y_arcsec, color='r', marker='x')
    # data_subplt.imshow(sum_data_inten, origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt.imshow(sum_data_inten1, norm=colors.PowerNorm(gamma=0.2) , origin='lower', extent = gauss_cutout1.meta['extent_arcsec'])
    data_subplt.set_title('Selected Cutout\n FeXVI')
    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    
    
    data_subplt = fig.add_subplot(plot_grid[1,0])
    data_subplt.scatter(x_arcsec, y_arcsec, color='r', marker='x')
    # data_subplt.imshow(sum_data_inten, origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt.imshow(sum_data_fullregion2, norm=colors.PowerNorm(gamma=0.2) ,origin='lower', extent = Gaussdata_cube2.meta['extent_arcsec'])
    data_subplt.set_title('Region Data FeXXIII\n'+' ' +gauss_cutout2.meta['mod_index']['date_obs'])
    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[1,1])
    data_subplt.scatter(x_arcsec, y_arcsec, color='r', marker='x')
    # data_subplt.imshow(sum_data_inten, origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt.imshow(sum_data_inten2, norm=colors.PowerNorm(gamma=0.2) , origin='lower', extent = gauss_cutout2.meta['extent_arcsec'])
    data_subplt.set_title('Selected Cutout\n FeXXIII')
    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')

    
    #####2nd Figure window for intensity, vel, width, and spectra#####


        #2nd Figure window for intensity, vel, width, and spectra
    from matplotlib.gridspec import SubplotSpec

    fig = plt.figure(figsize=[20,10])
    plot_grid = fig.add_gridspec(nrows=2, ncols=4, wspace=0.75, hspace = 0.25)
    


    #Subplots for Fe16

    data_subplt = fig.add_subplot(plot_grid[0,0])
    inten_map1 = fit_gauss1.get_map(component=0, measurement='intensity')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(inten_map1.data, norm=colors.PowerNorm(gamma=0.4, vmax  = 200000, vmin = 0) , origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap='PuBu_r')
    #define position for subplt colourbar
    pos1 = data_subplt.imshow(inten_map1.data, norm=colors.PowerNorm(gamma=0.4, vmax  = 200000, vmin = 0) , origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap='PuBu_r')
    
    data_subplt.scatter(x_arcsec, y_arcsec, color='r', marker='x')
    data_subplt.set_title('Selected Cutout FeXVI \n [SunPy Intensity]')
#Insert fraction=0.046, pad=0.04 to scale the colorbars appropriately
    fig.colorbar(pos1, fraction=0.046, pad=0.04, label = 'Intensity [erg / (s sr cm2)]')

    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[0,1])
    vel_map1 = fit_gauss1.get_map(component=0, measurement='vel')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(vel_map1.data,origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap = 'RdBu', vmax = 50, vmin = -50)
    #define position for subplt colourbar
    pos2 = data_subplt.imshow(vel_map1.data,origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap = 'RdBu', vmax = 50, vmin = -50)
    data_subplt.set_title('Selected Cutout FeXVI \n [SunPy Velocity]')
    fig.colorbar(pos2, fraction=0.046, pad=0.04, label = 'km/s')

    data_subplt.set_xlabel('Solar-X [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[0,2])
    width_map1 = fit_gauss1.get_map(component=0, measurement='width')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(width_map1.data,origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap = 'viridis', vmax = 0.05, vmin = 0.02)
    #define position for subplt colourbar
    pos3 = data_subplt.imshow(width_map1.data,origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap = 'viridis', vmax = 0.05, vmin = 0.02)
    data_subplt.set_title('Selected Cutout FeXVI \n [SunPy Width]')
    fig.colorbar(pos3,fraction=0.046, pad=0.04, label = 'Wavelength [$\AA$]')

    data_subplt.set_xlabel('Solar-X [arcsec]')



    #Subplots for Fe23

    data_subplt = fig.add_subplot(plot_grid[1,0])
    inten_map2 = fit_gauss2.get_map(component=0, measurement='intensity')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(inten_map2.data, norm=colors.PowerNorm(gamma=0.4, vmax  = 250000, vmin = 0) , origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap='PuBu_r')
    #define position for subplt colourbar
    pos4 = data_subplt.imshow(inten_map2.data, norm=colors.PowerNorm(gamma=0.4, vmax  = 250000, vmin = 0) , origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap='PuBu_r')
    data_subplt.scatter(x_arcsec, y_arcsec, color='r', marker='x')
    data_subplt.set_title('Selected Cutout FeXXIII \n [SunPy Intensity]')
#Insert fraction=0.046, pad=0.04 to scale the colorbars appropriately
    fig.colorbar(pos4, fraction=0.046, pad=0.04, label = 'Intensity [erg / (s sr cm2)]')

    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[1,1])
    vel_map2 = fit_gauss2.get_map(component=0, measurement='vel')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(vel_map2.data,origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap = 'RdBu', vmax = 50, vmin = -50)
    #define position for subplt colourbar
    pos5 = data_subplt.imshow(vel_map2.data,origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap = 'RdBu', vmax = 50, vmin = -50)
    data_subplt.set_title('Selected Cutout FeXXIII \n [SunPy Velocity]')
    fig.colorbar(pos5, fraction=0.046, pad=0.04, label = 'km/s')

    data_subplt.set_xlabel('Solar-X [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[1,2])
    width_map2 = fit_gauss2.get_map(component=0, measurement='width')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(width_map2.data,origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap = 'viridis', vmax = 0.08, vmin = 0.02)
    #define position for subplt colourbar
    pos6 = data_subplt.imshow(width_map2.data,origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap = 'viridis', vmax = 0.08, vmin = 0.02)
    data_subplt.set_title('Selected Cutout FeXXIII \n [SunPy Width]')
    fig.colorbar(pos6,fraction=0.046, pad=0.04, label = 'Wavelength [$\AA$]')


    data_subplt.set_xlabel('Solar-X [arcsec]')

#plotting spectra FEXVI
    profile_subplt = fig.add_subplot(plot_grid[0,3])
    profile_subplt.errorbar(data_x1, data_y1, yerr=data_err1, ls='', marker='o', color='k')
    profile_subplt.plot(data_x1, data_y1, ls='', marker='o', color='k')
    profile_subplt.plot(fit_x1, fit_y1, color='b', label='Combined profile FEXVI')
    # profile_subplt.plot(c0_x1, c0_y1, color='r', label=fit_res1.fit['line_ids'][0])
#plotting spectra FEXXIII
    profile_subplt.errorbar(data_x2, data_y2, yerr=data_err2, ls='', marker='o', color='k')
    profile_subplt.plot(data_x2, data_y2, ls='', marker='o', color='k')
    profile_subplt.plot(fit_x2, fit_y2, color='m', label='Combined profile FEXXIII')
    # profile_subplt.plot(c0_x2, c0_y2, color='g', label=fit_res2.fit['line_ids'][0])

    # profile_subplt.plot(c1_x, c1_y, color='r', ls='--', label=fit_res.fit['line_ids'][1])
    # profile_subplt.plot(c2_x, c2_y, color='g', label='Background')
    # profile_subplt.set_title(f'Cutout indices: iy = {iy}, ix = {ix}')
    profile_subplt.set_title('Spectra for Heliocentric Coordinates \n' f'{x_arcsec}, {y_arcsec}')
    profile_subplt.set_xlabel('Wavelength [$\AA$]')
    profile_subplt.set_ylabel('Intensity\n ['+gauss_cutout1.unit.to_string()+']')
    profile_subplt.legend(loc='upper left', frameon=False)
    profile_subplt.set_ylim(0)
    profile_subplt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    data_subplt.set_xlabel('Solar-X [arcsec]')

#Residuals plotting
    fit_x_vals1, fit_y_vals1 = fit_gauss1.get_fit_profile(coords=[iy,ix],use_mask=True) #use len(data_y) to get compatible array length for residuals, if necessary
    fit_x_vals2, fit_y_vals2 = fit_gauss2.get_fit_profile(coords=[iy,ix],use_mask=True) #use len(data_y) to get compatible array length for residuals, if necessary
   

# print(fit_x_vals)
    resids1 = fit_y_vals1 - data_y1
    resids2 = fit_y_vals2 - data_y2
    profile_subplt2 = fig.add_subplot(plot_grid[1,3], sharex= profile_subplt)
    
    plt.errorbar(x=fit_x_vals1, y=resids1, yerr=data_err1, fmt='.', color ='b', label=fit_gauss1.fit['line_ids'][0])
    plt.errorbar(x=fit_x_vals2, y=resids2, yerr=data_err2, fmt='.', color ='m', label=fit_gauss2.fit['line_ids'][0])
    plt.axhline(0.0, linestyle='--', color= 'black', linewidth=1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc='upper left')
    plt.ylabel('Residual Intensity \n [erg / (s sr cm2)]')
    plt.xlabel('Wavelength [$\AA$]')

    
    #####3rd Fig Window for spectra and residuals only#####
     #3rd Fig Window for spectra and residuals only
    fig = plt.figure(figsize=[20,10])
    
    
    
    #Define values for residuals and plot spectra


    ###Gauss Fits###

    #Define values for residuals and plot spectra
    
    fit_x_vals1, fit_y_vals1 = fit_gauss1.get_fit_profile(coords=[iy,ix],use_mask=True) #use len(data_y) to get compatible array length for residuals, if necessary
    fit_x_vals2, fit_y_vals2 = fit_gauss2.get_fit_profile(coords=[iy,ix],use_mask=True) #use len(data_y) to get compatible array length for residuals, if necessary
    
    plot_grid = fig.add_gridspec(nrows=2, ncols=2, wspace=0.05, hspace = 0.10)
    
    profile_subplt = fig.add_subplot(plot_grid[0,0])
    profile_subplt.errorbar(data_x1, data_y1, yerr=data_err1, ls='', marker='o', color='k')
    profile_subplt.plot(data_x1, data_y1, ls='', marker='o', color='k')
    profile_subplt.plot(fit_x1, fit_y1, color='b', label='Combined profile FEXVI')
    # profile_subplt.plot(c0_x1, c0_y1, color='r', label=fit_res1.fit['line_ids'][0])
    #plotting spectra FEXXIII
    profile_subplt.errorbar(data_x2, data_y2, yerr=data_err2, ls='', marker='o', color='k')
    profile_subplt.plot(data_x2, data_y2, ls='', marker='o', color='k')
    profile_subplt.plot(fit_x2, fit_y2, color='m', label='Combined profile FEXXIII')
    # profile_subplt.plot(c0_x2, c0_y2, color='g', label=fit_res2.fit['line_ids'][0])
    
    # profile_subplt.plot(c1_x, c1_y, color='r', ls='--', label=fit_res.fit['line_ids'][1])
    # profile_subplt.plot(c2_x, c2_y, color='g', label='Background')
    # profile_subplt.set_title(f'Cutout indices: iy = {iy}, ix = {ix}')
    profile_subplt.set_title('Gaussian Fit Spectra for Helio-centric Coordinates \n' f'{x_arcsec}, {y_arcsec}')
    profile_subplt.set_xlabel('Wavelength [$\AA$]')
    profile_subplt.set_ylabel('Intensity \n ['+raster_cutout1.unit.to_string()+']')
    profile_subplt.legend(loc='upper left', frameon=False)
    
    # data_subplt.set_xlabel('Solar-X [arcsec]')
    # data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    #Remove tick labels from upper spectra to 'blend' with residual plot
    profile_subplt.label_outer()
    
    
    
    #Set up and plot residuals
    resids1 = fit_y_vals1 - data_y1
    resids2 = fit_y_vals2 - data_y2
    # print(resids1)
    # print(np.sum(resids1**2/fit_y_vals1))

# ###chi2 testing
# #data[data.mask == False]
#     from scipy import stats
#     print(fit_y_vals1)
#     fit_y_unmasked = fit_y_vals1[fit_y_vals1.mask == False]
#     print(fit_y_unmasked)
#     data_ymask = data_y1[fit_y_vals1.mask == True]
#     print(data_y1)
#     print(data_y1[3:][:14])
#     print(sum(data_y1[3:][:14]))
# ####    
          
    profile_subplt2 = fig.add_subplot(plot_grid[1,0], sharex= profile_subplt)
    
    plt.errorbar(x=fit_x_vals1, y=resids1, yerr=data_err1, fmt='.', color ='b', label=fit_gauss1.fit['line_ids'][0] + '\n '+ '$\chi^{2}$:' + ' ' + "{:.3f}".format(chi2fit1))
    plt.errorbar(x=fit_x_vals2, y=resids2, yerr=data_err2, fmt='.', color ='m', label=fit_gauss2.fit['line_ids'][0] + '\n '+ '$\chi^{2}$:' + ' ' + "{:.3f}".format(chi2fit2))
    plt.axhline(0.0, linestyle='--', color= 'black', linewidth=1)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc='upper left')
    plt.ylabel('Residual Intensity \n [erg / (s sr cm2)]')
    plt.xlabel('Wavelength [$\AA$]')
    #define limits for residual axis
    plt.ylim(-4000,3300)



    ######Kappa Fits #######
    ########################
    #Define values for residuals and plot spectra


   
    kfit_x_vals1, kfit_y_vals1 = fit_res1.get_fit_profile(coords=[iy,ix],use_mask=True) #use len(data_y) to get compatible array length for residuals, if necessary
    kfit_x_vals2, kfit_y_vals2 = fit_res2.get_fit_profile(coords=[iy,ix],use_mask=True) #use len(data_y) to get compatible array length for residuals, if necessary
    
    profile_subplt3 = fig.add_subplot(plot_grid[0,1], sharey= profile_subplt)
    profile_subplt3.errorbar(KAPdata_x1, KAPdata_y1, yerr=KAPdata_err1, ls='', marker='o', color='k')
    profile_subplt3.plot(KAPdata_x1, KAPdata_y1, ls='', marker='o', color='k')
    profile_subplt3.plot(KAPfit_x1, KAPfit_y1, color='b', label='Combined profile FEXVI' + '\n '+ '$\kappa$:' + ' ' + "{:.3f}".format(kappaparam1))
    # profile_subplt.plot(c0_x1, c0_y1, color='r', label=fit_res1.fit['line_ids'][0])
    #plotting spectra FEXXIII
    profile_subplt3.errorbar(KAPdata_x2, KAPdata_y2, yerr=KAPdata_err2, ls='', marker='o', color='k')
    profile_subplt3.plot(KAPdata_x2, KAPdata_y2, ls='', marker='o', color='k')
    profile_subplt3.plot(KAPfit_x2, KAPfit_y2, color='m', label='Combined profile FEXXIII' + '\n '+ '$\kappa$:' + ' ' + "{:.3f}".format(kappaparam2))

    profile_subplt3.set_title('Kappa Fit Spectra for Helio-centric Coordinates \n' f'{x_arcsec}, {y_arcsec}')
    profile_subplt3.set_xlabel('Wavelength [$\AA$]')
    profile_subplt3.set_ylabel('Intensity \n ['+raster_cutout1.unit.to_string()+']')
    profile_subplt3.legend(loc='upper left', frameon=False)
    
    #Remove tick labels from upper spectra to 'blend' with residual plot
    profile_subplt3.label_outer()
    plt.ylim(0)




    #Set up and plot residuals
    kresids1 = kfit_y_vals1 - KAPdata_y1
    kresids2 = kfit_y_vals2 - KAPdata_y2

    profile_subplt4 = fig.add_subplot(plot_grid[1,1], sharex= profile_subplt3, sharey= profile_subplt2)
    
    profile_subplt4.errorbar(x=kfit_x_vals1, y=kresids1, yerr=KAPdata_err1, fmt='.', color ='b', label=fit_res1.fit['line_ids'][0] + '\n '+ r'$\chi_{\kappa}^{2}$:' + ' ' + "{:.3f}".format(KAPchi2fit1))
    profile_subplt4.errorbar(x=kfit_x_vals2, y=kresids2, yerr=KAPdata_err2, fmt='.', color ='m', label=fit_res2.fit['line_ids'][0] + '\n '+ r'$\chi_{\kappa}^{2}$:' + ' ' + "{:.3f}".format(KAPchi2fit2))
    profile_subplt4.axhline(0.0, linestyle='--', color= 'black', linewidth=1)
    # profile_subplt4.yaxis.set_visible(False)
    # profile_subplt4.get_yticklabels(False)
    ####Hide ticklabels for y axis
    plt.setp(profile_subplt4.get_yticklabels(), visible=False)
    # profile_subplt4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    profile_subplt4.legend(loc='upper left')
    profile_subplt4.set_xlabel('Wavelength [$\AA$]')

    
    
    ####Making plots for Intensities ONLY###
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=[10,8]) 
    plt.suptitle(gauss_cutout1.meta['mod_index']['date_obs'][0:10] + '\n' +gauss_cutout1.meta['mod_index']['date_obs'][10:])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=None)
    # plot_grid = fig.add_gridspec(nrows=1, ncols=2, wspace=0.4)
        
    ####FE16####

    ax1.set_xlabel('Solar-X [arcsec]')
    ax1.set_ylabel('Solar-Y [arcsec]')
    
    ax1.set_title('FeXVI Intensity')
    fe16intensity = ax1.imshow(inten_map1.data, norm=colors.PowerNorm(gamma=0.4, vmax  = 200000, vmin = 0) , origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap='PuBu_r')

    
    ####FE23#####
    
    ax2.set_xlabel('Solar-X [arcsec]')
    # ax2.set_ylabel('Solar-Y [arcsec]')
    
    ax2.set_title('FeXXIII Intensity')
    fe23intensity = ax2.imshow(inten_map2.data, norm=colors.PowerNorm(gamma=0.4, vmax  = 200000, vmin = 0) , origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap='PuBu_r')

    fmt = plt.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    cbar = fig.colorbar(fe16intensity, ax=(ax1,ax2), orientation= 'horizontal', fraction=0.146 ,pad=0.1,label = "Intensity [erg / (s sr $cm^{2}$)]", format = fmt)
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    
    plt.savefig('intenplots/'+ gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), dpi = 300)
    
    ###Saving intensity arrays###
    
    np.save('datainten/1/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), inten_map1.data)
    np.save('datainten/2/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), inten_map2.data)
    
    ###Saving vel arrays###
    
    np.save('datavel/1/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), vel_map1.data)
    np.save('datavel/2/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), vel_map2.data)
    
    ###Saving width arrays###
    
    np.save('datawidth/1/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), width_map1.data)
    np.save('datawidth/2/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), width_map2.data)
    
    
    
    #####PLOTTING CHI2 AND KAPPA MAPS###
    fig = plt.figure(figsize=[15,10])
        
    plt.suptitle(gauss_cutout1.meta['mod_index']['date_obs'][0:10] + '\n' +gauss_cutout1.meta['mod_index']['date_obs'][10:])
        
    plot_grid = fig.add_gridspec(nrows=2, ncols=3, hspace=0.25, wspace=0.25)
        
    ####FE16 KAPPA####
    subplt = fig.add_subplot(plot_grid[0,0])
    subplt.set_xlabel('Solar-X [arcsec]')
    subplt.set_ylabel('Solar-Y [arcsec]')
    
    
    #     print(fit_res1.fit['params'][iy, ix][4])
    
    
    # #slice the 3d array to select the last value (kappa) from each sub-array
    #     print(fit_res1.fit['params'][:,:,-1][iy,ix])
    
    ###saving np arrays for kappas, kappa chi2, gauss chi2###
    np.save('datakappa/1/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), fit_res1.fit['params'][:,:,-1])
    np.save('datakappa/2/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), fit_res2.fit['params'][:,:,-1])
    
    np.save('datachi2/1/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), fit_res1.fit['chi2'])
    np.save('datachi2/2/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), fit_res2.fit['chi2'])
    
    np.save('datachi2gauss/1/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), fit_gauss1.fit['chi2'])
    np.save('datachi2gauss/2/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), fit_gauss2.fit['chi2'])
    
    
    subplt.set_title('FEXVI $\kappa$ Map')
    subplt.imshow(fit_res1.fit['params'][:,:,-1], origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'CMRmap', vmax = 20, vmin = 0.0)
    pos = subplt.imshow(fit_res1.fit['params'][:,:,-1], origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'CMRmap', vmax = 20, vmin = 0.0)
    
    fig.colorbar(pos,fraction=0.046, pad=0.04, label = "$\kappa$")
    
    ####FE23 KAPPA#####
    subplt = fig.add_subplot(plot_grid[1,0])
    subplt.set_xlabel('Solar-X [arcsec]')
    subplt.set_ylabel('Solar-Y [arcsec]')
    
    subplt.set_title('FEXXIII $\kappa$ Map')
    subplt.imshow(fit_res2.fit['params'][:,:,-1], origin='lower', extent = raster_cutout2.meta['extent_arcsec'], cmap = 'CMRmap', vmax = 20, vmin = 0.0)
    pos = subplt.imshow(fit_res2.fit['params'][:,:,-1], origin='lower', extent = raster_cutout2.meta['extent_arcsec'], cmap = 'CMRmap', vmax = 20, vmin = 0.0)
    
    fig.colorbar(pos,fraction=0.046, pad=0.04, label = "$\kappa$")
    
    
    ####FE16 CHI2####
    
    subplt = fig.add_subplot(plot_grid[0,1])
    subplt.set_xlabel('Solar-X [arcsec]')
    
    subplt.set_title('FEXVI' + ' ' + r'$\chi_{\kappa}^{2}$' + ' ' +'Map')
    subplt.imshow(fit_res1.fit['chi2'], norm=colors.PowerNorm(gamma=0.4, vmax=1, vmin=0), origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'plasma')
    pos = subplt.imshow(fit_res1.fit['chi2'], norm=colors.PowerNorm(gamma=0.4, vmax=1, vmin=0), origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'plasma')
    plt.setp(subplt.get_yticklabels(), visible=False)
    fig.colorbar(pos,fraction=0.046, pad=0.04, label = r'$\chi_{\kappa}^{2}$')
    
    
    ####FE23 CHI2####
    subplt = fig.add_subplot(plot_grid[1,1])
    subplt.set_xlabel('Solar-X [arcsec]')
    
    subplt.set_title('FEXXIII' + ' ' + r'$\chi_{\kappa}^{2}$' + ' ' +'Map' )
    subplt.imshow(fit_res2.fit['chi2'], norm=colors.PowerNorm(gamma=0.4, vmax=1, vmin=0), origin='lower', extent = raster_cutout2.meta['extent_arcsec'], cmap = 'plasma')
    pos = subplt.imshow(fit_res2.fit['chi2'], norm=colors.PowerNorm(gamma=0.4, vmax=1, vmin=0), origin='lower', extent = raster_cutout2.meta['extent_arcsec'], cmap = 'plasma')
    plt.setp(subplt.get_yticklabels(), visible=False)
    fig.colorbar(pos,fraction=0.046, pad=0.04, label = r'$\chi_{\kappa}^{2}$')
    
    ####FE16 Gauss CHI2####

    subplt = fig.add_subplot(plot_grid[0,2])
    subplt.set_xlabel('Solar-X [arcsec]')
    
    subplt.set_title('FEXVI' + ' ' + r'$\chi_{G}^{2}$' + ' ' +'Map (Gaussian)')
    subplt.imshow(fit_gauss1.fit['chi2'], norm=colors.PowerNorm(gamma=0.4, vmax=1, vmin=0), origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'plasma')
    pos = subplt.imshow(fit_gauss1.fit['chi2'], norm=colors.PowerNorm(gamma=0.4, vmax=1, vmin=0), origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'plasma')
    plt.setp(subplt.get_yticklabels(), visible=False)
    fig.colorbar(pos,fraction=0.046, pad=0.04, label = r'$\chi_{G}^{2}$')
    
    
    ####FE23 Gauss CHI2####
    subplt = fig.add_subplot(plot_grid[1,2])
    subplt.set_xlabel('Solar-X [arcsec]')
    
    subplt.set_title('FEXXIII' + ' ' + r'$\chi_{G}^{2}$' + ' ' +'Map (Gaussian)' )
    subplt.imshow(fit_gauss2.fit['chi2'], norm=colors.PowerNorm(gamma=0.4, vmax=1, vmin=0), origin='lower', extent = raster_cutout2.meta['extent_arcsec'], cmap = 'plasma')
    pos = subplt.imshow(fit_gauss2.fit['chi2'], norm=colors.PowerNorm(gamma=0.4, vmax=1, vmin=0), origin='lower', extent = raster_cutout2.meta['extent_arcsec'], cmap = 'plasma')
    plt.setp(subplt.get_yticklabels(), visible=False)
    fig.colorbar(pos,fraction=0.046, pad=0.04, label = r'$\chi_{G}^{2}$')


    
    
    plt.savefig('kappachiplots/'+ gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), dpi = 300)
    
    # plt.show()
    
    
    ###Flattened Intensity Array 
    flat_inten1 = np.ravel(inten_map1.data)
    flat_inten2 = np.ravel(inten_map2.data)
    
    flat_reschi1 = np.ravel(fit_res1.fit['chi2'])
    flat_gausschi1 = np.ravel(fit_gauss1.fit['chi2'])
    flat_kap1 = np.ravel(fit_res1.fit['params'][:,:,-1])
    
    flat_reschi2 = np.ravel(fit_res2.fit['chi2'])
    flat_gausschi2 = np.ravel(fit_gauss2.fit['chi2'])
    flat_kap2 = np.ravel(fit_res2.fit['params'][:,:,-1])
    
    
    
    #Filtering Function, pass flattened arrays into the function (fit_res is not flattened), outputs resized filtered kappa 
    def kappamap_filt(gausschi, reschi, kap, inten, fit_res):
        kapcontain = []
        for i in range(len(kap)):
            #change kappa values to higher vals if gaussian fit is better
            if gausschi[i] < reschi[i]:
                kapcontain.append(50)
            #change kappa values to higher vals if chi2 is 0 (chi2 = 0 is not reasonable)
            elif reschi[i] == 0:
                kapcontain.append(50)
            #change kappa values to higher vals if chi2 is >= 1
            elif reschi[i] >= 1:
                kapcontain.append(50)
            #change kappa values to higher vals if kappa is 0 (kappa = 0 is not reasonable)
            elif kap[i] == 0:
                kapcontain.append(50)
            #change kappa values to higher vals if intensities are lower than 1000 inten
            elif inten[i] < 1000:
                kapcontain.append(50)
            #preserve 'good' kappa values in the array
            else:
                kapcontain.append(kap[i])
    
        return np.resize(np.array(kapcontain), fit_res.fit['params'][:,:,-1].shape)
    
    
    
    
    kap_filt1 = kappamap_filt(flat_gausschi1, flat_reschi1, flat_kap1, flat_inten1, fit_res1)
    kap_filt2 = kappamap_filt(flat_gausschi2, flat_reschi2, flat_kap2, flat_inten2, fit_res2)
    
    ###Save Numpy Arrays###
    np.save('datakappafilt/1/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), kap_filt1)
    np.save('datakappafilt/2/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), kap_filt2)
    
    #function to provide the chi2 vals for the filtered kappa values
    
    def chi2_filt(reschi, kapfilt):
        chifiltcontain = []
        for i in range(len(reschi)):
            if np.ravel(kapfilt)[i] <= 20:
                chifiltcontain.append(reschi[i])
        #Ensures that if no corresponding chi2 vals are found for good kappas, container becomes 0 so np.mean() works with it
        #NOTE: WHEN GOING THROUGH CHI2 FOR GOOD VALS, EXCLUDE 0 CHI2 VALUES WHEN AVERAGING!!!!
        if chifiltcontain == []:
            chifiltcontain = 0
        return chifiltcontain
    
    
    #Function to produce container of all vals where kappa <= 20:
    def goodkapvals(kapfilt):
        goodkapcontain = []
        for i in range(len(np.ravel(kapfilt))):
            if np.ravel(kapfilt)[i] <= 20:
                goodkapcontain.append(np.ravel(kapfilt)[i])
                
        #NOTE: WHEN GOING THROUGH KAPPAS FOR GOOD VALS, EXCLUDE 0 KAPPA VALUES WHEN AVERAGING!!!!
        if goodkapcontain == []:
            goodkapcontain = 0
        return goodkapcontain
    
        
    goodchi1 = chi2_filt(flat_reschi1, kap_filt1)
    goodchi2 = chi2_filt(flat_reschi2, kap_filt2)
    
    ###Counting Number of 'good' kappa fits in the array
    goodkapct1 = np.count_nonzero(kap_filt1 <= 20)
    goodkapct2 = np.count_nonzero(kap_filt2 <= 20)
    
    
    ###Assigning vars to goodkappa value container
    goodkapval1 = goodkapvals(kap_filt1)
    goodkapval2 = goodkapvals(kap_filt2)
    
    # print("goodkapval1", goodkapval1)
    # print("goodkapval2", goodkapval2)
    
    #goodvalue array recording number of good kappas for Fe16, Fe23, Mean Kappas for Fe16, Fe23, Mean Chi2 for Fe16, Fe23, and the recorded observation time
    GoodVals = np.array([[goodkapct1, goodkapct2, np.mean(goodkapval1), np.mean(goodkapval2), np.mean(goodchi1), np.mean(goodchi2), gauss_cutout1.meta['mod_index']['date_obs'][11:19]]])
    
    print("GoodVals: ", GoodVals)
    
    np.savetxt('goodkappavals/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':','')+'.csv', GoodVals, fmt='%s', delimiter=',')
    
    print(kap_filt1.shape)
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=[10,8]) 
    plt.suptitle(gauss_cutout1.meta['mod_index']['date_obs'][0:10] + '\n' +gauss_cutout1.meta['mod_index']['date_obs'][10:])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=None)

        
    ####FE16 KAPPA####
    ax1.set_xlabel('Solar-X [arcsec]')
    ax1.set_ylabel('Solar-Y [arcsec]')
    
    ax1.set_title('FeXVI $\kappa$ Map')
    fe16kappa = ax1.imshow(kap_filt1, origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'CMRmap', vmax = 20, vmin = 0.0)

    
    ####FE23 KAPPA#####
    
    ax2.set_xlabel('Solar-X [arcsec]')

    ax2.set_title('FeXXIII $\kappa$ Map')
    ax2.imshow(kap_filt2, origin='lower', extent = raster_cutout2.meta['extent_arcsec'], cmap = 'CMRmap', vmax = 20, vmin = 0.0)

    
    cbar = fig.colorbar(fe16kappa, ax=(ax1,ax2), orientation= 'horizontal', fraction=0.146, pad=0.1,label = "$\kappa$")
    
    plt.setp(ax2.get_yticklabels(), visible=False)
        
    plt.savefig('kfiltplots/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), dpi = 300)
    # plt.show()


####CREATE KAPPA INTENSITY OVERLAY####

    
    from numpy.ma import masked_array
    
    
    
    
    
    
    
    flat_kapfilt1, flat_kapfilt2 = np.ravel(kap_filt1), np.ravel(kap_filt2)
    
    flat_inten1filt, flat_inten2filt = np.where(flat_inten1 <= 20, 21, flat_inten1), np.where(flat_inten2 <= 20, 21, flat_inten2) 
    
    kapoverlay1 = np.resize(np.where(flat_kapfilt1 > 20, flat_inten1filt, flat_kapfilt1), fit_res1.fit['params'][:,:,-1].shape)
    kapoverlay2 = np.resize(np.where(flat_kapfilt2 > 20, flat_inten2filt, flat_kapfilt2), fit_res2.fit['params'][:,:,-1].shape)
    ###Save Numpy Arrays###
    np.save('datakappaoverlay/1/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), kapoverlay1)
    np.save('datakappaoverlay/2/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), kapoverlay2)
    
    intenmask1 = masked_array(kapoverlay1, kapoverlay1 <= 20)
    kapmask1 = masked_array(kapoverlay1, kapoverlay1 > 20)
    
    intenmask2 = masked_array(kapoverlay2, kapoverlay2 <= 20)
    kapmask2 = masked_array(kapoverlay2, kapoverlay2 > 20)
    
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=[10,8]) 
    plt.suptitle(gauss_cutout1.meta['mod_index']['date_obs'][0:10] + '\n' +gauss_cutout1.meta['mod_index']['date_obs'][10:])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05, hspace=None)
        
    ####FE16 KAPPA####

    ax1.set_xlabel('Solar-X [arcsec]')
    ax1.set_ylabel('Solar-Y [arcsec]')
    
    ax1.set_title('FeXVI Intensity $\kappa$ Overlay Map')
    intenmap1 = ax1.imshow(intenmask1, norm=colors.LogNorm(vmax  = 250000) , origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap='PuBu_r')
    kapmap1 = ax1.imshow(kapmask1, origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'RdYlGn', vmax = 20, vmin = 0.0)

    ####FE23 KAPPA#####
    
    ax2.set_xlabel('Solar-X [arcsec]')
    # ax2.set_ylabel('Solar-Y [arcsec]')
    
    ax2.set_title('FeXXIII Intensity $\kappa$ Overlay Map')
    intenmap2 = ax2.imshow(intenmask2, norm=colors.LogNorm(vmax  = 250000) , origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap='PuBu_r')
    kapmap2 = ax2.imshow(kapmask2, origin='lower', extent = raster_cutout2.meta['extent_arcsec'], cmap = 'RdYlGn', vmax = 20, vmin = 0.0)

    
    cbar = fig.colorbar(kapmap1, ax=(ax1,ax2), orientation= 'horizontal', fraction=0.146, pad=0.1,label = "$\kappa$")
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.savefig('koverlayplots/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), dpi = 300)
    # plt.show()
    
    
    ###############################################################
    ###############################################################
    ###############################################################
    #####Figure window for intensity kappa overlay, vel, width#####


    
    #####Figure window for intensity kappa overlay, vel, width#####

    fig = plt.figure(figsize=[15,10])
    plot_grid = fig.add_gridspec(nrows=2, ncols=3, wspace=0.25, hspace = 0.3)
    plt.suptitle(gauss_cutout1.meta['mod_index']['date_obs'][0:10] + '\n' +gauss_cutout1.meta['mod_index']['date_obs'][10:])


    #Subplots for Fe16

    data_subplt = fig.add_subplot(plot_grid[0,0])
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    intenmap1 = data_subplt.imshow(intenmask1, norm=colors.LogNorm(vmax  = 250000) , origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap='PuBu_r')
    kapmap1 = data_subplt.imshow(kapmask1, origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'RdYlGn', vmax = 20, vmin = 0.0)
    #define position for subplt colourbar
    data_subplt.set_title('FeXVI Intensity $\kappa$ Overlay')
    pos1 = kapmap1

    fig.colorbar(pos1,fraction=0.046, pad=0.04, label = "$\kappa$")
    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[0,1])
    vel_map = fit_gauss1.get_map(component=0, measurement='vel')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(vel_map.data,origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap = 'RdBu', vmax = 50, vmin = -50)
    #define position for subplt colourbar
    pos2 = data_subplt.imshow(vel_map.data,origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap = 'RdBu', vmax = 50, vmin = -50)
    data_subplt.set_title('FeXVI Velocity')
    fig.colorbar(pos2, fraction=0.046, pad=0.04, label = 'km/s')
    plt.setp(data_subplt.get_yticklabels(), visible=False)
    data_subplt.set_xlabel('Solar-X [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[0,2])
    width_map = fit_gauss1.get_map(component=0, measurement='width')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(width_map.data,origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap = 'viridis', vmax = 0.05, vmin = 0.02)
    #define position for subplt colourbar
    pos3 = data_subplt.imshow(width_map.data,origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap = 'viridis', vmax = 0.05, vmin = 0.02)
    data_subplt.set_title('FeXVI Line-Width')
    fig.colorbar(pos3,fraction=0.046, pad=0.04, label = 'Wavelength [$\AA$]')
    plt.setp(data_subplt.get_yticklabels(), visible=False)
    data_subplt.set_xlabel('Solar-X [arcsec]')



    #Subplots for Fe23

    data_subplt = fig.add_subplot(plot_grid[1,0])
    data_subplt.set_title('FeXXIII Intensity $\kappa$ Overlay')
    intenmap2 = data_subplt.imshow(intenmask2, norm=colors.LogNorm(vmax  = 250000) , origin='lower', extent = gauss_cutout1.meta['extent_arcsec'], cmap='PuBu_r')
    kapmap2 = data_subplt.imshow(kapmask2, origin='lower', extent = raster_cutout1.meta['extent_arcsec'], cmap = 'RdYlGn', vmax = 20, vmin = 0.0)
    pos4 = kapmap2

    fig.colorbar(pos4,fraction=0.046, pad=0.04, label = "$\kappa$")

    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[1,1])
    vel_map = fit_gauss2.get_map(component=0, measurement='vel')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(vel_map.data,origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap = 'RdBu', vmax = 50, vmin = -50)
    #define position for subplt colourbar
    pos5 = data_subplt.imshow(vel_map.data,origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap = 'RdBu', vmax = 50, vmin = -50)
    data_subplt.set_title('FeXXIII Velocity')
    fig.colorbar(pos5, fraction=0.046, pad=0.04, label = 'km/s')
    plt.setp(data_subplt.get_yticklabels(), visible=False)
    data_subplt.set_xlabel('Solar-X [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[1,2])
    width_map = fit_gauss2.get_map(component=0, measurement='width')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(width_map.data,origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap = 'viridis', vmax = 0.08, vmin = 0.02)
    #define position for subplt colourbar
    pos6 = data_subplt.imshow(width_map.data,origin='lower', extent = gauss_cutout2.meta['extent_arcsec'], cmap = 'viridis', vmax = 0.08, vmin = 0.02)
    data_subplt.set_title('FeXXIII Line-Width')
    fig.colorbar(pos6,fraction=0.046, pad=0.04, label = 'Wavelength [$\AA$]')
    plt.setp(data_subplt.get_yticklabels(), visible=False)

    data_subplt.set_xlabel('Solar-X [arcsec]')
    plt.savefig('koverlay_vel_width/'+gauss_cutout1.meta['mod_index']['date_obs'][11:19].replace(':',''), dpi = 300)

    plt.show()

