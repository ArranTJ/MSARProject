# -*- coding: utf-8 -*-
"""
Created on Mon May 20 10:44:39 2024

@author: Arran
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.cbook as cbook
import matplotlib.colors as colors
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import wcs_to_celestial_frame
import eispac

from astropy.visualization import astropy_mpl_style
from astropy.io import fits

plt.style.use(astropy_mpl_style)

if __name__ == '__main__':
    # Read in the fit template and EIS observation
    data_filepath = './eis_20170906_120233.data.h5'
    
    #testing template fitting for both lines
    
    tmplt = eispac.EISFitTemplate(value=[ 19253.7287, 262.9824, 0.0341, 118.9812, 1601.5393,263.7539,0.0649, 342.9053], line_ids=['Fe XVI 262.984', 'Fe XXIII 263.760' ], wmin=262.78399658203125, wmax= 263.96002197265625)

    ##original tmplt fit
    # template_filepath = './fe_16_262_984.1c.template.h5'
    # tmplt = eispac.read_template(template_filepath)
    ##
    data_cube = eispac.read_cube(data_filepath, tmplt.central_wave)
    
    
    
    # Select a cutout of the raster
    eis_frame = wcs_to_celestial_frame(data_cube.wcs)
    lower_left = [None, SkyCoord(Tx=500, Ty=-260, unit=u.arcsec, frame=eis_frame)]
    upper_right = [None, SkyCoord(Tx=600, Ty=-150, unit=u.arcsec, frame=eis_frame)]
    raster_cutout = data_cube.crop(lower_left, upper_right)
    
    # # Fit the data and save it to disk
    fit_res = eispac.fit_spectra(raster_cutout, tmplt, ncpu='max')
    save_filepaths = eispac.save_fit(fit_res, save_dir='cwd')

    # # Find indices and world coordinates of max intensity
    sum_data_inten = raster_cutout.sum_spectra().data
    sum_data_fullregion = data_cube.sum_spectra().data  #use for full region image on subplot?
    iy, ix = np.unravel_index(sum_data_inten.argmax(), sum_data_inten.shape)
    ex_world_coords = raster_cutout.wcs.array_index_to_world(iy, ix, 0)[1]
    y_arcsec, x_arcsec = ex_world_coords.Ty.value, ex_world_coords.Tx.value
    
    print(y_arcsec)
    print(x_arcsec)
    print(ex_world_coords)
    # Extract data profile and interpolate fit at higher spectral resolution
    data_x = raster_cutout.wavelength[iy, ix, :]
    data_y = raster_cutout.data[iy, ix, :]
    data_err = raster_cutout.uncertainty.array[iy, ix, :]
    fit_x, fit_y = fit_res.get_fit_profile(coords=[iy,ix], num_wavelengths=100)
    c0_x, c0_y = fit_res.get_fit_profile(0, coords=[iy,ix], num_wavelengths=100)
    c1_x, c1_y = fit_res.get_fit_profile(1, coords=[iy,ix], num_wavelengths=100)
    # c2_x, c2_y = fit_res.get_fit_profile(2, coords=[iy,ix], num_wavelengths=100)
    
    
    
    
    
    
    # Make a multi-panel figure with the cutout and example profile
    
    
    fig = plt.figure(figsize=[10,5])
    plot_grid = fig.add_gridspec(nrows=2, ncols=3, wspace=0.3)
    
    # data_subplt = fig.add_subplot(plot_grid[0,0])
    # data_subplt.imshow(sum_data_inten, origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt = fig.add_subplot(plot_grid[0,0])
    data_subplt.scatter(x_arcsec, y_arcsec, color='r', marker='x')
    # data_subplt.imshow(sum_data_inten, origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt.imshow(sum_data_fullregion, norm=colors.PowerNorm(gamma=0.2) ,origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt.set_title('Region Data\n'+' ' +raster_cutout.meta['mod_index']['date_obs'])
    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[0,1])
    data_subplt.scatter(x_arcsec, y_arcsec, color='r', marker='x')
    # data_subplt.imshow(sum_data_inten, origin='lower', extent = data_cube.meta['extent_arcsec'])
    data_subplt.imshow(sum_data_inten, norm=colors.PowerNorm(gamma=0.2) , origin='lower', extent = raster_cutout.meta['extent_arcsec'])
    data_subplt.set_title('Selected Cutout\n')
    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
        
    
    profile_subplt = fig.add_subplot(plot_grid[0,2])
    # profile_subplt.errorbar(data_x, data_y, yerr=data_err, ls='', marker='o', color='k')
    profile_subplt.plot(data_x, data_y, ls='', marker='o', color='k')
    profile_subplt.plot(fit_x, fit_y, color='b', label='Combined profile')
    profile_subplt.plot(c0_x, c0_y, color='r', label=fit_res.fit['line_ids'][0])
    profile_subplt.plot(c1_x, c1_y, color='r', ls='--', label=fit_res.fit['line_ids'][1])
    # profile_subplt.plot(c2_x, c2_y, color='g', label='Background')
    profile_subplt.set_title(f'Cutout indices: iy = {iy}, ix = {ix}')
    profile_subplt.set_xlabel('Wavelength [$\AA$]')
    profile_subplt.set_ylabel('Intensity ['+raster_cutout.unit.to_string()+']')
    profile_subplt.legend(loc='upper left', frameon=False)
    
    data_subplt = fig.add_subplot(plot_grid[1,0])
    inten_map = fit_res.get_map(component=0, measurement='intensity')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(inten_map.data, norm=colors.PowerNorm(gamma=0.5) , origin='lower', extent = raster_cutout.meta['extent_arcsec'], cmap='PuBu_r')
    #define position for subplt colourbar
    pos1 = data_subplt.imshow(inten_map.data, norm=colors.PowerNorm(gamma=0.5) , origin='lower', extent = raster_cutout.meta['extent_arcsec'], cmap='PuBu_r')
    data_subplt.set_title('Selected Cutout\n [SunPy Intensity]')
    fig.colorbar(pos1, label = 'Intensity [erg / (s sr cm2)]')

    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[1,1])
    vel_map = fit_res.get_map(component=0, measurement='vel')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(vel_map.data,origin='lower', extent = raster_cutout.meta['extent_arcsec'], cmap = 'RdBu')
    #define position for subplt colourbar
    pos2 = data_subplt.imshow(vel_map.data,origin='lower', extent = raster_cutout.meta['extent_arcsec'], cmap = 'RdBu')
    data_subplt.set_title('Selected Cutout\n [SunPy Velocity]')
    fig.colorbar(pos2, label = 'km/s')

    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    data_subplt = fig.add_subplot(plot_grid[1,2])
    width_map = fit_res.get_map(component=0, measurement='width')
    #sunpy measurement map can be used in matplotlib subplt using .imshow by including .data at the end of our maps (e.g inten_map.data)
    data_subplt.imshow(width_map.data,origin='lower', extent = raster_cutout.meta['extent_arcsec'], cmap = 'viridis')
    #define position for subplt colourbar
    pos3 = data_subplt.imshow(width_map.data,origin='lower', extent = raster_cutout.meta['extent_arcsec'], cmap = 'viridis')
    data_subplt.set_title('Selected Cutout\n [SunPy Width]')
    fig.colorbar(pos3, label = 'Angstrom')

    data_subplt.set_xlabel('Solar-X [arcsec]')
    data_subplt.set_ylabel('Solar-Y [arcsec]')
    
    
    # inten_map = fit_res.get_map(component=0, measurement='intensity')
    # inten_map.peek()
    
    # vel_map = fit_res.get_map(0, 'vel')
    # vel_map.peek()
    
    linewidth_map = fit_res.get_map(0, 'width')
    linewidth_map.peek()
    
    plt.show()
    

    