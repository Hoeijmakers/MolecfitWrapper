# MolecfitWrapper
A python wrapper for running molecfit on a time-series of high-resolution spectra.


This repo contains the python source code to call Molecfit on a time-series of high-resolution spectra, and the ability to examine the result using a GUI interface. Call as `list_of_wls,list_of_t_spectra = do_molecfit(headers,spectra,mode='HARPS')`, where `headers` is a list containing astropy header objects of the 1D spectra contained in the list `spectra`.

The `mode` keyword currently only accepts `HARPS` and `HARPSN`.

The output of `do_molecfit` is a list of telluric transmission spectra corresponding to each spectrum in the time-series, and a list of the corresponding wavelength axes.

# Warning:
This code is experimental and therefore not well documented or tested; and you may need to tweak it more to get it to run on your system. For example, you will need to alter the `molecfit_input_folder` and `molecfit_prog_folder` variables on the first few lines of `do_molecfit`, which are currently pointing to the correct folders on my system.

