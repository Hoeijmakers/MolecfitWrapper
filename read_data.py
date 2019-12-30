#################################
#ESPRESSO
#################################

def read_ESPRESSO_S2D(inpath,outname,air=True,molecfit=True):
    """
    reads in the ESPRESSO files and prepares them for use in molecfit
    this functions then calls do_mmolecfit from the molecfit module
    input:
        inpath: type: string, path to the s2d ESPRESSO files
        outpath: type: string, path to where the telluric correction should be saved
        air: type:boolean, is the wavelength in air or vacuum
        molecfit: type: boolean, function can be run with or without starting molecfit
    note for Jens: The wavelength is much easier in ESPRESSO compared to HARPS, so I kicked out
    all the parts needed for that in the HARPS function
    Since Romain has for some reason only given me the S2D files and not the S1D files (God knows why),
    ans also only fiber A (again, God knows why), this function only does molecfit on the sodium orders
    """
    import os
    import pdb
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import utils as ut
    import molecfit as mol
    import pyfits
    import copy
    import scipy.interpolate as interp
    import pickle

    #First check the input:
    ut.typetest('inpath in read_ESPRESSO_S2D ',inpath,str)
    ut.typetest('outname in read_ESPRESSO_S2D ',outname,str)
    ut.typetest('air in read_ESPRESSO_S2D ',air,bool)
    if os.path.exists(inpath) != True:
        print("ERROR in read_ESPRESSO_S2D: Data input path (%s) does not exist." % inpath)
        sys.exit()

    filelist=os.listdir(inpath)
    N=len(filelist)

    if len(filelist) == 0:
        print("ERROR in read_ESPRESSO_S2D: input folder (%s) is empty." % inpath)
        sys.exit()

    #The following variables define the lists in which all the necessary data will be stored.
    framename=[]
    header=[]
    type=[]
    texp=np.array([])
    date=[]
    npx=np.array([])
    norders=np.array([])
    mjd=np.array([])
    ccfmjd=np.array([])
    s1dmjd=np.array([])
    s2d=[]
    airmass=np.array([])
    berv=np.array([])
    wave=[]

    outpath = ut.path(outname)
    if os.path.exists(outpath) != True:
        os.makedirs(outpath)

    #ccftotal = 0 #This will hold the sum of the CCFs
    s2d_count = 0
    sci_count = 0



    for i in range(N):
        if filelist[i].endswith('S2D_A.fits'):
            s2d_count += 1
            print(filelist[i])
            #data,hdr=fits.getdata(inpath+filelist[i],header=True)

            hdul = fits.open(inpath+filelist[i])
            data = copy.deepcopy(hdul[1].data)
            hdr1 = hdul[0].header
            hdr2 = hdul[1].header
            wavedata=copy.deepcopy(hdul[5].data)
            wave.append(wavedata)
            hdul.close()
            del hdul[0].data
            if hdr2['EXTNAME'] == 'SCIDATA':
                framename.append(filelist[i])
                header.append(hdr1)
                type.append(hdr2['EXTNAME'])
                texp=np.append(texp,hdr1['EXPTIME'])
                date.append(hdr1['DATE-OBS'])
                mjd=np.append(mjd,hdr1['MJD-OBS'])
                npx=np.append(npx,hdr2['NAXIS1'])
                norders=np.append(norders,hdr2['NAXIS2'])
                s2d.append(data)
                s2d_count += 1

                sci_count += 1
                berv=np.append(berv,hdr1['HIERARCH ESO QC BERV'])
                airmass=np.append(airmass,0.5*(hdr1['HIERARCH ESO TEL3 AIRM START']+hdr1['HIERARCH ESO TEL3 AIRM END']))



         

       
    #Now we catch some errors:
    #-The above should have read a certain number of s2d files.
    #-A certain number of these should be SCIENCE frames.
    #-There should be at least one WAVE file.
    #-All exposures should have the same number of spectral orders.
    #-All orders should have the same number of pixels (this is true for ESPRESSO).
    #-The wave frame should have the same dimensions as the order frames.
    #-If nowave is set, test that all frames used the same wave_A calibrator.
    #-The blaze file needs to have the same shape as the s2d files.
    #-The number of s1d files should be the same as the number of s2d files.



    # if s2d_count != s1d_count:
#    #     print('ERROR in read_ESPRESSO_s2d: The numbers of 1ds and s2d files are different.')
#    #     print("These are the files and their types:")
#    #     for i in range(len(type)):
#    #         print('   '+framename[i]+'  %s' % type[i])
#    #     sys.exit()
#    if s2d_count == 0:
#        print("ERROR in read_ESPRESSO_s2d: The input folder (%s) does not contain files ending in s2d.fits." % inpath)
#        sys.exit()
#    if sci_count == 0:
#        print("ERROR in read_ESPRESSO_s2d: The input folder (%2) contains s2d files, but none of them are classified as SCIENCE frames with the HIERARCH ESO DPR CATG keyword.")
#        print("These are the files and their types:")
#        for i in range(len(type)):
#            print('   '+framename[i]+'  %s' % type[i])
#        sys.exit()
#    if np.max(np.abs(norders-norders[0])) == 0:
#        norders=int(norders[0])
#    else:
#        print("ERROR in read_ESPRESSO_s2d: Not all files have the same number of orders.")
#        print("These are the files and their number of orders:")
#        for i in range(len(type)):
#            print('   '+framename[i]+'  %s' % norders[i])
#        sys.exit()
#    if np.max(np.abs(npx-npx[0])) == 0:
#        npx=int(npx[0])
#    else:
#        print("ERROR IN read_ESPRESSO_s2d: Not all files have the same number of pixels.")
#        print("These are the files and their number of pixels:")
#        for i in range(len(type)):
#            print('   '+framename[i]+'  %s' % npx[i])
#        sys.exit()
#    if np.max(np.abs(nrv-nrv[0])) == 0:
#        nrv=int(nrv[0])
#    else:
#        print("ERROR IN read_ESPRESSO_s2d: Not all files have the same number of pixels.")
#        print("These are the files and their number of pixels:")
#        for i in range(len(type)):
#            print('   '+framename[i]+'  %s' % npx[i])
#        sys.exit()

  
   
#
#    if len(s1dhdr) != len(s2d) and molecfit == True:
#        print('ERROR in read_ESPRESSO_s2d: The number of s1d SCIENCE files and s2d SCIENCE files is not the same. (%s vs %s)' % (len(s1dhdr),len(s2d)))
#        print('Switching off the molecfit option will suppress this error.')



    #Ok, so now we should have ended up with a number of lists that contain all
    #the relevant information of our science frames.
    #We determine how to sort the resulting lists in time:
#    sorting = np.argsort(mjd)
    s2dsorting = np.argsort(mjd)


    #First sort the s1d files for application of molecfit.
    if molecfit == True:
        s2dhdr_sorted=[]
        s2d_sorted=[]
        wave_sorted=[]
        for i in s2dsorting:
            s2dhdr_sorted.append(header[i])
            s2d_sorted.append(s2d[i])
            wave_sorted.append(wave[i])

        # print('Molecfit will be executed onto the files in this order:')
        # for x in s1dhdr_sorted:
        #     print(x['DATE-OBS'])
        list_of_wls,list_of_trans = mol.do_molecfit(s2dhdr_sorted,s2d_sorted,mode='ESPRESSO',load_previous=False,order=116,wave=wave_sorted)
        mol.write_telluric_transmission_to_file(list_of_wls,list_of_trans,outpath+'telluric_transmission_spectra.pkl')




#    ccftotal = 0.0
    #Now we loop over all exposures and collect the i-th order from each exposure,
    #put these into a new matrix and save them to FITS images:
    f=open(outpath+'obs_times','w',newline='\n')
    headerline = 'MJD'+'\t'+'DATE'+'\t'+'EXPTIME'+'\t'+'MEAN AIRMASS'+'\t'+'BERV (km/s)'+'\t'+'FILE NAME'
    
    for i in range(int(norders[0])):
        order = np.zeros((sci_count,int(npx[0])))
        
        wave_axis = wave[0][i]/10.0#Convert to nm.
#        ccforder = np.zeros((ccf_count,nrv))
        print('CONSTRUCTING ORDER %s' % i)
        c = 0#To count the number of science frames that have passed. The counter
        # c is not equal to j because the list of files contains not only SCIENCE
        # frames.
#        cc = 0#Same for ccfs
#        for j in range(len(ccfsorting)):
#            ccf=ccfs[ccfsorting[j]]
#            ccforder[cc,:] = ccf[i,:]
#            cc+=1
        for j in range(len(s2dsorting)):#Loop over exposures
            if i ==0:
                print('---'+type[s2dsorting[j]]+'  '+date[s2dsorting[j]])
           
            exposure = s2d[s2dsorting[j]]

            
            order[c,:] = exposure[i,:]
            #T_i = interp.interp1d(list_of_wls[j],list_of_trans[j])#This should be time-sorted, just as the s2d files.
            #Do a manual check here that the MJDs are identical.
            #Also, determiine what to do with airtovac.
            #tel_order[c,:] = T_i[wave_axis]
            #Now I also need to write it to file.
            if i ==0:#Only do it the first time, not for every order.
                line = str(mjd[s2dsorting[j]])+'\t'+date[s2dsorting[j]]+'\t'+str(texp[s2dsorting[j]])+'\t'+str(airmass[s2dsorting[j]])+'\t'+str(berv[s2dsorting[j]])+'\t'+framename[s2dsorting[j]]+'\n'
                f.write(line)
            c+=1
#        ccftotal+=ccforder



#        fits.writeto(outpath+'ccf_'+str(i)+'.fits',ccforder,overwrite=True)
        fits.writeto(outpath+'order_'+str(i)+'.fits',order,overwrite=True)
        fits.writeto(outpath+'wave_'+str(i)+'.fits',wave_axis,overwrite=True)
#    fits.writeto(outpath+'ccftotal.fits',ccftotal,overwrite=True)
    f.close()
    print('Time-table written to '+outpath+'obs_times')


    #ADD OUTPUT OF THE TELLURIC TRANSMISSION SPECTRUM HERE
    #Also output the list of trans spectra to the datafolder while checking that the
    #dates are the same as in the just exported config file.






#######################################
#HARPS
#######################################

def read_HARPS_e2ds(inpath,outname,air=True,nowave=False,molecfit=True):
    """THIS IS A PYTHON TRANSLATION OF READ_DATA (BELOW). IT SHOULD NOT WORK
    WITH A PATHS FILE, JUST A FOLDER THAT CONTAINS ONLY FITS FILES AND THEN
    IT WORKS FROM THE KEYWORDS TO DO EVERYTHING AUTOMATICALLY.
   if __name__ == "__main__":
    read_HARPS_e2ds('/Users/julias/Documents/E2DS_pipeline/data/stars/WASP-127/2017-02-28','/Users/julias/Documents/E2DS_pipeline/data/stars/WASP-127/tell_corr/night1',True,False,False)
    WRITE GOOD TESTS AND DOCUMENTATION.
    #inpath
    /Users/julias/Documents/E2DS_pipeline/data/stars/WASP-127/2017-02-28
    #outname
    /Users/julias/Documents/E2DS_pipeline/data/stars/WASP-127/tell_corr/night1
    Set the nowave keyword to True if the dataset has no wave files associated with it.
    This may happen if you downloaded ESO Advanced Data Products, which include
    reduced science e2ds's but not reduced wave e2ds's. The wavelength solution
    is still encoded in the fits header however, so we take it from there, instead."""
    import os
    import pdb
    from astropy.io import fits
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import utils as ut
    import molecfit as mol
    import pyfits
    import copy
    import scipy.interpolate as interp
    import pickle

    #First check the input:
    ut.typetest('inpath in read_HARPS_e2ds ',inpath,str)
    ut.typetest('outname in read_HARPS_e2ds ',outname,str)
    ut.typetest('air in read_HARPS_e2ds ',air,bool)
    if os.path.exists(inpath) != True:
        print("ERROR in read_HARPS_e2ds: Data input path (%s) does not exist." % inpath)
        sys.exit()

    filelist=os.listdir(inpath)
    N=len(filelist)

    if len(filelist) == 0:
        print("ERROR in read_HARPS_e2ds: input folder (%s) is empty." % inpath)
        sys.exit()

    #The following variables define the lists in which all the necessary data will be stored.
    framename=[]
    header=[]
    s1dhdr=[]
    type=[]
    texp=np.array([])
    date=[]
    mjd=np.array([])
    ccfmjd=np.array([])
    s1dmjd=np.array([])
    npx=np.array([])
    nrv=np.array([])
    norders=np.array([])
    e2ds=[]
    s1d=[]
    wave1d=[]
    airmass=np.array([])
    berv=np.array([])
    wave=[]
    blaze=[]
    ccfs=[]
    wavefile_used = []
    outpath = ut.path(outname)
    if os.path.exists(outpath) != True:
        os.makedirs(outpath)

    #ccftotal = 0 #This will hold the sum of the CCFs
    e2ds_count = 0
    sci_count = 0
    wave_count = 0
    ccf_count = 0
    blaze_count = 0
    s1d_count = 0

    for i in range(N):
        if filelist[i].endswith('e2ds_A.fits'):
            e2ds_count += 1
            print(filelist[i])
            #data,hdr=fits.getdata(inpath+filelist[i],header=True)

            hdul = fits.open(inpath+filelist[i])
            data = copy.deepcopy(hdul[0].data)
            hdr = hdul[0].header
            hdul.close()
            del hdul[0].data
            if hdr['HIERARCH ESO DPR CATG'] == 'SCIENCE':
                framename.append(filelist[i])
                header.append(hdr)
                type.append(hdr['HIERARCH ESO DPR CATG'])
                texp=np.append(texp,hdr['EXPTIME'])
                date.append(hdr['DATE-OBS'])
                mjd=np.append(mjd,hdr['MJD-OBS'])
                npx=np.append(npx,hdr['NAXIS1'])
                norders=np.append(norders,hdr['NAXIS2'])
                e2ds.append(data)

                sci_count += 1
                berv=np.append(berv,hdr['HIERARCH ESO DRS BERV'])
                airmass=np.append(airmass,0.5*(hdr['HIERARCH ESO TEL AIRM START']+hdr['HIERARCH ESO TEL AIRM END']))

                if nowave == True:
                    #Record which wavefile was used by the pipeline to
                        #create the wavelength solution.
                    wavefile_used.append(hdr['HIERARCH ESO DRS CAL TH FILE'])

                    wavedata=read_wave_from_HARPS_header(hdr)
                    wave.append(wavedata)
            # else:
                # berv=np.append(berv,np.nan)
                # airmass=np.append(airmass,np.nan)
        if filelist[i].endswith('wave_A.fits'):
            print(filelist[i]+' (wave)')
            if nowave == True:
                print("WARNING in read_HARPS_e2ds: nowave was set to True but a wave_A file")
                print("was detected. This wave file is now ignored in favor of the header.")
            else:
                wavedata=fits.getdata(inpath+filelist[i])
                wave.append(wavedata)
                wave_count += 1
        if filelist[i].endswith('ccf_G2_A.fits'):
            #ccf,hdr=fits.getdata(inpath+filelist[i],header=True)
            hdul = fits.open(inpath+filelist[i])
            ccf = copy.deepcopy(hdul[0].data)
            hdr = hdul[0].header
            hdul.close()
            del hdul[0].data

            if hdr['HIERARCH ESO DPR CATG'] == 'SCIENCE':
                #ccftotal+=ccf
                ccfs.append(ccf)
                ccfmjd=np.append(ccfmjd,hdr['MJD-OBS'])
                nrv=np.append(nrv,hdr['NAXIS1'])
                ccf_count += 1

        if filelist[i].endswith('blaze_A.fits'):
                print(filelist[i]+' (blaze)')
                blazedata=fits.getdata(inpath+filelist[i])
                blaze.append(blazedata)
                blaze_count += 1
        if filelist[i].endswith('s1d_A.fits'):
            hdul = fits.open(inpath+filelist[i])
            data_1d = copy.deepcopy(hdul[0].data)
            hdr = hdul[0].header
            hdul.close()
            del hdul[0].data
            if hdr['HIERARCH ESO DPR CATG'] == 'SCIENCE':
                s1d.append(data_1d)
                s1dhdr.append(hdr)
                s1dmjd=np.append(s1dmjd,hdr['MJD-OBS'])
                s1d_count += 1
    #Now we catch some errors:
    #-The above should have read a certain number of e2ds files.
    #-A certain number of these should be SCIENCE frames.
    #-There should be at least one WAVE file.
    #-All exposures should have the same number of spectral orders.
    #-All orders should have the same number of pixels (this is true for HARPS).
    #-The wave frame should have the same dimensions as the order frames.
    #-If nowave is set, test that all frames used the same wave_A calibrator.
    #-The blaze file needs to have the same shape as the e2ds files.
    #-The number of s1d files should be the same as the number of e2ds files.


    if ccf_count != sci_count:
        print("ERROR in read_HARPS_e2ds: There is a different number of science CCFs as there is science frames.")
        sys.exit()
    # if e2ds_count != s1d_count:
    #     print('ERROR in read_HARPS_e2ds: The numbers of 1ds and e2ds files are different.')
    #     print("These are the files and their types:")
    #     for i in range(len(type)):
    #         print('   '+framename[i]+'  %s' % type[i])
    #     sys.exit()
    if e2ds_count == 0:
        print("ERROR in read_HARPS_e2ds: The input folder (%s) does not contain files ending in e2ds.fits." % inpath)
        sys.exit()
    if sci_count == 0:
        print("ERROR in read_HARPS_e2ds: The input folder (%2) contains e2ds files, but none of them are classified as SCIENCE frames with the HIERARCH ESO DPR CATG keyword.")
        print("These are the files and their types:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % type[i])
        sys.exit()
    if np.max(np.abs(norders-norders[0])) == 0:
        norders=int(norders[0])
    else:
        print("ERROR in read_HARPS_e2ds: Not all files have the same number of orders.")
        print("These are the files and their number of orders:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % norders[i])
        sys.exit()
    if np.max(np.abs(npx-npx[0])) == 0:
        npx=int(npx[0])
    else:
        print("ERROR IN read_HARPS_e2ds: Not all files have the same number of pixels.")
        print("These are the files and their number of pixels:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % npx[i])
        sys.exit()
    if np.max(np.abs(nrv-nrv[0])) == 0:
        nrv=int(nrv[0])
    else:
        print("ERROR IN read_HARPS_e2ds: Not all files have the same number of pixels.")
        print("These are the files and their number of pixels:")
        for i in range(len(type)):
            print('   '+framename[i]+'  %s' % npx[i])
        sys.exit()
    if wave_count >= 1:
        wave=wave[0]#SELECT ONLY THE FIRST WAVE FRAME. The rest is ignored.
    else:
        if nowave == False:
            print("ERROR in read_HARPS_e2ds: No wave_A.fits file was detected.")
            print("These are the files in the folder:")
            for i in range(N):
                print(filelist[i])
            print("This may have happened if you downloaded the HARPS data from the")
            print("ADP query form, which doesn't include wave_A files (as far as I")
            print("have seen). Set the /nowave keyword in your call to read_HARPS_e2ds")
            print("if you indeed do not expect a wave_A file to be present.")
    if nowave == True:
        if all(x == wavefile_used[0] for x in wavefile_used):
            print("Nowave is set, and simple wavelength calibration extraction")
            print("works, as all files in the dataset used the same wave_A file.")
            wave=wave[0]
        else:
            print("ERROR IN read_HARPS_e2ds: Nowave is set, but not all files")
            print("in the dataset used the same wave_A file when the pipeline was")
            print("run. Catching this requres an interpolation step that is currently")
            print("not yet implemented. Exiting. These are the filenames and their")
            print("wave_A file used:")
#            for i in range(N):
#                print('   '+framename[i]+'  %s' % wavefile_used[0])
            wave=wave[0]
            print("I ALLOW YOU TO CONTINUE BUT USING ONLY THE FIRST WAVELENGTH")
            print("SOLUTION. A PART OF THE DATA MAY BE AFFECTED BY HAVING ASSUMED")
            print("THE WRONG SOLUTION. If you are doing transits, you don't need")
            print("this kind of precision.")
    if blaze_count >= 1:
        blaze=blaze[0]#SELECT ONLY THE FIRST WAVE FRAME. The rest is ignored.
    if np.shape(wave) != np.shape(e2ds[0]):
        print("ERROR in read_HARPS_e2ds: A wave file was detected but its shape (%s,%s) does not match that of the orders (%s,%s)" % (np.shape(wave)[0],np.shape(wave)[1],np.shape(e2ds[0])[0],np.shape(e2ds[0])[1]))
    if np.shape(blaze) != np.shape(e2ds[0]) and blaze_count > 0:
        print("ERROR in read_HARPS_e2ds: A blaze file was detected but its shape (%s,%s) does not match that of the orders (%s,%s)" % (np.shape(blaze)[0],np.shape(wave)[1],np.shape(e2ds[0])[0],np.shape(e2ds[0])[1]))
    if len(s1dhdr) != len(e2ds) and molecfit == True:
        print('ERROR in read_HARPS_e2ds: The number of s1d SCIENCE files and e2ds SCIENCE files is not the same. (%s vs %s)' % (len(s1dhdr),len(e2ds)))
        print('Switching off the molecfit option will suppress this error.')



    #Ok, so now we should have ended up with a number of lists that contain all
    #the relevant information of our science frames.
    #We determine how to sort the resulting lists in time:
    sorting = np.argsort(mjd)
    ccfsorting = np.argsort(ccfmjd)
    s1dsorting = np.argsort(s1dmjd)


    #First sort the s1d files for application of molecfit.
    if molecfit == True:
        s1dhdr_sorted=[]
        s1d_sorted=[]
        for i in range(len(s1dsorting)):
            s1dhdr_sorted.append(s1dhdr[s1dsorting[i]])
            s1d_sorted.append(s1d[s1dsorting[i]])

        # print('Molecfit will be executed onto the files in this order:')
        # for x in s1dhdr_sorted:
        #     print(x['DATE-OBS'])
        list_of_wls,list_of_trans = mol.do_molecfit(s1dhdr_sorted,s1d_sorted,load_previous=False)
        mol.write_telluric_transmission_to_file(list_of_wls,list_of_trans,outpath+'telluric_transmission_spectra.pkl')




    ccftotal = 0.0
    #Now we loop over all exposures and collect the i-th order from each exposure,
    #put these into a new matrix and save them to FITS images:
    f=open(outpath+'obs_times','w',newline='\n')
    headerline = 'MJD'+'\t'+'DATE'+'\t'+'EXPTIME'+'\t'+'MEAN AIRMASS'+'\t'+'BERV (km/s)'+'\t'+'FILE NAME'
    for i in range(norders):
        order = np.zeros((sci_count,npx))
        ccforder = np.zeros((ccf_count,nrv))
        wave_axis = wave[i,:]/10.0#Convert to nm.
        print('CONSTRUCTING ORDER %s' % i)
        c = 0#To count the number of science frames that have passed. The counter
        # c is not equal to j because the list of files contains not only SCIENCE
        # frames.
        cc = 0#Same for ccfs
        for j in range(len(ccfsorting)):
            ccf=ccfs[ccfsorting[j]]
            ccforder[cc,:] = ccf[i,:]
            cc+=1
        for j in range(len(sorting)):#Loop over exposures
            if i ==0:
                print('---'+type[sorting[j]]+'  '+date[sorting[j]])
            if type[sorting[j]] == 'SCIENCE':#This check may be redundant.
                exposure = e2ds[sorting[j]]
                order[c,:] = exposure[i,:]
                #T_i = interp.interp1d(list_of_wls[j],list_of_trans[j])#This should be time-sorted, just as the e2ds files.
                #Do a manual check here that the MJDs are identical.
                #Also, determiine what to do with airtovac.
                #tel_order[c,:] = T_i[wave_axis]
                #Now I also need to write it to file.
                if i ==0:#Only do it the first time, not for every order.
                    line = str(mjd[sorting[j]])+'\t'+date[sorting[j]]+'\t'+str(texp[sorting[j]])+'\t'+str(airmass[sorting[j]])+'\t'+str(berv[sorting[j]])+'\t'+framename[sorting[j]]+'\n'
                    f.write(line)
                c+=1
        ccftotal+=ccforder



        fits.writeto(outpath+'ccf_'+str(i)+'.fits',ccforder,overwrite=True)
        fits.writeto(outpath+'order_'+str(i)+'.fits',order,overwrite=True)
        fits.writeto(outpath+'wave_'+str(i)+'.fits',wave_axis,overwrite=True)
    fits.writeto(outpath+'ccftotal.fits',ccftotal,overwrite=True)
    f.close()
    print('Time-table written to '+outpath+'obs_times')


    #ADD OUTPUT OF THE TELLURIC TRANSMISSION SPECTRUM HERE
    #Also output the list of trans spectra to the datafolder while checking that the
    #dates are the same as in the just exported config file.



def read_wave_from_HARPS_header(h):
    """This reads the wavelength solution from the HARPS header keywords that
    encode the coefficients as a 4-th order polynomial."""
    import numpy as np
    import utils as ut
    npx = h['NAXIS1']
    no = h['NAXIS2']
    x = ut.findgen(npx) #np.linspace(0,npx-1,npx)
    wave=np.zeros((npx,no))

    key_counter = 0
    for i in range(no):
        l = x*0.0
        for j in range(4):
            l += h['ESO DRS CAL TH COEFF LL%s' %key_counter]*x**j
            key_counter +=1
        wave[:,i] = l
    wave = wave.T
    return(wave)


if __name__ == "__main__":
#    read_HARPS_e2ds('/Users/julias/Documents/E2DS_pipeline/data/stars/WASP-166/2017-03-14/', '/Users/julias/Documents/E2DS_pipeline/data/stars/WASP-166/tell_corr/night3/',True,True,False)
    
    read_ESPRESSO_S2D( '/Users/julias/Documents/E2DS_pipeline/data/stars/WASP-76/ESPRESSO/2018-09-03/',  '/Users/julias/Documents/E2DS_pipeline/data/stars/WASP-76/ESPRESSO/tell_corr/night1/', True,True)
 

