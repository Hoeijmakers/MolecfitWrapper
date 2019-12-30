def do_molecfit(headers,spectra,mode='HARPS',load_previous=False,order=-1,wave=0.):
    """This is a function that pipes a list of s1d spectra into molecfit, and
    executes it. It first launches the molecfit gui on the middle spectrum of the
    sequence, and then loops through the entire list, returning the transmission
    spectra of the Earths atmosphere in the same order as the list provided.
    These can then be used to correct the s1d spectra or the e2ds spectra.
    Note that the s1d spectra are assumed to be in the barycentric frame in vaccuum,
    but that the output transmission spectrum is in the observers frame, and e2ds files
    are in air wavelengths by default.

    If you have run do_molecfit before, and want to reuse the output of the previous run
    for whatever reason, set the load_previous keyword to True. This will reload the
    list of transmission spectra created last time, if available.
    """


    import pdb
    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    import os.path
    import utils as ut
    import pickle
    import copy
    molecfit_input_folder='/Users/julias/Documents/E2DS_pipeline/MOLECfit/MOLECfit_new/share/molecfit/'
    molecfit_prog_folder='/Users/julias/Documents/E2DS_pipeline/MOLECfit/MOLECfit_new/bin/'
    temp_specname = copy.deepcopy(mode)#The name of the temporary file used (without extension).
    #The spectrum will be named like this.fits There should be a this.par file as well,
    #that contains a line pointing molecfit to this.fits:
    parname=temp_specname+'.par'





    #====== ||  START OF PROGRAM   ||======#
    N = len(headers)
    if N != len(spectra):
        print('ERROR in prep_for_molecfit: Length of list of headers is not equal to length of list of spectra (%s , %s)' % (N,len(spectra)))
        sys.exit()

    #Test that the input root and molecfit roots exist; that the molecfit root contains the molecfit executables.
    #that the input root contains the desired parfile and later fitsfile.
    molecfit_input_root=ut.path(molecfit_input_folder)
    molecfit_prog_root=ut.path(molecfit_prog_folder)
    if os.path.isdir(molecfit_input_root) != True:
        print('ERROR in prep_for_molecfit: '+molecfit_input_root+' does not exist!')
        sys.exit()
    if os.path.isdir(molecfit_prog_root) != True:
        print('ERROR in prep_for_molecfit: '+molecfit_prog_root+' does not exist!')
        sys.exit()
    if os.path.isfile(molecfit_input_root+temp_specname+'.par') != True:
        print('ERROR in prep_for_molecfit: '+molecfit_input_root+temp_specname+'.par does not exist!')
        sys.exit()
    if os.path.isfile(molecfit_prog_root+'molecfit') != True:
        print('ERROR in prep_for_molecfit: '+molecfit_prog_root+'molecfit does not exist!')
        sys.exit()
    if os.path.isfile(molecfit_prog_root+'molecfit_gui') != True:
        print('ERROR in do_molecfit: '+molecfit_prog_root+'molecfit_gui does not exist!')
        sys.exit()

    pickle_outpath = molecfit_input_root+'previous_run_of_do_molecfit.pkl'


    if load_previous == True:
        if os.path.isfile(pickle_outpath) ==  False:
            print('WARNING in do_molecfit: Previously saved run is not available.')
            print('The user will have to re-fit.')
            print('That run will then be saved.')
            load_previous = False
        else:
            pickle_in = open(pickle_outpath,"rb")
            list_of_wls,list_of_fxc,list_of_trans = pickle.load(pickle_in)

    if load_previous == False:
        list_of_wls = []
        list_of_fxc = []
        list_of_trans = []

        if mode=='HARPS':
            middle_i = int(round(0.5*N))#We initialize molecfit on the middle spectrum of the time series.
            write_file_to_molecfit(molecfit_input_root,temp_specname+'.fits',headers,spectra,middle_i)
        elif mode=='ESPRESSO':
            middle_i = int(round(0.5*N))#We initialize molecfit on the middle spectrum of the time series.
            write_file_to_molecfit_ESPRESSO(molecfit_input_root,temp_specname+'.fits',headers,spectra,wave,middle_i,order)
        else:
            print("Mode does not exist in do_molecfit.")
            sys.exit()
        execute_molecfit(molecfit_prog_root,molecfit_input_root+parname,gui=True)
        wl,fx,trans = retrieve_output_molecfit(molecfit_input_root+temp_specname)
        remove_output_molecfit(molecfit_input_root,temp_specname)

        for i in range(N):#range(len(spectra)):
            print('Fitting spectrum %s from %s' % (i+1,len(spectra)))
            t1=ut.start()
            write_file_to_molecfit(molecfit_input_root,temp_specname+'.fits',headers,spectra,i)
            execute_molecfit(molecfit_prog_root,molecfit_input_root+parname,gui=False)
            wl,fx,trans = retrieve_output_molecfit(molecfit_input_root+temp_specname)
            remove_output_molecfit(molecfit_input_root,temp_specname)
            list_of_wls.append(wl*1000.0)#Convert to nm.
            list_of_fxc.append(fx/trans)
            list_of_trans.append(trans)
            ut.end(t1)

        pickle_outpath = molecfit_input_root+'previous_run_of_do_molecfit.pkl'
        with open(pickle_outpath, 'wb') as f: pickle.dump((list_of_wls,list_of_fxc,list_of_trans),f)

    to_do_manually = check_fit_gui(list_of_wls,list_of_fxc,list_of_trans)
    if len(to_do_manually) > 0:
        print('The following spectra were selected to redo manually:')
        print(to_do_manually)
        #CHECK THAT THIS FUNCIONALITY WORKS:
        for i in to_do_manually:
            i = int(i)
            write_file_to_molecfit(molecfit_input_root,temp_specname+'.fits',headers,spectra,i)
            execute_molecfit(molecfit_prog_root,molecfit_input_root+parname,gui=True)
            wl,fx,trans = retrieve_output_molecfit(molecfit_input_root+temp_specname)
            list_of_wls[i] = wl*1000.0#Convert to nm.
            list_of_fxc[i] = fx
            list_of_trans[i] = trans
    return(list_of_wls,list_of_trans)


class molecfit_gui(object):
    """This class defines most of the behaviour of the molecfit GUI."""
    def __init__(self,wls,fxc,trans):
        import matplotlib.pyplot as plt
        import math
        self.wls=wls
        self.fxc=fxc
        self.trans=trans
        self.i = 0#The current spectrum.
        self.N = len(wls)#total number of spectra.
        self.fig,self.ax = plt.subplots(2,1,sharex=True,figsize=(14,6))
        self.maxboxes = 20.0
        self.nrows = math.ceil(self.N/self.maxboxes)#number of rows
        plt.subplots_adjust(left=0.05)#Create space for the interface (#rhyme comments ftw <3)
        plt.subplots_adjust(right=0.75)
        plt.subplots_adjust(bottom=0.1+0.05*self.nrows)
        plt.subplots_adjust(top=0.95)
        self.set_spectrum(self.i)
        self.ax[0].set_title('Corrected / uncorrected s1d #%s/%s' % (self.i,self.N-1))#start counting at 0.
        self.ax[1].set_title('Transmission spectrum')
        self.img1a=self.ax[0].plot(wls[self.i],fxc[self.i]*trans[self.i],alpha=0.5)
        self.img1b=self.ax[0].plot(wls[self.i],fxc[self.i],alpha=0.5)
        self.img2 =self.ax[1].plot(wls[self.i],trans[self.i])
        self.ax[0].set_ylim(0,self.img_max)
        self.selected=[]#list of selected spectra. Is empty when we start.
        self.crosses=[]#List of checkboxes that have crosses in them. Starts empty, too.

    def set_spectrum(self,i):
        #This modifies the currently active spectrum to be plotted.
        import numpy as np
        import utils as ut
        import matplotlib.pyplot as plt

        ut.typetest('i in molecfit_gui/set_order',i,int)
        self.wl = self.wls[i]
        self.spectrum = self.fxc[i]
        self.Tspectrum= self.trans[i]
        self.img_max = np.nanmean(self.spectrum[ut.selmax(self.spectrum,0.02,s=0.02)])*1.3

    def update_plots(self):
        #This redraws the plot planels, taking care to reselect axis ranges and such.
        import matplotlib.pyplot as plt
        import pdb
        import numpy as np
        import copy
        self.img1a[0].set_xdata(self.wl)
        self.img1a[0].set_ydata(self.spectrum*self.Tspectrum)
        self.img1b[0].set_xdata(self.wl)
        self.img1b[0].set_ydata(self.spectrum)
        self.img2[0].set_xdata(self.wl)
        self.img2[0].set_ydata(self.Tspectrum)
        self.ax[0].set_title('Corrected / uncorrected s1d #%s/%s' % (self.i,self.N-1))
        self.ax[0].set_ylim(0,self.img_max)
        self.fig.canvas.draw_idle()

    def slide_spectrum(self,event):
        #This handles the spectrum-selection slide bar.
        self.i = int(self.spectrum_slider.val)
        self.set_spectrum(self.i)
        self.update_plots()

    def previous(self,event):
        #The button to go to the previous spectrum.
        self.i -= 1
        if self.i <0:#If order tries to become less than zero, loop to the highest spectrum.
            self.i = self.N-1
        self.set_spectrum(self.i)
        self.spectrum_slider.set_val(self.i)#Update the slider value.
        self.update_plots()#Redraw everything.

    def next(self,event):
        #The button to go to the next spectrum. Similar to previous().
        self.i += 1
        if self.i > self.N-1:
            self.i = 0
        self.set_spectrum(self.i)
        self.spectrum_slider.set_val(self.i)
        self.update_plots()

    def cancel(self,event):
        #This is a way to crash hard out of the python interpreter.
        import sys
        print('------Canceled by user')
        sys.exit()

    def save(self,event):
        #The "save" button is actually only closing the plot. Actual saving of things
        #happens after plt.show() below.
        import matplotlib.pyplot as plt
        print('------Closing GUI')
        plt.close(self.fig)


    def draw_crosses(self):
        import math
        import pdb
        for c in self.crosses:
            self.selec.lines.remove(c)
        self.crosses=[]#Delete the references and start again.
        self.fig.canvas.draw_idle()
        for s in self.selected:
            cx = s % (self.maxboxes)
            cy = self.nrows-math.floor(s/self.maxboxes)-0.5
            x = [float(cx)-0.5,float(cx)+0.5]
            y1 = [float(cy)-0.5,float(cy)+0.5]
            y2 = y1[::-1]
            self.crosses.append(self.selec.plot(x,y1,color='red')[0])
            self.crosses.append(self.selec.plot(x,y2,color='red')[0])
        self.fig.canvas.draw_idle()





def check_fit_gui(wls,fxc,trans):
    """This code initializes the GUI that plots the telluric-corrected spectra
    from Molecfit. The user may select spectra to be re-fit manually via the Molecfit GUI. Note that
    since molecfit takes between a few and 10 minutes to run on a single spectrum,
    this becomes arduous when more than a few spectra are selected in this way.
    It quickly becomes worthwile to redo the entire sequence with different inclusion
    regions overnight. The code returns the list of spectra that need to be done
    manually via the Molecfit GUI.

    Input: The list of wl axes, each of which was returned by a call to molecfit;
    and similarly the corrected spectra fxc and the transmission spectra.

    """


    import sys
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button, RadioButtons, CheckButtons
    import utils as ut
    import numpy as np

    M = molecfit_gui(wls,fxc,trans)

    #The slider to cycle through orders:
    rax_slider = plt.axes([0.8, 0.2, 0.1, 0.02])
    rax_slider.set_title('Order')
    M.spectrum_slider = Slider(rax_slider,'', 0,M.N-1,valinit=0,valstep=1)#Store the slider in the model class
    M.spectrum_slider.on_changed(M.slide_spectrum)

    #The Previous order button:
    rax_prev = plt.axes([0.8, 0.1, 0.04, 0.05])
    bprev = Button(rax_prev, ' <<< ')
    bprev.on_clicked(M.previous)

    #The Next order button:
    rax_next = plt.axes([0.86, 0.1, 0.04, 0.05])
    bnext = Button(rax_next, ' >>> ')
    bnext.on_clicked(M.next)

    #The save button:
    rax_save = plt.axes([0.92, 0.1, 0.07, 0.05])
    bsave = Button(rax_save, 'Continue')
    bsave.on_clicked(M.save)

    #The cancel button:
    rax_cancel = plt.axes([0.92, 0.025, 0.07, 0.05])
    bcancel = Button(rax_cancel, 'Cancel')
    bcancel.on_clicked(M.cancel)

    #This is to rescale the x-size of the checkboxes so that they are squares.
    bbox = M.fig.get_window_extent().transformed(M.fig.dpi_scale_trans.inverted())
    width, height = bbox.width*M.fig.dpi, bbox.height*M.fig.dpi


    M.selec=plt.axes([0.05,0.03,0.7,0.05*M.nrows])
    M.selec.spines['bottom'].set_color('white')
    M.selec.spines['top'].set_color('white')
    M.selec.spines['left'].set_color('white')
    M.selec.spines['right'].set_color('white')
    vlines = ut.findgen(M.N-1)+0.5

    row = M.nrows
    offset = 0
    for i in range(M.N):
        #print(i,float(i)-offset)

        if float(i)-offset > M.maxboxes-1.0:
            row -= 1
            offset += M.maxboxes
        M.selec.plot(float(i)-offset+np.array([-0.5,-0.5,0.5,0.5,-0.5]),[row,row-1,row-1,row,row],color='black')
        M.selec.text(float(i)-offset,row-0.5,'%s' % i,color='black',horizontalalignment='center',verticalalignment='center')



    M.selec.set_xlim(-0.55,M.maxboxes-1.0+0.55)#A little margin to make sure that the line thickness is included.
    M.selec.set_ylim(-0.05,1.0*M.nrows+0.05)
    #M.selec.set_yticklabels([])
    M.selec.xaxis.set_tick_params(labelsize=8)
    M.selec.yaxis.set_tick_params(labelsize=8)



    def select_spectrum_box(event):

            #This handles with a mouseclick in either of the three plots while in add mode.
        if event.inaxes in [M.selec]:#Check that it occurs in one of the subplots.
            cc = event.xdata*1.0#xdata is the column that is selected.
            cr = event.ydata*1.0
            spectrum = np.round(cc)+np.round((M.nrows-cr-0.5))*M.maxboxes
            if spectrum < M.N:
                if spectrum in M.selected:
                    M.selected.remove(spectrum)
                    print('---Removed spectrum %s from manual' % spectrum)
                else:
                    M.selected.append(spectrum)
                    print('---Added spectrum %s to manual' % spectrum)
            M.draw_crosses()
    M.click_connector = M.fig.canvas.mpl_connect('button_press_event',select_spectrum_box)#This is the connector that registers clicks

    plt.show()
    print('Closed GUI, returning.')
    return(M.selected)




def remove_output_molecfit(path,name):
    import os
    try:
        os.remove(path+name+'_out.fits')
    except:
        pass
    try:
        os.remove(path+name+'_out_gui.fits')
    except:
        pass
    try:
        os.remove(path+name+'_out_molecfit_out.txt')
    except:
        pass
    try:
        os.remove(path+name+'_out_fit.par')
    except:
        pass
    try:
        os.remove(path+name+'_out_fit.atm')
    except:
        pass
    try:
        os.remove(path+name+'_out_apply_out.txt')
    except:
        pass
    try:
        os.remove(path+name+'_out_fit.atm.fits')
    except:
        pass
    try:
        os.remove(path+name+'_out_fit.res')
    except:
        pass
    try:
        os.remove(path+name+'_out_tac.asc')
    except:
        pass
    try:
        os.remove(path+name+'_out_tac.fits')
    except:
        pass
    try:
        os.remove(path+name+'_TAC.fits')
    except:
        pass
    try:
        os.remove(path+name+'_out_fit.res.fits')
    except:
        pass
    # os.remove(path+name+'_out_fit_res.fits')
    try:
        os.remove(path+name+'_out_fit.asc')
    except:
        pass
    try:
        os.remove(path+name+'_out_fit.fits')
    except:
        pass
    return


def retrieve_output_molecfit(path):
    import astropy.io.fits as fits
    import os.path
    import sys
    file = path+'_out_tac.fits'
    if os.path.isfile(file) != True:
        print('ERROR in retrieve_output_molecfit: '+file+' does not exist!')
        sys.exit()


    with fits.open(path+'_out_tac.fits') as hdul:
        wl=hdul[1].data['lambda']
        fx=hdul[1].data['flux']
        trans=hdul[1].data['mtrans']
    return(wl,fx,trans)

def execute_molecfit(molecfit_prog_root,molecfit_input_file,gui=False):
    import os
    if gui == False:
        command = molecfit_prog_root+'molecfit '+molecfit_input_file
        command2 = molecfit_prog_root+'calctrans '+molecfit_input_file
        os.system(command)
        os.system(command2)
    if gui == True:
        command = 'python3 '+molecfit_prog_root+'molecfit_gui '+molecfit_input_file
        os.system(command)
    #python3 /Users/hoeijmakers/Molecfit/bin/molecfit_gui /Users/hoeijmakers/Molecfit/share/molecfit/spectra/cross_cor/test.par

def write_file_to_molecfit(molecfit_file_root,name,headers,spectra,ii):
    """This is a wrapper for writing a spectrum from a list to molecfit format.
    name is the filename of the fits file that is the output.
    headers is the list of astropy header objects associated with the list of spectra
    in the spectra variable. ii is the number from that list that needs to be written.
    """
    import astropy.io.fits as fits
    from scipy import stats
    import copy
    import utils as ut
    import numpy as np
    import scipy.constants as const
    ii = int(ii)
    spectrum = spectra[ii]
    npx = len(spectrum)
    berv = headers[ii]['HIERARCH ESO DRS BERV']*1000.0#Need to un-correct the s1d spectra to go back to the frame of the Earths atmosphere.
    wave = (headers[ii]['CDELT1']*ut.findgen(len(spectra[ii]))+headers[ii]['CRVAL1'])*(1.0-berv/const.c)
    #at the end, when the transmission spectrum is corrected, we stay in the barycentric frame because these will be used to
    #correct the e2ds spectra which are not yet berv-corrected.
    err = np.sqrt(spectrum)

    #Write out the s1d spectrum in a format that molecfit eats.
    #This is a fits file with an empty primary extension that contains the header of the original s1d file.
    #Plus an extension that contains a binary table with 3 columns.
    #The names of these columns need to be indicated in the molecfit parameter file,
    #as well as the name of the file itself. This is currently hardcoded.
    col1 = fits.Column(name = 'wavelength', format = '1D', array = wave)
    col2 = fits.Column(name = 'flux', format       = '1D', array = spectrum)
    col3 = fits.Column(name = 'err_flux', format   = '1D', array = err)
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = fits.Header()
    prihdr = copy.deepcopy(headers[ii])
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(molecfit_file_root+name,overwrite=True)
    print('Spectrum %s written' % ii)
    return(0)
    
def write_file_to_molecfit_ESPRESSO(molecfit_file_root,name,headers,spectra, waves, ii,order):
    """This is a wrapper for writing a spectrum from a list to molecfit format.
    name is the filename of the fits file that is the output.
    headers is the list of astropy header objects associated with the list of spectra
    in the spectra variable. ii is the number from that list that needs to be written.
    """
    import astropy.io.fits as fits
    from scipy import stats
    import copy
    import utils as ut
    import numpy as np
    import scipy.constants as const
    ii = int(ii)
    spectrum = spectra[ii][order]
    npx = len(spectrum)
    berv = headers[ii]['HIERARCH ESO QC BERV']*1000.0#Need to un-correct the s1d spectra to go back to the frame of the Earths atmosphere.
    wave = waves[ii][order]*(1.0-berv/const.c)
    
    #at the end, when the transmission spectrum is corrected, we stay in the barycentric frame because these will be used to
    #correct the e2ds spectra which are not yet berv-corrected.
    err = np.sqrt(spectrum)

    #Write out the s1d spectrum in a format that molecfit eats.
    #This is a fits file with an empty primary extension that contains the header of the original s1d file.
    #Plus an extension that contains a binary table with 3 columns.
    #The names of these columns need to be indicated in the molecfit parameter file,
    #as well as the name of the file itself. This is currently hardcoded.
    col1 = fits.Column(name = 'wavelength', format = '1D', array = wave)
    col2 = fits.Column(name = 'flux', format       = '1D', array = spectrum)
    col3 = fits.Column(name = 'err_flux', format   = '1D', array = err)
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = fits.Header()
    prihdr = copy.deepcopy(headers[ii])
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(molecfit_file_root+name,overwrite=True)
    print('Spectrum %s written' % ii)
    return(0)




#The following are not strictly molecfit-specific
def write_telluric_transmission_to_file(wls,T,outpath):
    """This saves a list of wl arrays and a corresponding list of transmission-spectra
    to a pickle file, to be read by the function below."""
    import pickle
    print('------Saving telluric transmission to '+outpath)
    with open(outpath, 'wb') as f: pickle.dump((wls,T),f)

def read_telluric_transmission_from_file(inpath):
    import pickle
    print('------Reading telluric transmission from '+inpath)
    pickle_in = open(inpath,"rb")
    return(pickle.load(pickle_in))#This is a tuple that can be unpacked into 2 elements.







def apply_telluric_correction(inpath,list_of_wls,list_of_orders):
    import scipy.interpolate as interp
    import sys
    import pdb
    import utils as ut
    """This applies a telluric correction (done by molecfit) that was saved in the
    telluric"""
    wlT,fxT = read_telluric_transmission_from_file(inpath)

    No = len(list_of_wls)
    x = ut.findgen(No)

    if No != len(list_of_orders):
        print('ERROR in telluric correction: List of data wls and List of orders do not have the same length.')
        sys.exit()

    Nexp = len(wlT)

    if Nexp != len(fxT):
        print('ERROR in telluric correction: List of telluric wls and telluric spectra read from file do not have the same length.')
        sys.exit()


    list_of_orders_cor = []
    # ut.save_stack('test.fits',list_of_orders)
    # pdb.set_trace()
    for i in range(No):#Do the correction order by order.
        order = list_of_orders[i]
        order_cor = order*0.0
        wl = list_of_wls[i]

        for j in range(Nexp):
            #interpolates over the s1d file
            T_i = interp.interp1d(wlT[j],fxT[j],fill_value="extrapolate")
            #calculates the interpolation on the wwavelength grid of the order
            #then divides out the order
            order_cor[j]=order[j]/T_i(wl)
        list_of_orders_cor.append(order_cor)
        ut.statusbar(i,x)
    return(list_of_orders_cor)


