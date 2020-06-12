# =============================================================================
import numpy as np
import h5py
import datetime
import os
import PQFit as FT
# =============================================================================
class Data:
    """this is a container for the data collected in multiple or single run and information regarding that measurement
        Attributes
        ----------
        sample: str
            The name given to the sample during the data acquisition (default is empty string)
        type: str
            If the measurement is a qens spectra or a fws (default is empty string)
        instrument: str
            String characterizing the instrument
        run_num: array
            the run numbers of the files (default is empty string)
        temperature: float
            this is the temperature recorded during the acquisition (default is Not a Number)
        time: date value
           The time at which the measurement started (default the current date and time, meaning the date and time at which the initialisation happens )
        end_time: date value
           The time at which the measurement ended (default the current date and time, meaning the date and time at which the initialisation happens )
        time_delta: float
            this is the time in second since the time at which the set temperature was changed (default is Not A Number)
        time_initial: date value
            The time at which the set temperature was changed (default the current date and time, meaning the date and time at which the initialisation happens )
        temperature_initial: float
            The temperature at which the set temperature was changed (default is Not a Number)
        temperature_final: float
            The temperature of the set temperature(default is Not a Number)
        x:?? depending on what analysis is performed the value of x and y can change, it is safer to not use it
        y:??depending on what analysis is performed the value of x and y can change, it is safer to not use it
        sqw: array
            The dynamic structure factor (default is Not a Number)
        dsqw: array
            The error of dynamic structure factor (default is Not a Number)
        q: array
            the value of the scattering vector at which the data were collected
        hw: the energy transfer for this measurement

        Fit:list of fits

        Methods
        -------
        -load data
        -clean_nonsense
        -get_values_qens
        """

    # initialising the class. written so you can also put values in it if you have already loaded the data in another variable by writing:

    # Measurement(sample:'Ig with whatever',type:'fws', temperature:300)
    # with this an instance of the class Measurement is created and has default values for every attribute except  sample,type and temperature
#%% initialize the class
    def __init__(self,
                 ProtEntry=None,
                 sqw=None,
                 dsqw=None,
                 q=None,
                 hw=None):
        # initialising all the variables
        self.sqw = [] if sqw is None else sqw
        self.dsqw = [] if dsqw is None else dsqw
        self.q = [] if q is None else q
        self.hw = [] if hw is None else hw
        self.T=[]
        self.FWS_q=[]
        self.FWS_hw=[]
        self.FWS_sqw=[]
        self.FWS_dsqw=[]
        self.FWS_T=[]
        self.ProtEntry=ProtEntry
        self.Fits=[]
#%% functions to load data
    def load_mantid(self,filename):
        '''
        This function loads data for different instruments from reduced mantid files

        Returns
        -------
        None.

        '''
        if os.path.exists(self.ProtEntry['filepath_analysed'] + filename + 'QENS'):
            self.hw, self.q, self.sqw,self.dsqw,self.T=self.load_MANTID_QENS(self.ProtEntry['filepath_analysed'] + filename + 'QENS')
        if os.path.exists(self.ProtEntry['filepath_analysed'] + filename + 'FWS'):
            self.FWS_hw, self.FWS_q, self.FWS_sqw,self.FWS_dsqw,self.FWS_T=self.load_MANTID_FWS(self.ProtEntry['filepath_analysed'] + filename + 'FWS')
        
    def load_ascii(self,filename):
        '''
        This function reads data from ascii which are organized as follows
        """""""""""""""""""""""
        #T
        value
        #q
        values
        #hw
        values
        #sqw
        values as matrice
        #dsqw
        errors as matrice
        Returns
        -------
        None.

        '''
        totalFileName=self.ProtEntry['filepath_analysed'] +filename
        T=[]
        q=[]
        hw=[]
        sqw=[]
        dsqw=[]
        self.hw=hw
        self.sqw=sqw
        self.dsqw=dsqw
        self.q=q
        self.T=T
        #%% load mantid functions
    def load_MANTID_FWS(self, filename ):   # read Mantid-reduced IN16B FWS data
        f      = h5py.File( filename, 'r' )           # open hdf file for reading
        k      = 1                                    # stored Mantid workspace index (one for each energy)
        h5path = '/mantid_workspace_'+str(k)+'/'      # internal hdf node name for that workspace
        sqw=[]
        dsqw=[]
        T=[]
        q=[]
        hw=[]
        while h5path+'workspace/axis1/' in f.keys():  # loop over all energy offsets
            #x= np.transpose( np.array( f[ h5path+'workspace/axis1/' ]  ) )
            y= np.transpose( np.array( f[ h5path+'workspace/axis2/' ]  ) )
            sqw.append([np.transpose( np.array( f[ h5path+'workspace/values/' ] ) )])
            dsqw.append([np.transpose( np.array( f[ h5path+'workspace/errors/' ] ) )])
            T.append(np.transpose( np.array( f[ h5path+ 'logs/sample.temperature/value/'])))
            q.append([ 4 * np.pi / 6.271 * np.sin( y * np.pi / 360.0 )])
            hw.append( np.array( f[ h5path+'logs/Doppler.maximum_delta_energy/value' ] ))
            k += 1
            h5path   = '/mantid_workspace_'+str(k)+'/'
            del y
        return sqw,dsqw,T,q,hw
    # read Mantid-reduced IN16B QENS data
    def load_MANTID_QENS(self, filename ):
        f = h5py.File( filename, 'r' )
        x =np.transpose( np.array( f['/mantid_workspace_1/workspace/axis1/']  )[0:-1]*1000),
        y=np.transpose( np.array( f['/mantid_workspace_1/workspace/axis2/']  ) ),
        sqw=np.transpose( np.array( f['/mantid_workspace_1/workspace/values/'] ) ),
        dsqw=np.transpose( np.array( f['/mantid_workspace_1/workspace/errors/'] ) ),
        T=np.transpose( np.array( f['/mantid_workspace_1/logs/sample.temperature/value/'] ) )
        q=4 * np.pi / 6.271 * np.sin( np.asarray(y) * np.pi / 360.0 )
        return np.squeeze(x), np.squeeze(q), np.squeeze(sqw), np.squeeze(dsqw),T
    #%% other functions
    def clean_nonsense(self):
        '''
        This function is setting all errors of non-reasonable values to infinity and the corresponding counts to zero.

        Returns
        -------
        None.

        '''
        self.sqw[np.isnan(self.sqw)] = 0
        self.sqw[np.isinf(self.sqw)] = 0
        self.sqw[np.isnan(self.dsqw)] = 0
        self.dsqw[np.isnan(self.sqw)] = np.inf
        self.dsqw[np.isinf(self.sqw)] = np.inf
        self.sqw[self.sqw < 5e-7 * self.sqw.max()] = 0  # reasonable signal-to-noise range
        self.dsqw[np.isnan(self.dsqw)] = np.inf
        self.dsqw[self.sqw == 0] = np.inf  # no signal at all should have infinite error
        self.dsqw[self.dsqw == 0] = np.inf  # "a measurement without error is nonsense"

    def Fit(self,fitname,Vanadium):
        self.Fits.append(FT.Fit(fitname))
        self.clean_nonsense()
        self.Fits[-1].fit_QENS(self,fitname,Vanadium[0])
    
    def Subtract(self,subdat,scaling):
        if len(scaling)==1:
            scaling=scaling*np.ones(len(self.q))
        if not len(scaling)==len(self.q):
            raise NameError('Provided scaling dimensions for subtraction is not correct!')
        for hiq,_ in enumerate(self.q):
            self.sqw[:,hiq]=self.sqw[:,hiq]-scaling[hiq]*subdat[0].sqw[:,hiq]
            self.dsqw[:,hiq]=(self.dsqw[:,hiq]**2+scaling[hiq]*subdat[0].dsqw[:,hiq]**2)**0.5