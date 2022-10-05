
# import libraries 
# ------------------------------------------------------------------------
# A library to use OpenSees via Python
import EzGM

from time import time
# A Library to visualize data from Python
import matplotlib.pyplot as plt
# A library provides high-performance vector, matrix and higher-dimensional data structures for Python
import numpy as np
import openseespy.opensees as ops
# A visulization package for OpenSees in Python
#import openseespy.postprocessing.Get_Rendering as opsplt  # another visulization package for opensees
## A visulization package for OpenSees in Python
#import openseespy.postprocessing.ops_vis as opsv

# A command required to print the figures on notebook
#%matplotlib inline 
class Unit():
        # Define units
        # ------------------------------------------------------------------------
        # Basic Units
        m = 1.0
        kN = 1.0
        sec = 1.0

        # Length
        mm = m / 1000.0
        cm = m / 100.0
        inch = 25.4 * mm
        ft = 12.0 * inch

        # Force
        N = kN / 1000.0
        kips = kN * 4.448221615
        lb = kips / 1.0e3

        # Stress (kN/m2 or kPa)
        Pa = N / (m ** 2)
        kPa = Pa * 1.0e3
        MPa = Pa * 1.0e6
        GPa = Pa * 1.0e9
        ksi = 6.8947573 * MPa
        psi = 1e-3 * ksi

        # Mass - Weight
        tonne = kN * sec ** 2 / m
        kg = N * sec ** 2 / m
        lb = psi*inch**2

        # Gravitational acceleration
        g = 9.81*m/sec**2

        # Time
        min = 60*sec
        hr = 60*min

class GM:
    def __init__(self):
        pass
    def define_GM_pattern(self,filepath="Records"):
            with open(f'{filepath}//GMR_names.txt') as file:    # read the record names
                gm_names = [line.rstrip() for line in file]
            print(f"gm_names length = {len(gm_names),}{gm_names}")
            dts = np.loadtxt(f'{filepath}//GMR_dts.txt')        # load the time steps for records
            gm_idx = 0                                      # index for record being applied
            tsTag = 1                                       # tag for time series to use
            pTag = 1                                        # tag for load pattern to use
            for item in dts:
                A_g = -np.loadtxt(f'{filepath}//'+gm_names[gm_idx]) # load the record file as an array
                print(f"{gm_names[gm_idx]} depreminin zaman serisi oluşturuluyor...")
                dt = dts[gm_idx]                                # time step of record
                values = list(-1 * A_g)                         # should be negative
                ops.wipe()
                ops.timeSeries('Path', tsTag, '-dt', dt, '-values', *A_g, '-factor', Unit.g) # time series object
                ops.pattern('UniformExcitation', pTag, 1, '-accel', tsTag)              # pattern object
                print(f"{gm_names[gm_idx]} depreminin zaman serisi opensees domaininde tanımlandı...")
                gm_idx += 1
                tsTag += 1
                pTag += 1

    def tbdy_load_records(self,file='Records',database_name='NGA_W2',SD1=1.073, SDS=2.333, PGA=0.913, 
    nGM=11, selection=1, Tp=1, Mw_lim=[5.5, 8], Vs30_lim=[180, 360], Rjb_lim=[0, 20], fault_lim=2, opt=1, maxScale=2, weights=[1, 1]):
                        
            
                """
                file : Txt formatındaki deprem kayıtlarının bulunduğu klasörün dosya yolu
                Rule 1: Mean of selected records should remain above the lower bound target spectra.
                    For selection = 1: Sa_rec = (Sa_1 or Sa_2) - lower bound = 1.0 * SaTarget(0.2Tp-1.5Tp)
                    For Selection = 2: Sa_rec = (Sa_1**2+Sa_2**2)**0.5 - lower bound = 1.3 * SaTarget(0.2Tp-1.5Tp)

                Rule 2:
                    No more than 3 records can be selected from the same event! In other words, rec_eqID cannot be the same for 
                    more than 3 of the selected records.

                Rule 3:
                    At least 11 records (or pairs) must be selected.


                    Parameters
            |      ----------

                SD1 : float, optional, the default is 1.073.
                    Short period design spectral acceleration coefficient.
                SDS : float, optional, the default is 2.333.
                    Design spectral acceleration coefficient for a period of 1.0 seconds.
                PGA : float, optional, the default is 0.913.
                    Peak ground acceleration.
                nGM : int, optional, the default is 11.
                    Number of records to be selected.
                selection : int, optional, the default is 1.
                    Number of ground motion components to select.
                Tp : float, optional, the default is 1.
                    Predominant period of the structure.
                Mw_lim : list, optional, the default is None.
                    The limiting values on magnitude.
                Vs30_lim : list, optional, the default is None.
                    The limiting values on Vs30.
                Rjb_lim : list, optional, the default is None.
                    The limiting values on Rjb.
                fault_lim : int, optional, the default is None.
                    The limiting fault mechanism. 0 for unspecified fault 1 for strike-slip fault 2 for normal fault 3 for reverse fault
                opt : int, optional, the default is 1.
                    If equal to 0, the record set is selected using method of “least squares”.
                     If equal to 1, the record set selected such that scaling factor is closer to 1.
                      If equal to 2, the record set selected such that both scaling factor and standard deviation is lowered.
                maxScale : float, optional, the default is 2.
                    Maximum allowed scaling factor, used with opt=2 case.
                weights = list, optional, the default is [1,1].
                    Error weights (mean,std), used with opt=2 case.

                """
                
                
                startTime = time()
                # 1.) Initialize the tbdy_2018 object for record selection
                spec = EzGM.Selection.tbdy_2018(database=database_name, outdir=file)

                # 2.) Select the ground motions
                """
                spec.select(Lat=41.0582, Long=29.00951, DD=2, Soil='ZC', nGM=11, selection=1, Tp=1,
                            Mw_lim=[6.5, 8], Vs30_lim=[200, 700], Rjb_lim=[0, 20], fault_lim=None, opt=0, 
                            maxScale=2)

                # selected records can be plotted at this stage
                spec.plot(save=1, show=1)
                
                spec.select(SD1, SDS, PGA, nGM, selection, Tp, Mw_lim, Vs30_lim, Rjb_lim, fault_lim, opt, maxScale, weights)

                spec.plot(save=1, show=1)
                """

                spec.select(SD1, SDS, PGA, nGM, selection, Tp, Mw_lim, Vs30_lim,
                 Rjb_lim, fault_lim, opt, maxScale, weights)

                spec.plot(save=1, show=1)
                # 3.) If database == 'NGA_W2' you can first download the records via nga_download method
                # from NGA-West2 Database [http://ngawest2.berkeley.edu/] and then use write method

                spec.ngaw2_download(username = 'muhammedsural@gmail.com', pwd = 'sural6177')
                # 4.) If you have records already inside recs_f\database.zip\database or
                # downloaded records for database = NGA_W2 case, write whatever you want,
                # the object itself, selected and scaled time histories
                spec.write(obj=1, recs=1, recs_f='')
                EzGM.utility.RunTime(startTime)

                # Let's also save design response spectrum
                import numpy as np
                Periods = np.array([spec.T]).T
                Sa = np.array([spec.target]).T
                SaT = np.concatenate((Periods,Sa),axis=1)
                np.savetxt(f'{file}//SaT.txt', SaT, fmt = '%.5f')

                