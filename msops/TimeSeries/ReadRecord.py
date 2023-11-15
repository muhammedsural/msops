import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def calc_timegap(Td,damp_ratio,T) -> float:
        """
        INFO
            "Required time gap between mainshock and aftershock for dynamic analysis of 
structures" makalesinde iki ardışık deprem arasında doğal hareketi sıfırlamak için gerekli olan zaman aralığının hesaplanmasındaki  önerilen formülasyon kullanılmıştır.
        INPUT
            Td          : Kuvvetli yer hareketi süresi
            damp_ratio  : Yapının sönüm oranı
            T           : Yapının doğal titreşim periyodu
        OUTPUT
            R_timegap = durgun geçmesi gereken zaman (sn)
        """
        R_rest = Td*(0.05/damp_ratio)*(((21.8559*T)+0.0258)*(Td**(-0.9982))+0.0214)
        R_timegap = round(R_rest,0)
        return R_timegap
    
def ReadRecord (filePath:str,gap:float,g=9.81,plot=1) -> pd.DataFrame:
    acceleration = []

    with open(filePath,"r") as file:
        for count,line in enumerate(file):
            if count >= 4:
                #newfile = filePath.split('/')[-1].rsplit(".VT2")[0].rsplit(".AT2")[0].rsplit(".DT2")[0]
                #dosya = open(newfile.join(".txt"),"w",encoding="utf-8")
                #for satir in line:
                #    dosya.write(line)
                #dosya.close()
                for i in line.strip().split():
                    acceleration.append(float(i)*g)
            if count == 3:
                npts = float(line.replace(",","").replace("SEC","").split()[1])
                dt = float(line.replace(",","").replace("SEC","").split()[3])
    time = np.arange(0,npts*dt,dt)
    
    if gap is not None:
        timegap = np.arange(time[-1],time[-1]+gap+dt,dt)
        accgap  = np.arange(0,len(timegap))*0
    
        time = np.append(time,timegap)
        acceleration = np.append(acceleration,accgap)
    
    if plot == 1:
        import matplotlib.pyplot as plt 
        fig, ax = plt.subplots(figsize=(20,10))
        fig.subplots_adjust(bottom=0.15, left=0.2)
        ax.grid()
        ax.plot(time,acceleration)
        ax.set_xlabel('Time [Sec]')
        ax.set_ylabel('Acceleration [cm/sn2]')
        #ax.axhline(0, color='black', lw=2)
        if g == 1:
            ax.set_ylabel('Acceleration [g]')
    TimeSeries = pd.DataFrame(columns=["Time","Acceleration"])
    
    TimeSeries["Time"] = time
    TimeSeries["Acceleration"] = acceleration

    return TimeSeries
    
def load_PEERNGA_record(filepath):

    '''
        Load record in .at2 format (PEER NGA Databases)

        Input:
            filepath : file path for the file to be load
            
        Returns:
        
            acc : vector wit the acceleration time series
            dt : time step
            npts : number of points in record
            eqname : string with year_name_station_component info

    '''

    import numpy as np

    with open(filepath) as fp:
        line = next(fp)
        line = next(fp).split(',')
        year = (line[1].split('/'))[2]
        eqname = (year + '_' + line[0].strip() + '_' + 
                  line[2].strip() + '_comp_' + line[3].strip())
        line = next(fp)
        line = next(fp).split(',')
        npts = int(line[0].split('=')[1])
        dt = float(line[1].split('=')[1].split()[0])
        acc = np.array([p for l in fp for p in l.split()]).astype(float)
    
    return acc,dt,npts,eqname

# Kullanılmaya fonksiyon
def main_after_record(eventname="Mammoth_Lakes",plot=False):
    """
    INPUT :
            eventname: Deprem kaydının ismi bu klasör içerisinde ilgili depremin farklı kayıtları bulunmakta
            plot     : Default False, Birleştirilmiş ivme kaydını çizer
    OUTPUT :
            mainafterTime : Zaman array i
            mainafterAcc  : İvme  array i
            dt            : Kayıtların zaman adımları
    """
    
    eventlist = [i for i in os.listdir(f"{eventname}") if i.endswith('.AT2')] #verilen deprem kayıtlarının bulunduğu klasördeki .AT2 uzantılı dosyaların listelenmesi
    print(eventlist)
    mainafterTime = [] #deprem kayıtlarına ait zaman serilerinin eklendiği liste 
    mainafterAcc = [] # deprem kayıtlarına ait ivme serilerinin eklendiği liste
    for count,event in enumerate(eventlist):
        with open(f"{eventname}//{event}") as fp:
            line = next(fp)
            line = next(fp).split(',')
            year = (line[1].split('/'))[2]
            eqname = (year + '_' + line[0].strip() + '_' + 
                    line[2].strip() + '_comp_' + line[3].strip())
            line = next(fp)
            line = next(fp).split(',')
            npts = int(line[0].split('=')[1]) #kayıt sayısı
            dt = float(line[1].split('=')[1].split()[0]) # zaman aralığı
            acc = np.array([p for l in fp for p in l.split()]).astype(float) #ivme değerlerinin eklendiği array
            time = np.arange(0,npts*dt,dt) # zaman serisinin oluşturduğu array

            timegap = np.arange(time[-1],time[-1]+100+dt,dt) # zaman boşluğu serisi
            accgap = np.arange(0,len(timegap))*0 # ivme değerlerinin 0 olduğu seri


            mainafterTime.append(np.append(time,timegap)) #birleştirilmiş zaman serisi
            mainafterAcc.append(np.append(acc,accgap)) # birleştirilmiş ivme serisi
    if plot :
        for i in range(len(mainafterAcc)):
            fig, ax = plt.subplots(figsize=(20,5))
            fig.subplots_adjust(bottom=0.15, left=0.2)
            ax.grid()
            ax.plot(mainafterTime[i],mainafterAcc[i])
            ax.set_xlabel('Time [Sec]')
            ax.set_ylabel('Acceleration [cm/sn]')

    return mainafterTime,mainafterAcc,dt
          #plt.figure(figsize=[60,30])
            #plt.plot(mainafterTime[i],mainafterAcc[i]) 

def ChangeTimeSeriesForPGA(TimeSeries : pd.DataFrame,targetPGA : float) -> pd.DataFrame:
    """ Zaman serisini PGA değerine göre verilen ivme değeri ile orantılayarak yeniden oluşturur."""
    realPGA = TimeSeries.Acceleration.max()
    coef = targetPGA/realPGA
    newTimeSeries = [acc*coef for acc in TimeSeries.Acceleration]
    TimeSeries["ChangedAcc"] = newTimeSeries
    del newTimeSeries
    return TimeSeries


