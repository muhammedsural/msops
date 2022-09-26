
from dataclasses import dataclass,field,asdict
from typing import Optional
from ..Units.Unit import Unit as un
from ..Material.material import defaultopsMat
from ..Material.material import opsmaterial
from numpy import array

class rcsectionprop():
    def __init__(self,sectionname:str,width:float,height:float,concmaterial:None,rebarmaterial:None,numrebars : None,dialongrebars:None,density=24.99):
        """
            sectionname:str,width:float,height:float,
            concmaterial: Default None must be string,
            rebarmaterial: Default None must be string,
            numrebars : Default None must be array [number of top of section rebars,number of intermediate of section rebars,number of bot of section rebars ],
            dialongrebars: Default None rebars diameter [diameter,area]
        """     
        
        """
        f_sy    : akma anındaki gerilme
        eps_sy  : akma başlangıcı şekildeğiştirme değeri
        eps_sh  : akma sonu pekleşme başlangıcı şekildeğiştirme değeri 
        eps_su  : kopma şekildeğiştirme değeri
        f_su    : kopma gerilmesi/akma gerilmesi oranı minimumlar alınmıştır TBDY-bölüm5 Tablo5A.1
        Es      : 2*10**5
        """
        
        #                      f_sy    ,eps_sy,eps_sh,eps_su, Kres, Es
        steel = {
                    "S220" :[220*un.MPa,0.0011, 0.011, 0.12 , 1.20,2*10**5*un.MPa],
                    "S420" :[420*un.MPa,0.0021, 0.008, 0.08 , 1.15,2*10**5*un.MPa],
                    "B420C":[420*un.MPa,0.0021, 0.008, 0.08 , 1.15,2*10**5*un.MPa],
                    "B500C":[500*un.MPa,0.0025, 0.008, 0.08 , 1.15,2*10**5*un.MPa]
                }
        
        #                     fck    ,    fctk   ,     Ec
        conc = {
                    "C16": [16*un.MPa, 1.4*un.MPa, 27000*un.MPa],
                    "C18": [18*un.MPa, 1.5*un.MPa, 27500*un.MPa],
                    "C20": [20*un.MPa, 1.6*un.MPa, 28000*un.MPa],
                    "C25": [25*un.MPa, 1.8*un.MPa, 30000*un.MPa],
                    "C30": [30*un.MPa, 1.9*un.MPa, 32000*un.MPa],
                    "C35": [35*un.MPa, 2.1*un.MPa, 33000*un.MPa],
                    "C40": [40*un.MPa, 2.2*un.MPa, 34000*un.MPa],
                    "C45": [45*un.MPa, 2.3*un.MPa, 36000*un.MPa],
                    "C50": [50*un.MPa, 2.5*un.MPa, 37000*un.MPa]
               }
        
        
        if concmaterial not in conc.keys():
            print(f"beton malzemesi girilmemiş default olarak C20 atandı...")
            concmaterial = "C20"
        if rebarmaterial not in steel.keys():
            print(f"çelik malzemesi girilmemiş default olarak B500C atandı...")
            rebarmaterial = "B500C"
        if numrebars is None:
            print(f"Kesitte bulunan donatı sayıları girilmemiş default değerler atandı...")
            numrebars = [3,2,3]
        if dialongrebars is None:
            print(f"Kesitteki boyuna donatı çapı girilmemiş default değer atandı...")
            dialongrebars = 22 *un.mm
        self.properties = {
            "name"               : sectionname,
            "b"                  : width,
            "h"                  : height,
            "area"               : width*height,
            "density"            : density,
            "I33"                : width*height**3/12,
            "I22"                : height*width**3/12,
            "I23"                : 0,
            "concmat"            : conc[concmaterial],
            "rebarmat"           : steel[rebarmaterial],
            "numrebars"          : numrebars,
            "longnitudediarebars": [dialongrebars,3.14*dialongrebars**2/4]
        }


@dataclass
class RecSection(object):
    
    name               :  Optional[str]
    b                  :  Optional[float]
    h                  :  Optional[float]
    area               :  Optional[float]        = None
    I33                :  Optional[float]        = None
    I22                :  Optional[float]        = None
    I23                :  Optional[float]        = 0. 
    coreConc           :  Optional[opsmaterial]  = None
    coverConc          :  Optional[opsmaterial]  = None
    matRebars          :  Optional[opsmaterial]  = None
    numrebars          :  Optional[list[int]]    = field(default_factory=list,metadata={'info': ['top','int','Bot']})
    dia_area_rebars    :  Optional[list[float]]  = field(default_factory=list,metadata={'info': ['diameter','rebar_area']}) 
    fiberData          :  Optional[list[float]]  = None       
    
    
    def calcArea(self):
        return self.b*self.h
    
    def calcI33(self):
        return self.b*self.h**3/12
    
    def calcI22(self):
        return self.h*self.b**3/12
    
    def __post_init__(self):
        self.area = self.calcArea()
        self.I33  = self.calcI33() 
        self.I22  = self.calcI22()
        
        if len(self.numrebars) == 0:
            self.numrebars = [3,2,3]
            
        if len(self.dia_area_rebars) == 0:
            fi = 22 * un.mm
            self.dia_area_rebars = [fi , 3.14*fi**2/4]
            
        if self.coreConc is None or self.coverConc is None or self.matRebars is None:
            y = defaultopsMat()
            if self.coreConc is None:
                print('Çekirdek beton malzemesinin bilgileri girilmemiş default değerler atandı...')
                defaultcoreconc = opsmaterial('Concrete02',1,[-26922.92,-0.00546,-24433.71,-0.01390],y)
                self.coreConc   = defaultcoreconc
            if self.coverConc is None:
                print('Kabuk beton malzemesinin bilgileri girilmemiş default değerler atandı...')
                defaultcoverConc = opsmaterial('Concrete02',2,[-20000.0, -0.002,-17800.05,-0.00349],y)
                self.coverConc   = defaultcoverConc
            if self.matRebars is None:
                print('Donatı malzemesinin bilgileri girilmemiş default değerler atandı...')
                fsy = 500*un.MPa;     # Yield stress
                Es  = 2*10**5*un.MPa;     # Young's modulus
                bs  = 0.01           # strain-hardening ratio
                R0  = 18             # control the transition from elastic to plastic branches
                cR1 = 0.925         # control the transition from elastic to plastic branches
                cR2 = 0.15          # control the transition from elastic to plastic branches
                defaultmatRebars = opsmaterial('Steel02',3,[ fsy, Es, bs,  R0, cR1, cR2],y)
                self.matRebars   = defaultmatRebars
    
    def asdict(self):
        return asdict()
