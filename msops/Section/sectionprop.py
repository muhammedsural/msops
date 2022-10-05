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
    """
        INPUT
            name           = Kesit adı           
            b              = Kesit genişliği
            h              = Kesit yüksekliği
            coreConc       = Çekirdek malzemesi
            coverConc      = Kabuk malzemesi
            matRebars      = Çelik malzemesi
            numrebars      = Donatı adetleri [üst başlık, gövde , alt başlık]
            dia_rebars     = Donatı çapı [üst başlık, gövde , alt başlık] => [üst başlık, gövde , alt başlık,total donatı alanı,% pursantaj]
            fiberData      =
    """
    name               :  Optional[str]
    b                  :  Optional[float]
    h                  :  Optional[float]
    coreConc           :  Optional[opsmaterial]  = opsmaterial(3,1,[-26922.92,-0.00546,-24433.71,-0.01390],defaultopsMat(),stress_strain_test=False)
    coverConc          :  Optional[opsmaterial]  = opsmaterial(3,2,[-20000.0, -0.002,-17800.05,-0.00349],defaultopsMat(),stress_strain_test=False)
    matRebars          :  Optional[opsmaterial]  = opsmaterial(22,3,[ 500*un.MPa, 2*10**5*un.MPa, 0.01,18, 0.925, 0.15],defaultopsMat(),stress_strain_test=False)
    numrebars          :  Optional[list[int]]    = field(default_factory=list,metadata={'info': ['top','int','Bot']})
    dia_rebars         :  Optional[list[float]]  = field(default_factory=list,metadata={'info': ['diameter','rebar_area']}) 
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
        self.I23 = 0
        
        if len(self.numrebars) == 0:
            self.numrebars = [3,2,3]
            
        if len(self.dia_rebars) <= 3:
            fi = 22 * un.mm
            rebars = [fi,fi,fi]
            top = 3.14*rebars[0]**2/4
            int = 3.14*rebars[1]**2/4
            bot = 3.14*rebars[2]**2/4
            total = top*self.numrebars[0]+int*self.numrebars[1]+bot*self.numrebars[2]
            pursantaj = (total/self.area)*100
            self.dia_rebars = [top,int,bot,total,pursantaj]
        else:
            top = 3.14*self.dia_rebars[0]**2/4
            int = 3.14*self.dia_rebars[1]**2/4
            bot = 3.14*self.dia_rebars[2]**2/4
            total = top*self.numrebars[0]+int*self.numrebars[1]+bot*self.numrebars[2]
            pursantaj = (total/self.area)*100
            self.dia_rebars.append(total) 
            self.dia_rebars.append(pursantaj)
            
        if self.coreConc.stress_strain_test:
            self.coreConc.testMaterial(materialName=self.coreConc.MaterialType)
        if self.coverConc.stress_strain_test:
            self.coverConc.testMaterial(materialName=self.coverConc.MaterialType)
        if self.matRebars.stress_strain_test:
            self.matRebars.testMaterial(materialName=self.matRebars.MaterialType)
        
    def asdict(self):
        return asdict()
