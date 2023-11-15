from dataclasses import dataclass,field,asdict
from typing import Dict, Optional
from msops.Units.Unit import Unit as un
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt


# conc = {
#                     "C16": [16*un.MPa, 1.4*un.MPa, 27000*un.MPa],
#                     "C18": [18*un.MPa, 1.5*un.MPa, 27500*un.MPa],
#                     "C20": [20*un.MPa, 1.6*un.MPa, 28000*un.MPa],
#                     "C25": [25*un.MPa, 1.8*un.MPa, 30000*un.MPa],
#                     "C30": [30*un.MPa, 1.9*un.MPa, 32000*un.MPa],
#                     "C35": [35*un.MPa, 2.1*un.MPa, 33000*un.MPa],
#                     "C40": [40*un.MPa, 2.2*un.MPa, 34000*un.MPa],
#                     "C45": [45*un.MPa, 2.3*un.MPa, 36000*un.MPa],
#                     "C50": [50*un.MPa, 2.5*un.MPa, 37000*un.MPa]
#                }

class MaterialPropManager:
    
    def calc_shear_modules(Ec : float, poisson : float) -> float:
        return round(Ec/(2*(1+poisson)),4)
        
    def calc_young_modules(fck : int) -> float:
        return round(57000*un.MPa*(fck)**0.5,4)
    
    def defineStrainHistory(peaksArray,scaleFactor,nSteps,nCycles) -> list:
        strain = []
        for thisPeak in peaksArray:
            for i in range(nCycles):
                strain = np.append(strain,np.linspace(0,thisPeak*scaleFactor,nSteps))
                strain = np.append(strain,np.linspace(thisPeak*scaleFactor,-thisPeak*scaleFactor,nSteps))
                strain = np.append(strain,np.linspace(-thisPeak*scaleFactor,0,nSteps))

        return strain
    
    def default_Openseespy_Material() -> dict:
        OpenSeesMaterialDefaultValues = {}
    
        OpenSeesMaterialDefaultValues['Bond_SP01'         ]      =[60.0,0.01,75.0,0.1,0.4,0.75]
        OpenSeesMaterialDefaultValues['Cast'              ]      =[10.0,1.0,0.1,60.0,29000.0,1.0,0.05,18.0,0.925,0.15,0.0,1.0,0.0,1.0]
        OpenSeesMaterialDefaultValues['Concrete01'        ]      =[-4.4,-0.002,-4.576,-0.04]
        OpenSeesMaterialDefaultValues['Concrete02'        ]      =[-4.4,-0.002,-4.576,-0.04,0.1,0.572,286.0]
        OpenSeesMaterialDefaultValues['Concrete04'        ]      =[-4.4,-0.002,-0.2,3700,0.5,0.001,0.1]
        OpenSeesMaterialDefaultValues['Concrete06'        ]      =[-4.4,-0.002,2.0,1.0,0.32,0.44,0.0002,4.0,0.08]
        OpenSeesMaterialDefaultValues['Concrete07'        ]      =[-4.4,-0.002,3700,0.44,0.0002,2.0,2.3,3.97]
        OpenSeesMaterialDefaultValues['Elastic'           ]      =[29000.0,0.0,29000.0]
        OpenSeesMaterialDefaultValues['ElasticPP'         ]      =[29000.0,0.0020689655,-0.0020689655,0.0]
        OpenSeesMaterialDefaultValues['ElasticPPGap'      ]      =[29000.0,60.0,0.001,0.0]
        OpenSeesMaterialDefaultValues['ENT'               ]      =[29000.0]
        OpenSeesMaterialDefaultValues['Hysteretic'        ]      =[60.0,0.003,78.0,0.024,61.2,0.1,-60.0,-0.003,-78.0,-0.024,-61.2,-0.1,1.0,1.0,0.0,0.0,0.0]
        OpenSeesMaterialDefaultValues['Bilin'             ]      =[29000.0,0.01,0.01,60.0,-60.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0,0.0]
        OpenSeesMaterialDefaultValues['ModIMKPeakOriented']      =[29000.0,0.01,0.01,60.0,-60.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0]
        OpenSeesMaterialDefaultValues['ModIMKPinching'    ]      =[29000.0,0.01,0.01,60.0,-60.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0]
        OpenSeesMaterialDefaultValues['Pinching4'         ]      =[60.0,0.003,78.0,0.024,61.2,0.1,61.0,1.0,-60.0,-0.003,-78.0,-0.024,-61.2,-0.1,-61.0,-1.0,0.5,0.25,0.05,0.5,0.25,0.05,1.0,0.2,0.3,0.2,0.9,0.5,0.5,2.0,2.0,0.5,1.0,0.0,1.0,1.0,0.9,10.0,'energy']
        OpenSeesMaterialDefaultValues['PySimple1'         ]      =[1,40.0,0.01,200.0,0.0]
        OpenSeesMaterialDefaultValues['QzSimple1'         ]      =[1,40.0,0.01,0.0,0.0]
        OpenSeesMaterialDefaultValues['ReinforcingSteel'  ]      =[60.0,66.0,29000.0,2900.0,0.008,0.02,'-GABuck',6,1,0.4,0.5,'-DMBuck',6,1,'-CMFatigue',0.26,0.506,0.389,'-IsoHard',4.3,0.01,'-MPCurveParams',0.33,18,4]
        OpenSeesMaterialDefaultValues['SAWS'              ]      =[15.799848,0.545094768764215,1.04095,159.83706,0.1022018,-0.0324118361176701,1.0,0.0692552,0.8,1.1]
        OpenSeesMaterialDefaultValues['SelfCentering'     ]      =[29000.0,2900.0,60.0,0.1,0,0,1]
        OpenSeesMaterialDefaultValues['Steel01'           ]      =[60.0,29000.0,0.05,0.0,1.0,0.0,1.0]
        OpenSeesMaterialDefaultValues['Steel02'           ]      =[60.0,29000.0,0.05,18.0,0.925,0.15,0.0,1.0,0.0,1.0,0.0]
        OpenSeesMaterialDefaultValues['SteelMPF'          ]      =[60.0,40.0,29000.0,0.05,0.01,20.0,0.925,0.15,0.0,1.0,0.0,1.0]
        OpenSeesMaterialDefaultValues['TzSimple1'         ]      =[1,40.0,0.01,0.0]
        OpenSeesMaterialDefaultValues['UVCuniaxial'       ]      =[29000.0,60.0,122.63,19.74,143.49,248.14,2,31638.0,277.32,1548.6,9.04]
        OpenSeesMaterialDefaultValues['ViscousDamper'     ]      =[29000.0,200.0,0.3,0.0,1,1e-6,1e-10,15.0]
        
        return OpenSeesMaterialDefaultValues
    
    def plot_stressstrain(self,strain : list, stress : list ,figSizeH = 6, figSizeV = 4, DPI = 200 ,Name = None,**kwargs):
        
        figModel = plt.figure(f'Material Response {Name}',figsize=(figSizeH,figSizeV), dpi=DPI, facecolor='w', edgecolor='k' )
        axModel = figModel.add_subplot(1,1,1)    
        axModel.plot(strain, stress,**kwargs)
        axModel.grid()
        axModel.set_xlabel('Strain')
        axModel.set_ylabel('Stress')
        axModel.set_title(Name + ' Material Response')
        plt.show()

@dataclass
class MechanicalProperty:
    """
    Material mechanical properties 
    E : young modulus
    U : poission ratio
    A : thermal coef
    G : shear modulus
    """
    E : float = field(default_factory=float)
    U : float = field(default_factory=float)
    A : float = field(default_factory=float)
    G : float = field(default_factory=float)

    def __repr__(self) -> str:
        return f'Young Modulus : {self.E}, poisson : {self.U}, Thermal coef : {self.A}, shear_modules : {self.G} '

    def asdict():
        return asdict()
    
@dataclass        
class concrete(MaterialPropManager):
    Name          : str  
    Fck           : float 
    Fctk          : Optional[float] = field(default_factory=float)
    Density       : float = 24.99
    PropMech      : Optional[MechanicalProperty] = field(default_factory=MechanicalProperty)

    def __repr__(self) -> str:
        return f'Name : {self.Name}, fck : {self.Fck}, Ec : {self.PropMech.E}'
    
    def asdict():
        return asdict()
    
    def __post_init__(self):
        if self.PropMech.E == 0.0:
            poission = 0.2
            young_mod = MaterialPropManager.calc_young_modules(fck= self.Fck)
            shear_mod = MaterialPropManager.calc_shear_modules(Ec= young_mod,poisson=poission)
            self.PropMech = MechanicalProperty(E= young_mod, U= poission, A= 0.0, G= shear_mod)
        self.Fctk = self.Fck * 0.1

class defaultopsMat():
    
    def __init__(self) -> None:
        
        self.OpenSeesMaterialDefaultValues = {}
        
        self.OpenSeesMaterialDefaultValues['Bond_SP01']          =[60.0,0.01,75.0,0.1,0.4,0.75]
        self.OpenSeesMaterialDefaultValues['Cast']               =[10.0,1.0,0.1,60.0,29000.0,1.0,0.05,18.0,0.925,0.15,0.0,1.0,0.0,1.0]
        self.OpenSeesMaterialDefaultValues['Concrete01']         =[-4.4,-0.002,-4.576,-0.04]
        self.OpenSeesMaterialDefaultValues['Concrete02']         =[-4.4,-0.002,-4.576,-0.04,0.1,0.572,286.0]
        self.OpenSeesMaterialDefaultValues['Concrete04']         =[-4.4,-0.002,-0.2,3700,0.5,0.001,0.1]
        self.OpenSeesMaterialDefaultValues['Concrete06']         =[-4.4,-0.002,2.0,1.0,0.32,0.44,0.0002,4.0,0.08]
        self.OpenSeesMaterialDefaultValues['Concrete07']         =[-4.4,-0.002,3700,0.44,0.0002,2.0,2.3,3.97]
        self.OpenSeesMaterialDefaultValues['Elastic']            =[29000.0,0.0,29000.0]
        self.OpenSeesMaterialDefaultValues['ElasticPP']          =[29000.0,0.0020689655,-0.0020689655,0.0]
        self.OpenSeesMaterialDefaultValues['ElasticPPGap']       =[29000.0,60.0,0.001,0.0]
        self.OpenSeesMaterialDefaultValues['ENT']                =[29000.0]
        self.OpenSeesMaterialDefaultValues['Hysteretic']         =[60.0,0.003,78.0,0.024,61.2,0.1,-60.0,-0.003,-78.0,-0.024,-61.2,-0.1,1.0,1.0,0.0,0.0,0.0]
        self.OpenSeesMaterialDefaultValues['Bilin']              =[29000.0,0.01,0.01,60.0,-60.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0,0.0]
        self.OpenSeesMaterialDefaultValues['ModIMKPeakOriented'] =[29000.0,0.01,0.01,60.0,-60.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0]
        self.OpenSeesMaterialDefaultValues['ModIMKPinching']     =[29000.0,0.01,0.01,60.0,-60.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0]
        self.OpenSeesMaterialDefaultValues['Pinching4']          =[60.0,0.003,78.0,0.024,61.2,0.1,61.0,1.0,-60.0,-0.003,-78.0,-0.024,-61.2,-0.1,-61.0,-1.0,0.5,0.25,0.05,0.5,0.25,0.05,1.0,0.2,0.3,0.2,0.9,0.5,0.5,2.0,2.0,0.5,1.0,0.0,1.0,1.0,0.9,10.0,'energy']
        self.OpenSeesMaterialDefaultValues['PySimple1']          =[1,40.0,0.01,200.0,0.0]
        self.OpenSeesMaterialDefaultValues['QzSimple1']          =[1,40.0,0.01,0.0,0.0]
        self.OpenSeesMaterialDefaultValues['ReinforcingSteel']   =[60.0,66.0,29000.0,2900.0,0.008,0.02,'-GABuck',6,1,0.4,0.5,'-DMBuck',6,1,'-CMFatigue',0.26,0.506,0.389,'-IsoHard',4.3,0.01,'-MPCurveParams',0.33,18,4]
        self.OpenSeesMaterialDefaultValues['SAWS']               =[15.799848,0.545094768764215,1.04095,159.83706,0.1022018,-0.0324118361176701,1.0,0.0692552,0.8,1.1]
        self.OpenSeesMaterialDefaultValues['SelfCentering']      =[29000.0,2900.0,60.0,0.1,0,0,1]
        self.OpenSeesMaterialDefaultValues['Steel01']            =[60.0,29000.0,0.05,0.0,1.0,0.0,1.0]
        self.OpenSeesMaterialDefaultValues['Steel02']            =[60.0,29000.0,0.05,18.0,0.925,0.15,0.0,1.0,0.0,1.0,0.0]
        self.OpenSeesMaterialDefaultValues['SteelMPF']           =[60.0,40.0,29000.0,0.05,0.01,20.0,0.925,0.15,0.0,1.0,0.0,1.0]
        self.OpenSeesMaterialDefaultValues['TzSimple1']          =[1,40.0,0.01,0.0]
        self.OpenSeesMaterialDefaultValues['UVCuniaxial']        =[29000.0,60.0,122.63,19.74,143.49,248.14,2,31638.0,277.32,1548.6,9.04]
        self.OpenSeesMaterialDefaultValues['ViscousDamper']      =[29000.0,200.0,0.3,0.0,1,1e-6,1e-10,15.0]
        pass
    
    def defineStrainHistory(self,peaksArray,scaleFactor,nSteps,nCycles):
        strain = []
        for thisPeak in peaksArray:
            for i in range(nCycles):
                strain = np.append(strain,np.linspace(0,thisPeak*scaleFactor,nSteps))
                strain = np.append(strain,np.linspace(thisPeak*scaleFactor,-thisPeak*scaleFactor,nSteps))
                strain = np.append(strain,np.linspace(-thisPeak*scaleFactor,0,nSteps))

        return strain
    
    def testMaterial(self,materialName=None,scaleFactor = 0.01):
        AllStressStrain = {}

        peaksArray=[1,2,4,8,10]
        nSteps = 100
        nCycles = 3
        strain=self.defineStrainHistory(peaksArray,scaleFactor,nSteps,nCycles)

        thisCount = 0
        for thisMaterial in self.OpenSeesMaterialDefaultValues.keys():
            if materialName == None:
                return None
            if materialName != None and materialName == thisMaterial:
                print(thisMaterial)
                counter = thisCount + 1
                ops.wipe()
                materialTag = 99
                
                inputArray = self.OpenSeesMaterialDefaultValues[thisMaterial]
                
                print(inputArray)
                ops.uniaxialMaterial(thisMaterial,materialTag,*inputArray)
                ops.testUniaxialMaterial(materialTag)
                stress = []
                for eps in strain:
                    ops.setStrain(eps)
                    stress.append(ops.getStress())
                    tangent = ops.getTangent() # Not used
                #global AllStressStrain
                thisCount = len(list(AllStressStrain.keys()))
                thisKey = 'Run' + str(thisCount+1) + ' ' + thisMaterial
                AllStressStrain[thisKey] = {}
                AllStressStrain[thisKey]['strain'] = strain
                AllStressStrain[thisKey]['stress'] = stress
                self.plot_stressstrain(strain=AllStressStrain[thisKey]['strain'],stress=AllStressStrain[thisKey]['stress'],Name=thisKey,linewidth='1',label=thisMaterial,marker = '')
                
    def plot_stressstrain(self,strain : list, stress : list ,figSizeH = 6, figSizeV = 4, DPI = 200 ,Name = None,**kwargs):
        
        figModel = plt.figure(f'Material Response {Name}',figsize=(figSizeH,figSizeV), dpi=DPI, facecolor='w', edgecolor='k' )
        axModel = figModel.add_subplot(1,1,1)    
        axModel.plot(strain, stress,**kwargs)
        axModel.grid()
        axModel.set_xlabel('Strain')
        axModel.set_ylabel('Stress')
        axModel.set_title(Name + ' Material Response')
        plt.show()

@dataclass
class Concrete():
    name          : str   = None
    fck           : float = None
    fctk          : float = None
    poisson       : float = 0.2
    density       : float = 24.99
    Ec            : float = field(default_factory=float)
    shear_modules : float = field(default_factory=float)
    
    
    def __repr__(self) -> str:
        return f'Name : {self.name}, fck : {self.fck}, Ec : {self.Ec}, poisson : {self.poisson}, shear_modules : {self.shear_modules} '
    
    def __calc_shear_modules(self):
        Gc = self.Ec/(2*(1+self.poisson))
        self.shear_modules = Gc
        
    def __calc_young_modules(self):
        Ec = 57000*un.MPa*(self.fck/un.MPa)**0.5
        self.Ec = Ec
    
    def __calc_fctk(self):
        self.fctk = self.fck * 0.1

    def __post_init__(self):
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
        if self.name is not None:
            if self.name not in conc.keys():
                raise KeyError(future_index = conc.keys(), target_index=self.name, message=f"Girilen beton malzemesi bulunamadi.... malzeme listesi => {conc.keys()} ")
            self.fck  = conc[self.name][0]
            self.fctk = conc[self.name][1]
            self.Ec   = conc[self.name][2]

            self.__calc_shear_modules()

        if self.Ec is None:
            self.__calc_young_modules()

        if self.shear_modules is None:
            self.__calc_shear_modules()
        
    def asdict():
        return asdict()
@dataclass
class Steel():
    name    : str   = 'B500C'
    density : float = 80.
    f_sy    : float = None
    eps_sy  : float = None
    eps_sh  : float = None
    eps_su  : float = None
    Kres    : float = None
    Es      : float = 2*10**5*un.MPa
    
    def __repr__(self) -> str:
        return f'Name : {self.name}, fsy : {self.f_sy}, eps_sy : {self.eps_sy}, eps_sh : {self.eps_sh}, eps_su : {self.eps_su}, {self.eps_su}, Es : {self.Es} '
    
    def asdict():
        return asdict()
    
    def __post_init__(self):
        steel = {
                    "S220" :[220,0.0011, 0.011, 0.12 , 1.20],
                    "S420" :[420,0.0021, 0.008, 0.08 , 1.15],
                    "B420C":[420,0.0021, 0.008, 0.08 , 1.15],
                    "B500C":[500,0.0025, 0.008, 0.08 , 1.15]
                }
        if self.name in steel.keys():
            self.f_sy   = steel[self.name][0]
            self.eps_sy = steel[self.name][1]
            self.eps_sh = steel[self.name][2]
            self.eps_su = steel[self.name][3]
            self.Kres   = steel[self.name][4]

@dataclass
class opsmaterial:
    """
    MaterialTypeIndex =>    0:'Bond_SP01'        ,
                            1:'Cast'             ,
                            2:'Concrete01'       ,
                            3:'Concrete02'       ,
                            4:'Concrete04'       ,
                            5:'Concrete06'       ,
                            6:'Concrete07'       ,
                            7:'Elastic'          ,
                            8:'ElasticPP'        ,
                            9:'ElasticPPGap'     ,
                            10:'ENT'              ,
                            11:'Hysteretic'       ,
                            12:'Bilin'            ,
                            13:'ModIMKPeakOriented'
                            14:'ModIMKPinching'   ,
                            15:'Pinching4'        ,
                            16:'PySimple1'        ,
                            17:'QzSimple1'        ,
                            18:'ReinforcingSteel' ,
                            19:'SAWS'             ,
                            20:'SelfCentering'    ,
                            21:'Steel01'          ,
                            22:'Steel02'          ,
                            23:'SteelMPF'         ,
                            24:'TzSimple1'        ,
                            25:'UVCuniaxial'      ,
                            26:'ViscousDamper'    ,
        matTag =
        inputArray =
        young_Module =
        stress_strain_test =
    """
    MaterialTypeIndex  : Optional[int]           = None
    matTag             : Optional[int]           = None
    inputArray         : Optional[np.array]      = None
    young_Module       : Optional[float]         = None
    stress_strain_test : Optional[bool]          = False
    
    def __post_init__(self) -> None: 
        material=self.get_MaterialType()
        if self.young_Module is None:
            self.young_Module = 0
            
        if self.inputArray is None:
            defaultmat = self.defaultOpsMaterial()
            self.inputArray = defaultmat[material]
            
        if self.stress_strain_test is True:
            self.testMaterial(materialName=material)
    
    def get_MaterialType(self) -> str:
        matType= {
            0:'Bond_SP01'        ,
            1:'Cast'             ,
            2:'Concrete01'       ,
            3:'Concrete02'       ,
            4:'Concrete04'       ,
            5:'Concrete06'       ,
            6:'Concrete07'       ,
            7:'Elastic'          ,
            8:'ElasticPP'        ,
            9:'ElasticPPGap'     ,
            10:'ENT'              ,
            11:'Hysteretic'       ,
            12:'Bilin'            ,
            13:'ModIMKPeakOriented',
            14:'ModIMKPinching'   ,
            15:'Pinching4'        ,
            16:'PySimple1'        ,
            17:'QzSimple1'        ,
            18:'ReinforcingSteel' ,
            19:'SAWS'             ,
            20:'SelfCentering'    ,
            21:'Steel01'          ,
            22:'Steel02'          ,
            23:'SteelMPF'         ,
            24:'TzSimple1'        ,
            25:'UVCuniaxial'      ,
            26:'ViscousDamper'    ,
        }
        return matType[self.MaterialTypeIndex]
    
    def asdict(self) -> dict:
        return asdict(self)
    
    def calc_young_modules(self):
        pass

    def defaultOpsMaterial(self):
        OpenSeesMaterialDefaultValues = {}
        OpenSeesMaterialDefaultValues['Bond_SP01'         ]      =[60.0,0.01,75.0,0.1,0.4,0.75]
        OpenSeesMaterialDefaultValues['Cast'              ]      =[10.0,1.0,0.1,60.0,29000.0,1.0,0.05,18.0,0.925,0.15,0.0,1.0,0.0,1.0]
        OpenSeesMaterialDefaultValues['Concrete01'        ]      =[-4.4,-0.002,-4.576,-0.04]
        OpenSeesMaterialDefaultValues['Concrete02'        ]      =[-4.4,-0.002,-4.576,-0.04,0.1,0.572,286.0]
        OpenSeesMaterialDefaultValues['Concrete04'        ]      =[-4.4,-0.002,-0.2,3700,0.5,0.001,0.1]
        OpenSeesMaterialDefaultValues['Concrete06'        ]      =[-4.4,-0.002,2.0,1.0,0.32,0.44,0.0002,4.0,0.08]
        OpenSeesMaterialDefaultValues['Concrete07'        ]      =[-4.4,-0.002,3700,0.44,0.0002,2.0,2.3,3.97]
        OpenSeesMaterialDefaultValues['Elastic'           ]      =[29000.0,0.0,29000.0]
        OpenSeesMaterialDefaultValues['ElasticPP'         ]      =[29000.0,0.0020689655,-0.0020689655,0.0]
        OpenSeesMaterialDefaultValues['ElasticPPGap'      ]      =[29000.0,60.0,0.001,0.0]
        OpenSeesMaterialDefaultValues['ENT'               ]      =[29000.0]
        OpenSeesMaterialDefaultValues['Hysteretic'        ]      =[60.0,0.003,78.0,0.024,61.2,0.1,-60.0,-0.003,-78.0,-0.024,-61.2,-0.1,1.0,1.0,0.0,0.0,0.0]
        OpenSeesMaterialDefaultValues['Bilin'             ]      =[29000.0,0.01,0.01,60.0,-60.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0,0.0]
        OpenSeesMaterialDefaultValues['ModIMKPeakOriented']      =[29000.0,0.01,0.01,60.0,-60.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0]
        OpenSeesMaterialDefaultValues['ModIMKPinching'    ]      =[29000.0,0.01,0.01,60.0,-60.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,0.04,0.04,0.01,0.01,0.33,0.33,0.225349487,0.225349487,0.0,0.0]
        OpenSeesMaterialDefaultValues['Pinching4'         ]      =[60.0,0.003,78.0,0.024,61.2,0.1,61.0,1.0,-60.0,-0.003,-78.0,-0.024,-61.2,-0.1,-61.0,-1.0,0.5,0.25,0.05,0.5,0.25,0.05,1.0,0.2,0.3,0.2,0.9,0.5,0.5,2.0,2.0,0.5,1.0,0.0,1.0,1.0,0.9,10.0,'energy']
        OpenSeesMaterialDefaultValues['PySimple1'         ]      =[1,40.0,0.01,200.0,0.0]
        OpenSeesMaterialDefaultValues['QzSimple1'         ]      =[1,40.0,0.01,0.0,0.0]
        OpenSeesMaterialDefaultValues['ReinforcingSteel'  ]      =[60.0,66.0,29000.0,2900.0,0.008,0.02,'-GABuck',6,1,0.4,0.5,'-DMBuck',6,1,'-CMFatigue',0.26,0.506,0.389,'-IsoHard',4.3,0.01,'-MPCurveParams',0.33,18,4]
        OpenSeesMaterialDefaultValues['SAWS'              ]      =[15.799848,0.545094768764215,1.04095,159.83706,0.1022018,-0.0324118361176701,1.0,0.0692552,0.8,1.1]
        OpenSeesMaterialDefaultValues['SelfCentering'     ]      =[29000.0,2900.0,60.0,0.1,0,0,1]
        OpenSeesMaterialDefaultValues['Steel01'           ]      =[60.0,29000.0,0.05,0.0,1.0,0.0,1.0]
        OpenSeesMaterialDefaultValues['Steel02'           ]      =[60.0,29000.0,0.05,18.0,0.925,0.15,0.0,1.0,0.0,1.0,0.0]
        OpenSeesMaterialDefaultValues['SteelMPF'          ]      =[60.0,40.0,29000.0,0.05,0.01,20.0,0.925,0.15,0.0,1.0,0.0,1.0]
        OpenSeesMaterialDefaultValues['TzSimple1'         ]      =[1,40.0,0.01,0.0]
        OpenSeesMaterialDefaultValues['UVCuniaxial'       ]      =[29000.0,60.0,122.63,19.74,143.49,248.14,2,31638.0,277.32,1548.6,9.04]
        OpenSeesMaterialDefaultValues['ViscousDamper'     ]      =[29000.0,200.0,0.3,0.0,1,1e-6,1e-10,15.0]
        return OpenSeesMaterialDefaultValues

    def defineStrainHistory(self,peaksArray,scaleFactor,nSteps,nCycles):
        strain = []
        for thisPeak in peaksArray:
            for i in range(nCycles):
                strain = np.append(strain,np.linspace(0,thisPeak*scaleFactor,nSteps))
                strain = np.append(strain,np.linspace(thisPeak*scaleFactor,-thisPeak*scaleFactor,nSteps))
                strain = np.append(strain,np.linspace(-thisPeak*scaleFactor,0,nSteps))

        return strain
    
    def testMaterial(self,materialName=None,nCycles=3,nSteps=100,scaleFactor = 0.01,peaksArray =[1,2,4,8,10]):
        AllStressStrain = {}
        
        strain=self.defineStrainHistory(peaksArray,scaleFactor,nSteps,nCycles)
        thisCount = 0
        OpenSeesMaterialDefaultValues = self.defaultOpsMaterial()
        OpenSeesMaterialDefaultValues[self.get_MaterialType()]=self.inputArray

        for thisMaterial in OpenSeesMaterialDefaultValues.keys():
            if materialName == None:
                print(f"Material name ={materialName}")
                return None
            if materialName != None and materialName == thisMaterial:
                print(thisMaterial)
                counter = thisCount + 1
                ops.wipe()
                materialTag = 99
                
                inputArray = OpenSeesMaterialDefaultValues[thisMaterial]
                
                print(inputArray)
                ops.uniaxialMaterial(thisMaterial,materialTag,*inputArray)
                ops.testUniaxialMaterial(materialTag)
                stress = []
                for eps in strain:
                    ops.setStrain(eps)
                    stress.append(ops.getStress())
                    tangent = ops.getTangent() # Not used
                #global AllStressStrain
                thisCount = len(list(AllStressStrain.keys()))
                thisKey = 'Run' + str(thisCount+1) + ' ' + thisMaterial
                AllStressStrain[thisKey] = {}
                AllStressStrain[thisKey]['strain'] = strain
                AllStressStrain[thisKey]['stress'] = stress
                
                self.plot_stressstrain(strain=AllStressStrain[thisKey]['strain'],stress= AllStressStrain[thisKey]['stress'],Name=thisMaterial)
                
    def set_OpenSeesMaterialDefaultValues(self):
        self.OpenSeesMaterialDefaultValues[self.MaterialType]=self.inputArray
        self.testMaterial()
        
    def plot_stressstrain(self,strain : list, stress : list ,figSizeH = 6, figSizeV = 4, DPI = 200 ,Name = None,**kwargs):
        
        figModel = plt.figure(f'Material Response {Name}',figsize=(figSizeH,figSizeV), dpi=DPI, facecolor='w', edgecolor='k' )
        axModel = figModel.add_subplot(1,1,1)    
        axModel.plot(strain, stress,**kwargs)
        axModel.grid()
        axModel.set_xlabel('Strain')
        axModel.set_ylabel('Stress')
        axModel.set_title(Name + ' Material Response')
        plt.show()
        

"""def main() -> None:
   '''Main function'''
   mat = concrete(Name="C25",Fck=25)
   print(mat)

if __name__ == "__main__":
   main()"""