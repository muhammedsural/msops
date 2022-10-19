from ast import Return
from dataclasses import dataclass,field,asdict
from typing import Dict, Optional
from ..Units.Unit import Unit as un
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
# https://github.com/jupyter-widgets/ipywidgets/issues/1853
from ipywidgets.widgets.interaction import show_inline_matplotlib_plots



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
        show_inline_matplotlib_plots()

@dataclass
class Concrete():
    name          : str   = None
    density       : float = 24.99
    fck           : float = None
    Ec            : float = field(default_factory=float)
    poisson       : float = 0.2
    shear_modules : float = field(default_factory=float)
    
    def __calc_shear_modules(self):
        Gc = self.Ec/(2*(1+self.poisson))
        self.shear_modules = Gc
        
    def __calc_young_modules(self):
        Ec = 57000*un.MPa*(self.fck/un.MPa)**0.5
        self.Ec = Ec
        
    def __post_init__(self):
        self.__calc_young_modules()
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
    
    def asdict():
        return asdict()
    def __post_init__(self):
        steel = {
                    "S220" :[220*un.MPa,0.0011, 0.011, 0.12 , 1.20],
                    "S420" :[420*un.MPa,0.0021, 0.008, 0.08 , 1.15],
                    "B420C":[420*un.MPa,0.0021, 0.008, 0.08 , 1.15],
                    "B500C":[500*un.MPa,0.0025, 0.008, 0.08 , 1.15]
                }
        if self.name is not None:
            self.f_sy   = steel[self.name][0]
            self.eps_sy = steel[self.name][1]
            self.eps_sh = steel[self.name][2]
            self.eps_su = steel[self.name][3]
            self.Kres   = steel[self.name][4]
            
@dataclass
class opsmaterial(object):
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
                            26:'ViscousDamper'    
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
    
    def __post_init__(self): 
        material=self.get_MaterialType() 
        if self.young_Module is None:
            self.young_Module = 0      
        if self.stress_strain_test is True:
            
            self.testMaterial(materialName=material,scaleFactor=self.inputArray[3])
            
        self.OpenSeesMaterialDefaultValues = self.defaultOpsMaterial()
        self.OpenSeesMaterialDefaultValues[material]=self.inputArray
        pass
    
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
            26:'ViscousDamper'    
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
        OpenSeesMaterialDefaultValues[self.MaterialType]=self.inputArray

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
        show_inline_matplotlib_plots()