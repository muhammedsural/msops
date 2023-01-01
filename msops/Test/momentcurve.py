from pandas import DataFrame
from msops.Analysis.builder import opsbuild as bops
from msops.Material.material import Concrete
from msops.Section.sectionprop import RecSection
from msops.Units.Unit import Unit
from msops.Plotter.msplotter import plotter as pltr
import openseespy.opensees as ops

# MOMENT - CURVATURE ANALYSİS
# Muhammed Sural
#=================================================================================================================================
# Beton malzeme özellikleri
conc = Concrete(name='C25/30',fck=25)
#=================================================================================================================================
#Kesit özellikleri ile ilgili bilgiler
#covers =  opsm(2,2,[-25000.0, -0.002,-21569.71,-0.00349],stress_strain_test=False)
#covers.testMaterial(materialName=covers.get_MaterialType(),scaleFactor=covers.inputArray[3])
#cores =   opsm(3,2,[-28876.095, -0.00355, -21886.619, -0.01223, 2.5, 1250.0, 0.1],stress_strain_test=False)
#cores.testMaterial(materialName=cores.get_MaterialType(),scaleFactor=cores.inputArray[3])

sec = RecSection(Id=1,name='B3050',b=300,h=500,cover=25,k=0.35,numrebars=[4,0,3],dia_rebars=[14,14,14])
sec.set_reinforcement_conc(fck = conc.fck,s = 100, etriye_capi = 10, x_koladeti=2,y_koladeti=2,tension=True,plot=False)
coreMaterial = sec.coreConc.testMaterial(materialName=sec.coreConc.get_MaterialType())
coverMaterial = sec.coverConc.testMaterial(materialName=sec.coverConc.get_MaterialType())
coverMatName,coreMatName = sec.coverConc.get_MaterialType(),sec.coreConc.get_MaterialType()
#=================================================================================================================================
#Concrete and rebar material defined
# REINFORCING STEEL
#=============================================================================
fsy = 420*Unit.MPa;     # Yield stress
Es = 2*10**5*Unit.MPa;  # Young's modulus
bs = 0.01           # strain-hardening ratio
R0 = 18             # control the transition from elastic to plastic branches
cR1 = 0.925         # control the transition from elastic to plastic branches
cR2 = 0.15          # control the transition from elastic to plastic branches
minStrain = -0.1    # minimum steel strain in the fibers (steel buckling)
maxStrain = 0.1     # maximum steel strain in the fibers (steel rupture)
IDconcCover = 1   # Tag for unconfined concrete material
IDconcCore = 2    # Tag for confined concrete material
IDSteel = 3       # Tag for steel material without min-max properties
IDMinMaxSteel = 4 # Tag for steel material with min-max properties

coverInputarray = sec.coverConc.inputArray
coreInputarray  = sec.coreConc.inputArray
print(f"cover conc inputs => {sec.coverConc.inputArray}")
print(f"core conc inputs  => {sec.coreConc.inputArray}")
ops.wipe()
ops.model('basic','-ndm',2,'-ndf',3)

ops.uniaxialMaterial(coverMatName,IDconcCover,*coverInputarray) 
ops.uniaxialMaterial(coreMatName,IDconcCore , *coreInputarray)
ops.uniaxialMaterial('Steel02',IDSteel, fsy, Es, bs,  R0, cR1, cR2)       
#Fiber section defined
#=============================================================================================================================
Fiber = 1
bararea = sec.dia_rebars[3]*Unit.mm*10**-3
# print(bararea)
# print(Fiber,sec.h,sec.b,sec.cover*Unit.mm,sec.cover*Unit.mm,sec.numrebars[0],bararea,sec.numrebars[2],bararea,sec.numrebars[1],bararea,)
fiber_sec = bops.BuildRCrectSection(Fiber,
                            sec.h*Unit.mm,
                            sec.b*Unit.mm,
                            sec.cover*Unit.mm,
                            sec.cover*Unit.mm,
                            IDconcCore,
                            IDconcCover,
                            IDSteel,
                            sec.numrebars[0],
                            bararea,
                            sec.numrebars[2],
                            bararea,
                            sec.numrebars[1],
                            bararea,
                            nfCoreY=20,
                            nfCoreZ=20,
                            nfCoverY=20,
                            nfCoverZ=20,
                            pflag=1
                            )
pltr.plot_fiber_section(SecID=Fiber,fiber_sec=fiber_sec)

#Moment-Curvature Analysis
#================================================================================================================================
axial = 0
targetCurvature = 0.1
Load,Curvature = bops.MomentCurvature(secTag=Fiber,axialLoad=axial,maxK=targetCurvature,HSec=sec.h*Unit.mm,BSec=sec.b*Unit.mm,numIncr=20)
pltr.plot_Moment_Curvature(Curvature,Load,label=f'Axial={axial}')
ops.wipe('all')

MomCurv = DataFrame({"Moment":Load,"Curvature":Curvature})
print(MomCurv)
MomCurv.to_excel("Beam33050.xlsx",sheet_name='769-462')