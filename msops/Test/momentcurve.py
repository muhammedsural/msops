from cProfile import label
from msops.Analysis.builder import opsbuild as bops
from msops.Material.material import Concrete
from msops.Section.sectionprop import RecSection
from msops.Units.Unit import Unit
from msops.Material.tbdymaterials import tbdy_mander
from msops.Plotter.msplotter import plotter as pltr
import matplotlib.pyplot as plt
import openseespy.opensees as ops

# MOMENT - CURVATURE ANALYSİS
# Muhammed Sural
for axial in range(0,10,5):
    #=================================================================================================================================
    # Beton malzeme özellikleri
    conc = Concrete(name='C25/30',fck=25)
    #=================================================================================================================================
    #Kesit özellikleri ile ilgili bilgiler
    sec = RecSection(Id=1,name='C4040',b=700,h=300,cover=25,k=1,numrebars=[6,2,6],dia_rebars=[14,14,14])
    sec.set_reinforcement_conc(plot=False)
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
    ops.wipe()
    ops.model('basic','-ndm',2,'-ndf',3)
    ops.uniaxialMaterial('Concrete01',IDconcCover,*coverInputarray) 
    ops.uniaxialMaterial('Concrete01',IDconcCore , *coreInputarray)
    ops.uniaxialMaterial('Steel02',IDSteel, fsy, Es, bs,  R0, cR1, cR2)       

    #Fiber section defined
    #=============================================================================================================================
    Fiber = 1
    bararea = sec.dia_rebars[3]*Unit.mm*10**-3
    print(Fiber,sec.h,sec.b,sec.cover,sec.cover,sec.numrebars[0],bararea,sec.numrebars[2],bararea,sec.numrebars[1],bararea,)
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
    Load,Curvature = bops.MomentCurvature(secTag=Fiber,axialLoad=axial,maxK=15,HSec=sec.h*Unit.mm,BSec=sec.b,numIncr=100)
    pltr.plot_Moment_Curvature(Curvature,Load,label=f'Axial={axial}')
    ops.wipe('all')