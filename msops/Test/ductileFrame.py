from msops.Analysis.builder import opsbuild as op
from msops.CreateModel.model import createPortalFrame
from msops.Codes.Tbdy import FundemantelParameters
from msops.Material import material,tbdymaterials
from msops.Plotter import msplotter
from msops.TimeSeries import ReadRecord
from msops.Units.Unit import Unit
from msops.Section import sectionprop as sp
from msops.Elements.Frame import ElasticBeamColumn,Column,Beam,FrameType,FrameDatas
from msops.Elements.Node import Node,NodeDatas
from msops.Elements.BeamIntegration import Integration,IntegrationType
from msops.Section.sectionprop import RecSection
import re
import pandas as pd
import openseespy.opensees as ops


path = "C:\\Users\\Lenovo\\Desktop\\github\\prj3386\\Tez\\EtabsModelleri\\S-5k-7-DTS1-ZB-Tasarim_overwrite.$et"
Frame_section = pd.DataFrame()

with open(path,"r") as et:
    data = et.read()

"""
frame_section = re.match("FRAMESECTION",data)  
print(frame_section)  """

"""spl = re.split("\$",data)
frame_section = re.split("\n",spl[10])
frame_section = [line.strip().split() for line in frame_section]
print(frame_section[0])
"""

#==========================================================================================================================================================================
#Create Auto Portal Frame
numBay=3
numFloor=5
firstFloor=3
bayWidht=6
storyHeight=3
pf = createPortalFrame(numBay,numFloor,firstFloor,bayWidht,storyHeight)
print(f'beam dict =>{pf.beam_dict}\n column dict => {pf.column_dict}\n nodes dict => {pf.node_dict} ')

#==========================================================================================================================================================================
#Create Definitions
"""interiorColumn_H = 70 * Unit.cm 
interiorColumn_B = 30 * Unit.cm
exteriorColumn_H = 30 *Unit.cm
exteriorColumn_B = 70 *Unit.cm
beam_H = 50*Unit.cm
beam_B = 30*Unit.cm"""

conc = material.Concrete(name='C25/30',fck=25*Unit.MPa)
steel = material.Steel(name='S420')
wBeam = 6.5 * Unit.kN/Unit.m  # dist. load on beam


Col7030 = RecSection(Id=1,name='C7030',b=700,h=300,cover=25,k=0.7,numrebars=[3,4,3],dia_rebars=[14,14,14])
Col7030.set_reinforcement_conc(plot=True,s=100)
Col3070 = RecSection(Id=2,name='C7030',b=300,h=700,cover=25,k=0.7,numrebars=[3,4,3],dia_rebars=[14,14,14])
Col3070.set_reinforcement_conc(plot=True,s=500)



"""real_rebar_area = [462,615,769,923]
num_rebar_beam  = [round(area/Col3070.dia_rebars[3],0) for area in real_rebar_area ]
print(num_rebar_beam)
Beam1FL_I = RecSection(Id=3,name='Beam1FL',b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)
Beam1FL_J = RecSection(Id=3,name='Beam1FL',b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)
Beam1M_I = RecSection(Id=4,name='Beam1M' ,b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)
Beam1M_J = RecSection(Id=4,name='Beam1M' ,b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)

Beam5F_I = RecSection(Id=5,name='Beam5F' ,b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)
Beam5F_J = RecSection(Id=5,name='Beam5F' ,b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)
Beam5L_I = RecSection(Id=6,name='Beam5L' ,b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)
Beam5L_J = RecSection(Id=6,name='Beam5L' ,b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)
Beam5M_I = RecSection(Id=7,name='Beam5M' ,b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)
Beam5M_J = RecSection(Id=7,name='Beam5M' ,b=300,h=500,cover=25,k=0.35,numrebars=[3,4,3],dia_rebars=[14,14,14]).set_reinforcement_conc(plot=True)"""



#"""kBeam = 0.35
#kCol =0.7
#
#Area_Int_Col = interiorColumn_B*interiorColumn_H
#I_Int_Col    = kCol * 1/12 * interiorColumn_B * interiorColumn_H**3
#
#Area_Ext_Col = exteriorColumn_B*exteriorColumn_H
#I_Ext_Col    = kCol * 1/12 * exteriorColumn_B * exteriorColumn_H**3
#
#Area_Beam    = beam_B*beam_H
#I_Beam       = kBeam * 1/12 * beam_B * beam_H**3
#
#bardiameterbeam        = 14*Unit.mm
#bardiametercol         = 14*Unit.mm
#cover                  = 25 * Unit.mm
#barareabeam            = 3.14*bardiameterbeam**2/4
#barareacol             = 3.14*bardiametercol**2/4
#   
#ColNumBarsTop          = 3
#ColNumBarsBot          = 3
#ColNumBarsInterior     = 4"""
#
#op.modelbuild()
#
#"""B3050_mass = beam_B*beam_H*concdensity*bayWidht/Unit.g
#C3070_mass = interiorColumn_B*interiorColumn_H*concdensity*storyHeight/Unit.g
#C7030_mass = interiorColumn_B*interiorColumn_H*concdensity*storyHeight/Unit.g"""
#
##==========================================================================================================================================================================
##Create Nodes
#for node in pf.node_dict.keys():
#    if node in pf.floorNodes[0]:
#        ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',C3070_mass/2, C3070_mass/2, 0.0)
#    
#    if node in pf.floorNodes[1]:
#        if node in [pf.floorNodes[1][0],pf.floorNodes[1][-1]]:
#            ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',(C3070_mass+B3050_mass+C3070_mass)/2, (C3070_mass+B3050_mass+C3070_mass)/2, 0.0)
#        else:
#            ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',C3070_mass+B3050_mass, C3070_mass+B3050_mass, 0.0)
#            
#    if node in pf.floorNodes[2]:
#        if node in [pf.floorNodes[2][0],pf.floorNodes[2][-1]]:
#            ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',(C3070_mass+B3050_mass+C3070_mass)/2, (C3070_mass+B3050_mass+C3070_mass)/2, 0.0)
#        else:
#            ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',(C3070_mass+C3070_mass)/2+B3050_mass,(C3070_mass+C3070_mass)/2+B3050_mass, 0.0)
#    
#    if node in pf.floorNodes[3]:
#        if node in [pf.floorNodes[3][0],pf.floorNodes[3][-1]]:
#            ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',(C3070_mass+B3050_mass+C3070_mass)/2,(C3070_mass+B3050_mass+C3070_mass)/2, 0.0)
#        else:
#            ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',C3070_mass+B3050_mass,C3070_mass+B3050_mass, 0.0)
#    
#    if node in pf.floorNodes[4]:
#        if node in [pf.floorNodes[4][0],pf.floorNodes[4][-1]]:
#            ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',(C3070_mass+B3050_mass+C3070_mass)/2,(C3070_mass+B3050_mass+C3070_mass)/2, 0.0)
#        else:
#            ops.node(node,pf.node_dict[node][0],pf.node_dict[node][1],'-mass',C3070_mass+B3050_mass,C3070_mass+B3050_mass, 0.0)
#    
#for i in pf.floorNodes[0]:
#    ops.fix(i, 1, 1, 1)
#
##==========================================================================================================================================================================
##Create Transformation
#ColTransfTag = 1
#BeamTransfTag = 2
#ops.geomTransf('Linear', BeamTransfTag) #Beam Tranformation Tag
#ops.geomTransf('PDelta', ColTransfTag)  #Column Tranformation Tag
#
##==========================================================================================================================================================================
##Create Integration
#ElBeamSec_Elastic  = 30       # Tag for elastic sbeam section 
#ElColSec_Elastic1  = 31       # Tag for elastic col section 
#ElColSec_Elastic2  = 32       # Tag for elastic col section 
#ops.section('Elastic', ElBeamSec_Elastic , Ec, Area_Ext_Col, I_Ext_Col) #[P,Mz] davranışını sonuç olarak veriyor. B60/30
#ops.section('Elastic', ElColSec_Elastic1 , Ec, Area_Ext_Col, I_Ext_Col) #[P,Mz] davranışını sonuç olarak veriyor. B60/30
#ops.section('Elastic', ElColSec_Elastic1 , Ec, Area_Ext_Col, I_Ext_Col) #[P,Mz] davranışını sonuç olarak veriyor. B60/30
#
## Moment - Curvature lardan alınan akma momentleri ve akma birim şekildeğiştirmeleri girilmiştir
#BeamMatTagFlex1 = 1333
#BeamMatTagFlex2 = 1334
#BeamMatTagFlex3 = 1335
#BeamMatTagFlex4 = 1336
#BeamMatTagFlex5 = 1337
#
#MyBeam     = 40         		                                                # yield moment kNm
#PhiYCol   = 1.02*10**-5			                                                # yield curvature
#EIBeamCrack= MyBeam/PhiYCol	                                                    # cracked section inertia
#b         = 0.01                                                                # strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent)
#ops.uniaxialMaterial ('Steel01' ,BeamMatTagFlex1 ,MyBeam ,EIBeamCrack ,b) 		# bilinear behavior for flexure
#