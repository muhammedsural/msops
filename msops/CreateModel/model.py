from dataclasses import dataclass,field
from unicodedata import name
import numpy as np
import pandas as pd
from msops.Elements.BeamIntegration import Integration, IntegrationType
from msops.Elements.Frame import Column, FrameDatas
from msops.Elements.Node import NodeDatas,Node
from msops.Material.material import Concrete,Steel
from msops.Units.Unit import Unit
from msops.Section.sectionprop import RecSection,SectionDatas

class createPortalFrame():
    """
    numBay,
    numFloor,
    firstFloor,
    bayWidht,
    storyHeight
    """
    
    def __init__(self,numBay,numFloor,firstFloor,bayWidht,storyHeight):
        self.numBay     = numBay
        self.numFloor   = numFloor
        self.bayWidht   = bayWidht
        self.storyHeight= storyHeight
        self.firstFloor = firstFloor
        self.eleTag = 1

        self.axis_h,self.axis_v        = self.createaxis()
        
        self.floorNodes,self.node_dict = self.createnodes()
        self.column_dict               = self.createcolumns()
        self.beam_dict                 = self.createbeams()
        #self.sectiondict               = zip(self.beam_dict,self.column_dict)
        #self.Section_dict = self.CalcSectionLength()
        self.floorFrames = self.Floor_Frames()

    def createaxis(self):
        axis_h = np.linspace( 0 , (self.numBay * self.bayWidht) ,self.numBay+1)
        #axis_v = np.linspace( 0 , (self.numFloor * self.storyHeight) ,self.numFloor+1)
        axis_v = [0,self.firstFloor]
        for i in range(2,self.numFloor+1):
            axis_v.append((axis_v[1]-self.storyHeight)+i*self.storyHeight)
        
        return axis_h,axis_v
    
    def createnodes(self):
        node_no = 1
        node_dict = {} # node_dict={node no:(x,y)}
        floorNo = 0
        floorNodes = {} #master node için floorNodes={floorno : [all nodes in floor]}
        nodes = []
        for i in self.axis_v:
            floorNodes[floorNo] = [] #floor nodelar oluşturuldu
            for j in self.axis_h:
                nodes.append(node_no) #node no lar yerleştirildi
                
                # OPENSEES COMMAND
                #op.node( node_no , float(j) , float(i))

                # hangi düğüm noktaları aynı kattadır tanımlanır
                node_dict[node_no] = [float(j) , float(i)]
                floorNodes[floorNo].append(node_no)
                node_no +=1
            floorNo +=1
        return floorNodes,node_dict

    def createcolumns(self):
        
        column_dict = {}
        
        for j in range(0 , self.numBay + 1):
            end1 = j + 1
            end2 = end1 + self.numBay + 1
            for i in range( 0 ,self.numFloor):
                column_dict[self.eleTag] = [end1, end2]
                end1 =end2
                end2 += self.numBay + 1
                self.eleTag += 1
        
        colLength = 0
        for key in column_dict.keys():
            if key != "secprop":
                collenx = self.node_dict[column_dict[key][0]][0] - self.node_dict[column_dict[key][1]][0]
                colleny = self.node_dict[column_dict[key][0]][1] - self.node_dict[column_dict[key][1]][1]
                colLength = (collenx**2 + colleny**2)**0.5
                column_dict[key].append(colLength)
        return column_dict
    
    def createbeams(self):
        beam_dict = {}
        for j in range( 1 , self.numFloor +1):
            end1 = (self.numBay + 1) * j + 1
            end2 = end1+1
            for i in range(0, self.numBay):
                #ElasticBeamColum( eleTag, end1 , end2 , secType ,E*1000000, 1, M , massType)
                beam_dict[self.eleTag] = [end1, end2]
                end1 = end2
                end2 = end1 + 1
                self.eleTag += 1
        
        beamLength = 0
        for key in beam_dict.keys():
            if key != "secprop":
                beamlenx = self.node_dict[beam_dict[key][0]][0] - self.node_dict[beam_dict[key][1]][0]
                beamleny = self.node_dict[beam_dict[key][0]][1] - self.node_dict[beam_dict[key][1]][1]
                beamLength = (beamlenx**2 + beamleny**2)**0.5
                beam_dict[key].append(beamLength) 
                
        return beam_dict

    def CalcSectionLength(self):
        SectionLength={}
        sectlength = 0
        for key in self.sectiondict.keys():
            if key != "secprop":
                lenx = self.node_dict[self.sectiondict[key][0]][0] - self.node_dict[self.sectiondict[key][1]][0]
                leny = self.node_dict[self.sectiondict[key][0]][1] - self.node_dict[self.sectiondict[key][1]][1]
                SectionLength[key] = (lenx**2 + leny**2)**0.5
        return SectionLength

    def Floor_Frames(self):
        
        floorframes = {"EleId":[],"EleType":[],"Floor":[]}
        for colId,nodes in self.column_dict.items():    
            for floor in self.floorNodes.keys():
                if floor != max(self.floorNodes.keys()):
                    if nodes[0] in self.floorNodes[floor] and nodes[1] in self.floorNodes[floor+1]:
                        floorframes["EleId"].append(colId)
                        floorframes["EleType"].append("Column")
                        floorframes["Floor"].append(floor+1)

        for beamId,nodes in self.beam_dict.items():
            for floor in self.floorNodes.keys():
                if floor != max(self.floorNodes.keys()):
                    if nodes[0] in self.floorNodes[floor] and nodes[1] in self.floorNodes[floor]:
                        floorframes["EleId"].append(beamId)
                        floorframes["EleType"].append("Beam")
                        floorframes["Floor"].append(floor)
                else:
                    if nodes[0] in self.floorNodes[floor] and nodes[1] in self.floorNodes[floor]:
                        floorframes["EleId"].append(beamId)
                        floorframes["EleType"].append("Beam")
                        floorframes["Floor"].append(floor)
        floorFrames = pd.DataFrame(floorframes,index = floorframes["EleId"])
        del floorframes
        return floorFrames

@dataclass
class createPortalFrame2():
    numBay              : int
    numFloor            : int 
    firstFloor          : int
    bayWidht            : int
    storyHeight         : int
    concrete_material   : Concrete      = Concrete(name='C25/30',fck=30*Unit.MPa)
    rebar_material      : Steel         = Steel(name='S420')
    dbAxis              : list          = field(default_factory=list)
    dbFloorNodes        : dict          = field(default_factory=dict)
    dbNodes             : NodeDatas     = NodeDatas()
    dbSections          : SectionDatas  = SectionDatas()
    dbFrames            : FrameDatas    = FrameDatas()
    
    def __post_init__(self) -> None:
        pass
    
    def create_axis(self) -> None:
        """Auto create axis vertical and horizontal"""
        
        axis_h = np.linspace( 0 , (self.numBay * self.bayWidht) ,self.numBay+1)
        #axis_v = np.linspace( 0 , (self.numFloor * self.storyHeight) ,self.numFloor+1)
        axis_v = [0,self.firstFloor]
        for i in range(2,self.numFloor+1):
            axis_v.append((axis_v[1]-self.storyHeight)+i*self.storyHeight)
        self.dbAxis = [axis_h,axis_v]
        
    def create_nodes(self) -> None:
        """Create automatic portal frame nodes and define floor nodes"""
        
        node_no = 1
        floorNo = 0
        for i in self.dbAxis[1]:
            self.dbFloorNodes[floorNo] = []
            for j in self.dbAxis[0]:
                node = Node(Id=node_no, Coords=[float(j),0.,float(i)])
                self.dbNodes.add_node(node)
                self.dbFloorNodes[floorNo].append(node.Id)
                node_no +=1
            floorNo +=1
    
    def set_concrete_material(self,**kwargs) -> None:
        """
            Set Concrete material
            
            Attributes:
                name : Concrete material name
                fck  : Concrete maximum stress Unit MPa
        """
        self.concrete_material = Concrete(**kwargs)
    
    def set_rebar_material(self,**kwargs) -> None:
        """
            Set Concrete material
            
            Attributes:
                name    : str   = 'B500C'          rebar name
                density : float = 80.              rebar density
                f_sy    : float = None             rebar yield stress
                eps_sy  : float = None             rebar yield strain
                eps_sh  : float = None             rebar
                eps_su  : float = None             rebar ultimate strain
                Kres    : float = None             
                Es      : float = 2*10**5*un.MPa   young modules of rebar
                //Name: str= 'B500C',Density : float = 80.,fsy: float = None,epssy: float = None,epssh: float = None,epssu: float = None,Kres: float = None,Es: float = 2*10**5*Unit.MPa
        
        """
        self.rebar_material = Steel(**kwargs)
    
    def create_section(self,**kwargs) -> None:
        """
            Create section and added dbsection
            
            Attributes :
                name       Section name
                b          Section width
                h          Section height
                coreConc   Section core concrete
                coverConc  Section cover concrete
                matRebars  Section rebar material
                numrebars  Section rebars count    Donatı adetleri [üst başlık, gövde , alt başlık]
                dia_rebars Section rebars diameter Donatı çapı [üst başlık, gövde , alt başlık] => [üst başlık, gövde , alt başlık,total donatı alanı,% pursantaj]
                fiberData  Section fiber coords values
        """
        self.dbSections.add_section(RecSection(**kwargs))
        
    def create_column(self) -> None:
        """Create column element (Elasticbeamcolumn -force based)"""
        frameNo = 1
        ##Create Section
        #coreConc=None,coverConc=None,matRebars=None
        self.create_section(Id = 1,name='C30/70',b=0.3,h=0.7,numrebars=[3,8,3],dia_rebars=[14,14,14])
        self.create_section(Id=2 ,name='C70/30',b=0.7,h=0.3,numrebars=[3,8,3],dia_rebars=[14,14,14])
        #Material
            #Core Concrete
            #Cover Concrete
            #Rebar
        
        #BeamColumnIntegration
        #Intg  = Integration(Id=1,IntegrationType=IntegrationType.HingeMidpoint,Args=[section1, 0.5*section1.h, section1, 0.5*section1.h, section1])
        #Intg  = Integration(Id=2,IntegrationType=IntegrationType.HingeMidpoint,Args=[section2, 0.5*section2.h, section2, 0.5*section2.h, section2])
        
        for floor in self.dbFloorNodes.keys():
            if floor != max(self.dbFloorNodes):
                for iNode,jNode in zip(self.dbFloorNodes[floor],self.dbFloorNodes[floor+1]):
                    #Nodes
                    nodes = self.dbNodes.find_nodes([iNode,jNode])
                    ##Column
                    #Column(Id=frameNo,EleNodes=nodes,TransfTag=1,Integration=colIntg)
    
    def create_beam(self) -> None:
        pass
    
    def get_total_mass_buildins(self) -> float:
        pass
    
    """    
    def createFrames(self):
        frameNo = 1        
        framelist = []
        for floor in self.dbFloorNodes.keys():
            if floor != max(self.dbFloorNodes):
                for inode,jnode in zip(self.dbFloorNodes[floor],self.dbFloorNodes[floor+1]):
                    
                    frame = opsFrame(id=frameNo,FrameType='forceBeamColumn',
                                     iNode=inode,jNode=jnode,
                                     isection= self.dbSection[1], jsection=self.dbSection[1])
                    framelist.append(frame)
                    frameNo += 1
                
                if floor != 0:
                    for index,node in enumerate(self.dbFloorNodes[floor]):
                        if index != len(self.dbFloorNodes[floor])-1:
                            inode = node
                            jnode = self.dbFloorNodes[floor][index+1]
                            frame = opsFrame(id=frameNo,FrameType='forceBeamColumn',
                                            iNode=inode,jNode=jnode,
                                            isection= self.dbSection[1], jsection=self.dbSection[1])
                            framelist.append(frame)
                            frameNo += 1
            if floor == max(self.dbFloorNodes) :
                for index,node in enumerate(self.dbFloorNodes[floor]):
                        if index != len(self.dbFloorNodes[floor])-1:
                            inode = node
                            jnode = self.dbFloorNodes[floor][index+1]
                            frame = opsFrame(id=frameNo,FrameType='forceBeamColumn',
                                            iNode=inode,jNode=jnode,
                                            isection= self.dbSection[1], jsection=self.dbSection[1])
                            framelist.append(frame)
                            frameNo += 1
        self.dbFrame = framelist
        pass
    
    def __post_init__(self):
        self.createAxis()
        self.createNodes()
        self.createSection()
        self.createFrames()
        """

def main():
    model = createPortalFrame2(numBay=2,numFloor=2,firstFloor=4,bayWidht=3,storyHeight=3)
    model.create_axis()
    model.create_nodes()
    model.set_concrete_material(name = 'C30',fck = 30*Unit.MPa)
    model.set_rebar_material(name='S420')
    #print(model.dbSections.find_sections(SectionIdList=[1,2]))
    model.create_column()
    """Kolon yaratırken oluşturduğu betonarme kesitlerde eğer kabuk ve çekirdek için bir opsmaterial nesnesi yaratmazsam default değer atar. Bunun test grafiğinide aşağıdaki kodlar ile elde ederim
    bir nesne ataması yaparsam ona uygun olarak test çizimi gelir.
    
    C1 = model.dbSections.find_sections([1])
    sec1 = C1[0].coverConc
    sec1.testMaterial(materialName=sec1.get_MaterialType())
    """
if __name__ == "__main__":
    main()