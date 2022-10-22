from abc import ABC
from dataclasses import dataclass, field
from enum import Enum,auto
import math
from msops.Elements.BeamIntegration import Integration
from msops.Units.Unit import Unit 
from Node import Node

@dataclass
class FrameType(Enum):
    Column = auto()
    Beam   = auto()
    Wall   = auto()

@dataclass
class ElasticBeamColumn(ABC):
    """
            eleTag (int)	tag of the element
            eleNodes (list (int))	a list of two element nodes
        EleArgs for forcebeamcolumn must be contains ;
            transfTag (int)	tag of transformation
            integrationTag (int)	tag of beamIntegration()
        and optional parameters ;
            maxIter (int)	maximum number of iterations to undertake to satisfy element compatibility (optional)
            tol (float)	tolerance for satisfaction of element compatibility (optional)
            mass (float)	element mass density (per unit length), from which a lumped-mass matrix is formed (optional)
    """
    Id          : int
    EleNodes    : list[Node] = field(default_factory=list)
    TransfTag   : str = field(default_factory=str)
    Integration : Integration = field(default_factory=Integration)
    #Section     : RecSection = field(default_factory=RecSection)
    MaxIter     : int = field(default_factory=int)
    Tolerance   : int = field(default_factory=int)
    Mass        : int = field(default_factory=int)
    Length      : float = 0
    
    def frame_length(self) -> None:
        """Calculation frame length function"""
        lengthx = self.EleNodes[0].Coords[0]  - self.EleNodes[1].Coords[0] 
        lengthy = self.EleNodes[0].Coords[1]  - self.EleNodes[1].Coords[1] 
        
        if len(self.EleNodes[0].Coords) == 3:
            lengthz = self.EleNodes[0].Coords[2]  - self.EleNodes[1].Coords[2] 
            length  = math.sqrt(lengthx**2+lengthy**2+lengthz**2)
        else:
            length  = math.sqrt(lengthx**2+lengthy**2)
            
        self.Length = length
        
    def calc_mass(self) -> None:
        """Calculation frame mass function"""
        self.Mass = self.Integration.Args[0].area * 24.99 * self.Length / Unit.g

@dataclass
class Column(ElasticBeamColumn):
    """ ElasticBeamcolumn for column elements"""
    frameType : FrameType = FrameType.Column
    
    """def __post__init__(self):
        super().frame_length()
        super().calc_mass()"""
    
    def __post_init__(self):
        self.frame_length()
        self.calc_mass()
    
    def __repr__(self) -> str:
        return f'EleId = {self.Id}, EleNodes = {self.EleNodes}, transfTag = {self.TransfTag}, integrationTag = {self.Integration} length = {self.Length} mass = {self.Mass} '

@dataclass
class Beam(ElasticBeamColumn):
    """ ElasticBeamcolumn for column elements"""
    frameType : FrameType = FrameType.Beam
    mass : float = 0.
    
    def __post_init__(self):
        self.frame_length()
        self.calc_mass()
    
    def __repr__(self) -> str:
        return f'EleId = {self.Id}, EleNodes = {self.EleNodes}, transfTag = {self.TransfTag}, integrationTag = {self.Integration} length = {self.Length} mass = {self.Mass} '

"""@dataclass
class FrameDatas:
    Frames : List[ElasticBeamColumn] = []
    
    def add_frame(self, frame: ElasticBeamColumn) -> None:
        Add an ElasticBeamColumn to the list of Frames.
        self.Frames.append(frame)

    def find_frame(self, frameType: FrameType) -> List[ElasticBeamColumn]:
        Find all frames with a particular role in the employee list
        return [frame for frame in self.Frames if frame.FrameType is FrameType]"""

"""def main() -> None:

    #dbFrame = FrameDatas()
    node1 = Node(Id=1,Coords=[0.,0.],BoundryCondition=[0,0,0])
    print(node1)
    node2 = Node(Id=2,Coords=[0.,5.],BoundryCondition=[0,0,0])
    print(node2)
    #Create Material
    conc = Concrete(name='C20',fck=20)
    
    #Create Section
    #,coreConc=None,coverConc=None,matRebars=None
    section = RecSection(name='C4040',b=40,h=40,numrebars=[3,2,3],dia_rebars=[1,1],fiberData=None)
    
    Intg  = Integration(Id=1,IntegrationType=IntegrationType.HingeMidpoint,Args=[section,0.5,section,0.5,section])
    print(Intg)
    print(Intg.Args[0].area)
    
    #,Mass=0.,Length=0.
    frame1 = Column(Id=1,EleNodes=[node1,node2],TransfTag="PDelta",Integration=Intg,MaxIter=300,Tolerance=0.01)
    print(frame1)
    frame2 = Beam(Id=1,EleNodes=[node1,node2],TransfTag="PDelta",Integration=Intg,MaxIter=300,Tolerance=0.01)
    print(frame2)
    


if __name__ == "__main__":
    main()"""
