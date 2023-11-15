from abc import ABC
from dataclasses import dataclass, field
from enum import Enum,auto
import math
from msops.Elements.BeamIntegration import Integration
from msops.Units.Unit import Unit 
from msops.Elements.Node import Node

@dataclass
class FrameType(Enum):
    Column = auto()
    Beam   = auto()
    Wall   = auto()

@dataclass
class ElasticBeamColumn(ABC):
    """
    Abstract class for frame elements
        Attributes:
            Id
            EleNodes
            TransfTag
            Integration
            MaxIter
            Tolerance
            Mass
            Length
    """
    Id          : int
    EleNodes    : list[Node] = field(default_factory=list)
    TransfTag   : str = field(default_factory=str)
    Integration : Integration = field(default_factory=Integration)
    MaxIter     : int = field(default_factory=int)
    Tolerance   : int = field(default_factory=int)
    Mass        : float = field(default_factory=float)
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
        self.Mass = round(self.Integration.Args[0].area * 24.99 * self.Length / Unit.g,3)

    def addNodalMass(self):
        nodalMass = self.Mass/2
        self.EleNodes[0].add_mass(mass=nodalMass)
        self.EleNodes[1].add_mass(mass=nodalMass)

@dataclass
class Column(ElasticBeamColumn):
    """ ElasticBeamcolumn for column elements"""
    frameType : FrameType = FrameType.Column
    #iSection  : RecSection = field(default_factory=RecSection)
    #jSection  : RecSection = field(default_factory=RecSection)
    
    def __post_init__(self):
        self.frame_length()
        self.calc_mass()
        self.addNodalMass()
    
    def __repr__(self) -> str:
        return f'EleId = {self.Id}, EleNodes = {self.EleNodes}, transfTag = {self.TransfTag}, integrationTag = {self.Integration} length = {self.Length} mass = {self.Mass} '

@dataclass
class Beam(ElasticBeamColumn):
    """ ElasticBeamcolumn for column elements"""
    frameType : FrameType = FrameType.Beam
    
    def __post_init__(self):
        self.frame_length()
        self.calc_mass()
        self.addNodalMass()
    
    def __repr__(self) -> str:
        return f'EleId = {self.Id}, EleNodes = {self.EleNodes}, transfTag = {self.TransfTag}, integrationTag = {self.Integration} length = {self.Length} mass = {self.Mass} '

@dataclass
class FrameDatas:
    Frames : list[ElasticBeamColumn] = field(default_factory=list)
    
    def add_frame(self, frame: ElasticBeamColumn) -> None:
        """Add an ElasticBeamColumn to the list of Frames."""
        self.Frames.append(frame)

    def find_frame(self, frameType: FrameType) -> list[ElasticBeamColumn]:
        """Find all frames with a particular role in the employee list"""
        return [frame for frame in self.Frames if frame.frameType is frameType]
    
    


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
    section = RecSection(Id=1,name='C4040',b=40,h=40,cover=2.5,k=1,numrebars=[3,2,3],dia_rebars=[1,1],fiberData=None)
    
    Intg  = Integration(Id=1,IntegrationType=IntegrationType.HingeMidpoint,Args=[section,0.5,section,0.5,section])
    print(Intg)
    print(Intg.Args[0].area)
    
    #,Mass=0.,Length=0.
    frame1 = Column(Id=1,EleNodes=[node1,node2],TransfTag="PDelta",Integration=Intg,MaxIter=300,Tolerance=0.01)
    print(frame1)
    frame2 = Beam(Id=2,EleNodes=[node1,node2],TransfTag="PDelta",Integration=Intg,MaxIter=300,Tolerance=0.01)
    print(frame2)
    
if __name__ == "__main__":
    main()"""

