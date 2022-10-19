from dataclasses import dataclass, field
from enum import Enum, auto


class IntegrationType(Enum):
    """
        Integration type in openseespy. distributed plasticy [1-10] and plastic hinge [11-15]
    """
    Lobatto           = auto()
    Legendre          = auto()
    NewtonCotes       = auto()
    Radau             = auto()
    Trapezoidal       = auto()
    CompositeSimpson  = auto()
    UserDefined       = auto()
    FixedLocation     = auto()
    LowOrder          = auto()
    MidDistance       = auto()
    
    UserHinge         = auto()
    HingeMidpoint     = auto()
    HingeRadau        = auto()
    HingeRadauTwo     = auto()
    HingeEndpoint     = auto()

@dataclass
class Integration:
    """
        BeamIntegration object for  openseespy beamintegration commands
        
        Attributes :
            Id : int
            IntegrationType : IntegrationType 
            Args : list  -> for example if you are create HingeMidpoint integration type args must include [sectionI,lpl,sectionj,lpj,sectionElastic]
    """
    Id : int
    IntegrationType : IntegrationType 
    Args : list = field(default_factory=list)
    
    def __repr__(self) -> str:
        return f'IntegrationId : {self.Id}'