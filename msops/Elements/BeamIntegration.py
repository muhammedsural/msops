from dataclasses import dataclass, field
from enum import Enum, auto


class IntegrationType(Enum):
    """
        Integration type in openseespy. distributed plasticy [1-10] and plastic hinge [11-15]
    """
    Lobatto           = "Lobatto"
    Legendre          = "Legendre"
    NewtonCotes       = "NewtonCotes"
    Radau             = "Radau"
    Trapezoidal       = "Trapezoidal     "
    CompositeSimpson  = "CompositeSimpson"
    UserDefined       = "UserDefined     "
    FixedLocation     = "FixedLocation   "
    LowOrder          = "LowOrder        "
    MidDistance       = "MidDistance     "
    
    UserHinge         = "UserHinge    "
    HingeMidpoint     = "HingeMidpoint"
    HingeRadau        = "HingeRadau   "
    HingeRadauTwo     = "HingeRadauTwo"
    HingeEndpoint     = "HingeEndpoint"

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