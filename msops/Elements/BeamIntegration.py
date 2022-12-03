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
            Args : list  -> for example if you are create HingeMidpoint integration type, args must include [sectionI_id,lpl,sectionj_id,lpj,sectionElastic_id]
    """
    Id : int
    IntegrationType : IntegrationType 
    Args : list = field(default_factory=list)
    
    def __repr__(self) -> str:
        return f'IntegrationId : {self.Id}, Integration type : {self.IntegrationType}, Integration arguments : {self.Args}'

@dataclass
class IntegrationDatas:
    Integrations : list[Integration] = field(default_factory=list)
    
    def add_integration(self, integration: Integration) -> None:
        """Add an ElasticBeamColumn to the list of Frames."""
        self.Integrations.append(integration)

    def delete_integration(self):
        """Deleted integration"""
        pass

    def update_integration(self,**kwargs):
        """Update integration"""
        pass

    def Get_integration(self,IntegrationId : int) -> None:
        """Get integration"""
        for integ in self.Integrations:
            if IntegrationId == integ.Id:
                return integ
            else:
                continue

    def Get_all_integration(self):
        """Get all defined integration"""
        pass


