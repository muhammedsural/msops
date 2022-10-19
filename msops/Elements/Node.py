from dataclasses import dataclass, field
from typing import Optional

class CreateNodeError(Exception):
    """
    Exception raised for errors in create node element
    
    Attributes : 
        Coords -- input node coords which caused the error
        message -- explanation of the error
    """
    def __init__(self, *args: object) -> None:
        super().__init__(*args)
    
    def __str__(self) -> str:
        return f'Node Tag :{self.args[0]} --> Wrong coordinates input, None or node not have x,y coords => Coords = {self.args[1]}'


    
@dataclass
class Node:
    Id : int
    Coords : list
    BoundryCondition : Optional[list] = field(default_factory=list)
    mass : Optional[float] = field(default_factory=float)
    
    def __post__init(self) -> None:
        if self.Coords is None or len(self.node.Coords) < 2:
            raise CreateNodeError(self.node.Id,self.node.Coords)
    
    def __repr__(self) -> str:
        return f'NodeId = {self.Id}'
    
    def fully_fix(self) -> None:
        self.BoundryCondition = [1,1,1]
    
"""def main() -> None:
    node = Node(Id=1,Coords=[0.,1.],BoundryCondition=[0,0,0])
    if node.Coords is None or len(node.Coords) < 2:
        raise CreateNodeError(node.Id,node.Coords)
    print(node)
    
    node.fully_fix()
    print(node)
    
if __name__ == "__main__":
    main()"""
    
    