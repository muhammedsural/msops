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
    Mass : Optional[float] = field(default_factory=float)
    
    def __post_init__(self) -> None:
        if self.Coords is None or len(self.Coords) < 2:
            raise CreateNodeError(self.node.Id,self.node.Coords)
    
    def __repr__(self) -> str:
        return f'NodeId = {self.Id}, NodeCoords = {self.Coords} '
    
    def fully_fix(self) -> None: #"UX UY UZ RX RY RZ"
        self.BoundryCondition = [1,1,1]
    
    def add_mass(self,mass):
        self.Mass += mass

@dataclass
class NodeDatas:
    Nodes : list[Node] = field(default_factory=list)
    
    def add_node(self, node: Node) -> None:
        """Add an Node to the list of Nodes."""
        self.Nodes.append(node)

    def find_nodes(self, nodeIdList : list) -> list[Node]:
        """Find all frames with a particular role in the employee list"""
        return [node for node in self.Nodes if node.Id in nodeIdList]

"""def main() -> None:
    node = Node(Id=1,Coords=[0.,1.],BoundryCondition=[0,0,0])
    if node.Coords is None or len(node.Coords) < 2:
        raise CreateNodeError(node.Id,node.Coords)
    print(node)
    
    node.fully_fix()
    print(node)
    
if __name__ == "__main__":
    main()"""
    
    