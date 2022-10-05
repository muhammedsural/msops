import numpy as np
import pandas as pd

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

