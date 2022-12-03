

Rebars = {
    6  : 3.14*0.6 **2/4,
    8  : 3.14*0.8 **2/4,
    10 : 3.14*0.10**2/4,
    12 : 3.14*0.12**2/4,
    14 : 3.14*0.14**2/4,
    16 : 3.14*0.16**2/4,
    18 : 3.14*0.18**2/4,
    20 : 3.14*0.20**2/4,
    22 : 3.14*0.22**2/4,
    24 : 3.14*0.24**2/4,
    25 : 3.14*0.25**2/4,
    26 : 3.14*0.26**2/4,
    28 : 3.14*0.28**2/4,
    30 : 3.14*0.30**2/4,
    32 : 3.14*0.32**2/4,
    34 : 3.14*0.34**2/4,
    36 : 3.14*0.36**2/4,
    38 : 3.14*0.38**2/4,
    40 : 3.14*0.40**2/4,
}

def CalculateRebarArea(diameter:int) -> float:
    """
        Calculate rebar area unit cm2

        diameter : Rebar size unit mm

        Output
            Rebar Area : Area of rebar unit cm2
    
    """
    return 3.14*(diameter/10)**2/4

def CalculateRebarNumbers(NeedRebarArea : float,UseRebarDiameter : int) -> int:
    """
        Calculate rebar number from given needed rebar area. if rebar number equal 1 function fix it 2. Because we can not construct the frame 1 rebar. İts montage rebars
        
        NeedRebarArea    : Need rebar area, Units cm2
        UseRebarDiameter : Use rebar diameter from needed rebar area, Units mm

        Output
            NumberRebar : Number of rebar 

    """
    RebarArea = CalculateRebarArea(diameter=UseRebarDiameter)
    needRebars = NeedRebarArea/RebarArea
    rounded = round(needRebars)
    if needRebars-rounded == 0:
        NumberRebar = 0
    if needRebars-rounded < 0 :
        NumberRebar = needRebars + (-1 * (needRebars-rounded))
    if needRebars-rounded > 0:
        NumberRebar = round(needRebars + 1-(round(needRebars-rounded)))
    
    if NumberRebar == 1 :
        NumberRebar = 2
    

    return int(NumberRebar)

"""if __name__ == '__main__':
    for rebar in Rebars.keys():
        for area in [i for i in range(1,26)]:
            numRebar = CalculateRebarNumbers(area, rebar)
            if numRebar == 1 :
                print("fucked up")
                break
            print(numRebar)"""


    

