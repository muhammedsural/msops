from msops.Analysis.builder import opsbuild as op
from msops.CreateModel.model import createPortalFrame
from msops.Codes.Tbdy import FundemantelParameters
from msops.Material import material,tbdymaterials
from msops.Plotter import msplotter
from msops.TimeSeries import ReadRecord
from msops.Units.Unit import Unit
from msops.Section import sectionprop as sp



numBay,numFloor,firstFloor,bayWidht,storyHeight=3,8,4,5,3
pf = createPortalFrame(2,3,3,3,4)

"""
INPUT
    name = Kesit adı 
    b = Kesit genişliği 
    h = Kesit yüksekliği 
    coreConc = Çekirdek malzemesi 
    coverConc = Kabuk malzemesi 
    matRebars = Çelik malzemesi 
    numrebars = Donatı adetleri [üst başlık, gövde , alt başlık]
    dia_rebars = Donatı çapı [üst başlık, gövde , alt başlık] => [üst başlık, gövde , alt başlık,total donatı alanı,% pursantaj] 
    fiberData =
"""
Col1numrebars = [3,2,3]
Col1diabars   = [22*Unit.mm,22*Unit.mm,22*Unit.mm]
fc = 20 * Unit.MPa
BCol1=0.4
HCol1=0.4
s=150 #mm
etriye_çapı = 10 #mm
boyuna_donatı_çapı = 22 #mm
pas_payı = 50 #mm
baslık_donatı_top = 3 #adet
baslık_donatı_bot = 3 #adet
gövde_donatı_adeti =2 #adet
x_koladeti = 3 #kesitin x eksenini kesen kol sayısı
y_koladeti = 3 #kesitin y eksenini kesen kol sayısı
unconfined,confined=tbdymaterials.tbdy_mander("B500C",fc/Unit.MPa,BCol1*10**3,HCol1*10**3,s,etriye_çapı,Col1diabars[0]*10**3,pas_payı,Col1numrebars[0],Col1numrebars[2],Col1numrebars[1],x_koladeti,y_koladeti)
msplotter.plotter.plot_mander(confined['strain_stress'][0],confined['strain_stress'][1],label="Confined model")
msplotter.plotter.plot_mander(unconfined['strain_stress'][0],unconfined['strain_stress'][1],label="Unonfined model")

#Col1core      = material.opsmaterial(3,)
#Col1 = rc(name='C4040',b=0.4,h=0.4,coreConc=None,coverConc=None,matRebars=None,numrebars=Col1numrebars,dia_rebars=Col1diabars,fiberData=None)


