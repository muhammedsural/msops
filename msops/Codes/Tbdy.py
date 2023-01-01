import pandas as pd
import math as mt
import openseespy.opensees as ops
import numpy as np
from scipy.interpolate import interp1d
from statistics import median
import matplotlib.pyplot as plt

class FundemantelParameters:

    def StoryDrift(columndict,nodalOutput):
        """
            INFO
                Zaman tanım alanında analizden sonra çalıştırılması gerekir.

            INPUT
                columndict   : Kolon bilgilerinin tutulduğu sözlük
                nodalOutput  : Düğüm noktalarının deplasman sonuçlarının tutulduğu DataFrame

            OUTPUT
                drift        : Düşey elemanların drift kontrollerinin yapıldığı DataFrame
            
        """
        storyDrift = {}
        for colId in columndict.keys():
            storyDrift[colId] = {"xDrift":[],
                                 "yDrift":[],
                                 "xDriftmax":0,
                                 "yDriftmax":0,
                                 "Length"   :0,
                                 "xStoryDrift":0,
                                 "xStoryDriftCheck":"",
                                 "yStoryDrift":0,
                                 "yStoryDriftCheck":""}
        for colId,nodes in columndict.items():
            for step in range(1,len(nodalOutput["nodeX_disp"][nodes[0]])):
                ix= nodalOutput["nodeX_disp"][nodes[0]][step]
                jx= nodalOutput["nodeX_disp"][nodes[1]][step]
                iy= nodalOutput["nodeY_disp"][nodes[0]][step]
                jy= nodalOutput["nodeY_disp"][nodes[1]][step]
                driftx = jx-ix
                drifty = jy-iy
                storyDrift[colId]["xDrift"].append(driftx)       
                storyDrift[colId]["yDrift"].append(drifty)    
                storyDrift[colId]["Length"] = nodes[2]
                
        drift = pd.DataFrame(storyDrift).T
        for i in drift.index:
            drift["xDriftmax"][i]= abs(max(drift["xDrift"][i]))
            drift["yDriftmax"][i] =abs(max(drift["yDrift"][i]))
            #print(max(drift["xDrift"][i]),max(drift["yDrift"][i]))
            drift["xStoryDrift"][i]= drift["xDriftmax"][i]/drift["Length"][i]
            drift["yStoryDrift"][i]= drift["yDriftmax"][i]/drift["Length"][i]
            
            if drift["xStoryDrift"][i]>0.008 :
                drift["xStoryDriftCheck"][i] =f"{round(drift['xStoryDrift'][i],5)} > 0.008 ×" 
            else : 
                drift["xStoryDriftCheck"][i] =f"{round(drift['xStoryDrift'][i],5)} < 0.008 ✓" 
            
            
            if drift["yStoryDrift"][i]>0.008 :
                drift["yStoryDriftCheck"][i] =f"{round(drift['yStoryDrift'][i],5)} > 0.008 ×" 
            else : 
                drift["yStoryDriftCheck"][i] =f"{round(drift['yStoryDrift'][i],5)} < 0.008 ✓"

        return drift

    def StoryDrift2(NodalDisplacement : pd.DataFrame,column_dict : dict,floorFrames : pd.DataFrame) -> pd.DataFrame:
        """
            INFO
                Zaman tanım alanında analizde kolonların göreli kat ötelemesi sonuçları getirilir.

            INPUT
                NodalDisplacement  : Düğüm noktalarının deplasman sonuçlarının tutulduğu DataFrame
                columndict         : Kolon bilgilerinin tutulduğu sözlük
                floorFrames        : Katta bulunan elemanların bilgilerinin tutulduğu dataframe
            OUTPUT
                ColumnDrift        : Düşey elemanların göreli kat öteleme sonuçları
            
        """
        ColumnProp = pd.DataFrame(column_dict).T
        ColumnProp.columns = ["iNode","jNode","Length","bw","h","cover","matcovtag","matcoretag","matsteeltag","Lpl","Location"]

        ColumnDrift = pd.DataFrame(columns=["EleTags","Floor","DeltaXMax","DeltaYMax","Length","XDriftCheck","YDriftCheck"],index=column_dict.keys())
        ColumnDrift["EleTags"] = column_dict.keys()
        ColumnDrift["Floor"] = floorFrames[(floorFrames.EleType == "Column")]["Floor"]
        ColumnDrift["Length"] = ColumnProp["Length"]
        for colId,nodes in column_dict.items():
            tempix = NodalDisplacement.query(f"Nodetags == {nodes[0]}")["NodeDispX"].reset_index() #Series
            tempiy = NodalDisplacement.query(f"Nodetags == {nodes[0]}")["NodeDispY"].reset_index() #Series
            tempjx = NodalDisplacement.query(f"Nodetags == {nodes[1]}")["NodeDispX"].reset_index() #Series
            tempjy = NodalDisplacement.query(f"Nodetags == {nodes[1]}")["NodeDispY"].reset_index() #Series
            DeltaX = tempjx["NodeDispX"]-tempix["NodeDispX"]
            DeltaY = tempjy["NodeDispY"]-tempiy["NodeDispY"]
            ColumnDrift["DeltaXMax"][colId] = DeltaX.max()
            ColumnDrift["DeltaYMax"][colId] = DeltaY.max()
            ColumnDrift["XDriftCheck"] = ColumnDrift["DeltaXMax"]/ColumnDrift["Length"]
            ColumnDrift["YDriftCheck"] = ColumnDrift["DeltaYMax"]/ColumnDrift["Length"]
        del tempix,tempiy,tempjx,tempjy,DeltaX,DeltaY,ColumnProp

        return ColumnDrift
class Performance:
        
    def Steelstrain_Perform_Level(self,eps_su):
        """
        INPUT
            eps_su : Ultimate strain value 
        OUTPUT
            eps_perf = [eps_sgö,eps_skh,eps_ssh]
        """
        eps_sgö = eps_su*0.4
        eps_skh = 0.75 * eps_sgö
        eps_ssh = 0.0075
        eps_perf = [eps_sgö,eps_skh,eps_ssh]
        return eps_perf

    def Concstrain_Perform_Level(self,h,bw,s,f_sy,f_co,pas_payı,etriye_çapı,boyuna_donatı_çapı,numBarsTop,numBarsBot,gövde_donatı_adeti,x_koladeti,y_koladeti):

        """
        INPUT:
            
            h                   : Kesitin yüksekliği
            bw                  : Kesitin genişliği [mm]
            s                   : Etriye aralığı
            f_sy                : Çelik akma dayanımı
            f_co                : Kabuk beton basınç dayanımı
            pas_payı            : Beton pas payı (mm)
            etriye_çapı         : Etriye donatı çapı (mm)
            boyuna_donatı_çapı  : Boyuna donatı çapı (mm)
            numBarsTop          : Kesit başlık bölgesindeki donatı sayısı 2 başlıkta bulunan toplam adet
            gövde_donatı_adeti  : Kesit gövde bölgesindeki donatı sayısı  2 tarafta bulunan toplam adet
            x_koladeti          : x eksenini kesen sargı kol adeti
            y_koladeti          : y eksenini kesen sargı kol adeti


        OUTPUT:
            verilen kesit bilgilerine göre performans levellerinin strain değerlerinin listesini döndürür:
            perform_Level = [eps_cgö,eps_ckh,eps_csh]

        """
        
        bar_area = 3.14*boyuna_donatı_çapı**2/4
        top_bar_area = numBarsTop * bar_area
        int_bar_area = gövde_donatı_adeti * bar_area
        bot_bar_area = numBarsBot * bar_area
        A_s = top_bar_area + bot_bar_area + int_bar_area

        #ke değerinin bulunması
        b_0 = bw-(pas_payı+etriye_çapı/2)*2 #core_x
        h_0 = h-(pas_payı+etriye_çapı/2)*2  #core_y
        #birim_x = (bw-2*pas_payı-2*etriye_çapı-boyuna_donatı_çapı)/(baslık_donatı_adeti-1) #birim aralık x
        birim_x_top= (bw-2*pas_payı-2*etriye_çapı-boyuna_donatı_çapı)/(numBarsTop-1)
        birim_x_bot= (bw-2*pas_payı-2*etriye_çapı-boyuna_donatı_çapı)/(numBarsBot-1)
        birim_y =(h-2*pas_payı-2*etriye_çapı-boyuna_donatı_çapı)/(gövde_donatı_adeti+1) #birim aralık y

        #ai_x = 2*(baslık_donatı_adeti-1)*birim_x**2
        ai_x_top= (numBarsTop-1)*birim_x_top**2
        ai_x_bot= (numBarsBot-1)*birim_x_bot**2
        ai_x_total= ai_x_top+ai_x_bot
        ai_y = 2*(gövde_donatı_adeti+1)*birim_y**2

        #toplam_ai2 =ai_x+ai_y
        toplam_ai2_tot = ai_x_total+ai_y

        #a = 1-(toplam_ai2/(6*b_0*h_0))
        a = 1-(toplam_ai2_tot/(6*b_0*h_0))
        b = 1-(s/(2*b_0))
        c = 1-(s/(2*h_0))
        d = (1-(A_s/(b_0*h_0)))**-1
        k_e =round((a*b*c*d),3)
        #ro_x = A_sx/(s*b_0)
        #ro_y = A_sx/(s*h_0)
        #Hacimsel oranların bulunması
        #check kol adetleri
        x_kol_max = max(numBarsBot,numBarsTop)
        y_kol_max = gövde_donatı_adeti+2 #Çift sıra başlık donatısı göz ardı edilmiştir.

        if x_koladeti > x_kol_max:
            x_koladeti = x_kol_max
            print(f"x_koladeti {x_kol_max} olarak değiştirildi maksimum {x_kol_max} kadar atılabiliyor")
        if y_koladeti > y_kol_max:
            y_koladeti = y_kol_max
            print(f"y_koladeti {y_kol_max} olarak değiştirildi maksimum {y_kol_max} kadar atılabiliyor")

        ro_x = round((x_koladeti*3.14*etriye_çapı**2/4)/(s*b_0),5) #hacimsel oran x
        ro_y = round((y_koladeti*3.14*etriye_çapı**2/4)/(s*h_0),5) #hacimsel oran y

        alfa_se = a*b*c
        ro_sh_min = min(ro_x,ro_y)
        omega_we =alfa_se*ro_sh_min*f_sy/f_co

        eps_cgö = 0.0035+0.04*mt.sqrt(omega_we) #<= 0.018 
        eps_ckh=eps_cgö*0.75
        eps_csh = 0.0025
        perform_Level = [eps_cgö,eps_ckh,eps_csh]
        return perform_Level

    def Rotation_Perform_Level(self,ultimate_curvature,yield_curvature,Lp,Ls,db) -> pd.DataFrame:
        """
        INPUT
            ultimate_curvature  : Kesitin maksimum eğrilik değeri. Kesitin moment-curvature eğrisinin idealleştirilmesinden tespit edilebilir.
            yield_curvature     : Kesitin akma eğrilik değeri. Kesitin moment-curvature eğrisinin idealleştirilmesinden tespit edilebilir.
            Lp                  : Plastik mafsal boyu çerçeve sistemler için 0.5* etkin doğrultudaki yüz
            Ls                  : Kesme açıklığı. Kesit yüksekliğinin 2 katı alınabilir.
            db                  : Boyuna donatı çapı
        OUTPUT
            rotation_performs = [göçme_bölgesi,kontollu_hasar,sınırlı_hasar]
        """
        a = (ultimate_curvature-yield_curvature)
        b = 1-0.5*(Lp/Ls)
        c = 4.5*ultimate_curvature*db
        collapse_Rotation = (2/3)*(a*Lp*(b+c))
        kontollu_hasar    = 0.75*collapse_Rotation
        sınırlı_hasar     = 0

        rotation_performs = [collapse_Rotation,kontollu_hasar,sınırlı_hasar]
        rotation_performs_limits = pd.DataFrame(rotation_performs).T
        rotation_performs_limits.columns = ["GÖ","KH","SH"]
        return rotation_performs_limits
    
    """def Frame_Performance_Check(self,steelStrainLevel,concStrainLevel,rotationLevel):
         
            BİRİM ŞEKİLDEĞİŞTİRMELER İÇİN 
                                0 => Göçme  ;
                                1 => Göçmenin önlenmesi ;
                                2 => Kontrollü hasar ;
                                3 => Sınırlı hasar
            DÖNMELER İÇİN 
                                0 => Göçme  ;
                                1 => Göçmenin önlenmesi ;
                                2 => Kontrollü hasar ;
                                3 => Sınırlı hasar
            KESİT İÇİN 
                        0 => Göçme Bölgesi
                        1 => İleri Hasar bölgesi
                        2 => Belirgin hasar bölgesi
                        3 => Sınırlı hasar bölgesi
        
        
        performance = pd.DataFrame(columns=["TopCoverConcPerform","TopSteelPerform","TopCoreConcPerform",
                                    "BotCoreConcPerform","BotCoverConcPerform","BotSteelPerform",
                                    "iRotationPerform","jRotationPerform","SectionPerform"],index=[i for i in range(1,len(ops.getEleTags())+1)])

        Col1exterior_rotation_perform =mam.rotation_Perform_Level(ultimate_curvature=0.01,yield_curvature=0.002,Lp=Lpl1,Ls=Lpl1,db=bardiameterexteriorcol)
        print(Col1exterior_rotation_perform)


        for ele in ops.getEleTags():
            # Top cover conc performance only one fiber from given location 
            for strain in sectionStrainStress["top_cover"][ele]["strain"]:
                if strain > Col1exterior_conc_Performs_Level[0]:
                    #print("Göçme")
                    performance["TopCoverConcPerform"][ele]=["Collapse"]
                if strain < Col1exterior_conc_Performs_Level[0] and strain > Col1exterior_conc_Performs_Level[1]:
                    performance["TopCoverConcPerform"][ele]=["Before Collapse"]
                    
                    #print("Göçmenin önlenmesi")
                if strain < Col1exterior_conc_Performs_Level[1] and strain > Col1exterior_conc_Performs_Level[2]:
                    performance["TopCoverConcPerform"][ele]=["Kontrollü hasar"]
                    
                    #print("Kontrollü hasar")
                if strain < Col1exterior_conc_Performs_Level[2]:
                    performance["TopCoverConcPerform"][ele]=["Sınırlı hasar"]
            
            # Top steel performance
            for strain in sectionStrainStress["steel_top"][ele]["strain"]:
                if strain > Col1_TopSteel_performs[0]:
                    #print("Göçme")
                    performance["TopSteelPerform"][ele]=["Collapse"]
                if strain < Col1_TopSteel_performs[0] and strain > Col1_TopSteel_performs[1]:
                    performance["TopSteelPerform"][ele]=["Before Collapse"]
                    
                    #print("Göçmenin önlenmesi")
                if strain < Col1_TopSteel_performs[1] and strain > Col1_TopSteel_performs[2]:
                    performance["TopSteelPerform"][ele]=["Kontrollü hasar"]
                    
                    #print("Kontrollü hasar")
                if strain < Col1_TopSteel_performs[2]:
                    performance["TopSteelPerform"][ele]=["Sınırlı hasar"]
            
            # Top cover conc performance only one fiber from given location 
            for strain in sectionStrainStress["top_core"][ele]["strain"]:
                if strain > Col1exterior_conc_Performs_Level[0]:
                    #print("Göçme")
                    performance["TopCoreConcPerform"][ele]=["Collapse"]
                if strain < Col1exterior_conc_Performs_Level[0] and strain > Col1exterior_conc_Performs_Level[1]:
                    performance["TopCoreConcPerform"][ele]=["Before Collapse"]
                    
                    #print("Göçmenin önlenmesi")
                if strain < Col1exterior_conc_Performs_Level[1] and strain > Col1exterior_conc_Performs_Level[2]:
                    performance["TopCoreConcPerform"][ele]=["Kontrollü hasar"]
                    
                    #print("Kontrollü hasar")
                if strain < Col1exterior_conc_Performs_Level[2]:
                    performance["TopCoreConcPerform"][ele]=["Sınırlı hasar"]
            
            # Top steel performance
            for strain in sectionStrainStress["bot_core"][ele]["strain"]:
                if strain > Col1_TopSteel_performs[0]:
                    #print("Göçme")
                    performance["BotCoreConcPerform"][ele]=["Collapse"]
                if strain < Col1_TopSteel_performs[0] and strain > Col1_TopSteel_performs[1]:
                    performance["BotCoreConcPerform"][ele]=["Before Collapse"]
                    
                    #print("Göçmenin önlenmesi")
                if strain < Col1_TopSteel_performs[1] and strain > Col1_TopSteel_performs[2]:
                    performance["BotCoreConcPerform"][ele]=["Kontrollü hasar"]
                    
                    #print("Kontrollü hasar")
                if strain < Col1_TopSteel_performs[2]:
                    performance["BotCoreConcPerform"][ele]=["Sınırlı hasar"]
            
            # Top cover conc performance only one fiber from given location 
            for strain in sectionStrainStress["bot_cover"][ele]["strain"]:
                if strain > Col1exterior_conc_Performs_Level[0]:
                    #print("Göçme")
                    performance["BotCoverConcPerform"][ele]=["Collapse"]
                if strain < Col1exterior_conc_Performs_Level[0] and strain > Col1exterior_conc_Performs_Level[1]:
                    performance["BotCoverConcPerform"][ele]=["Before Collapse"]
                    
                    #print("Göçmenin önlenmesi")
                if strain < Col1exterior_conc_Performs_Level[1] and strain > Col1exterior_conc_Performs_Level[2]:
                    performance["BotCoverConcPerform"][ele]=["Kontrollü hasar"]
                    
                    #print("Kontrollü hasar")
                if strain < Col1exterior_conc_Performs_Level[2]:
                    performance["BotCoverConcPerform"][ele]=["Sınırlı hasar"]
            
            # Top steel performance
            for strain in sectionStrainStress["steel_bot"][ele]["strain"]:
                if strain > Col1_TopSteel_performs[0]:
                    #print("Göçme")
                    performance["BotSteelPerform"][ele]=["Collapse"]
                if strain < Col1_TopSteel_performs[0] and strain > Col1_TopSteel_performs[1]:
                    performance["BotSteelPerform"][ele]=["Before Collapse"]
                    
                    #print("Göçmenin önlenmesi")
                if strain < Col1_TopSteel_performs[1] and strain > Col1_TopSteel_performs[2]:
                    performance["BotSteelPerform"][ele]=["Kontrollü hasar"]
                    
                    #print("Kontrollü hasar")
                if strain < Col1_TopSteel_performs[2]:
                    performance["BotSteelPerform"][ele]=["Sınırlı hasar"]
                    
            for rotationi,rotationj in zip(sectionOutput["itotalRot"][ele],sectionOutput["jtotalRot"][ele]):
                if rotationi > Col1exterior_rotation_perform[0]:
                    #print("Göçme")
                    performance["iRotationPerform"][ele]=["Collapse"]
                if rotationi < Col1exterior_rotation_perform[0] and rotationi > Col1exterior_rotation_perform[1]:
                    performance["iRotationPerform"][ele]=["Before Collapse"]
                    #print("Göçmenin önlenmesi")
                if rotationi < Col1exterior_rotation_perform[1] and rotationi > Col1exterior_rotation_perform[2]:
                    performance["iRotationPerform"][ele]=["Life Safety"]
                    #print("Kontrollü hasar")
                if rotationi < Col1exterior_rotation_perform[2]:
                    performance["iRotationPerform"][ele]=["Immediate Occupuancy"]
                    
                if rotationj > Col1exterior_rotation_perform[0]:
                    #print("Göçme")
                    performance["jRotationPerform"][ele]=["Collapse"]
                if rotationj < Col1exterior_rotation_perform[0] and rotationj > Col1exterior_rotation_perform[1]:
                    performance["jRotationPerform"][ele]=["Before Collapse"]
                    #print("Göçmenin önlenmesi")
                if rotationj < Col1exterior_rotation_perform[1] and rotationj > Col1exterior_rotation_perform[2]:
                    performance["jRotationPerform"][ele]=["Life Safety"]
                    #print("Kontrollü hasar")
                if rotationj < Col1exterior_rotation_perform[2]:
                    performance["jRotationPerform"][ele]=["Immediate Occupuancy"]
                    #print("Sınırlı hasar")

        pass"""


    def FramePerformanceBoundries(self,column_dict : dict,beam_dict : dict ,floorFrames : pd.DataFrame,important_points_ext : dict,important_points_int : dict) -> pd.DataFrame:
        """Calculate rotation and strain performance limits for all frame elements"""
        steelstrainlimits = self.Steelstrain_Perform_Level(eps_su=0.08)

        ColumnProp = pd.DataFrame(column_dict).T
        BeamProp = pd.DataFrame(beam_dict).T
        ColumnProp.columns = ["iNode","jNode","Length","bw","h","cover","matcovtag","matcoretag","matsteeltag","Lpl","Location"]
        BeamProp.columns = ["iNode","jNode","Length","bw","h","cover","Lpl","ult_curv","yield_curv"]

        # Rotation and strain performance limits according to TBDY
        #==================================================================================
        PerfLimit = floorFrames.copy()
        PerfLimit["H"] = ColumnProp["h"]
        PerfLimit["Lpl"] = ColumnProp["Lpl"]
        PerfLimit["H"][ColumnProp.last_valid_index():PerfLimit.last_valid_index()] = BeamProp["h"]
        PerfLimit["Lpl"][ColumnProp.last_valid_index():PerfLimit.last_valid_index()] = BeamProp["Lpl"]

        #PerfLimit.H.fillna(value= 0.5,inplace=True)
        #PerfLimit.Lpl.fillna(value=Lpl,inplace=True)
        limitsPerform = pd.DataFrame(columns=[  "ConcStrainGÖ",
                                                "ConcStrainKH",
                                                "ConcStrainSH",
                                                "SteelStrainGÖ",
                                                "SteelStrainKH",
                                                "SteelStrainSH",
                                                "RotationGÖ",
                                                "RotationKH",
                                                "RotationSH"],index=ops.getEleTags())
        for eleid,eletype in zip(PerfLimit["EleId"],PerfLimit["EleType"]):

            if eletype == "Column":
                if ColumnProp["Location"][eleid] == "Exterior":
                    colext_perf =self.Rotation_Perform_Level(ultimate_curvature=0.009,yield_curvature=0.004,Lp=PerfLimit["Lpl"][eleid],Ls=2*PerfLimit["H"][eleid],db=14)
                    
                    limitsPerform.loc[eleid] = [
                                                    round(important_points_ext['performance'][0][0],4),
                                                    round(important_points_ext['performance'][1][0],4),
                                                    round(important_points_ext['performance'][2][0],4),
                                                    steelstrainlimits[0],
                                                    steelstrainlimits[1],
                                                    steelstrainlimits[2],
                                                    round(colext_perf["GÖ"][0],4),
                                                    round(colext_perf["KH"][0],4),
                                                    round(colext_perf["SH"][0],4)
                                               ]         
                else:
                    colint_perf =self.Rotation_Perform_Level(ultimate_curvature=0.009,yield_curvature=0.004,Lp=PerfLimit["Lpl"][eleid],Ls=2*PerfLimit["H"][eleid],db=14)
                    
                    limitsPerform.loc[eleid] = [
                                                    round(important_points_int['performance'][0][0],4),
                                                    round(important_points_int['performance'][1][0],4),
                                                    round(important_points_int['performance'][2][0],4),
                                                    steelstrainlimits[0],
                                                    steelstrainlimits[1],
                                                    steelstrainlimits[2],
                                                    round(colint_perf["GÖ"][0],4),
                                                    round(colint_perf["KH"][0],4),
                                                    round(colint_perf["SH"][0],4)
                                                ]
            
            #beam only rotation
            else:
                beam_perf =self.Rotation_Perform_Level(ultimate_curvature=BeamProp["ult_curv"][eleid],yield_curvature=BeamProp["yield_curv"][eleid],Lp=BeamProp["Lpl"][eleid],Ls=2*BeamProp["Lpl"][eleid],db=14)
                limitsPerform.loc[eleid] = [
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                round(beam_perf["GÖ"][0],4),
                                                round(beam_perf["KH"][0],4),
                                                round(beam_perf["SH"][0],4)
                                            ]
            
        PerfLimit = pd.concat([PerfLimit,limitsPerform],axis=1)

        return PerfLimit

    def MaxCoreFiberStrain(self,StressStrain : pd.DataFrame) -> pd.DataFrame:

        # Fiber strain max values core con and steel materials
        #==========================================================================================================
        FiberStressStrainMax = StressStrain.copy()
        FiberStressStrainMax.drop(columns=["TopCoverStress","TopCoverStrain","TopCoreStress","TopSteelStress","BotCoverStress","BotCoverStrain","BotCoreStress","BotSteelStress"], axis=1,inplace=True)
        corestrain = [max(abs(top),abs(bot)) for top,bot in zip(FiberStressStrainMax["TopCoreStrain"],FiberStressStrainMax["BotCoreStrain"])]
        steelstrain = [max(abs(top),abs(bot)) for top,bot in zip(FiberStressStrainMax["TopSteelStrain"],FiberStressStrainMax["BotSteelStrain"])]
        FiberStressStrainMax["CoreStrainMax"]  = corestrain
        FiberStressStrainMax["SteelStrainMax"] = steelstrain
        FiberStressStrainMax.drop(columns=['TopCoreStrain', 'TopSteelStrain','BotCoreStrain', 'BotSteelStrain'], axis=1,inplace=True)
        FiberStressStrainMax.columns = ["Eleid","CoreStrainMax","SteelStrainMax"]
        #==========================================================================================================
        for i in range(FiberStressStrainMax.last_valid_index()+1,max(ops.getEleTags())):
            FiberStressStrainMax.loc[i] = [0,0,0]
            
        return FiberStressStrainMax

    def FramePerformanceCheck(self,column_dict,floorFrames,Lpl,beam_H,important_points_ext,important_points_int,FiberStressStrain,MomentRotation) -> pd.DataFrame:
        """All frame performance level calculate"""
        
        steelstrainlimits = self.Steelstrain_Perform_Level(eps_su=0.08)

        ColumnProp = pd.DataFrame(column_dict).T
        ColumnProp.columns = ["iNode","jNode","Length","bw","h","cover","matcovtag","matcoretag","matsteeltag","Lpl"]

        # Rotation and strain performance limits according to TBDY
        #==================================================================================
        PerfLimit = floorFrames.copy()
        PerfLimit["H"] = ColumnProp["h"]
        PerfLimit["Lpl"] = ColumnProp["Lpl"]
        PerfLimit.H.fillna(value= 0.5,inplace=True)
        PerfLimit.Lpl.fillna(value=Lpl,inplace=True)
        limitsPerform = pd.DataFrame(columns=["ConcStrainGÖ",
                                                "ConcStrainKH",
                                                "ConcStrainSH",
                                                "SteelStrainGÖ",
                                                "SteelStrainKH",
                                                "SteelStrainSH",
                                                "RotationGÖ",
                                                "RotationKH",
                                                "RotationSH"],index=ops.getEleTags())
        for eleid,eletype in zip(PerfLimit["EleId"],PerfLimit["EleType"]):

            if eletype == "Column":

                if PerfLimit["H"][eleid] == 0.3:
                    colext_perf =self.Rotation_Perform_Level(ultimate_curvature=0.009,yield_curvature=0.004,Lp=Lpl,Ls=2*0.3,db=14)
                    
                    limitsPerform.loc[eleid] = [round(important_points_ext['performance'][0][0],4),round(important_points_ext['performance'][1][0],4),round(important_points_ext['performance'][2][0],4),steelstrainlimits[0], steelstrainlimits[1],steelstrainlimits[2], round(colext_perf["GÖ"][0],4),round(colext_perf["KH"][0],4),round(colext_perf["SH"][0],4)]
                
                else:
                    colint_perf =self.Rotation_Perform_Level(ultimate_curvature=0.009,yield_curvature=0.004,Lp=Lpl,Ls=2*0.7,db=14)
                    
                    limitsPerform.loc[eleid] = [
                                                    round(important_points_int['performance'][0][0],4),
                                                    round(important_points_int['performance'][1][0],4),
                                                    round(important_points_int['performance'][2][0],4),
                                                    steelstrainlimits[0],
                                                    steelstrainlimits[1],
                                                    steelstrainlimits[2],
                                                    round(colint_perf["GÖ"][0],4),
                                                    round(colint_perf["KH"][0],4),
                                                    round(colint_perf["SH"][0],4)
                                                ]
            #beam only rotation
            else:
                
                beam_perf =self.Rotation_Perform_Level(ultimate_curvature=0.009,yield_curvature=0.004,Lp=Lpl,Ls=2*beam_H,db=16)
                limitsPerform.loc[eleid] = [
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                0,
                                                round(beam_perf["GÖ"][0],4),
                                                round(beam_perf["KH"][0],4),
                                                round(beam_perf["SH"][0],4)
                                            ]
            
        PerfLimit = pd.concat([PerfLimit,limitsPerform],axis=1)

        # Fiber strain max values core con and steel materials
        #==========================================================================================================
        FiberStressStrainMax = FiberStressStrain.copy()
        FiberStressStrainMax.drop(columns=["TopCoverStress","TopCoverStrain","TopCoreStress","TopSteelStress","BotCoverStress","BotCoverStrain","BotCoreStress","BotSteelStress"], axis=1,inplace=True)
        corestrain = [max(abs(top),abs(bot)) for top,bot in zip(FiberStressStrainMax["TopCoreStrain"],FiberStressStrainMax["BotCoreStrain"])]
        steelstrain = [max(abs(top),abs(bot)) for top,bot in zip(FiberStressStrainMax["TopSteelStrain"],FiberStressStrainMax["BotSteelStrain"])]
        FiberStressStrainMax["CoreStrainMax"]  = corestrain
        FiberStressStrainMax["SteelStrainMax"] = steelstrain
        FiberStressStrainMax.drop(columns=['TopCoreStrain', 'TopSteelStrain','BotCoreStrain', 'BotSteelStrain'], axis=1,inplace=True)
        FiberStressStrainMax.columns = ["Eleid","CoreStrainMax","SteelStrainMax"]
        #==========================================================================================================
        FramePerfmCheck = MomentRotation.copy()
        FramePerfmCheck.drop(columns=['iMoment', 'jMoment'], axis=1,inplace=True)
        #==========================================================================================================
        for i in range(FiberStressStrainMax.last_valid_index()+1,FramePerfmCheck.last_valid_index()+1):
            FiberStressStrainMax.loc[i] = [0,0,0]
            
        FramePerfmCheck = pd.concat([FramePerfmCheck,FiberStressStrainMax],axis=1)

        # Strain performance check
        #==========================================================================================================
        # 0 -> Sınırlı Hasar
        # 1 -> Belirgi Hasar
        # 2 -> İleri Hasar
        # 3 -> Göçme Durumu
        # 4 -> Not Fiber
        core_strainlevel = []
        steel_strainslevel = []
        for index,(core_strain,steel_strains) in zip(FramePerfmCheck.index,zip(FramePerfmCheck["CoreStrainMax"],FramePerfmCheck["SteelStrainMax"])):
            if FramePerfmCheck["Eletags"][index] == FramePerfmCheck["Eleid"][index]:
                if core_strain > PerfLimit["ConcStrainGÖ"][FramePerfmCheck["Eletags"][index]]:
                    core_strainlevel.append(3)
                elif core_strain < PerfLimit["ConcStrainGÖ"][FramePerfmCheck["Eletags"][index]] and core_strain >= PerfLimit["ConcStrainKH"][FramePerfmCheck["Eletags"][index]]:
                    core_strainlevel.append(2)
                elif core_strain < PerfLimit["ConcStrainKH"][FramePerfmCheck["Eletags"][index]] and core_strain >= PerfLimit["ConcStrainSH"][FramePerfmCheck["Eletags"][index]]:
                    core_strainlevel.append(1)
                elif core_strain < PerfLimit["ConcStrainSH"][FramePerfmCheck["Eletags"][index]]:
                    core_strainlevel.append(0)
                
                if steel_strains > PerfLimit["SteelStrainGÖ"][FramePerfmCheck["Eletags"][index]]:
                    steel_strainslevel.append(3)
                elif steel_strains < PerfLimit["SteelStrainGÖ"][FramePerfmCheck["Eletags"][index]] and steel_strains >= PerfLimit["SteelStrainKH"][FramePerfmCheck["Eletags"][index]]:
                    steel_strainslevel.append(2)
                elif steel_strains < PerfLimit["SteelStrainKH"][FramePerfmCheck["Eletags"][index]] and steel_strains >= PerfLimit["SteelStrainSH"][FramePerfmCheck["Eletags"][index]]:
                    steel_strainslevel.append(1)
                elif steel_strains < PerfLimit["SteelStrainSH"][FramePerfmCheck["Eletags"][index]]:
                    steel_strainslevel.append("0")
            else:
                core_strainlevel.append(4)
                steel_strainslevel.append(4)

                    
        FramePerfmCheck["ConcretePerformLevel"] = core_strainlevel
        FramePerfmCheck["SteelPerformLevel"] = steel_strainslevel

        #rotation performance check
        #==========================================================================================================
        # 0 -> Sınırlı Hasar
        # 1 -> Belirgi Hasar
        # 2 -> İleri Hasar
        # 3 -> Göçme Durumu
        irotlevel = []
        jrotlevel = []
        for eletag,(irot,jrot) in zip(FramePerfmCheck["Eletags"],zip(FramePerfmCheck["iRotation"],FramePerfmCheck["jRotation"])):
            if irot > PerfLimit["RotationGÖ"][eletag]:
                irotlevel.append(3)
            elif irot < PerfLimit["RotationGÖ"][eletag] and irot >= PerfLimit["RotationKH"][eletag]:
                irotlevel.append(2)
            elif irot < PerfLimit["RotationKH"][eletag] and irot >= PerfLimit["RotationSH"][eletag]:
                irotlevel.append(1)
            elif irot < PerfLimit["RotationSH"][eletag]:
                irotlevel.append(0)
            
            if jrot > PerfLimit["RotationGÖ"][eletag]:
                jrotlevel.append(3)
            elif jrot < PerfLimit["RotationGÖ"][eletag] and jrot >= PerfLimit["RotationKH"][eletag]:
                jrotlevel.append(2)
            elif jrot < PerfLimit["RotationKH"][eletag] and jrot >= PerfLimit["RotationSH"][eletag]:
                jrotlevel.append(1)
            elif jrot < PerfLimit["RotationSH"][eletag]:
                jrotlevel.append(0)

        FramePerfmCheck["iRotPerformLevel"] = irotlevel
        FramePerfmCheck["jRotPerformLevel"] = jrotlevel
        FramePerfmCheck.drop(columns=['Eleid'], axis=1,inplace=True)
        print(" 0 -> Sınırlı Hasar; 1 -> Belirgi Hasar; 2 -> İleri Hasar; 3 -> Göçme Durumu; 4 -> Not Fiber")
        return FramePerfmCheck

class TargetSpectrum:

    def HorizontalElasticSpectrum(self,Ss : float, S1 : float, soil : str)-> pd.DataFrame:
        """
        Args:
        Ss: Spectral Acceleration Parameter at Short Periods
        S1: Spectral Acceleration Parameter at 1-sec
        soil: Soil Type for example "ZA" or "ZB"
        Output
            Targetspektra : TargetSpektra according to TBDY
        """

        Ss_range = [0.25 , 0.50 , 0.75, 1.00 , 1.25 , 1.50 ]
        FS_table = {"ZA": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                    "ZB": [0.9 , 0.9 , 0.9 , 0.9 , 0.9 , 0.9], 
                    "ZC": [1.3 , 1.3 , 1.2 , 1.2 , 1.2 , 1.2],
                    "ZD": [1.6 , 1.4 , 1.2 , 1.1 , 1.0 , 1.0],
                    "ZE": [2.4 , 1.7 , 1.3 , 1.1 , 0.9 , 0.8]}

        S1_range = [0.10 , 0.20 , 0.30, 0.40 , 0.50 , 0.60 ]
        F1_table = {"ZA": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                    "ZB": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                    "ZC": [1.5 , 1.5 , 1.5 , 1.5 , 1.5 , 1.4],
                    "ZD": [2.4 , 2.2 , 2.0 , 1.9 , 1.8 , 1.7],
                    "ZE": [4.2 , 3.3 , 2.8 , 2.4 , 2.2 , 2.0]}

        if Ss < Ss_range[0]:
            FS_satir = np.polyfit(Ss_range[0:2], list(FS_table[soil])[0:2], 1)
            FS_katsayisi = np.poly1d( FS_satir )
            Fs = float( format(FS_katsayisi(Ss) , '.2f') )
            SDs = Ss * Fs
        elif Ss > Ss_range[-1]:
            FS_satir = np.polyfit(Ss_range[-3:-1], list(FS_table[soil])[-3:-1], 1)
            FS_katsayisi = np.poly1d( FS_satir )
            Fs = float( format(FS_katsayisi(Ss) , '.2f') )
            SDs = Ss * Fs    
        else:
            FS_satir = interp1d(Ss_range, FS_table[soil], kind='linear')
            FS_katsayisi = FS_satir(Ss)
            Fs = round( float(FS_katsayisi) , 2) 
            SDs = Ss * Fs

        if S1 < S1_range[0] :
            F1_satir = np.polyfit(S1_range[0:2], list(F1_table[soil])[0:2], 1)
            F1_katsayisi = np.poly1d( F1_satir )
            F1 = float( format(F1_katsayisi(S1) , '.2f') )
            SD1 = S1 * F1
        elif S1 > S1_range[-1]:
            F1_satir = np.polyfit(S1_range[-3:-1], list(F1_table[soil])[-3:-1], 1)
            F1_katsayisi = np.poly1d( F1_satir )
            F1 = float( format(F1_katsayisi(S1) , '.2f') )
            SD1 = S1 * F1
        else:    
            F1_satir = interp1d(S1_range, F1_table[soil], kind='linear')
            F1_katsayisi = F1_satir(S1)
            F1 = round(float(F1_katsayisi) , 2)
            SD1 = S1 * F1
            
        TA = 0.2 * SD1 / SDs
        TB = SD1 / SDs
        TL = 6
        
        
        T_list = np.arange(0.0, TL,.005)
            
        Sa = []
        
        for i in T_list:
            
            if i <TA:
                Sa.append(round((0.4 + 0.6*(i/TA))*SDs, 4))
                
            elif i >= TA and i<=TB:
                Sa.append(round(SDs, 4))
                
            elif i>TB and i <=TL:
                Sa.append(round(SD1/i, 4))
                
            elif i>TL:
                Sa.append(round(SD1*TL/(i**2), 4))
                
        target_spec = {"T" : T_list,
                    "Sa" : Sa}

        target_spec_df = pd.DataFrame().from_dict(target_spec)
        
        return target_spec_df

    def VerticalElasticSpektrum(self):
        pass

    def HorizontalDisplacementSpectrum(self):
        pass
    
    def Calc_Ra(self,R : int, T:float, I : float, D : float, SD1:float, SDs:float):
        """Deprem yükü azaltma katsayisi"""
        TB = SD1 / SDs
        if T > TB:
            Ra = R/I
        else:
            Ra = D + ((R/I)-D)*(T/TB)
        return Ra

    def ReducedTargetSpectrum(self,TargetSpectrum : pd.DataFrame,R : int, I : float, D : float, SD1 : float, SDs : float) -> pd.DataFrame:
        """Reduced Target Spectra according to TBDY. TargetSpectrum -> HorizontalElasticSpectrum"""
        Tw = TargetSpectrum.T
        RaT = [ self.Calc_Ra(R = R, T = T, I = I, D = D, SD1 =SD1 ,SDs = SDs) for T in Tw ]
        SaR = [(Sa/Ra) for Sa,Ra in zip(TargetSpectrum["Sa"],RaT)]
        TargetSpectrum["RaT"] = RaT
        TargetSpectrum["SaR"] = SaR

        return TargetSpectrum
    
    def TimeSeriesSpectra(self,Acceleration , Time ):    
        sampling_interval = Time[1]-Time[0]
        damping_ratio = 0.05
        Sd = []
        Sv = []
        Sa = []
        
        T = np.arange(0.05, 6.0,.001)
        for i in T:
            omega = 2*np.pi/i 
            mass = 1 
            k = ((omega)**2)*mass
            c = 2*mass*omega*damping_ratio
            K = k+3*c/sampling_interval + 6*mass/(sampling_interval**2)
            a = 6*mass / sampling_interval + 3*c
            b = 3*mass + sampling_interval*c/2
            u= [0]
            v= [0]
            ac= [0] # INITIAL CONDITIONS
            for j in range(len(Acceleration)-1) :
                df = - ( Acceleration[j+1] - Acceleration[j])+ a*v[j] + b*ac[j] # delta force
                du = df / K
                dv = 3*du / sampling_interval - 3*v[j] - sampling_interval * ac[j] /2    
                dac = 6* (du - sampling_interval*v[j]) / (sampling_interval**2) - 3* ac[j]
                u.append(u[j] + du)
                v.append(v[j] + dv)
                ac.append(ac[j] + dac)
            Sd.append(max([abs(x) for x in u]))
            #Sv.append(max([abs(x) for x in v]))
            #Sa.append(max([abs(x) for x in ac]))
            Sv.append(Sd[-1]*omega)
            Sa.append(Sd[-1]*omega**2)

        # Gorselleştirmesi   
        # plt.figure(figsize=[10,5] );
        # plt.suptitle(' Response Spectra' )
        # plt.subplot(3,1,1),plt.plot(T,Sd) ; plt.ylabel('Sd (m)') ; plt.grid()
        # plt.subplot(3,1,2),plt.plot(T,Sv) ; plt.ylabel('Sv (m/s)'); plt.grid()
        # plt.subplot(3,1,3),plt.plot(T,Sa) ; plt.ylabel('Sa (m/s2)'); plt.grid()

        return T,Sa,Sv,Sd
    
    def LocationSeriesSpectra(self,T,Accelertions,Periods):
        targetSa = 0
        for index,T_series in enumerate(Periods):
            if round(T_series,3) == T:
                targetSa = Accelertions[index]
                break
        return round(targetSa,4)

    def LocationHorizontalSpectra(self,R:float,I:float,D:float,T:float,SD1:float,SDs:float) -> list:
        """Elastik spektrum ile Ra katsayısı ile azaltılmış spektrumda yapının doğal titreşim periyoduna denk gelen spektral ivme değerinin listesini döndürür.
            ilk değer elastik spektral ivmedir, ikinci değer azaltılmış spektral ivme değeridir.
        """
        TA = 0.2 * SD1 / SDs
        TB = SD1 / SDs
        TL = 6

        if T <TA:
            Sae = round((0.4 + 0.6*(T/TA))*SDs, 4)
            
        elif T >= TA and T<=TB:
           Sae=round(SDs, 4)
            
        elif T>TB and T <=TL:
            Sae=round(SD1/T, 4)
            
        elif T>TL:
            Sae=round(SD1*TL/(T**2), 4)
        
        Ra = self.Calc_Ra(R, T, I , D , SD1, SDs)

        SaR = round(Sae/Ra,4)
        coeff = [Sae,SaR]
        return coeff
    

