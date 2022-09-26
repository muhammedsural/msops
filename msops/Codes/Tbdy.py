import pandas as pd
import math as mt
import openseespy.opensees as ops

class FundemantelParameters:
    def __init__(self) -> None:
        pass
    def __init__(self,columndict,nodalOutput) -> None:
        #self.columndict   = columndict
        #self.nodalOutput  = nodalOutput
        self.drift        = self.StoryDrift(columndict,nodalOutput)

    def StoryDrift(self,columndict,nodalOutput):
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
            storyDrift[colId] = {"xDrift":[],"yDrift":[],"xDriftmax":0,"yDriftmax":0,"Length":0,"xStoryDrift":0,"xStoryDriftCheck":"","yStoryDrift":0,"yStoryDriftCheck":""}
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

class Performance:
    def __init__(self) -> None:
        pass
    def __init__(self,eps_su,
                 h,bw,s,f_sy,f_co,pas_payı,etriye_çapı,boyuna_donatı_çapı,numBarsTop,numBarsBot,gövde_donatı_adeti,x_koladeti,y_koladeti,
                 ultimate_curvature,yield_curvature,Lp,Ls) -> None:
        self.eps_su = eps_su
        self.h,self.bw,self.s,self.f_sy,self.f_co,self.pas_payı,self.etriye_çapı,self.boyuna_donatı_çapı,self.numBarsTop,self.numBarsBot,self.gövde_donatı_adeti,self.x_koladeti,self.y_koladeti = h,bw,s,f_sy,f_co,pas_payı,etriye_çapı,boyuna_donatı_çapı,numBarsTop,numBarsBot,gövde_donatı_adeti,x_koladeti,y_koladeti
        self.ultimate_curvature,self.yield_curvature,self.Lp,self.Ls = ultimate_curvature,yield_curvature,Lp,Ls
        
        self.steelstrain_Perform_Level = self.Steelstrain_Perform_Level(eps_su=self.eps_su)
        self.concstrain_Perform_Level  = self.Concstrain_Perform_Level(self.h,self.bw,self.s,self.f_sy,self.f_co,self.pas_payı,self.etriye_çapı,self.boyuna_donatı_çapı,self.numBarsTop,self.numBarsBot,self.gövde_donatı_adeti,self.x_koladeti,self.y_koladeti)
        self.rotation_Perform_Level    = self.Rotation_Perform_Level(self.ultimate_curvature,self.yield_curvature,self.Lp,self.Ls,db=self.boyuna_donatı_çapı)
        
        
        pass
    
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
        h_0 = h-(pas_payı+etriye_çapı/2)*2 #core_y
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

    def Rotation_Perform_Level(self,ultimate_curvature,yield_curvature,Lp,Ls,db):
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
        return rotation_performs

    def Frame_Perform_Level(self,steelStrainLevel,concStrainLevel,rotationLevel):
        frame_Level = pd.DataFrame(columns=["EleId","eps_CP","eps_LS","eps_IO","epc_CP","epc_LS","epc_IO","rot_CP","rot_LS","rot_IO"],index=[i for i in range(1,len(ops.getEleTags())+1)])
        for steel,(conc,rot) in zip(steelStrainLevel,zip(concStrainLevel,rotationLevel)):
            print(steel,conc,rot)
        pass
    
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
