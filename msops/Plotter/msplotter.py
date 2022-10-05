import os
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import opsvis as opsv
import numpy as np

class plotter:
    @staticmethod
    def top_disp_plot(nodalDispDict,time,dt):
        """
            nodalDispDict : Displacement values for all nodes. get, from run analysis function
            time          : Time values get from run analysis func.
            dt            : Time step
        """
        index = 0
        value = 0
        for key,val in enumerate(nodalDispDict[max(ops.getNodeTags())]):
            if abs(val) == max(abs(nodalDispDict[max(ops.getNodeTags())])):
                index = key 
                value = val

        fig, ax = plt.subplots(figsize=(25,10))
        fig.subplots_adjust(bottom=0.15, left=0.2)
        ax.grid()
        plt.scatter(index*dt,value, s=80, c="r", alpha=0.5)
        ax.plot(time,nodalDispDict[max(ops.getNodeTags())])

        meandispvalue=nodalDispDict[max(ops.getNodeTags())].mean()
        ax.axhline(0, color='black', lw=2)
        #ax.text(0, nodalDispDict[max(ops.getNodeTags())].mean()+meandispvalue, "Mean displacement line",color ="r", fontsize=12)
        #ax.axhline(nodalDispDict[max(ops.getNodeTags())].mean(), color='red', lw=2)
        
        
        ax.annotate(f'Max Displacement = {round(value,3)}m', xy =(index*dt, value),xytext =(index*dt+meandispvalue, value+meandispvalue),arrowprops = dict(facecolor ='green',shrink = 0.05),)
        ax.set_xlabel('Time [Sec]')
        ax.set_ylabel(f'Horizontal Displacement of node {max(ops.getNodeTags())}')
        plt.show()
    
    def plot_MomRot(itotalRot,ibasicForce,jtotalRot,jbasicForce,count=5):
        """
         Time History analizi yapıldıktan sonra elde edilen sectionDataframeindeki ilgili kolonlar input olarak verilmeli
        """
        for ele in ops.getEleTags():
            if ele <= count:
                fig , ax = plt.subplots( 1,2 , sharex = True , sharey = True  , figsize=(20,5) )
                ax[0].plot( itotalRot[ele], ibasicForce[ele] ) , ax[0].axhline(c = "k") , ax[0].axvline(c = "k") , ax[0].set_xlabel("Rotation (rad)"), ax[0].set_ylabel( "Moment (kNM)") , ax[0].legend( ["i node"]) , ax[0].set_frame_on(False) ; 
                ax[1].plot( jtotalRot[ele], jbasicForce[ele] ) , ax[1].axhline(c = "k") , ax[1].axvline(c = "k") , ax[1].set_xlabel("Rotation (rad)"), ax[1].set_ylabel( "Moment (kNM)") , ax[1].legend( ["j node"]) , ax[1].set_frame_on(False) ; plt.suptitle( f"Frame{ele} Hysteresis Graphs".upper(), fontsize = 20 );
    
    def plot_Rotation(totalRot,plasticRot,elasticRot,count=5):
        for ele in ops.getEleTags():
            if ele <= count:
                fig , ax = plt.subplots( 1,3 , sharex = True , sharey = True  , figsize=(20,5) )
                ax[0].plot( totalRot[ele]   ) , ax[0].axhline(c = "k") , ax[0].axvline(c = "k") , ax[0].set_xlabel("Step"), ax[0].set_ylabel( "Rotation (rad)") , ax[0].legend( ["Total rot"]) , ax[0].set_frame_on(False) ; 
                ax[1].plot( plasticRot[ele] ) , ax[1].axhline(c = "k") , ax[1].axvline(c = "k") , ax[1].set_xlabel("Step"), ax[1].set_ylabel( "Rotation (rad)") , ax[1].legend( ["Plastic"]) , ax[1].set_frame_on(False) ; 
                ax[2].plot( elasticRot[ele] ) , ax[2].axhline(c = "k") , ax[2].axvline(c = "k") , ax[2].set_xlabel("Step"), ax[2].set_ylabel( "Rotation (rad)") , ax[2].legend( ["Elastic"]) , ax[2].set_frame_on(False) ; 
                plt.suptitle( f"{ele} Rotations Graphs".upper(), fontsize = 20 );
    
    def plot_Energy(ibasicForce,jbasicForce,itotalrot,iplasticrot,ielasticrot,jtotalrot,jplasticrot,jelasticrot,count=5):
        from scipy.integrate import cumtrapz
        for ele in ops.getEleTags():
            if ele<= count:
                EH_i_total   = cumtrapz(jbasicForce[ele],jtotalrot[ele])
                EH_i_plastic = cumtrapz(jbasicForce[ele],jplasticrot[ele])
                EH_i_elastic = cumtrapz(jbasicForce[ele],jelasticrot[ele])
                EH_j_total   = cumtrapz(ibasicForce[ele],itotalrot[ele])
                EH_j_plastic = cumtrapz(ibasicForce[ele],iplasticrot[ele])
                EH_j_elastic = cumtrapz(ibasicForce[ele],ielasticrot[ele])
                """
                for index,i in enumerate(EH_i_elastic):
                    if i<0:
                        EH_i_elastic[index] *= -1
                for index,i in enumerate(EH_j_elastic):
                    if i<0:
                        EH_j_elastic[index] *= -1
                 """       
                for count,(e,p) in enumerate(zip(EH_i_elastic,EH_i_plastic)):
                    EH_i_total[count] = e+p
                for count,(e,p) in enumerate(zip(EH_j_elastic,EH_j_plastic)):
                    EH_j_total[count] = e+p

                fig , (ax1,ax2) = plt.subplots( 2,3 , sharex = True , sharey = True  , figsize=(30,10) )
                ax1[0].plot( EH_i_total   ) , ax1[0].axhline(c = "k") , ax1[0].axvline(c = "k") , ax1[0].set_xlabel("Step"), ax1[0].set_ylabel( "EH (kNM)") , ax1[0].legend( ["i node total  "]) , ax1[0].set_frame_on(False), ax1[0].set_xlim(left = 0 ) ;
                ax1[1].plot( EH_i_plastic ) , ax1[1].axhline(c = "k") , ax1[1].axvline(c = "k") , ax1[1].set_xlabel("Step"), ax1[1].set_ylabel( "EH (kNM)") , ax1[1].legend( ["i node plastic"]) , ax1[1].set_frame_on(False), ax1[1].set_xlim(left = 0 ) ;
                ax1[2].plot( EH_i_elastic ) , ax1[2].axhline(c = "k") , ax1[2].axvline(c = "k") , ax1[2].set_xlabel("Step"), ax1[2].set_ylabel( "EH (kNM)") , ax1[2].legend( ["i node elastic"]) , ax1[2].set_frame_on(False), ax1[2].set_xlim(left = 0 ) ;
                ax2[0].plot( EH_j_total   ) , ax2[0].axhline(c = "k") , ax2[0].axvline(c = "k") , ax2[0].set_xlabel("Step"), ax2[0].set_ylabel( "EH (kNM)") , ax2[0].legend( ["j node total  "]) , ax2[0].set_frame_on(False), ax2[0].set_xlim(left = 0 ) ;
                ax2[1].plot( EH_j_plastic ) , ax2[1].axhline(c = "k") , ax2[1].axvline(c = "k") , ax2[1].set_xlabel("Step"), ax2[1].set_ylabel( "EH (kNM)") , ax2[1].legend( ["j node plastic"]) , ax2[1].set_frame_on(False), ax2[1].set_xlim(left = 0 ) ;
                ax2[2].plot( EH_j_elastic ) , ax2[2].axhline(c = "k") , ax2[2].axvline(c = "k") , ax2[2].set_xlabel("Step"), ax2[2].set_ylabel( "EH (kNM)") , ax2[2].legend( ["j node elastic"]) , ax2[2].set_frame_on(False), ax2[2].set_xlim(left = 0 ) ;
                plt.suptitle( f" Frame{ele} Dissipated Energy Graphs".upper(), fontsize = 20 );fig.dpi=300
                
                path="./OutputImages/EnergyOutputs"
                if os.path.exists(path):
                    pass
                else:
                    os.mkdir(path)
                fig.savefig(path+f"/Frame{ele}.png",facecolor='white')
                #plt.savefig(path+f"/Frame{ele}.png")
                #plt.show()
    
    def plot_StressStrain(sectionRegion,count=5):
        for ele in ops.getEleTags():
            if ele<= count:
                fig , ax = plt.subplots( 1,1 , sharex = True , sharey = True  , figsize=(20,5) )
                ax.plot( sectionRegion[ele]["strain"], sectionRegion[ele]["stress"] ) , ax.axhline(c = "k") , ax.axvline(c = "k") , ax.set_xlabel("Strain (m)"), ax.set_ylabel( "Stress (kN/m2)") , ax.legend( ["stress-strain"]) , ax.set_frame_on(False) ; 
                plt.suptitle( f"Frame{ele} Stress-Strain Graph".upper(), fontsize = 20 );

    #daha hazır değil
    @staticmethod
    def ms_as_top_disp_plot(dt,*values):
        """
        dt      : zaman serileri için ortak artım zamanı
        *values : içerisinde analizlerden elde edilen tepe node noktalarından birinin deplasman listeleri
        """
        concat_disp = [i for i in values[0]]
        for key in values[1::]:
            for i in key:
                concat_disp.append(i)
        
        concat_time = [i*dt for i in range(len(concat_disp))]
        fig, ax = plt.subplots(figsize=(25,10))
        fig.subplots_adjust(bottom=0.15, left=0.2)
        ax.grid()
        ax.plot(concat_time,concat_disp)
        ax.axhline(0, color='black', lw=2)
        ax.set(xlabel="Time [Sec]", ylabel=f"Horizontal Displacement of node {max(ops.getNodeTags())}",title="main-after shock top disp/time")
        plt.show()

    @staticmethod
    def animation(Eds,time,pf=20.8,sfac_a=50.,tkt=8.):
        """
            Eds,
            time,
            pf,
            sfac_a,
            tkt
        """
        input_parameters = (20.8, 50., 8.)
        pf, sfac_a, tkt = input_parameters
        fmt_defo = {'color': 'blue', 'linestyle': 'solid', 'linewidth': 3.0,
                'marker': '', 'markersize': 6}
        # 1. animate the deformated shape
        anim = opsv.anim_defo(Eds, time, sfac_a,fmt_defo=fmt_defo, xlim=[-5, 20],
                                        ylim=[-10, 20], fig_wi_he=(30., 30.))
        plt.show()
        anim.save(filename=".\aa")

    @staticmethod
    def spyMatrix():
        """
            Rijitlik matrisinin görsel şeklini verir.
        """
        Neqn = ops.systemSize()
        # Build spy matrix
        SpyMatrix = np.identity(Neqn)
        for e in ops.getEleTags():
            dofs = []
            # Build list of DOFs for this element
            for nd in ops.eleNodes(e):
                dofs += ops.nodeDOFs(nd)

                for idof in dofs:
                    if idof < 0: # Constrained DOF
                        continue
                    for jdof in dofs:
                        if jdof < 0: # Constrained DOF
                            continue
                        SpyMatrix[idof,jdof] = 1.0

                # Determine bandwidth
        bw = 0
        for i in range(Neqn):
            bwi = 0
            # Find non-zero farthest from diagonal on this row
            for j in range(i,Neqn):
                kij = SpyMatrix[i,j]
                if kij != 0.0:
                    bwi = j-i+1
            if bwi > bw:
                bw = bwi
        # Plot       
        plt.spy(SpyMatrix,markersize=.5)
        plt.xlabel(f'Bandwidth={bw}')

    @staticmethod
    def capacityCurve(DispCtrlNode,LoadFactor):
        """
            DispCtrlNode: # List containing load factors used throughout the analysis
            LoadFactor  : # List containing displacement of control node throughout the analysis  
            
        """
        # Plot the capacity curve
                # ------------------------------------------------------------------------
        fig, ax = plt.subplots(figsize=(10,10))
        fig.subplots_adjust(bottom=0.15, left=0.2)
        ax.grid()
        ax.plot(DispCtrlNode,LoadFactor,markersize=16)
        ax.axhline(0, color='black', lw=2)
        ax.set_xlabel('Top displacement [m]')
        ax.set_ylabel('Base Shear Foce [kN]')
        ax.set_title("Capacity Curve")
        plt.show()
    
    @staticmethod
    def plot_fiber_section(SecID,fiber_sec):
        """
            SecID     :tag for the section that is generated by BuildRCrectSection procedure
            fiber_sec
        """
        matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
        opsv.plot_fiber_section(fiber_sec, matcolor=matcolor)
        plt.title('Section ID:%d' %SecID)
        plt.axis('equal')
        plt.show()
    
    @staticmethod
    def plot_model(sfac=None):
        opsv.plot_model(fig_wi_he=(30., 25.),node_supports=True)
        plt.title('plot_model after defining elements')
        
        fig_lbrt=False
        fmt_model_loads={'color': 'black', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '', 'markersize': .5}
        node_supports=True
        ax=False
        opsv.plot_loads_2d(fig_wi_he=(30., 25.),nep=10,fig_lbrt=fig_lbrt,fmt_model_loads=fmt_model_loads,node_supports=node_supports,ax=ax)
        #opsv.plot_loads_2d()
        #sfac = 80.

        #opsv.plot_defo()
        # opsv.plot_defo(sfac)
        #fmt_interp = {'color': 'blue', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '.', 'markersize': 6}
        plt.show()
    
    @staticmethod
    def plot_mander(eps_c,f_c,**kwargs):
        """
            INPUT:
                eps_c  : list of strain value for confined conc
                f_c    : list of stress value for confined conc
        """
        fc_max = max(f_c)
        fc_max_index = f_c.index(fc_max)
        epsc_ = eps_c[fc_max_index]

        fig, ax = plt.subplots(figsize=(10,10))
        fig.subplots_adjust(bottom=0.15, left=0.2)
        ax.grid()
        
        ax.plot(eps_c,f_c,**kwargs)
        ax.plot(epsc_,fc_max)
        ax.set_xlabel('Strain (mm)')  # Add an x-label to the axes.
        ax.set_ylabel('Stress (MPa)')  # Add a y-label to the axes.
        ax.set_title("Mander Model")  # Add a title to the axes.
        ax.legend(loc='upper left')
        plt.show()
    
    @staticmethod
    def plot_internal_forces(sfacN=5.e-5, sfacV=5.e-5, sfacM=5.e-5):

        opsv.plot_defo()

        # 4. plot N, V, M forces diagrams

        opsv.section_force_diagram_2d('N', sfacN)
        plt.title('Axial force distribution')

        opsv.section_force_diagram_2d('T', sfacV)
        plt.title('Shear force distribution')

        opsv.section_force_diagram_2d('M', sfacM)
        plt.title('Bending moment distribution')

        plt.show()