from msops.Material.material import opsmaterial


#p1=[s1p, e1p] #stress and strain (or force & deformation) at first point of the envelope in the positive direction
#p2=[s2p, e2p] #stress and strain (or force & deformation) at second point of the envelope in the positive direction
#p3=[s3p, e3p] #stress and strain (or force & deformation) at third point of the envelope in the positive direction
#n1=[s1n, e1n] #stress and strain (or force & deformation) at first point of the envelope in the negative direction
#n2=[s2n, e2n] #stress and strain (or force & deformation) at second point of the envelope in the negative direction
#n3=[s3n, e3n] #stress and strain (or force & deformation) at third point of the envelope in the negative direction
#pinchX = 0.2  #pinching factor for strain (or deformation) during reloading
#pinchY = 0.8  #pinching factor for stress (or force) during reloading
#damage1 = 0.  #damage due to ductility: D1(mu-1)
#damage2 = 0.  #damage due to energy: D2(Eii/Eult)
#beta = 0.     #power used to determine the degraded unloading stiffness based on ductility, mu-beta (optional, default=0.0)

inputs = [60.0,0.003,78.0,0.024,61.2,0.1,-60.0,-0.003,-78.0,-0.024,-61.2,-0.1,1.0,1.0,0.0,0.0,0.0]
opsmaterial(MaterialTypeIndex=11,matTag=1,inputArray=inputs,young_Module=30_000_000,stress_strain_test=True)