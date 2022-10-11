import os
import numpy as np

def parse_psi4(filename):
    
    f = open(filename,'r')
    lines = f.readlines()
    
    geometries = []
    geometry = []
    gradients = []
    gradient = []
    geo = False   
    grad = False
    for line in lines:
        line = line.strip().split(' ')
    
        while '' in line:
            line.remove('')
    
        if len(line) != 0 and line[0] == "Center":
            geo = True       
    
        if geo and len(line) == 0:
            geometries.append(geometry)
            geometry = []
            geo = False
    
        if geo and len(line[0]) == 1:
            for coord in line[1:4]:
               geometry.append(np.array(float(coord)))
    
        if len(line)!= 0 and line[0] == '-Total':
            grad = True
    
        if grad and len(line) == 0:
            gradients.append(gradient)
            gradient = []
            grad = False
    
    
        if grad and len(line[0]) == 1:
            for gr in line[1:4]:
               gradient.append(np.array(float(gr)))

    return np.array(geometries[::2]), np.array(gradients)

dirs = os.listdir('./')

geo_geometries = []
geo_gradients = []
for d in sorted(dirs):
    if d.split('_')[0] == 'struct':
        print(d)
        geo, grad = parse_psi4(os.path.join(d, 'output.dat'))
        geo_geometries.append(geo[0])
        geo_gradients.append(grad[0])

   
qcf_geometries, qcf_gradients = parse_psi4('qcf_psi4_2.out')


for i in range(len(qcf_geometries)):
    geo_diff = np.linalg.norm(qcf_geometries[i]*0.529177249-geo_geometries[i])    
    geo_sum = (np.linalg.norm(qcf_geometries[i]*0.529177249)+np.linalg.norm(geo_geometries[i]))/2 
    print('\ngeo diff', geo_diff/geo_sum*100)

    grad_diff = np.linalg.norm(qcf_gradients[i]-geo_gradients[i])    
    grad_sum = (np.linalg.norm(qcf_gradients[i])+np.linalg.norm(geo_gradients[i]))/2 
    print('grad diff', grad_diff/grad_sum*100)

    #grad_diff = qcf_gradients[i]-geo_gradients[i]
    #grad_sum = qcf_gradients[i]+geo_gradients[i]
    #print('graddiff', np.mean(grad_diff/(grad_sum/2)*100))


