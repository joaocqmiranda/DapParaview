from tkinter import filedialog
from tkinter import Tk
from tkinter import *
import os
import numpy as np
import time


def root(what,tipos):
    root = Tk()
    # root.iconbitmap("c:/...") # Defenição Icon
    root.fileName = filedialog.askopenfilename(title=what,filetypes=tipos)

    return root.fileName

def Matrix0(i,j): # Vestor in Vector
    m=[0]*i
    for ii in range(i):
        m[ii]=[0]*j
    return m

#### TENTAR

File = root("Select .dat file",[("DAT", "*.dat")])
path, filename = os.path.split(File)

if os.path.isfile(os.path.splitext(File)[0]+str(".str")):
    Filedis = os.path.splitext(File)[0]+str(".str")
else:
    Filedis = root("File .str not found. Select .str file",[("STR", "*.str")])


start_time = time.time()
abrir_File = open(str(File), "r")

#####SEPARÇÃO ENTRE CAMINHO E NOME DE FICHEIRO ##################


vetor1 = [0]
it5 = 0
it51 = 0
it512 = 0

#def solve(lis): #Trans forma os numeros de string para numero
    #for x in lis:
        #try:
        #    yield float(x)
        #except ValueError:
        #    pass


for abrir1 in abrir_File.readlines():
    vetor1.insert(it5, abrir1.split())
    vetor1[it5][:] = [x for x in vetor1[it5] if x]  # REMOÇÃO DE CÉLULAS VAZIAS
    if len(vetor1[it5]) == 1:
        it51 = it5 - it512
        it512 = it512 + 1
    it5 = it5 + 1

nb = int(vetor1[0][0])
vindex = [0]

ii = -1
iii = -1
vindex = [0] * (nb)
#tsc=[0]*4
for i in range(0, len(vetor1) - 1):
    if ii <= nb - 2:
        if vetor1[i] == []:
            ii = ii + 1
            vindex[ii] = i
    if int(len(vetor1[i])) == 4:
        tsc = vetor1[i]

SAxes = Matrix0(nb,7)

Posieuler = Matrix0(nb,6)  # X	    Y	        Z	    e1		e2	    e3
Vel = Matrix0(nb,6)  # v0x         v0y         v0z         w0x         w0y         w0z
MassIner = Matrix0(nb,4)  # mass	    Ixx		Iyy	    Izz
Forces = Matrix0(nb,6)  # F0x         F0y         F0z         M0x         M0y         M0z
Bod = Matrix0(nb,1)
###### COM IF de MODO A SE mais flexivel a alterações nas linnhas

ii = 0

numbers = [0]*14  # F0x         F0y         F0z         M0x         M0y         M0z
gravity = [0]*3

for ii in range(14):
    numbers[ii] = int(vetor1[0][ii])
for ii in range(3):
    gravity[ii] = float(vetor1[1][ii])


for i in range(nb):
    for ii in range(6):
        Posieuler[i][ii] = float(vetor1[1 + vindex[i]][ii])  # X	    Y	        Z	    e1		e2	    e3
        Vel[i][ii] = float(vetor1[2 + vindex[i]][ii])  # v0x         v0y         v0z         w0x         w0y         w0z
        Forces[i][ii] = float(vetor1[4 + vindex[i]][ii])  # F0x         F0y         F0z         M0x         M0y         M0z
    Bod[i] = float(vetor1[5 + vindex[i]][0])

    for ii in range(4):
        MassIner[i][ii] = float(vetor1[3 + vindex[i]][ii])  # mass	    Ixx		Iyy	    Izz
    for ii in range(7):
        SAxes[i][ii] = float(vetor1[6 + vindex[i]][ii]) # a1 a2 a3 n #adiciona SAxes[i, :] = (SAxes[i, :]) + list(solve(vetor1[9 + vindex[i]][0:4]))

new_file = open(os.path.splitext(File)[0]+str(".py"), "w")

new_file.write("# Convertion from PC Crash to paraview 5.9.0 ")
new_file.write("\nimport time")
new_file.write('\nprint("Start Reset")')
new_file.write('\nstart_time = time.time()')
new_file.write('\nprint(time.strftime("%H:%M:%S", time.localtime()))')
new_file.write("\nResetSession()")
new_file.write('\nprint("Start Bodies Generation")')
new_file.write('\ntime.strftime("%H:%M:%S", time.localtime())')
new_file.write("\nfrom paraview.simple import *")
new_file.write("\nfrom paraview import servermanager")
new_file.write("\n#### disable automatic camera reset on 'Show'")
new_file.write("\nparaview.simple._DisableFirstRenderCameraReset()")
new_file.write('\nprint(time.strftime("%H:%M:%S", time.localtime()))')
new_file.write('\nimport numpy as np')
new_file.write('\nfilename="'+Filedis+'"')
new_file.write('\n                                                                                                                                                                                          ')
new_file.write('\ndef m(v):                                                                                                                                                                                 ')
new_file.write('\n  ma=[v[0],v[4],v[8],v[12],v[1],v[5],v[9],v[13],v[2],v[6],v[10],v[14],v[3],v[7],v[11],v[15]]                                                                                              ')
new_file.write('\n  return ma                                                                                                                                                                               ')
new_file.write('\n                                                                                                                                                                                          ')
new_file.write('\nanimdata = np.genfromtxt(filename, dtype=None, skip_header=1,                                                                                                                                 ')
new_file.write('\n    names=m(["M11","M12","M13","M14",                                                                                                                                             ')
new_file.write('\n             "M21","M22","M23","M24",                                                                                                                                             ')
new_file.write('\n             "M31","M32","M33","M34",                                                                                                                                             ')
new_file.write('\n             "M41","M42","M43","M44"]))                                                                                                                           ')
new_file.write('\n                                                                                                                                 ')
########################AminProgramableSource#################################
new_file.write("\nprogrammableSource1 = ProgrammableSource(registrationName='AnimProgrammableSource')                                                                                                        ")
new_file.write('\n# Properties modified on programmableSource1                                                                                                                                              ')
new_file.write("\nprogrammableSource1.OutputDataSetType = 'vtkTable'                                                                                                                                        ")
new_file.write('\nprogrammableSource1.Script = """import numpy as np                                                                                                                                        ')
new_file.write('\nimport math                                                                ')
new_file.write('\nfilename="'+Filedis+'"')
new_file.write('\n                                                                                                                                                                                          ')
new_file.write('\ndef m(v):                                                                                                                                                                                 ')
new_file.write('\n  ma=[v[0],v[4],v[8],v[12],v[1],v[5],v[9],v[13],v[2],v[6],v[10],v[14],v[3],v[7],v[11],v[15]]                                                                                              ')
new_file.write('\n  return ma                                                                                                                                                                               ')
new_file.write('\ndef rotationMatrixToEulerAngles(R):                                            ')
new_file.write('\n    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])                         ')
new_file.write('\n    singular = sy < 1e-6                                                       ')
new_file.write('\n    if  not singular :                                                         ')
new_file.write('\n        x = math.atan2(R[2,1] , R[2,2])                                        ')
new_file.write('\n        y = math.atan2(-R[2,0], sy)                                            ')
new_file.write('\n        z = math.atan2(R[1,0], R[0,0])                                         ')
new_file.write('\n    else :                                                                     ')
new_file.write('\n        x = math.atan2(-R[1,2], R[1,1])                                        ')
new_file.write('\n        y = math.atan2(-R[2,0], sy)                                            ')
new_file.write('\n        z = 0                                                                  ')
new_file.write('\n    return (x*180/np.pi, y*180/np.pi, z*180/np.pi)                        ')
new_file.write('\ndef rotationMatrixToEulerParameters(M1):                                                     ')
new_file.write('\n    #https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/')
new_file.write('\n    tr= float(M1[0, 0] + M1[1, 1] + M1[2, 2])                                               ')
new_file.write('\n    if tr > 0:                                                                              ')
new_file.write('\n        S = float(np.math.sqrt(tr+1.0)*2) # S=4*qw                                          ')
new_file.write('\n        qw = 0.25 * S                  #e0 r                                                ')
new_file.write('\n        qx = (M1[2, 1] - M1[1, 2]) / S #e1 i                                                ')
new_file.write('\n        qy = (M1[0, 2] - M1[2, 0]) / S #e2 j                                                ')
new_file.write('\n        qz = (M1[1, 0] - M1[0, 1]) / S #e3 k                                                ')
new_file.write('\n    elif M1[0, 0] > M1[1, 1] and M1[0, 0] > M1[2, 2]:                                       ')
new_file.write('\n        S = float(np.math.sqrt(1.0 + M1[0, 0] - M1[1, 1] - M1[2, 2])*2) # S=4*qw            ')
new_file.write('\n        qw = (M1[2, 1] - M1[1, 2]) / S #e0 r                                                ')
new_file.write('\n        qx = 0.25 * S                  #e1 i                                                ')
new_file.write('\n        qy = (M1[1, 0] + M1[0, 1]) / S #e2 j                                                ')
new_file.write('\n        qz = (M1[0, 2] + M1[2, 0]) / S #e3 k                                                ')
new_file.write('\n    elif M1[1, 1] > M1[2, 2]:                                                               ')
new_file.write('\n        S = float(np.math.sqrt(1.0 + M1[1, 1] - M1[0, 0] - M1[2, 2])*2) # S=4*qw            ')
new_file.write('\n        qw = (M1[0, 2] - M1[2, 0]) / S #e0 r                                                ')
new_file.write('\n        qx = (M1[1, 0] + M1[0, 1]) / S #e1 i                                                ')
new_file.write('\n        qy = 0.25 * S                  #e2 j                                                ')
new_file.write('\n        qz = (M1[2, 1] + M1[1, 2]) / S #e3 k                                                ')
new_file.write('\n    else :                                                                                   ')
new_file.write('\n        S = float(np.math.sqrt(1.0 + M1[2, 2] - M1[1, 1] - M1[0, 0])*2) # S=4*qw            ')
new_file.write('\n        qw = (M1[1, 0] - M1[0, 1]) / S #e0 r                                                ')
new_file.write('\n        qx = (M1[0, 2] + M1[2, 0]) / S #e1 i                                                ')
new_file.write('\n        qy = (M1[1, 2] + M1[2, 1]) / S #e2 j                                                ')
new_file.write('\n        qz = 0.25 * S                  #e3 k                                                ')
new_file.write('\n    return (qw,qx,qy,qz)                                                                    ')
new_file.write('\n                                                                                                                                                                                          ')
new_file.write('\nanimdata = np.genfromtxt(filename, dtype=None, skip_header=1,                                                                                                                                 ')
new_file.write('\n    names=m(["M11","M12","M13","M14",                                                                                                                                             ')
new_file.write('\n             "M21","M22","M23","M24",                                                                                                                                             ')
new_file.write('\n             "M31","M32","M33","M34",                                                                                                                                             ')
new_file.write('\n             "M41","M42","M43","M44"]))                                                                                                                           ')
new_file.write("\nrotc = [0]*int(''.join(map(str, animdata.shape)))                                                                                     ")
new_file.write("\nfor i in range(int(''.join(map(str, animdata.shape)))):                                                                                      ")
new_file.write('\n    rotc[i]=rotationMatrixToEulerAngles(np.matrix([[animdata["M11"][i],animdata["M12"][i],animdata["M13"][i]]  ')
new_file.write('\n                                                  ,[animdata["M21"][i],animdata["M22"][i],animdata["M23"][i]]  ')
new_file.write('\n                                                  ,[animdata["M31"][i],animdata["M32"][i],animdata["M33"][i]]]))')
new_file.write("\nrotc=np.array(rotc,dtype=[('RX', '<f8'), ('RY', '<f8'), ('RZ', '<f8')])                  ")
new_file.write('\n                                                                                                                                 ')
new_file.write('\nfor name in animdata.dtype.names:                                                                                                                                                             ')
new_file.write('\n  array = animdata[name]                                                                                                                                                                     ')
new_file.write('\n  output.RowData.append(array, name)    ')
new_file.write('\nfor name in rotc.dtype.names:                                                      ')
new_file.write('\n  array = rotc[name]                                                                              ')
new_file.write('\n  output.RowData.append(array, name)   ')
new_file.write("\nrotc1 = [0]*int(''.join(map(str, animdata.shape)))                                                                                     ")
new_file.write("\nfor i in range(int(''.join(map(str, animdata.shape)))):                                                                                      ")
new_file.write('\n    rotc1[i]=rotationMatrixToEulerParameters(np.matrix([[animdata["M11"][i],animdata["M12"][i],animdata["M13"][i]]  ')
new_file.write('\n                                                  ,[animdata["M21"][i],animdata["M22"][i],animdata["M23"][i]]  ')
new_file.write('\n                                                  ,[animdata["M31"][i],animdata["M32"][i],animdata["M33"][i]]]))')
new_file.write("\nrotc1=np.array(rotc1,dtype=[('e0', '<f8'), ('e1', '<f8'), ('e2', '<f8'),('e3', '<f8')])                  ")
new_file.write('\n                                                                                                                                 ')
new_file.write('\nfor name in animdata.dtype.names:                                                                                                                                                             ')
new_file.write('\n  array = animdata[name]                                                                                                                                                                     ')
new_file.write('\n  output.RowData.append(array, name)    ')
new_file.write('\nfor name in rotc1.dtype.names:                                                      ')
new_file.write('\n  array = rotc1[name]                                                                              ')
new_file.write('\n  output.RowData.append(array, name)"""    ')
new_file.write("\ntext1 = Text(registrationName='General Data')  ")
new_file.write('\ntext1.Text = """')
i1=0
nc=4
for i in range(14):
    if numbers[i]!= 0 and i==0 :
        new_file.write(str(numbers[0]) + ' Bodies/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==1:
        new_file.write(str(numbers[1]) + ' Spherical Joints/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==2:
        new_file.write(str(numbers[2]) + ' Revolute Joins/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==3:
        new_file.write(str(numbers[3]) + ' Cilindrical Joints/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==4:
        new_file.write(str(numbers[4]) + ' Ridig-Flex Joins/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==5:
        new_file.write(str(numbers[5]) + ' Grounded Bodies/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==6:
        new_file.write(str(numbers[6]) + ' Simple Constraints/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==7:
        new_file.write(str(numbers[7]) + ' Spring-Damper-Actuator Elements/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==8:
        new_file.write(str(numbers[8]) + ' Driver Constraints/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==9:
        new_file.write(str(numbers[9]) + ' Wheels/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==11:
        new_file.write(str(numbers[11]) + ' Steering Systems/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==12:
        new_file.write(str(numbers[12]) + ' Contacting Pars/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')
    if numbers[i] != 0 and i==13:
        new_file.write(str(numbers[13]) + ' Non Contacting Pars/')
        i1=i1+1
        if i1 == nc:
            new_file.write('\n')

new_file.write('\nGravity:')
if gravity[0]!=0:
    new_file.write(' X= '+str(gravity[0]))
if gravity[1]!=0:
    new_file.write(' Y= '+str(gravity[1]))
if gravity[2]!=0:
    new_file.write(' Z= '+str(gravity[2]))
new_file.write('\nType of system equation solver '+'use LU factorization.' if numbers[10]==(0) else 'use Augmented Lagrangian Formulation.' if numbers[10]==(1) else 'use Sparse matrix solver.')
new_file.write('\n"""')
new_file.write("\ntext1Display = Show(text1, GetActiveViewOrCreate('RenderView'), 'TextSourceRepresentation') ")
new_file.write("\ntext1Display.WindowLocation = 'Lower Center' ")
new_file.write("\ntext1Display.FontSize = 13")

##### FAZER FUNÇÃO DE PREVIEW COM O PYGRAPH #########################
print(Posieuler[1])

for i in range(nb):
    nbs=str(i)
    new_file.write("\nsuperquadric"+nbs+" = Superquadric(registrationName='Superquadric"+nbs+"')")
    #new_file.write("\n# Properties modified on superquadric")
    #new_file.write("\nsuperquadric"+nbs+".Center = [0.0,0.0,0.0]")
    new_file.write("\nsuperquadric"+nbs+".Scale = [" + str(SAxes[i][0])+","+str(SAxes[i][2])+","+str(SAxes[i][1])+"]")
    #new_file.write("\nsuperquadric"+nbs+".ThetaResolution = 16")
    #new_file.write("\nsuperquadric"+nbs+".PhiResolution = 16")
    #new_file.write("\nsuperquadric"+nbs+".Thickness = 0.3333")
    new_file.write("\nsuperquadric"+nbs+".ThetaRoundness = " + str(int(Bod[i])/SAxes[i][4]))
    new_file.write("\nsuperquadric"+nbs+".PhiRoundness = " + str(int(Bod[i])/SAxes[i][3]))
    new_file.write("\nsuperquadric"+nbs+".Size = 1")
    new_file.write("\nsuperquadric"+nbs+".Toroidal = 0") # Super elipsoid
    #new_file.write("\nrenderView1 = GetActiveViewOrCreate('RenderView')")
    # show data in view
    new_file.write("\ntransform" + nbs + " = Transform(registrationName='Transform"+nbs+"', Input=superquadric" + nbs + ")")
    new_file.write('\ntp' + nbs + ' = servermanager._getPyProxy(servermanager._getPyProxy(servermanager.CreateProxy("transforms", "Transform")))')
    new_file.write("\ntp" + nbs + ".Matrix=m(animdata[" + nbs + "])")
    new_file.write("\ntransform" + nbs + " .Transform = tp" + nbs + "")
    new_file.write("\ntransform" + nbs + "Display = Show(transform" + nbs + ")")
    new_file.write("\ntransform" + nbs + "Display.SetRepresentationType('Surface With Edges')")
    new_file.write("\ntext1 = Text(registrationName='Body " + nbs + " Data')  ")
    new_file.write('\ntext1.Text = """Body ' + nbs)
    new_file.write(('\nPosition [X,Y,Z] = '+str(Posieuler[i][0:3])) if sum(Posieuler[i][0:3])!=(0) else '')
    new_file.write(('\nEuler parameters [e1,e2,e3]= '+str(Posieuler[i][3:6])) if sum(Posieuler[i][3:6])!=(0) else '')
    new_file.write(('\nLinear Velocity [X,Y,Z] = '+str(Vel[i][0:3])+'  ') if sum(Vel[i][0:3])!=(0) else '')
    new_file.write(('\nAngular Velocity ['+chr(92)+'u03BE,'+chr(92)+'u03B7,'+chr(92)+'u03B6] = '+str(Vel[i][3:6])) if sum(Vel[i][3:6])!=(0) else '')
    new_file.write(('\nMass='+str(MassIner[i][0])+'  ')  if (MassIner[i][0])!=(0) else '')
    new_file.write(('\nMoments of Inertia ['+chr(92)+'u03BE,'+chr(92)+'u03B7,'+chr(92)+'u03B6] = '+str(MassIner[i][1:4])) if sum(MassIner[i][1:4])!=(0) else '')
    new_file.write(('\nContact Forces [X,Y,Z] = ' + str(Forces[i][0:3]) + '  ')  if sum(Forces[i][0:3])!=(0) else '')
    new_file.write(('\nContact Moments ['+chr(92)+'u03BE,'+chr(92)+'u03B7,'+chr(92)+'u03B6] = ' + str(Forces[i][3:6])) if sum(Forces[i][3:6])!=(0) else '')
    new_file.write(('\nSurface Properties [a1,a2,a3,epsilon1,epsilon2,friction,restitution]\n'+str(SAxes[i])) if sum(SAxes[i])!=(0) else '')
    new_file.write('\n*with [X,Y,Z] as coordinates of Global System')
    new_file.write('\nand $['+chr(92)+'u03BE,'+chr(92)+'u03B7,'+chr(92)+'u03B6]$ as coordinates of Local System')
    new_file.write('\n"""')
    #new_file.write('\nprogrammableSource' + nbs + '.Script = """')
    #new_file.write("\noutput.PointData.append(" + Mass[i] + ", 'Mass')")
    #new_file.write("\noutput.PointData.append(" + InertiaX[i] + ", 'InertiaX')")
    #new_file.write("\noutput.PointData.append(" + InertiaY[i] + ", 'InertiaY')")
    #new_file.write("\noutput.PointData.append(" + InertiaZ[i] + ", 'InertiaZ')")
    #new_file.write("\noutput.PointData.append(" + Stiff[i] + ", 'Stiff')")
    #new_file.write("\noutput.PointData.append(" + RestCoef[i] + ", 'RestCoef')")
    #new_file.write("\noutput.PointData.append(" + Friction1[i] + ", 'Friction1')")
    #new_file.write("\noutput.PointData.append(" + Friction2[i] + ", 'Friction2')")
    #new_file.write("\noutput.PointData.append(" + Color1[i] + ", 'Color1')")
    #new_file.write("\noutput.PointData.append(" + Color2[i] + ", 'Color2')")
    #new_file.write("\noutput.PointData.append(" + EulerAngles[i] + ", 'EulerAngles')")
    #new_file.write("\noutput.PointData.append(" + LinVel[i] + ", 'LinVel')")
    #new_file.write("\noutput.PointData.append(" + AngVel[i] + ", 'AngVel')")
    #new_file.write('\n"""')
    # trace defaults for the display properties.
    #new_file.write("\nsuperquadric" + nbs + "Display.Representation = 'Surface'") #Surface With Edges
    #new_file.write("\nsuperquadric" + nbs + "Display.ColorArrayName = [None, '']")
    #new_file.write("\nsuperquadric" + nbs + "Display.SelectTCoordArray = 'TextureCoords'")
    #new_file.write("\nsuperquadric" + nbs + "Display.SelectNormalArray = 'Normals'")
    #new_file.write("\nsuperquadric" + nbs + "Display.SelectTangentArray = 'None'")
    #new_file.write("\nsuperquadric" + nbs + "Display.OSPRayScaleArray = 'Normals'")
    #new_file.write("\nsuperquadric" + nbs + "Display.OSPRayScaleFunction = 'PiecewiseFunction'")
    #new_file.write("\nsuperquadric" + nbs + "Display.SelectOrientationVectors = 'None'")
    #new_file.write("\nsuperquadric" + nbs + "Display.ScaleFactor = 0.30000000000000004")
    #new_file.write("\nsuperquadric" + nbs + "Display.SelectScaleArray = 'None'")
    #new_file.write("\nsuperquadric" + nbs + "Display.GlyphType = 'Arrow'")
    #new_file.write("\nsuperquadric" + nbs + "Display.GlyphTableIndexArray = 'None'")
    #new_file.write("\nsuperquadric" + nbs + "Display.GaussianRadius = 0.015")
    #new_file.write("\nsuperquadric" + nbs + "Display.SetScaleArray = ['POINTS', 'Normals']")
    #new_file.write("\nsuperquadric" + nbs + "Display.ScaleTransferFunction = 'PiecewiseFunction'")
    #new_file.write("\nsuperquadric" + nbs + "Display.OpacityArray = ['POINTS', 'Normals']")
    #new_file.write("\nsuperquadric" + nbs + "Display.OpacityTransferFunction = 'PiecewiseFunction'")
    #new_file.write("\nsuperquadric" + nbs + "Display.Position = " + str(Position[i]))
    #new_file.write("\nsuperquadric" + nbs + "Display.DataAxesGrid.Position = " + str(Position[i]))
    #new_file.write("\nsuperquadric" + nbs + "Display.PolarAxes.Translation = " + str(Position[i]))
    ##new_file.write("\nsuperquadric" + nbs + "Display.Scale = " + str(SAxes[i][0:3]))
    #new_file.write("\nsuperquadric" + nbs + "Display.DataAxesGrid.Scale = " + str(SAxes[i][0:3]))
    #new_file.write("\nsuperquadric" + nbs + "Display.PolarAxes.Scale =" + str(SAxes[i][0:3]))
    # Properties modified on superquadric1Display
    #phi = EulerAngles[i][0]*180/np.pi
    #theta = EulerAngles[i][1]*180/np.pi
    #sigma = EulerAngles[i][2]*180/np.pi
    #new_file.write("\nsuperquadric" + nbs + "Display.Orientation = "+str([phi,theta,sigma]))
    # Properties modified on superquadric1Display.PolarAxes
    #new_file.write("\nsuperquadric" + nbs + "Display.PolarAxes.Orientation = "+str([phi,theta,sigma]))
    # change solid color
    #new_file.write("\nsuperquadric" + nbs + "Display.AmbientColor = "+str(Color1[i]))
    #new_file.write("\nsuperquadric" + nbs + "Display.DiffuseColor = "+str(Color1[i]))

    #new_file.write("\nsuperquadric" + nbs + "Display.EdgeColor = " + str([Color2[i][0]/255,Color2[i][1]/255,Color2[i][2]/255]))

#############DISTANCE##################################

abrir_File2 = open(str(Filedis), "r")

vetor2 = [0]
it5 = 0
it51 = 0
it512 = 0

for abrir2 in abrir_File2.readlines():
    vetor2.insert(it5, abrir2.split())
    it5 = it5 + 1
vector22 = vetor2[1:it5]

#############################SAVE MOTION DATA################################
timee = [0]*(it5-1)
Distance = [0]*(it5-1)

for i in range(0, it5-1): # String to float
    for ii in range(0, len(vector22[i])):
        vector22[i][ii] = float(vector22[i][ii])
new_file.write('\nprint("Writing Animation Scene")')
new_file.write('\nprint(time.strftime("%H:%M:%S", time.localtime()))')

new_file.write('\nscene = GetAnimationScene()                          ')

new_file.write('\nPythonAnimationCue1 = PythonAnimationCue()           ')
new_file.write('\nPythonAnimationCue1.Script= """                      ')

new_file.write('\nfrom paraview.simple import *')
new_file.write('\nfrom paraview import servermanager')
new_file.write('\nimport numpy as np                                                                                  ')
new_file.write('\n# assuming data.csv is a CSV file with the 1st row being the names names for                        ')
new_file.write('\n# the columns                                                                                       ')
new_file.write('\ndef m(v):                                                                                           ')
new_file.write('\n  ma=[v[0],v[4],v[8],v[12],v[1],v[5],v[9],v[13],v[2],v[6],v[10],v[14],v[3],v[7],v[11],v[15]]        ')
new_file.write('\n  return ma                                                                                         ')
new_file.write('\nanimdata = np.genfromtxt("'+Filedis+'", dtype=None, skip_header=1,                                   ')
new_file.write('\nnames=m(["M11","M12","M13","M14",                                                                   ')
new_file.write('\n         "M21","M22","M23","M24",                                                                   ')
new_file.write('\n         "M31","M32","M33","M34",                                                                   ')
new_file.write('\n         "M41","M42","M43","M44"]))                                                                 ')
new_file.write('\n                                                            ')
for i in range(0, nb):
    nbs=str(i)
    new_file.write('\ntransform'+ nbs +' = FindSource("Transform'+ nbs +'")                 ')
    new_file.write('\ntp'+ nbs +' = servermanager._getPyProxy(servermanager._getPyProxy(servermanager.CreateProxy("transforms", "Transform")))')

new_file.write('\n                                                            ')
new_file.write('\ndef start_cue(self): pass                                   ')
new_file.write('\n                                                            ')
new_file.write('\ndef end_cue(self): pass                                     ')
new_file.write('\n                                                            ')
#new_file.write('\nii=-'+nb-1)
new_file.write('\ndef tick(self):                                             ')
#new_file.write('\n  ii=ii+'+nb-1)
new_file.write('\n  ii = round(self.GetAnimationTime() * (GetAnimationScene().NumberOfFrames - 1))*'+str(nb))
new_file.write('\n  #print(ii/'+str(nb)+')')

for i in range(0, int(nb)):
    nbs=str(i)
    new_file.write("\n  tp"+ nbs +".Matrix = m(animdata[ii+" + str(i) + "])")
    new_file.write("\n  transform" + nbs + ".Transform = tp"+ nbs) #[Distance[ii]["+str(i)+"][0]-Distance[ii]["+str(i)+"][0],Distance[ii]["+str(i)+"][1]-Distance[ii]["+str(i)+"][1],Distance[ii]["+str(i)+"][2]-Distance[ii]["+str(i)+"][2]]")
new_file.write("\n  Render()"      )
new_file.write("\n  #ResetCamera()"      )
new_file.write('\n"""')
print(tsc)
new_file.write('\nscene.Cues.append(PythonAnimationCue1)               ')
new_file.write("\nscene.PlayMode = 'Sequence'"      )
new_file.write('\nscene.StartTime = '+ str(tsc[0]))
new_file.write('\nscene.EndTime = '+(tsc[1]))
new_file.write('\nscene.NumberOfFrames = int('+(tsc[1])+'/'+(tsc[2])+'/'+(tsc[3])+'+1)')


new_file.write("\ndef fp(sb):                                                                ")
new_file.write('\n  if str(sb).isdigit() == True and int(sb) in range(0,'+str(nb)+'):                                                                                                                          ')
new_file.write("\n	    plotData = FindSource('PlotData'+str(sb))")
new_file.write("\n	    lineChartView1 = GetActiveViewOrCreate('XYChartView')                  ")
new_file.write("\n	    # show data in view                                                    ")
new_file.write("\n	    plotData1Display = GetDisplayProperties(plotData, view=lineChartView1)")
new_file.write("\n	    plotData1Display.XArrayName = 'Time'                                   ")
new_file.write("\n	    lineChartView1.ChartTitle = 't={time}s'  #Body '+str(sb)+',           ")
new_file.write("\n	    lineChartView1.LeftAxisTitle = 'Position, Euler Angle("+chr(92)+"xb0)/Parameters'")
new_file.write("\n	    lineChartView1.BottomAxisTitle = 'Time'                                ")
new_file.write('\ndef g(sb):                                                                                                                                                                                  ')
new_file.write('\n	if str(sb).isdigit() == True and int(sb) in range(0,'+str(nb)+'):                                                                                                                          ')
new_file.write("\n		programmableSource1 = ProgrammableSource(registrationName='ProgrammableSource'+str(sb))                                                                                             ")
new_file.write('\n		# Properties modified on programmableSource1                                                                                                                                    ')
new_file.write("\n		programmableSource1.OutputDataSetType = 'vtkPolyData'                                                                                                                                 ")
new_file.write("\n		programmableSource1.OutputDataSetType = 'vtkTable'                                                                                                                                        ")
new_file.write('\n		programmableSource1.Script = """from paraview import servermanager                                                                                                              ')
new_file.write('\nfrom paraview.simple import GetDisplayProperties, FindSource, ResetCamera, Render, GetAnimationScene')
new_file.write('\nimport numpy as np')
new_file.write('\nfilename="'+Filedis+'"')
new_file.write('\nSelBody = """+str(sb)+"""                                                                                                                                                                     ')
new_file.write('\nTotBody = '+str(nb))
new_file.write('\n                                                                                                                                                                                          ')
new_file.write('\ndef ms(v):                                                                                                                                                                                 ')
new_file.write('\n  mas=[v[0],v[4],v[8],v[12],v[1],v[5],v[9],v[13],v[2],v[6],v[10],v[14],v[3],v[7],v[11],v[15],v[16],v[17],v[18],v[19],v[20],v[21],v[22]]                                                                                              ')
new_file.write('\n  return mas                                                                                                                                                                               ')
new_file.write('\n                                                                                                                                                                                          ')
new_file.write('\nits_data = servermanager.Fetch(FindSource(\'AnimProgrammableSource\'))     ')
new_file.write('\nnames=ms(["B"+str(SelBody)+"M11","B"+str(SelBody)+"M12","B"+str(SelBody)+"M13","B"+str(SelBody)+" X",                                                                                 ')
new_file.write('\n          "B"+str(SelBody)+"M21","B"+str(SelBody)+"M22","B"+str(SelBody)+"M23","B"+str(SelBody)+" Y",                                                                                 ')
new_file.write('\n          "B"+str(SelBody)+"M31","B"+str(SelBody)+"M32","B"+str(SelBody)+"M33","B"+str(SelBody)+" Z",                                                                                 ')
new_file.write('\n          "B"+str(SelBody)+"M41","B"+str(SelBody)+"M42","B"+str(SelBody)+"M43","B"+str(SelBody)+"M44",                 ')
new_file.write('\n          "B"+str(SelBody)+" $'+chr(92)+'u03B1$","B"+str(SelBody)+" $'+chr(92)+'u03B2$","B"+str(SelBody)+" $'+chr(92)+'u03B3$",                 ')
new_file.write('\n          "B"+str(SelBody)+" e0","B"+str(SelBody)+" e1","B"+str(SelBody)+" e2","B"+str(SelBody)+" e3"])                  ')
new_file.write('\n#Selection of lines to delete                                                                                                                                  ')
new_file.write('\ndata=int(its_data.GetRowData().GetNumberOfTuples())                                                                                                         ')
new_file.write('\ndelrow=[0]*(int((data)*(1-1/TotBody)))                                                                                                                                                 ')
new_file.write('\nii=0                                                                                                                                                                                      ')
new_file.write('\niii=0                                                                                                                                                                                     ')
new_file.write('\nfor i in range(0,(data)):                                                                                                                                                              ')
new_file.write('\n   if i == ii*TotBody + SelBody :                                                                                                                                                         ')
new_file.write('\n      ii=ii+1                                                                                                                                                                             ')
new_file.write('\n   else:                                                                                                                                                                                  ')
new_file.write('\n      delrow[iii]=int(i)                                                                                                                                                                  ')
new_file.write('\n      iii=iii+1                                                                                                                                                                           ')
new_file.write('\n                                                                                                                                                                                          ')
new_file.write('\nfim=(GetAnimationScene().EndTime)                       ')
new_file.write('\ninc=(fim/GetAnimationScene().NumberOfFrames)')
new_file.write('\noutput.RowData.append(np.arange(0,fim,inc),"Time") ')
new_file.write('\n#Append Data to Spreadsheet                                                                                                                                                               ')
new_file.write('\nfor i in range(0,len(names)):                                                                                                                                                             ')
new_file.write("\n   if names[i] in ['B'+str(SelBody)+' X', 'B'+str(SelBody)+' Y', 'B'+str(SelBody)+' Z','B'+str(SelBody)+' $"+chr(92)+"u03B1$','B'+str(SelBody)+' $"+chr(92)+"u03B2$','B'+str(SelBody)+' $"+chr(92)+"u03B3$','B'+str(SelBody)+' e0','B'+str(SelBody)+' e1','B'+str(SelBody)+' e2','B'+str(SelBody)+' e3']:      ")
new_file.write('\n      output.RowData.append(np.delete(np.array(its_data.GetRowData().GetArray(i)), delrow, 0), names[i])"""                                                                                                                            ')
new_file.write('\n	                                                                                                                            ')
new_file.write('\n		                                                                                                                                                                                          ')
new_file.write('\n		# create a new "Plot Data"                                                                                                                                                                ')
new_file.write("\n		plotData1 = PlotData(registrationName='PlotData'+str(sb), Input=programmableSource1)                                                                                                      ")
new_file.write('\n		                                                                                                                                                                                          ')
new_file.write('\n		# Create a new "Line Chart View"                                                                                                                                                          ')
new_file.write("\n		lineChartView1 = CreateView('XYChartView')                                                                                                                                                ")
new_file.write('\n		# show data in view                                                                                                                                                                       ')
new_file.write("\n		plotData1Display = Show(plotData1, lineChartView1, 'XYChartRepresentation')                                                                     ")
new_file.write("\n		plotData1Display.XArrayName = 'Time'                                                                     ")
new_file.write("\n		layout1 = GetLayout()")
new_file.write("\n		programmableSource1.OutputDataSetType = 'vtkPolyData'                                                                                                                                        ")
new_file.write("\n		programmableSource1.OutputDataSetType = 'vtkTable'                                                                                                                                        ")
new_file.write("\n		tableToPoints1 = TableToPoints(registrationName='TableToPoints'+str(sb), Input=programmableSource1)   ")
new_file.write("\n		tableToPoints1.XColumn = 'B'+str(sb)+' X'                                                              ")
new_file.write("\n		tableToPoints1.YColumn = 'B'+str(sb)+' Y'                                                              ")
new_file.write("\n		tableToPoints1.ZColumn = 'B'+str(sb)+' Z'                                                              ")
new_file.write("\n		programmableSource1.OutputDataSetType = 'vtkPolyData'                                                                                                                                        ")
new_file.write("\n		programmableSource1.OutputDataSetType = 'vtkTable'                                                                                                                                        ")
new_file.write("\n		fp(sb)                                                                                                                                  ")
new_file.write('\nfor i in range(0,'+str(nb)+'):                                                                                                                                                              ')
new_file.write('\n   g(i)                                                                                                                                                         ')
new_file.write("\nrenderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')")
new_file.write('\nrenderView1.AxesGrid.Visibility = 1')
new_file.write('\nrenderView1.AxesGrid.XTitleBold = 1')
new_file.write('\nrenderView1.AxesGrid.YTitleBold = 1')
new_file.write('\nrenderView1.AxesGrid.ZTitleBold = 1')
new_file.write("\nrenderView1.AxesGrid.XTitle = 'X'")
new_file.write("\nrenderView1.AxesGrid.YTitle = 'Y'")
new_file.write("\nrenderView1.AxesGrid.ZTitle = 'Z'")
new_file.write('\nrenderView1.AxesGrid.ShowGrid = 1')
new_file.write('\nrenderView1.AxesGrid.CustomBounds = [servermanager.Fetch(FindSource(\'AnimProgrammableSource\')).GetColumnByName("M14").GetRange()[0],')
new_file.write('servermanager.Fetch(FindSource(\'AnimProgrammableSource\')).GetColumnByName("M14").GetRange()[-1],')
new_file.write('servermanager.Fetch(FindSource(\'AnimProgrammableSource\')).GetColumnByName("M24").GetRange()[0],')
new_file.write('servermanager.Fetch(FindSource(\'AnimProgrammableSource\')).GetColumnByName("M24").GetRange()[-1], ')
new_file.write('servermanager.Fetch(FindSource(\'AnimProgrammableSource\')).GetColumnByName("M34").GetRange()[0],')
new_file.write('servermanager.Fetch(FindSource(\'AnimProgrammableSource\')).GetColumnByName("M34").GetRange()[-1]]')
new_file.write('\nrenderView1.AxesGrid.UseCustomBounds = 1')
new_file.write("\nResetCamera()"      )
new_file.write('\nprint("End")')
new_file.write('\nprint(time.strftime("%H:%M:%S", time.localtime()))')
new_file.write("\nexetime2=str(time.time() - start_time)")
exetime1=str(time.time() - start_time)
new_file.write("\nprint('Coneversion Time "+exetime1+" secunds, Exetcution Time '+exetime2+' secunds')")
new_file.write("\nprint('Type <g(Body Number)> to see the respective Plot')")
new_file.write("\nprint('Type <fp(Body Number)> to add Names to Axis (After Selecting an existing one or Just run to create a new chart)')")
print("Coneversion Time "+exetime1+" secunds")
