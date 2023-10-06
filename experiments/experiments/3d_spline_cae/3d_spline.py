from abaqus import *
from abaqusConstants import *

mdb = openMdb("test.cae")

rad = 5.0
zFactor = 0.25

p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D,
    type=DEFORMABLE_BODY)

xyz=[]

for i in range(100):
    x = rad*cos(i/5.0)
    y = rad*sin(i/5.0)
    z = i*zFactor
    xyz.append( (x,y,z) )
    
for j in range( len(xyz) ):
    dtm = p.DatumPointByCoordinate(coords=xyz[j])   

d = p.datums
p.WireSpline(points=d.values(), mergeType=IMPRINT, smoothClosedSpline=ON)

mdb.saveAs("post_test.cae")