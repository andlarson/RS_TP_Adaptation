
# Adapted from the "Using the Abaqus Scripting Interface with Abaqus/CAE".


from abaqus import *
from abaqusConstants import *
import part
import assembly
import step
import load
import interaction


myModel = mdb.models['Model-1']

myViewport = session.Viewport(name='Region syntax', origin=(20, 20), width=200, height=100)

mySketch = myModel.ConstrainedSketch(name='Sketch A', sheetSize=200)

mySketch.rectangle(point1=(-40, 30), point2=(-10,0))

mySketch.rectangle(point1=(10, 30), point2=(40, 0))

door = myModel.Part(name='Door', dimensionality=THREE_D, type=DEFORMABLE_BODY)

door.BaseSolidExtrude(sketch=mySketch, depth=20)

myAssembly = myModel.rootAssembly
doorInstance = myAssembly.Instance(name='Door-1', part=door)

pillarVertices = doorInstance.vertices.findAt(((-40, 30, 0),), ((40, 0, 0),))

myModel.StaticStep(name='impact', previous='Initial', initialInc=1, timePeriod=1)

myPillarLoad = myModel.ConcentratedForce(name='pillarForce', createStepName='impact', region=(pillarVertices,), cf1=12.5E4)

topFace = doorInstance.faces.findAt(((-25, 30, 10),))
bottomFace = doorInstance.faces.findAt(((-25, 0, 10),))

myFenderLoad = myModel.Pressure(name='pillarPressure', createStepName='impact', region=((topFace+bottomFace, SIDE1),), magnitude=10E4)

myEdge1 = doorInstance.edges.findAt(((10, 15, 20),))
myEdge2 = doorInstance.edges.findAt(((10, 15, 0),))

myDisplacementBc = myModel.DisplacementBC(name='xBC', createStepName='impact', region=(pillarVertices, myEdge1+myEdge2, topFace), u1=5.0)

faceRegion = doorInstance.faces.findAt(((-30, 15, 20), ), ((30, 15, 20), ))

mySurface = myModel.rootAssembly.Surface(name='exterior', side1Faces=faceRegion)

myFoundation = myModel.ElasticFoundation(name='elasticFloor', createStepName='Initial', surface=mySurface, stiffness=1500)

myViewport.setValues(displayedObject=door)
# myViewport.assemblyDisplay.setValues(step='impact', loads=ON, bcs=ON, fields=ON)
