import numpy as np
from math import *
import h5py
from pathlib import Path
import os
import json
import gpIO_c

class box:
    def __init__(self,limits):
        '''
        Takes into input a list of integers pairs  {(low,high) }_i , i=0:D with D the number of dimensions
        '''
        self.dimensions=len(limits)
        self.lower_edges=[ lim[0] for lim in limits]
        self.upper_edges=[ lim[1] for lim in limits]
        self.limits=limits
        self.shape=[self.upper_edges[d] - self.lower_edges[d] +1 for d in range(self.dimensions)] 


    def grow(self,nGrow):
        limits = [ [
            self.limits[d][0] - nGrow[d] ,
            self.limits[d][1] + nGrow[d] ,
              ]  for d in range(self.dimensions)       ]
        return box(limits)
    
    

    def toJson(self):
        boxRange=[[self.lower_edges[d],self.upper_edges[d] ] for d in range(self.dimensions) ]
        return {"range": boxRange }


class geometry:

    def __init__(self,shape,domain,symmetry="none",coordinates="cartesian",bc=None,grown_shape=None):
        '''
        Creates the geometry of a cell centered grid. Does not store fields.
        '''

        self.shape=shape
        self.lower_edges= np.array([float(bound[0]) for bound in domain])
        self.upper_edges= np.array([float(bound[1]) for bound in domain])

        self.step = [ (h - l)/n for l,h,n in zip(self.lower_edges,self.upper_edges,self.shape) ]
        self.dimensions=len(shape)
        self.symmetry=symmetry
        self.coordinates=coordinates
        self.bc=bc
        if grown_shape is None:
            grown_shape=shape
        self.grown_shape=grown_shape
    
    def domainBox(self):
        limits=[(0,self.shape[d] - 1 ) for d in range(self.dimensions)]
        selectionBox=box(limits)
        return selectionBox    

    def toJson(self):
        
        domain=[ [ self.lower_edges[d],self.upper_edges[d]] for d in range(self.dimensions) ]
        return {
                "domain" : domain,
                "shape" : list(self.shape),
                "coordinates" : "cartesian"
                }
                
               
        
    def positions(self,axis,selectionBox=None,nGhosts=None,nComponents=1):
        
        if selectionBox is None:
            selectionBox=self.domainBox()
        if nGhosts is not None:
            selectionBox=selectionBox.grow(nGhosts)
                    
        
        iRange=np.arange(selectionBox.lower_edges[axis],selectionBox.upper_edges[axis]+1,1)
        x=(iRange + 0.5 )*self.step[axis]  + self.lower_edges[axis]
        
        newShape=np.concatenate((selectionBox.shape,[nComponents]))


        if self.dimensions == 1 :
            return x

        if self.dimensions == 2:
            if axis==0:
                return np.outer(x,np.ones(shape=(box.shape[1])))
            else :
                    if axis==1:
                        return np.outer(np.ones(shape=box.shape[0]),x)

        if self.dimensions == 3 :
                if axis == 0:
                    return np.outer(x,np.ones(shape=(newShape[1],newShape[2],newShape[3]))).reshape(newShape)
                else:
                    if axis == 1:
                        tmp=np.outer(x,np.ones(shape=(newShape[2],newShape[3])))
                        return np.outer(np.ones(shape=(selectionBox.shape[0])),tmp).reshape(newShape)
                    else:
                        if axis == 2:
                            tmp=np.outer(np.ones(shape=(selectionBox.shape[0],selectionBox.shape[1])),x)
                            
                            return np.outer(tmp,np.ones(shape=(newShape[3]))).reshape(newShape)






class wavefunction:

    def __init__(self,geo,phis,boxes=None,nGhosts=None,name="name"):
        '''
        geo : geometry object
        phis : list of complex numpy arrays of rank 4 [incorrect argument will lead to crashes],
        boxes : lis of boxes. Use the domain box of geo if not specified
        name : name of the field
        '''

        if boxes is None:
            boxes=[ geo.domainBox() ]

        if nGhosts is None:
            nGhosts=[0 for d in range(geo.dimensions)]
        

        assert( len(boxes)==len(phis) )

        for cbox,phi in zip(boxes,phis):
            for d in range(geo.dimensions ):
                assert( (cbox.shape[d] + 2*nGhosts[d]) == phi.shape[d])
        

        self.geometry=geo
        self.phis=phis
        self.boxes=boxes
        self.name=name
        self.nGhosts=nGhosts


    def toJson(self, dirname ):

        components=0

        if ( len(self.phis) != 0 ):
            components=self.phis[0].shape[-1]

        j={
            "geometry" : self.geometry.toJson(),
            "boxes" : [ box.toJson() for box in self.boxes ],
            "name" :  self.name,
            "nGhosts" : self.nGhosts,
            "folder" : dirname,
            "components" : components
        }

        return j


    def save(self, dirname ):
        '''
        Saves the wavefunction in a format readable by a C++ code
        '''

        j=self.toJson(dirname)
        j["folder"]=dirname
        saveDir=Path(dirname)

        if not saveDir.exists():
            os.mkdir(saveDir)
        
        filename=os.path.join(dirname,"wave.json")
        with open ( filename,"w" ) as f:
            json.dump(j,f)

        gpIO_c.saveMultifab(self.phis,j)


def loadWavefunction(dirname):
        '''
        Loads a python wavefunction object from directory dirname
        '''

        with open(os.path.join(dirname,"wave.json")) as f:
            j=json.load(f)

        jGeom = j["geometry"]
        geo=geometry(shape=jGeom["shape"],domain=jGeom["domain"])
        boxes = [ box(jBox["range"]) for jBox in j["boxes"] ]

        nComponents=j["components"]
        nGhosts=j["nGhosts"]

        phis=gpIO_c.readMultifab(j)

        shapes = [ np.concatenate((currentBox.shape,[j["components"]]))    for currentBox in boxes ]

        phis = [ np.array(phi).reshape(np.concatenate((cbox.grow(nGhosts).shape,[nComponents])),order="F")  for cbox,phi in zip(boxes,phis)  ]




        return wavefunction(geo,phis,boxes=boxes,name=j["name"],nGhosts=j["nGhosts"])










