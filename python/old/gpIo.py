import numpy as np
import scipy as sp 
import matplotlib.pylab as plt
import pyarrow.parquet as pq
import pandas as pd 
import json
from os import listdir 
import re
import gp
from os.path import join,isfile



class gp_simulation:

    def __init__(self,dirname=None):
        if dirname is not None:
            self.load(dirname)

    def load(self,dirname):

        # load settings
        with open (dirname + "/description.json") as f:
            self.description=json.load(f)
        # load boxed data
        self.data=[None for i in range(len( self.description["boxes"]))]

        for f in listdir(dirname):
            name_match=re.match("box=(\d+)\.parquet",f)
            if name_match is not None:
                box_data=pq.read_table(dirname + "/" + f).to_pandas()
                index=int(name_match.group(1))
                self.data[ index ] = box_data  

    def getField( self,i , name ):

        box=self.description["boxes"][i]


        shape_grown= [ length_axis + ghosts_axis[0] + ghosts_axis[1] for length_axis , ghosts_axis in zip(box["shape"],box["ghosts"]) ]


        field_data=np.array( self.data[i][name] ).reshape(shape_grown)
        mask=[ slice(ghosts[0],length + ghosts[0]) for length, ghosts in zip(box["shape"],box["ghosts"])   ]
        

        return field_data[tuple(mask)]

    def getGeometry(self, i):
        box=self.description["boxes"][i]

        return gp.geometry(shape=box["shape"],domain=box["domain"],coordinates= self.description["coordinates"])
    @property
    def time(self):
        return self.description["time"]


class gp_simulations:

    def __init__(self,dirname=None):
        self.simulations=[]

        if dirname is not None:
            self.load(dirname)
    def __len__(self):
        return len(self.simulations)

    def load(self,dirname):
        
         sub_dirs=[f for f in listdir(dirname) if not isfile(join(dirname, f))]


         self.simulations = [gp_simulation(join(dirname,subdirname) ) for subdirname in sub_dirs]











        

