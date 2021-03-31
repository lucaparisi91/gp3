
import json
import numpy as np
import gp
import copy
import argparse

def createGaussianInitialCondition(shape,domain):
    geo=gp.geometry(shape=shape,domain=domain)
    r2=geo.positions(0)**2 + geo.positions(1)**2 + geo.positions(2)**2 
    y=np.exp(-r2)
    wave=gp.wavefunction(geo,[y] ,name="phi")
    return wave 


def createJsonInputFile(jTemplate,shape,domain,nProcessors):

    
    j=copy.deepcopy(jTemplate)

    
    j["geometry"]["shape"]=shape
    j["geometry"]["domain"]=domain

    maxGridSize=shape[0]/nProcessors

    assert( maxGridSize.is_integer() )

    j["maxGridSize"]=maxGridSize

    return j



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Create test input files for a run on n processors')
    parser.add_argument('nProcessors', type=int )
    args = parser.parse_args()
    nProcessors=args.nProcessors


    shape=[400,400,400]
    domain= [ [-5,5] , [-5,5] , [-5,5] ]

    wave=createGaussianInitialCondition(shape,domain)
    wave.save("gaussian_alpha1")

    with open("input-template.json") as f:
        jTemplate=json.load(f)
    
    j=createJsonInputFile(jTemplate,shape=shape,domain=domain,nProcessors=nProcessors)

    with open("input.json","w") as f:
        json.dump(j,f,indent=4)


    
    
