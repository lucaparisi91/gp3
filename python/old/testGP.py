import gp
import unittest 
import gp1D_c
import numpy as np
import matplotlib.pylab as plt 
from math import *
import gp3D_c


class testGeometry(unittest.TestCase):

    def test_initGeometry(self):
        shape=(128,128,128)
        domain=((-5,5) , (-5,5) ,(-5,5) )
        geo = gp.geometry(shape,domain)

        self.assertEqual(shape,tuple(geo.shape))
        self.assertEqual(geo.dimensions, len(shape))

        self.assertTrue(geo.positions(0).shape == shape )
        self.assertTrue(geo.positions(1).shape == shape )
        self.assertTrue(geo.positions(2).shape == shape )


class testModel(unittest.TestCase):

    def test_harmonic3D(self):

        settings = { "geometry" : 
            {
                "shape" : [128,128,128],
                "domain" : [[-5,5],[-5,5],[-5,5]],
                "coordinates" : "cartesian"
            },

            "functional" : 
	        {
		        "name" : "harmonic",
		        "omega" : 1.0 , 
                "laplacian" : {
                    "name" : "stencilLaplacian2",
                    "order" : 2
                }
	        }

        }

        geo = gp.geometry(**settings["geometry"])
        r2 =  geo.positions(0)**2 + geo.positions(1)**2 + geo.positions(2)**2 

        alpha=1.
        y= np.exp( - alpha * r2   ) + 0*1j

        hy= y * ( 3*alpha  + r2 * (-2*alpha*alpha + 0.5 ) ) 

        y2=gp3D_c.evaluate(y,settings)

        self.assertAlmostEqual( np.max( np.abs(np.real( y2 - hy))) , 0  , delta=1e-2 )


    def test_harmonicSphericalCart(self):

        settings = { "geometry" : 
            {
                "shape" : [1000],
                "domain" : [[0,5] ],
                "coordinates" : "cartesian",
                "bc" : [ "drichlet" ]
            },
            "functional" : 
	        {
		        "name" : "harmonic",
		        "omega" : 1.0 , 
                "laplacian" : {
                    "name" : "amrexLaplacian",
                    "order" : 2
                }
	        }
        }

        geo = gp.geometry(**settings["geometry"])
        r2 =  geo.positions(0)**2
        
        alpha=1.
        y= geo.positions(0)* np.exp( - alpha * r2   ) + 0*1j

        hy= y * ( 3*alpha  + r2 * (-2*alpha*alpha + 0.5 ) ) 

        y2=gp1D_c.evaluate(y,settings)

        plt.plot(geo.positions(0).flatten() , y2.flatten()   )
        plt.plot(geo.positions(0).flatten() , hy, label="exact")
        plt.legend()
        
        self.assertAlmostEqual( np.max( np.abs(np.real( y2 - hy))) , 0  , delta=1e-3 )
    def test_harmonic1D(self):

        settings = { "geometry" : 
            {
                "shape" : [1000],
                "domain" : [[-5,5] ],
                "coordinates" : "cartesian",
                "bc" : [ "periodic" ]
            },

            "functional" : 
	        {
		        "name" : "harmonic",
		        "omega" : 1.0 , 
                "laplacian" : {
                    "name" : "amrexLaplacian",
                    "order" : 2
                }
	        }
        }

        geo = gp.geometry(**settings["geometry"])
        r2 =  geo.positions(0)**2
        
        alpha=1.
        y= np.exp( - alpha * r2   ) + 0*1j

        hy=y*( alpha  + r2 * (-2*alpha*alpha + 0.5 ) ) 

        y2=gp1D_c.evaluate(y,settings)

        plt.plot(geo.positions(0).flatten() , y2.flatten()   )
        plt.plot(geo.positions(0).flatten() , hy, label="exact")
        plt.legend()
        
        self.assertAlmostEqual( np.max( np.abs(np.real( y2 - hy))) , 0  , delta=1e-3 )

    def test_harmonicSpherical(self):   

        settings = { "geometry" : 
            {
                "shape" : [1000],
                "domain" : [[0,5] ],
                "coordinates" : "spherical",
                "bc" : [ "neumann" ]
            },
            "functional" : 
	        {
		        "name" : "harmonic",
		        "omega" : 1.0 , 
                "laplacian" : {
                    "name" : "stencilLaplacian2",
                    "order" : 2
                }
	        }

        }

        geo = gp.geometry(**settings["geometry"])
        r2 =  geo.positions(0)**2
        
        alpha=1.
        y= np.exp( - alpha * r2   ) + 0*1j

        hy= y * ( 3*alpha  + r2 * (-2*alpha*alpha + 0.5 ) ) 

        y2=gp1D_c.evaluate(y,settings)

        plt.plot(geo.positions(0).flatten() , y2.flatten()   )
        plt.plot(geo.positions(0).flatten() , hy, label="exact")
        plt.legend()
        
        self.assertAlmostEqual( np.max( np.abs(np.real( y2 - hy))) , 0  , delta=1e-3 )



class testRun(unittest.TestCase):

    def test_harmonic3D(self):
        settings = {	
            "geometry" : 
	        {
		    "domain" : [ [-5.0,5.0] , [-5.0,5.0] , [-5.0,5.0] ],
		    "shape" : [128,128,128] , 
            "coordinates" : "cartesian",
            "bc" : [ "periodic", "periodic" , "periodic"]
            },
	        "run" : 
	            {
		        "label" : "gpGroundState",
		        "stepsPerBlock" : 10,
                "nBlocks" : 100,
		        "timeStep" : 1e-3,
		        "stepper" : "RK4",
		        "imaginaryTime" : True
	        },
	        "normalization" : 1.0,
	        "components" : 1 , 
	        "functional" : 
	        {
		        "name" : "harmonic",
		        "omega" : 1.0 , 
                "laplacian" : {
                    "name" : "stencilLaplacian2",
                    "order" : 2
                }
	        }

            }
        
        y=np.zeros(settings["geometry"]["shape"] , dtype=complex )
        geo=gp.geometry(**settings["geometry"])
        alpha = 1
        y+= np.exp(-alpha* ( geo.positions(0)**2 + geo.positions(1)**2 + geo.positions(2)**2) )

        gp3D_c.run(y, settings)

    def test_harmonic1D(self):
        
        settings = { "geometry" : 
            {
                "shape" : [180],
                "domain" : [[-5.0,5.0] ],
                "coordinates" : "cartesian",
                "bc" : [ "periodic" ]
            },
            "functional" : 
	        {
		        "name" : "harmonic",
		        "omega" : 1.0 , 
                "laplacian" : {
                    "name" : "stencilLaplacian2",
                    "order" : 2
                }
	        },
            "run" : 
	            {
		        "label" : "gpGroundState",
		        "stepsPerBlock" : 10000,
                "nBlocks" : 1000,
		        "timeStep" : 1e-4,
		        "stepper" : "RK4",
		        "imaginaryTime" : True
	        },
	        "normalization" : 1.0,
	        "components" : 1 , 
        }

        geo = gp.geometry(**settings["geometry"])
        r2 =  geo.positions(0)**2
        
        alpha=1.;
        y= np.exp( - alpha * r2   ) + 0*1j
        
        gp1D_c.run(y, settings)

    def test_harmonicSpherical_cartesian(self):
        
        settings = { "geometry" : 
            {
                "shape" : [1000],
                "domain" : [ [0,5.0] ],
                "coordinates" : "cartesian",
                "bc" : [ "drichlet" ]
            },
            "functional" : 
	        {
		        "name" : "harmonic",
		        "omega" : 1.0 , 
                "laplacian" : {
                    "name" : "stencilLaplacian2",
                    "order" : 2
                }
	        },
            "run" : 
	            {
		        "label" : "gpGroundState",
		        "stepsPerBlock" : 10000,
                "nBlocks" : 1000,
		        "timeStep" : 1e-5,
		        "stepper" : "RK4",
		        "imaginaryTime" : True
	        },
	        "normalization" : 1.0,
	        "components" : 1 , 
        }

        geo = gp.geometry(**settings["geometry"])
        r2 =  geo.positions(0)**2

        alpha=1.;
        y= geo.positions(0)*np.exp( - alpha * r2   ) + 0*1j

        gp1D_c.run(y, settings)


    def test_harmonicSpherical(self):
        
        settings = { "geometry" : 
            {
                "shape" : [1000],
                "domain" : [ [0,5.0] ],
                "coordinates" : "spherical",
                "bc" : [ "neumann" ]
            },
            "functional" : 
	        {
		        "name" : "harmonic",
		        "omega" : 1.0 , 
                "laplacian" : {
                    "name" : "stencilLaplacian2",
                    "order" : 2
                }
	        },
            "run" : 
	            {
		        "label" : "gpGroundState",
		        "stepsPerBlock" : 10000,
                "nBlocks" : 1000,
		        "timeStep" : 1e-5,
		        "stepper" : "RK4",
		        "imaginaryTime" : True
	        },
	        "normalization" : 1.0,
	        "components" : 1 , 
            "name" : "testHarmonicSpherical"

        }

        geo = gp.geometry(**settings["geometry"])
        r2 =  geo.positions(0)**2

        alpha=1.;
        y= np.exp( - alpha * r2   ) + 0*1j
        gp1D_c.run(y, settings)


    def test_gpDroplet3DSpherical(self):
        
        settings = { "geometry" : 
            {
                "shape" : [1000],
                "domain" : [ [0,50.0] ],
                "coordinates" : "spherical",
                "bc" : [ "neumann" ]
            },
            "functional" : 
	        {
		        "name" : "gpDroplet3D",
                "laplacian" : {
                    "name" : "stencilLaplacian2",
                    "order" : 2
                }
	        },
            "run" : 
	            {
		        "label" : "gpGroundState",
		        "stepsPerBlock" : 10000,
                "nBlocks" : 1000,
		        "timeStep" : 1e-4,
		        "stepper" : "RK4",
		        "imaginaryTime" : True
	        },
	        "normalization" : 19.,
	        "components" : 1 , 
            "name" : "testGPDropletSphericalBins"
        }


        geo = gp.geometry(**settings["geometry"])
        r2 =  geo.positions(0)**2

        alpha=0.1;
        y= np.exp( - alpha * r2   ) + 0*1j
        gp1D_c.run(y, settings)

if __name__ == "__main__":
    unittest.main()    
