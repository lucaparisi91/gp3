import gp
import unittest 
import gp_c
import numpy as np


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

    def test_evaluatePython(self):
        shape=(128,128,128)
        domain=((-5,5) , (-5,5) ,(-5,5) )
        geo = gp.geometry(shape,domain)

        y= np.exp( - geo.positions(0)**2 - geo.positions(1)**2 - geo.positions(2)**2 )

        gp_c.evaluate(y,y*0,geo)

        




    


if __name__ == "__main__":
    unittest.main()
    
