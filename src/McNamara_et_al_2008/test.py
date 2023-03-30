#!/usr/bin/env python3

import unittest
import mcnamara_et_al_2008

class TestMcNamara2008(unittest.TestCase):

    # test whether without any mutation you get 0 outcomes
    def test_init_zero(self):

        sim = mcnamara_et_al_2008.McNamara2008(muCoop = 0.0, muChoice = 0)

        self.assertEqual(sim.means(),(0.0,0.0))


    def test_init_half(self):

        sim2 = mcnamara_et_al_2008.McNamara2008(
                muCoop = 0.0, 
                muChoice = 0, 
                initvals = (0.5,0.5))

        self.assertEqual(sim2.means(),(0.5,0.5))
    
    def test_init_unequal(self):

        sim = mcnamara_et_al_2008.McNamara2008(
                muCoop = 0.0, 
                muChoice = 0, 
                initvals = (0,0.5))

        self.assertEqual(sim.means(),(0.0,0.5))


    # test for the pair formation function which is just u * u / sum_u
    def test_pair_formation(self):
        
        sim = mcnamara_et_al_2008.McNamara2008(
                muCoop = 0.0, 
                muChoice = 0, 
                initvals = (0,0))

        sim.pairing()



        self.assertEqual(sim.Q



if __name__ == '__main__':
    unittest.main()
