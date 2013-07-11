import unittest
import numpy as np
import trm.subs as subs

class TestFunctions(unittest.TestCase):

    def test_planck(self):
        # a particular value
        self.assertEqual(subs.planck(500.,6000.), 31749069741703.254)

        # Sending in two arrays should force an exception
        self.assertRaises(subs.SubsError, subs.planck, \
                              np.array([500.,600.]), np.array([5000.,6000.]))

"""
License
=======

Copyright 2011 Tom Marsh
All Rights Reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

if __name__ == "__main__":
    unittest.main()
