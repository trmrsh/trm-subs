import unittest
import numpy as np
import matplotlib.pyplot as plt
import trm.subs.dvect as dvect

class TestDvect(unittest.TestCase):

    def test_plot(self):
        x  = np.linspace(0.,10.,11)
        y  = x*x
        mx = abs(x-5) < 0.1
        e  = 5.*np.ones_like(y)
        dx = dvect.Dvect(x, name='X')
        dy = dvect.Dvect(y, err=e, mask=mx, name='Y')

        dvect.errorbar(plt, dx, dy)
        print ("The plot should show a parabola from 0 to 10. Click the exit if OK")
        plt.show()


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
