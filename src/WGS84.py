import math
import numpy as np


# Equator radius[m]
a = 6378137.0
# Reciprocal of oblateness
InversOblateness = 298.257223
# Eccentricity squared
E2 = 2.0 / InversOblateness - (1.0 / InversOblateness) ** 2.0
# Eccentricity
e = np.sqrt(E2)
# oblateness
f = 1.0 / InversOblateness
# Short radius[m]
b = a * (1.0 - f)
# Speed of light [m/s]
c = 2.99792458 * (10 ** 8)
# Earth's rotation rate (rad/s)
omega_e_dot = 7.2921151467E-05
# Earth's gravitational constant[m^3/s^2]
MUe = 3.986005e14
