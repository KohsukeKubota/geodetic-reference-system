import numpy as np
import coordinate as cd

print('=== Example 1 ===')
print('Station: lat = 30N, lon = 130E, height = 800km')
print('Receiver: lat = 35N, lon = 135E, height = 0km')
print('Convert BLH --> ECEF --> ENU --> find IPP at 300km')

# height in km
Sat_blh = cd.blh(30.0, 130.0, 800000.0)
RX_blh = cd.blh(35.0, 135.0, 0.0)

# Convert to ECEF
Sat_ecef = Sat_blh.to_ecef()
RX_ecef = RX_blh.to_ecef()

# Calculate ENU coordinate
# Specify origin by ECEF
Sat_enu = Sat_ecef.to_enu(RX_ecef)
# Specify origin by BLH
Sat_enu2 = Sat_ecef.to_enu(RX_blh)

print('=== Azimuth/deg, Elevation/deg, Distance/km ===')
print('from ECEF:{:.2f}'.format(Sat_enu.az_degree, Sat_enu.el_degree, Sat_enu.distance/1e3))
print('from BLH: {:.2f}'.format(Sat_enu2.az_degree, Sat_enu2.el_degree, Sat_enu2.distance/1e3))
