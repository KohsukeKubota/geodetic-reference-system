# geodetic-reference-system
We can calculate horizontal distance in meters between point1 and point2
along the Earth's ellipsoids.  Points are described by BLH coordinates(latitude, longtitude, and elipsoid height).
Height of the point1 is considered in the calculation.
Earth's parameter is WGS84.
Calculation is based on Vincenty.

# Geodetic reference system?
Geodetic reference system is the system which can representing the position on the earth by coordinates using longtitude, latitude and altitude.
This script use WGS84, which is famous for used in GPS.

# Requirement
It is assumed that users use numpy.
This script need Vincenty.We can calculate the distance between two GPS points.
We show how to install Vincenty in the below.
```bash
pip install vincenty
```

# Description
We can calculate horizontal distance in meters between point1 and point2
along the Earth's ellipsoids.  Points are described by BLH coordinates.
Height of the point1 is considered in the calculation.
Earth's parameter is WGS84.
Calculation is based on Vincenty.

We can get not only the distance between two points but also the azimuth and elevation of one point from the other point.
