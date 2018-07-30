# coding:utf-8
import numpy as np
import WGS84 as wgs84
import vincenty as vin
from copy import deepcopy


class h_distance:
    """
    Calculate horizontal distance in meters between point1 and point2
    along the Earth's ellipsoids.  Points are described by BLH coordinates.
    Height of the point1 is considered in the calculation.
    Earth's parameter is WGS84.
    Calculation is based on Vincenty.
    """

    def __init__(self, point1=None, point2=None):
        if isinstance(point1, blh) and isinstance(point2, blh):
            self.lat1 = point1.B
            self.lon1 = point1.L
            self.height = point1.H
            self.lat2 = point2.B
            self.lon2 = point2.L
            # Set height above wgs84 ellipsoid
            vin.a = (wgs84.a + self.height)
            vin.b = vin.a * (1.0 - wgs84.f)
            self.h_distance = vin.vincenty((self.lat1, self.lon1), (self.lat2, self.lon2)) * 1000.0
        else:
            self.lat1 = 0.0
            self.lon1 = 0.0
            self.height = 0.0
            self.lat2 = 0.0
            self.lon2 = 0.0
            self.h_distance = 0.0

    def __str__(self):
        return "{},{}".format(self.h_distance, self.height)


class enu:
    def __init__(self, e=0.0, n=0.0, u=0.0, origin=None):
        self.wgs84 = wgs84
        self.e = e
        self.n = n
        self.u = u
        self.origin = origin

    def __str__(self):
        return "{},{},{}".format(self.e, self.n, self.u)

    def AED2enu(self, AZ=0.0, EL=0.0, DIST=0.0):
        AZ_RAD = AZ * np.pi / 180.
        EL_RAD = EL * np.pi / 180.
        self.e = DIST * np.sin(AZ_RAD) * np.cos(EL_RAD)
        self.n = DIST * np.cos(AZ_RAD) * np.cos(EL_RAD)
        self.u = DIST * np.sin(EL_RAD)

    def to_ecef(self, shift=None):
        if isinstance(self.origin, blh):
            self.origin = self.origin.to_ecef()
        elif isinstance(self.origin, ecef):
            pass
        else:
            return None
        lat, lon = self.origin.to_blh().pos[:2] * np.pi / 180.0
        mx = np.array([
            -np.sin(lon), -np.sin(lat) * np.cos(lon), np.cos(lat) * np.cos(lon), np.cos(lon),
                          -np.sin(lat) * np.sin(lon), np.cos(lat) * np.sin(lon), 0, np.cos(lat), np.sin(lat)])
        mx = mx.reshape(3, 3)
        enu = self.pos
        geo = np.zeros(3)
        if shift:
            geo = self.origin.pos
            res = mx.dot(enu) + geo
        res = mx.dot(enu)
        return ecef(res[0], res[1], res[2])

    @property
    def el_rad(self):
        return np.arctan2(self.u, np.sqrt(self.e ** 2 + self.n ** 2))

    @property
    def el_degree(self):
        return np.rad2deg(self.el_rad)

    @property
    def az_rad(self):
        theta = np.arctan2(self.n, self.e)
        if isinstance(theta, np.ndarray):
            az = np.ones(theta.shape)
            az[theta < 0] = np.pi/2 - theta[theta < 0]
            az[theta < np.pi/2] = np.pi/2 - theta[theta < np.pi/2]
            az[theta == np.pi/2] = theta[theta == np.pi/2] - np.pi
            az[theta > np.pi/2] = 2*np.pi + (np.pi/2 - theta[theta > np.pi/2])
        else:
            if theta < 0:
                print(theta)
                az = np.pi/2 - theta
            elif theta < np.pi/2:
                az = np.pi/2 - theta
            elif theta == np.pi/2:
                az = theta - np.pi/2
            else:
                az = 2 * np.pi + (np.pi / 2 - theta)
        return az

    @property
    def az_degree(self):
        return np.rad2deg(self.az_rad)

    @property
    def distance(self):
        return np.sqrt(self.e ** 2 + self.n ** 2 + self.u ** 2)

    @property
    def pos(self):
        return np.array([self.e, self.n, self.u])


class ecef:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.wgs84 = wgs84
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "{},{},{}".format(self.x, self.y, self.z)

    def __add__(self, other):
        x, y, z = self.pos + other.pos
        return ecef(x, y, z)

    def __sub__(self, other):
        x, y, z = self.pos - other.pos
        return ecef(x, y, z)

    def __mul__(self, other):
        x, y, z = self.pos * other
        return ecef(x, y, z)

    def to_blh(self, unit="degree"):
        h = self.wgs84.a ** 2.0 - self.wgs84.b ** 2.0
        p = np.hypot(self.x, self.y)
        t = np.arctan2(self.z * self.wgs84.a, p * self.wgs84.b)
        sint = np.sin(t)
        cost = np.cos(t)

        lat = np.arctan2(self.z + h / self.wgs84.b * (sint ** 3.0), p - h / self.wgs84.a * (cost ** 3.0))  # 緯度[rad]を計算
        n = self.wgs84.a / np.sqrt(1.0 - self.wgs84.E2 * (np.sin(lat) ** 2.0))  # 卯酉線曲率半径
        lon = np.arctan2(self.y, self.x)
        height = (p / np.cos(lat)) - n

        lat = lat * 180.0 / np.pi
        lon = lon * 180.0 / np.pi
        return blh(lat, lon, height, unit)

    def to_enu(self, origin: object) -> object:
        ans = None
        if isinstance(origin, ecef) or isinstance(origin, blh):
            origin_blh = None
            if isinstance(origin, blh):
                origin_blh = blh(origin.B, origin.L, origin.H)
                origin_blh.change_unit_to_rad()
                origin = origin_blh.to_ecef()
            else:
                origin_blh = origin.to_blh()
                origin_blh.change_unit_to_rad()

            dx = self.x - origin.x
            dy = self.y - origin.y
            dz = self.z - origin.z
            sB = np.sin(origin_blh.B)
            cB = np.cos(origin_blh.B)
            sL = np.sin(origin_blh.L)
            cL = np.cos(origin_blh.L)
            e = -dx * sL + dy * cL
            n = -dx * cL * sB + -dy * sL * sB + dz * cB
            u = dx * cL * cB + dy * sL * cB + dz * sB
            ans = enu(e, n, u, origin)
        return ans

    @property
    def pos(self):
        return np.array([self.x, self.y, self.z])


class blh:
    def __init__(self, B=0.0, L=0.0, H=0.0, unit="degree"):
        self.wgs84 = wgs84
        # latitude
        self.B = B
        # longtitude
        self.L = L
        # ellipsoid height[m]
        self.H = H
        # unit
        self.unit = unit

    def __str__(self):
        return "{},{},{}".format(self.B, self.L, self.H)

    def __array__(self):
        return self.B, self.L, self.H

    def copy(self):
        return blh(self.B, self.L, self.H, self.unit)

    def change_unit_to_degree(self):

        if self.unit == "rad":
            self.B *= 180.0 / np.pi
            self.L *= 180.0 / np.pi
            self.unit = "degree"

    def change_unit_to_rad(self):
        if self.unit == "degree":
            self.B *= np.pi / 180.0
            self.L *= np.pi / 180.0
            self.unit = "rad"

    def to_ecef(self):
        copy = blh(deepcopy(self.B), deepcopy(self.L), deepcopy(self.H), deepcopy(self.unit))
        copy.change_unit_to_rad()
        n = self.wgs84.a / np.sqrt(1.0 - self.wgs84.E2 * np.sin(copy.B) ** 2.0)
        x = (n + copy.H) * np.cos(copy.B) * np.cos(copy.L)
        y = (n + copy.H) * np.cos(copy.B) * np.sin(copy.L)
        z = ((1.0 - self.wgs84.E2) * n + copy.H) * np.sin(copy.B)
        return ecef(x, y, z)

    @property
    def pos(self):
        return np.array([self.B, self.L, self.H])


def main():
    print("---self test---")
    hoge = ecef(-19163561.617619168, 5514390.2974610282, 17332227.130933248).to_blh()
    print((hoge.B, hoge.L, hoge.H))
    print(blh(24.3453, 124.1587, 20.0).to_ecef())


if __name__ == '__main__':
    main()
