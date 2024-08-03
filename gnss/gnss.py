import math
from typing import Tuple
import numpy as np
 
class GNSS:

    def __init__(self):
        """
        Initialize the GNSS class with the semi-major and semi-minor axis of the Earth ellipsoid.
        """
        # Semi-major axis of the Earth ellipsoid (equatorial radius) in meters
        self.a = 6378137
        # Semi-minor axis of the Earth ellipsoid (polar radius) in meters
        self.b = 6356752.3142
        # Flattening of the Earth ellipsoid
        self.f = (self.a - self.b) / self.a
        # First eccentricity of the Earth ellipsoid
        self.e = (self.f * (2 - self.f)) ** 2
    def geodetic_to_ecef(self: any, lat: float, lon: float, h: float) -> Tuple[float, float, float]:
        """
        Converts geodetic coordinates (latitude, longitude, altitude) to Earth-Centered Earth-Fixed (ECEF) coordinates.

        Args:
            lat (float): Latitude in degrees.
            lon (float): Longitude in degrees.
            h (float): Altitude in meters.

        Returns:
            Tuple[float, float, float]: ECEF coordinates (x, y, z) in meters.
        """
        # Convert (lat, lon) to radians
        lamb = math.radians(lat)
        phi = math.radians(lon)
        # Calculate sine and cosine of lamb and phi
        s = math.sin(lamb)
        N = self.a / math.sqrt(1 - self.e* s * s)
    
        sin_lambda = math.sin(lamb)
        cos_lambda = math.cos(lamb)
        sin_phi = math.sin(phi)
        cos_phi = math.cos(phi)
    
        # Calculate ECEF coordinates
        x = (h + N) * cos_lambda * cos_phi
        y = (h + N) * cos_lambda * sin_phi
        z = (h + (1 - self.e) * N) * sin_lambda
    
        return x, y, z
    
    def ecef_to_enu(self, x: float, y: float, z: float, lat0: float, lon0: float, h0: float) -> Tuple[float, float, float]:
        """
        Converts Earth-Centered Earth-Fixed (ECEF) coordinates to local tangent plane coordinates (ENU).

        :param x: ECEF x-coordinate in meters
        :param y: ECEF y-coordinate in meters
        :param z: ECEF z-coordinate in meters
        :param lat0: reference latitude in degrees
        :param lon0: reference longitude in degrees
        :param h0: reference altitude in meters
        :return: ENU coordinates (east, north, up) in meters
        """
        # Convert reference coordinates to radians
        lamb = math.radians(lat0)
        phi = math.radians(lon0)
        # Calculate sine and cosine of lamb and phi
        s = math.sin(lamb)
        N = self.a / math.sqrt(1 - self.e * s * s)

        sin_lambda = math.sin(lamb)
        cos_lambda = math.cos(lamb)
        sin_phi = math.sin(phi)
        cos_phi = math.cos(phi)

        # Calculate ECEF position of reference point
        x0 = (h0 + N) * cos_lambda * cos_phi
        y0 = (h0 + N) * cos_lambda * sin_phi
        z0 = (h0 + (1 - self.e) * N) * sin_lambda

        # Calculate ECEF coordinates relative to reference point
        xd = x - x0
        yd = y - y0
        zd = z - z0

        # Calculate ENU coordinates
        # Note: The ENU coordinate system is defined in the local tangent plane
        #       at the reference point, which is aligned with the local east,
        #       north, and up directions.
        xEast = -sin_phi * xd + cos_phi * yd  # east
        yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd  # north
        zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd  # up

        return xEast, yNorth, zUp
    
    def geodetic_to_enu(self, lat, lon, h, lat_ref, lon_ref, h_ref) -> Tuple[float, float, float]:
        """
        Converts geodetic coordinates to local tangent plane coordinates (ENU).

        :param lat: geodetic latitude in degrees
        :param lon: geodetic longitude in degrees
        :param h: geodetic altitude in meters
        :param lat_ref: reference latitude in degrees
        :param lon_ref: reference longitude in degrees
        :param h_ref: reference altitude in meters
        :return: tuple of ENU coordinates (x, y, z) in meters
        """
        # Convert geodetic to ECEF coordinates
        x, y, z = self.geodetic_to_ecef(lat, lon, h)

        # Convert ECEF to ENU coordinates
        return self.ecef_to_enu(x, y, z, lat_ref, lon_ref, h_ref)
    
    @staticmethod
    def Relative_Distance(x, y, z, x0, y0, z0):
        return math.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)