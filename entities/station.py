from shapely import Polygon
from skyfield.toposlib import GeographicPosition
from skyfield.api import wgs84
from geopy.distance import geodesic
from numpy import linspace

from .landmark import Landmark
from ..units import Distance, Angle

__all__ = ['Station']

class Station(Landmark):
    __coverage_area: Polygon
    __detectable: bool

    def __init__(self, title: str, latitude: Angle, longitude: Angle, radius: Distance):
        super().__init__(title, latitude, longitude)
        self.__coverage_area = self.__round_polygon(latitude.degrees, longitude.degrees, radius)
        self.__detectable = False


    def coverage_area(self) -> Polygon:
        return self.__coverage_area
    

    def set_detectable(self, detectable: bool):
        self.__detectable = detectable


    def detectable(self) -> bool:
        return self.__detectable


    def __round_polygon(self, latitude_degrees: float, longitude_degrees: float, radius: Distance):
        azimuths = linspace(0, 360, 100)
        boundaries: list[tuple[float, float]] = []
        for azimuth in azimuths:
            point = geodesic(meters=radius.m).destination(point=(latitude_degrees, longitude_degrees), bearing=azimuth)
            boundaries.append((point.longitude, point.latitude))

        return Polygon(boundaries)