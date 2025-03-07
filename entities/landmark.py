from skyfield.toposlib import GeographicPosition
from skyfield.api import wgs84
from ..units import Angle, Distance

__all__ = ['Landmark']

class Landmark:
    __title: str
    __position: GeographicPosition
    __bearing: Angle
    __distance: Distance

    def __init__(self, title: str, latitude: Angle, longitude: Angle):
        self.__title = title
        self.__position = wgs84.latlon(latitude.degrees, longitude.degrees)
        self.__bearing = Angle(degrees=0)
        self.__distance = Distance(m=0)


    def title(self) -> str:
        return self.__title


    def position(self) -> GeographicPosition:
        return self.__position
    
    
    def set_bearing(self, new_bearing: Angle):
        self.__bearing = new_bearing


    def bearing(self) -> Angle:
        return self.__bearing
    

    def set_distance(self, new_distance: Distance):
        self.__distance = new_distance


    def distance(self) -> Distance:
        return self.__distance