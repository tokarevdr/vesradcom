from skyfield.api import wgs84
from skyfield.toposlib import GeographicPosition
from datetime import datetime

from . import Satellite
from ..units import Angle

__all__ = ['Vessel']


class Vessel:
    __course: Angle
    __position: GeographicPosition

    def __init__(self, latitude: Angle = Angle(radians=0), longitude: Angle = Angle(radians=0), course: Angle = Angle(degrees=0)):
        self.__course = course
        self.__position = wgs84.latlon(latitude.degrees, longitude.degrees)


    def set_position(self, latitude: Angle, longitude: Angle):
        self.__position = wgs84.latlon(latitude.degrees, longitude.degrees)


    def position(self) -> GeographicPosition:
        return self.__position
    

    def set_course(self, new_course: Angle):
        self.__course = new_course


    def course(self) -> Angle:
        return self.__course