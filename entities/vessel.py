from skyfield.api import wgs84
from skyfield.toposlib import GeographicPosition

from ..units import Angle, Velocity

__all__ = ['Vessel']


class Vessel:
    __course: Angle
    __position: GeographicPosition
    __max_velocity: Velocity
    __velocity: Velocity

    def __init__(self, latitude: Angle = Angle(radians=0), longitude: Angle = Angle(radians=0),
                 course: Angle = Angle(degrees=0), max_velocity: Velocity = Velocity(km_per_s=0), velocity: Velocity = Velocity(km_per_s=0)):
        self.__course = course
        self.__position = wgs84.latlon(latitude.degrees, longitude.degrees)
        self.__max_velocity = max_velocity
        self.__velocity = velocity


    def set_position(self, latitude: Angle, longitude: Angle):
        self.__position = wgs84.latlon(latitude.degrees, longitude.degrees)


    def position(self) -> GeographicPosition:
        return self.__position
    

    def set_course(self, new_course: Angle):
        self.__course = new_course


    def course(self) -> Angle:
        return self.__course
    

    def max_velocity(self) -> Velocity:
        return self.__max_velocity
    

    def set_velocity(self, velocity: Velocity):
        self.__velocity = velocity


    def velocity(self) -> Velocity:
        return self.__velocity