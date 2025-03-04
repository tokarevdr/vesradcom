from typing import Optional
from datetime import datetime
from shapely import Polygon

from ..vessel import Vessel
from ..satellite import Satellite
from ..units import Angle, Velocity

__all__ = ['Simulation']


class Simulation:
    __current_datetime: datetime
    __area: Polygon
    __vessel: Optional[Vessel]
    __satellites: list[Satellite]

    def __init__(self):
        self.__current_datetime = datetime(2000, 1, 1)
        self.__area = Polygon()
        self.__vessel = None
        self.__satellites = list()


    def set_area(self, area: Polygon):
        self.__area = area


    def set_vessel(self, vessel: Vessel):
        self.__vessel = vessel


    def append_satellite(self, satellite: Satellite):
        self.__satellites.append(satellite)


    def remove_satellite(self, index: int):
        self.__satellites.pop(index)


    def satellites_count(self):
        return len(self.__satellites)
    

    def satellite_at(self, index: int) -> Satellite:
        return self.__satellites[index]


    def set_current_datetime(self, current_datetime: datetime):
        self.__current_datetime = current_datetime

        for satellite in self.__satellites:
            satellite.at(current_datetime)
            satellite.communication_possible = self.__vessel.can_communicate(satellite, self.__current_datetime)


    def optimize(self, start_time: datetime, end_time: datetime) -> tuple[Angle, Velocity, Angle, Angle]:
        return (Angle(degrees=0), Velocity(km_per_s=0), Angle(degrees=0), Angle(degrees=0))