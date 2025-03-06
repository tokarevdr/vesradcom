from skyfield.api import wgs84
from skyfield.toposlib import GeographicPosition
from datetime import datetime

from . import Satellite
from . import LinearAntenna
from ..units import Angle


__all__ = ['Vessel']


class Vessel:
    __course: Angle
    __position: GeographicPosition
    __antenna: LinearAntenna

    def __init__(self, antenna: LinearAntenna, latitude: Angle = Angle(radians=0), longitude: Angle = Angle(radians=0)):
        self.__antenna = antenna
        self.__position = wgs84.latlon(latitude.degrees, longitude.degrees)


    def set_position(self, latitude: Angle, longitude: Angle):
        self.__position = wgs84.latlon(latitude.degrees, longitude.degrees)


    def position(self) -> GeographicPosition:
        return self.__position
    

    def antenna(self) -> LinearAntenna:
        return self.__antenna


    def can_communicate(self, sat: Satellite, time: datetime):
        alt, az, dist = sat.altaz(self.position().latitude, self.position().longitude, time)

        print(alt, az, dist.km)

        visible = alt.degrees > 0

        if visible:
            print('visible')

        az_diff = Angle(radians=(az.radians - self.__antenna.azimuth().radians))
        alt_diff = Angle(radians=(alt.radians - self.__antenna.altitude().radians))
        detectable = self.__antenna.is_input_detectable(sat.power, sat.gain, az_diff, alt_diff, dist)

        if detectable:
            print('detectable')

        return visible and detectable