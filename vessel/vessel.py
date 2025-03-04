from skyfield.api import wgs84
from datetime import datetime

from ..satellite import Satellite
from ..antenna import LinearAntenna
from ..units import Angle


__all__ = ['Vessel']


class Vessel:
    def __init__(self, latitude: Angle, longitude: Angle):
        self.position = wgs84.latlon(latitude.degrees, longitude.degrees)


    def can_communicate(self, sat: Satellite, antenna: LinearAntenna, time: datetime):
        alt, az, dist = sat.altaz(self.position.latitude, self.position.longitude, time)

        print(alt, az, dist.km)

        visible = alt.degrees > 0

        if visible:
            print('visible')

        az_diff = Angle(radians=(az.radians - antenna.azimuth.radians))
        elev_diff = Angle(radians=(alt.radians - antenna.elevation.radians))
        detectable = antenna.is_input_detectable(sat.power, sat.gain, az_diff, elev_diff, dist)

        if detectable:
            print('detectable')

        return visible and detectable