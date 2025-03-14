from skyfield.api import EarthSatellite, load, wgs84, Timescale
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition
from shapely import Polygon
from datetime import datetime, timedelta
from numpy import cos, arccos, linspace
from geopy.distance import geodesic
from ..units import Angle, Distance, Power

__all__ = ['Satellite']


class Satellite:
    __ts: Timescale
    __skyfield_sat: EarthSatellite
    __title: str
    __position: Geocentric
    __subpoint_position: GeographicPosition
    __coverage_area: Polygon
    __visible: bool
    __detectable: bool
    __altitude: Angle
    __azimuth: Angle
    __distance: Distance
    __power: Power
    __gain: float

    def __init__(self, tle: str, power: Power, gain: float):
        lines = tle.split('\n')
        self.__ts = load.timescale()
        self.__skyfield_sat = EarthSatellite(lines[1].strip(), lines[2].strip(), lines[0].strip())
        self.__title = self.__skyfield_sat.name
        self.__position = Geocentric([0, 0, 0])
        self.__subpoint_position = wgs84.latlon(0, 0)
        self.__coverage_area = Polygon()
        self.__visible = False
        self.__detectable = False
        self.__altitude = Angle(degrees=0)
        self.__azimuth = Angle(degrees=0)
        self.__distance = Distance(m=0)
        self.__power = power
        self.__gain = gain
        

    def at(self, time: datetime):
        t = self.__ts.from_datetime(time)
        geocentric = self.__skyfield_sat.at(t)
        self.__position = geocentric
        self.__distance = geocentric.distance()
        self.__subpoint_position = wgs84.subpoint(geocentric)
        self.__coverage_area = self.__visible_area()


    def track(self, start_datetime: datetime, end_datetime: datetime, interval: timedelta) -> list[GeographicPosition]:
        current_datetime = start_datetime
        datetimes = []        
        while current_datetime <= end_datetime:
            datetimes.append(current_datetime)
            current_datetime += interval

        times = self.__ts.from_datetimes(datetimes)
        geocentric = self.__skyfield_sat.at(times)
        subpoint_positions = wgs84.subpoint(geocentric)

        return subpoint_positions


    def position(self) -> Geocentric:
        return self.__position
    

    def title(self) -> str:
        return self.__title
    

    def subpoint_position(self) -> GeographicPosition:
        return self.__subpoint_position
    

    def coverage_area(self) -> Polygon:
        return self.__coverage_area
    

    def set_visible(self, visible: bool):
        self.__visible = visible


    def visible(self) -> bool:
        return self.__visible
    

    def set_detectable(self, detectable: bool):
        self.__detectable = detectable


    def detectable(self) -> bool:
        return self.__detectable
    

    def set_altitude(self, new_altitude: Angle):
        self.__altitude = new_altitude


    def altitude(self) -> Angle:
        return self.__altitude
    

    def set_azimuth(self, new_azimuth: Angle):
        self.__azimuth = new_azimuth


    def azimuth(self) -> Angle:
        return self.__azimuth
    

    def set_distance(self, new_distance: Distance):
        self.__distance = new_distance


    def distance(self) -> Distance:
        return self.__distance
    

    def power(self) -> Power:
        return self.__power
    

    def gain(self) -> float:
        return self.__gain


    def __visible_area(self) -> Polygon:
        # Угол видимости с учетом минимального угла места
        lat = self.__subpoint_position.latitude.degrees
        lon = self.__subpoint_position.longitude.degrees
        azimuths = linspace(360, 0, 200, endpoint=False)
        eps = 0
        R = wgs84.radius.km
        h = wgs84.height_of(self.__position).km
        theta = arccos(R * cos(eps) / (R + h)) - eps
        
        # Радиус зоны покрытия
        d = R * theta

        boundaries: list[tuple[float, float]] = []
        for azimuth in azimuths:
            point = geodesic(kilometers=d).destination(point=(lat, lon), bearing=azimuth)
            # if point.latitude > 80:
            #     point.latitude = 80
            boundaries.append((point.longitude, point.latitude))

        return Polygon(boundaries)
    

    def altaz(self, observation_point: GeographicPosition, current_datetime: datetime) -> tuple[Angle, Angle, Distance]:
        difference = self.__skyfield_sat - observation_point
        t = self.__ts.from_datetime(current_datetime)
        topocentric = difference.at(t)

        return topocentric.altaz()