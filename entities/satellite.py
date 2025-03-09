from skyfield.api import EarthSatellite, load, wgs84, Timescale
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition
from shapely import Polygon
import datetime
from numpy import sin, cos, arccos, linspace, pi, zeros, stack, transpose, arctan, sqrt, array
from numpy.linalg import norm
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
        sub_lat, sub_lon = self.__subpoint_latlon(geocentric, time)
        
        self.__subpoint_position = wgs84.latlon(sub_lat.degrees, sub_lon.degrees)
        self.__coverage_area = self.__visible_area()


    def track(self, start_time: datetime, end_time: datetime, point_count: int) -> tuple[list[float], list[float]]:
        lats = []
        lons = []
        time = start_time
        dt = (end_time - start_time) / point_count

        while time < end_time:
            t = self.__ts.from_datetime(time)
            geocentric = self.__skyfield_sat.at(t)
            sub_lat, sub_lon = self.__subpoint_latlon(geocentric, time)
            lats.append(sub_lat.degrees)
            lons.append(sub_lon.degrees)
            time += dt

        lats.append(None)
        lons.append(None)

        return (lons, lats)


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


    # Satellite Orbits Models, Methods and Applications (Dr. Oliver Montenbruck), p. 33
    def __right_ascension_of_greenwich_meridian(self, t: datetime) -> Angle:
        J2000 = datetime.datetime(2000, 1, 2, 12, tzinfo=datetime.UTC)

        timedelta_since_epoch: datetime.timedelta = t - J2000
        days_since_epoch = timedelta_since_epoch.total_seconds() / 3600 / 24

        degs = (280.4606 + 360.9856473*days_since_epoch) % 2*pi

        return Angle(degrees=degs)


    def __visible_area(self) -> Polygon:
        # Угол видимости с учетом минимального угла места
        lat = self.__subpoint_position.latitude.degrees
        lon = self.__subpoint_position.longitude.degrees
        azimuths = linspace(0, 360, 100)
        eps = 0
        R = wgs84.radius.km
        h = wgs84.height_of(self.__position).km
        theta = arccos(R * cos(eps) / (R + h)) - eps
        
        # Радиус зоны покрытия
        d = R * theta

        boundaries: list[tuple[float, float]] = []
        for azimuth in azimuths:
            point = geodesic(kilometers=d).destination(point=(lat, lon), bearing=azimuth)
            boundaries.append((point.longitude, point.latitude))

        return Polygon(boundaries)
    

    def __subpoint_latlon(self, geocentric: Geocentric, time: datetime) -> tuple[Angle, Angle]:
        sub_lat, _ = wgs84.latlon_of(geocentric)
        (ra, _, _) = geocentric.radec()
        sub_lon = Angle(radians=(ra.radians - self.__right_ascension_of_greenwich_meridian(time).radians))

        return (sub_lat, sub_lon)
    

    def altaz(self, latitude: Angle, longitude: Angle, time: datetime) -> tuple[Angle, Angle, Distance]:
        THETA = self.__right_ascension_of_greenwich_meridian(time).radians

        Rz = zeros((3, 3))
        Rz[0][0] = cos(THETA)
        Rz[1][0] = -sin(THETA)
        Rz[0][1] = sin(THETA)
        Rz[1][1] = cos(THETA)
        Rz[2][2] = 1

        r_ef = Rz @ [[self.__position.position.km[0]], [self.__position.position.km[1]], [self.__position.position.km[2]]]

        R = wgs84.radius.km * array([[cos(latitude.radians) * cos(longitude.radians)], [cos(latitude.radians) * sin(longitude.radians)] ,[sin(latitude.radians)]])
    
        s_ef = r_ef - R

        e_E = [-sin(longitude.radians), cos(longitude.radians), 0]
        e_N = [-sin(latitude.radians) * cos(longitude.radians), -sin(latitude.radians) * sin(longitude.radians), cos(latitude.radians)]
        e_Z = [cos(latitude.radians) * cos(longitude.radians), cos(latitude.radians) * sin(longitude.radians), sin(latitude.radians)]

        E = transpose(stack((e_E, e_N, e_Z), axis=-1))

        s = E @ s_ef

        s_E = s[0][0]
        s_N = s[1][0]
        s_Z = s[2][0]

        A = arctan(s_E / s_N)
        El = arctan(s_Z / sqrt(s_E**2 + s_N**2))
        dist = norm(s)

        return (Angle(radians=El), Angle(radians=A), Distance(km=dist))