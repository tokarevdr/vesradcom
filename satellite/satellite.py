from skyfield.api import EarthSatellite, load, wgs84
import datetime
from numpy import sin, cos, arccos, linspace, pi, zeros, stack, transpose, arctan, sqrt, array
from numpy.linalg import norm
from geopy.distance import geodesic
from ..units import Angle, Distance

__all__ = ['Satellite']


class Satellite:
    def __init__(self, tle: str, power: float, gain: float):
        lines = tle.split('\n')
        self.ts = load.timescale()
        self.satellite = EarthSatellite(lines[1].strip(), lines[2].strip(), lines[0].strip())
        self.pos = 0
        self.dist = Distance(au=0)
        self.sub_pos = 0
        self.coverage_area = []
        self.power = power
        self.gain = gain
        self.communication_possible = False


    def at(self, time: datetime):
        t = self.ts.from_datetime(time)
        geocentric = self.satellite.at(t)
        self.pos = geocentric
        self.r = geocentric.distance()
        sub_lat, _ = wgs84.latlon_of(geocentric)
        (ra, _, _) = geocentric.radec()
        sub_lon = Angle(radians=(ra.radians - self.__right_ascension_of_greenwich_meridian(time).radians))
        
        self.sub_pos = wgs84.latlon(sub_lat.degrees, sub_lon.degrees)
        h = wgs84.height_of(geocentric).km
        self.coverage_area = self.__coverage_area(h, sub_lat, sub_lon, 100)


    # Satellite Orbits Models, Methods and Applications (Dr. Oliver Montenbruck), p. 33
    def __right_ascension_of_greenwich_meridian(self, t: datetime):
        J2000 = datetime.datetime(2000, 1, 2, 12, tzinfo=datetime.UTC)

        timedelta_since_epoch: datetime.timedelta = t - J2000
        days_since_epoch = timedelta_since_epoch.total_seconds() / 3600 / 24

        degs = (280.4606 + 360.9856473*days_since_epoch) % 2*pi

        return Angle(degrees=degs)


    def __coverage_area(self, h, lat, lon, m):
        # Угол видимости с учетом минимального угла места
        azimuths = linspace(0, 360, m)
        eps = 0
        R = wgs84.radius.km
        theta = arccos(R * cos(eps) / (R + h)) - eps

        # Радиус зоны покрытия
        d = R * theta

        boundary_lat = []
        boundary_lon = []
        for azimuth in azimuths:
            point = geodesic(kilometers=d).destination(point=(lat.degrees, lon.degrees), bearing=azimuth)
            boundary_lat.append(point.latitude)
            boundary_lon.append(point.longitude)

        coverage_area = wgs84.latlon(boundary_lat, boundary_lon)

        return coverage_area
    

    def altaz(self, latitude: Angle, longitude: Angle, time: datetime):
        THETA = self.__right_ascension_of_greenwich_meridian(time).radians

        Rz = zeros((3, 3))
        Rz[0][0] = cos(THETA)
        Rz[1][0] = -sin(THETA)
        Rz[0][1] = sin(THETA)
        Rz[1][1] = cos(THETA)
        Rz[2][2] = 1

        r_ef = Rz @ [[self.pos.position.km[0]], [self.pos.position.km[1]], [self.pos.position.km[2]]]

        R = wgs84.radius.km * array([[cos(latitude.radians) * cos(longitude.radians)], [cos(latitude.radians) * sin(longitude.radians)] ,[sin(latitude.radians)]])
    
        s_ef = r_ef - R

        e_E = [-sin(longitude.radians), cos(longitude.radians), 0]
        e_N = [-sin(latitude.radians) * cos(longitude.radians), -sin(latitude.radians) * sin(longitude.radians), cos(latitude.radians)]
        e_Z = [cos(latitude.radians) * cos(longitude.radians), cos(latitude.radians) * sin(longitude.radians), sin(latitude.radians)]

        E = transpose(stack((e_E, e_N, e_Z), axis=-1))

        s = E @ s_ef

        print(s)

        s_E = s[0][0]
        s_N = s[1][0]
        s_Z = s[2][0]

        A = arctan(s_E / s_N)
        El = arctan(s_Z / sqrt(s_E**2 + s_N**2))
        dist = norm(s)

        return (Angle(radians=El), Angle(radians=A), Distance(km=dist))

    def is_visible(self, pos, time: datetime):
        diff = self.pos - pos
        t = self.ts.from_datetime(time)
        topocentric = diff.at(t)

        alt, _, _ = topocentric.altaz()

        return alt.degrees > 0