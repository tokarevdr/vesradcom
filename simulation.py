from typing import Optional
from datetime import datetime
from shapely import Polygon, Point
from skyfield.toposlib import GeographicPosition
from geopy.distance import geodesic
import math

from .entities import Landmark, Station, Antenna, Vessel, Satellite
from .units import Angle, Velocity, Distance

__all__ = ['Simulation']


class Simulation:
    __current_datetime: datetime
    __area: Optional[Polygon]
    __antenna: Optional[Antenna]
    __vessel: Optional[Vessel]
    __landmarks: list[Landmark]
    __stations: list[Station]
    __satellites: list[Satellite]

    def __init__(self, area: Optional[Polygon] = None, antenna: Optional[Antenna] = None, vessel: Optional[Vessel] = None):
        self.__current_datetime = datetime(2000, 1, 1)
        self.__area = Polygon(area)
        self.__antenna = antenna
        self.__vessel = vessel
        self.__landmarks = []
        self.__stations = []
        self.__satellites = []


    def area(self) -> Polygon:
        return self.__area
    

    def antenna(self) -> Antenna:
        return self.__antenna


    def vessel(self) -> Vessel:
        return self.__vessel
    

    def set_vessel_course(self, new_course: Angle):
        if (self.__vessel != None):
            self.__vessel.set_course(new_course)


    def set_vessel_position(self, latitude: Angle, longitude: Angle):
        if (self.__vessel != None):
            self.__vessel.set_position(latitude, longitude)


    def set_current_datetime(self, new_datetime: datetime):
        self.__current_datetime = new_datetime


    def current_datetime(self) -> datetime:
        return self.__current_datetime


    def append_landmark(self, landmark: Landmark):
        self.__landmarks.append(landmark)

    def landmark_at(self, index: int) -> Landmark:
        return self.__landmarks[index]
    
    def set_landmark_at(self, landmark: Landmark, index: int):
        self.__landmarks[index] = landmark

    def remove_landmark(self, index: int):
        self.__landmarks.pop(index)

    def landmark_count(self):
        return len(self.__landmarks)
    

    def append_station(self, station: Station):
        self.__stations.append(station)

    def station_at(self, index: int) -> Station:
        return self.__stations[index]
    
    def set_station_at(self, station: Station, index: int):
        self.__stations[index] = station

    def remove_station(self, index: int):
        self.__stations.pop(index)

    def station_count(self):
        return len(self.__stations)


    def append_satellite(self, satellite: Satellite):
        self.__satellites.append(satellite)

    def satellite_at(self, index: int) -> Satellite:
        return self.__satellites[index]
    
    def set_satellite_at(self, satellite: Satellite, index: int):
        self.__satellites[index] = satellite

    def remove_satellite(self, index: int):
        self.__satellites.pop(index)

    def satellite_count(self):
        return len(self.__satellites)


    def update(self):
        for satellite in self.__satellites:
            satellite.at(self.__current_datetime)
            alt, az, dist = satellite.altaz(self.__vessel.position().latitude, self.__vessel.position().longitude, self.__current_datetime)
            satellite.set_altitude(alt)
            satellite.set_azimuth(az)
            satellite.set_distance(dist)
            satellite.set_visible(self.__is_satellite_visible(satellite))
            satellite.set_detectable(self.__is_satellite_detectable(satellite))
        
        for landmark in self.__landmarks:
            landmark.set_bearing(self.__landmark_bearing(landmark))
            landmark.set_distance(self.__landmark_distance(landmark))

        for station in self.__stations:
            station.set_bearing(self.__landmark_bearing(station))
            station.set_distance(self.__landmark_distance(station))
            station.set_detectable(self.__is_station_detectable(station))


    def optimize(self, start_time: datetime, end_time: datetime) -> tuple[Angle, Velocity, Angle, Angle]:
        return (Angle(degrees=0), Velocity(km_per_s=0), Angle(degrees=0), Angle(degrees=0))
    

    def __is_satellite_visible(self, satellite: Satellite) -> bool:
        alt = satellite.altitude()

        return alt.degrees > 0
    

    def __is_satellite_detectable(self, satellite: Satellite) -> bool:
        alt = satellite.altitude()
        az = satellite.azimuth()
        dist = satellite.distance()

        az_diff = Angle(radians=(az.radians - self.__vessel.course().radians))
        detectable = self.__antenna.is_input_detectable(satellite.power().w, satellite.gain(), az_diff, alt, dist)

        return detectable
    

    def __is_station_detectable(self, station: Station) -> bool:
        point: Point = Point(self.__vessel.position().longitude.degrees, self.__vessel.position().latitude.degrees)
        return station.coverage_area().contains(point)
    

    def __landmark_bearing(self, landmark: Landmark) -> Angle:
        lat1 = self.__vessel.position().latitude.radians
        lon1 = self.__vessel.position().longitude.radians
        lat2 = landmark.position().latitude.radians
        lon2 = landmark.position().longitude.radians

        dLon = lon2 - lon1
        y = math.sin(dLon) * math.cos(lat2)
        x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dLon)
        bearing = math.atan2(y, x)

        if bearing < 0: bearing += 2*math.pi
        elif bearing > 2*math.pi: bearing -= 2*math.pi

        return Angle(radians=bearing)
    

    def __landmark_distance(self, landmark: Landmark) -> Distance:
        ves_pos = (self.__vessel.position().latitude.degrees, self.__vessel.position().longitude.degrees)
        landmark_pos = (landmark.position().latitude.degrees, landmark.position().longitude.degrees)

        return Distance(m = geodesic(ves_pos, landmark_pos).m)