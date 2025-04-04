from typing import Optional
from datetime import datetime, timedelta, UTC
from shapely import Polygon, Point
from shapely.ops import transform
import pyproj
from geopy.distance import geodesic
import math
from copy import copy
import numpy as np
import random

from .entities import Landmark, Station, Antenna, Vessel, Satellite
from .units import Angle, Velocity, Distance, Power

__all__ = ['Simulation', 'SimulationResult', 'from_simulation_result']


class SimulationResult:
    time: datetime
    sat_num: int
    point: tuple[Angle, Angle]
    velocity: Velocity
    power: Power
    course: Angle
    intersection: Polygon

    def __init__(self):
        self.time = 0
        self.sat_num = 0
        self.point = (Angle(degrees=0), Angle(degrees=0))
        self.velocity = Velocity(km_per_s=0)
        self.power = Power(w=0)
        self.course = Angle(degrees=0)

    
    def __repr__(self):
        return f'SimulationResult({self.time}, {self.sat_num}, {self.point}, {self.velocity.km_per_s} km/s, {self.power.w} W, {self.course})'
    

    def __str__(self):
        return f'({self.time}, {self.sat_num}, {self.point}, {self.velocity.km_per_s} km/s, {self.power.w} W, {self.course})'
    
    
class Simulation:
    __current_datetime: datetime
    __area: Optional[Polygon]
    __antenna: Optional[Antenna]
    __vessel: Optional[Vessel]
    __landmarks: list[Landmark]
    __stations: list[Station]
    __satellites: list[Satellite]

    def __init__(self, area: Optional[Polygon] = None, antenna: Optional[Antenna] = None, vessel: Optional[Vessel] = None):
        self.__current_datetime = datetime(2000, 1, 1, tzinfo=UTC)
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


    def set_antenna_power(self, power: Power):
        self.__antenna.set_power(power)


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
            alt, az, dist = satellite.altaz(self.__vessel.position(), self.__current_datetime)
            satellite.set_altitude(alt)
            satellite.set_azimuth(az)
            satellite.set_distance(dist)
            satellite.set_visible(self.__is_satellite_visible(satellite))
            satellite.set_detectable(self.received_power_from_satellite(satellite, self.__vessel.course()) >= self.antenna().min_detectable_power().w)
        
        for landmark in self.__landmarks:
            landmark.set_bearing(self.__landmark_bearing(landmark))
            landmark.set_distance(self.__landmark_distance(landmark))

        for station in self.__stations:
            station.set_bearing(self.__landmark_bearing(station))
            station.set_distance(self.__landmark_distance(station))
            station.set_detectable(self.__is_station_detectable(station))

    
    def random_result(self, start_time: datetime, end_time: datetime) -> SimulationResult:
        result = SimulationResult()

        result.time = start_time + random.random() * (end_time - start_time)
        result.sat_num = random.randint(0, self.satellite_count() - 1)

        sat = copy(self.satellite_at(result.sat_num))
        sat.at(result.time)

        points, _ = self.__available_points(sat.coverage_area())

        for point in points:
            start_point = (self.__vessel.position().latitude.degrees, self.__vessel.position().longitude.degrees)
            end_point = (point[1], point[0])
            distance = geodesic(start_point, end_point).km
            time_diff = result.time - self.__current_datetime
            total_secs = time_diff.total_seconds()

            velocity = Velocity(km_per_s = (distance / total_secs))

            if velocity.km_per_s > self.__vessel.max_velocity().km_per_s:
                points.remove(point)

        if not points:
            return SimulationResult()
        
        result.point = random.choice(points)
        result.course = Angle(degrees=(random.random() * 360))

        return result


    def optimize(self, start_time: datetime, end_time: datetime) -> tuple[Angle, Velocity, Angle, Angle]:
        return (Angle(degrees=0), Velocity(km_per_s=0), Angle(degrees=0), Angle(degrees=0))
    

    def bruteforce(self, start_time: datetime, end_time: datetime):
        results: list[SimulationResult] = []
        dt = timedelta(minutes=5)
        dcourse = Angle(degrees=1)
        courses = np.arange(0, 360, dcourse.degrees)

        time = start_time

        while time < end_time:
            for sat_num in range(self.satellite_count()):
                sat = copy(self.satellite_at(sat_num))
                sat.at(time)

                points, intersection = self.__available_points(sat.coverage_area())

                for point in points:
                    start_point = (self.__vessel.position().latitude.degrees, self.__vessel.position().longitude.degrees)
                    end_point = (point[1], point[0])
                    distance = geodesic(start_point, end_point).km
                    time_diff = time - self.__current_datetime
                    total_secs = time_diff.total_seconds()
                    velocity = Velocity(km_per_s=0)

                    if total_secs == 0:
                        is_close_to_current_lon = math.isclose(point[0], self.__vessel.position().longitude.degrees, rel_tol=1e-2)
                        is_close_to_current_lat = math.isclose(point[1], self.__vessel.position().latitude.degrees, rel_tol=1e-2)
                        is_close_to_current_pos = is_close_to_current_lon and is_close_to_current_lat
                        
                        if is_close_to_current_pos:
                            point = (self.__vessel.position().longitude.degrees, self.__vessel.position().latitude.degrees)
                        else:
                            continue
                    else:
                        velocity = Velocity(km_per_s = (distance / total_secs))

                    if velocity.km_per_s <= self.__vessel.max_velocity().km_per_s:
                        for course in courses:
                            result = SimulationResult()
                            result.time = time
                            result.sat_num = sat_num
                            result.point = point
                            result.velocity = velocity

                            detectable = self.received_power_from_satellite(sat, Angle(degrees=course))

                            if detectable:
                                result.course = Angle(degrees=course)
                                result.intersection = intersection
                                results.append(result)

            time += dt

        return results
    

    def __available_points(self, coverage_area: Polygon):
        wgs84 = 'EPSG:4326'
        mollweide = 'ESRI:54030'

        transformer_to_mollweide = pyproj.Transformer.from_crs(wgs84, mollweide, always_xy=True)
        transformer_to_wgs84 = pyproj.Transformer.from_crs(mollweide, wgs84, always_xy=True)

        polygon1 = transform(transformer_to_mollweide.transform, coverage_area)
        polygon2 = transform(transformer_to_mollweide.transform, self.__area)

        polygon3 = polygon1.intersection(polygon2)

        # intersect = coverage_area.intersection(self.__area)

        intersect = transform(transformer_to_wgs84.transform, polygon3)

        if not intersect:
            return ([], Polygon())

        lonmin, latmin, lonmax, latmax = intersect.bounds

        dlon = 0.5
        dlat = 0.5
        X, Y = np.meshgrid(np.arange(lonmin, lonmax, dlon), np.arange(latmin, latmax, dlat))

        points = zip(X.flatten(), Y.flatten())

        return ([i for i in points if intersect.contains(Point(i[0], i[1]))], intersect)
    

    def __is_satellite_visible(self, satellite: Satellite) -> bool:
        alt = satellite.altitude()

        return alt.degrees > 0
    

    def received_power_from_satellite(self, satellite: Satellite, course: Angle) -> bool:
        alt = satellite.altitude()
        az = satellite.azimuth()
        dist = satellite.distance()

        az_diff = Angle(degrees=((az.degrees - course.degrees) % 360))
        power = self.__antenna.received_power(satellite.power().w, satellite.gain(), az_diff, alt, dist)

        return power
    

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
    

def from_simulation_result(original_simulation: Simulation, result: SimulationResult) -> Simulation:
    antenna = copy(original_simulation.antenna())
    antenna.set_power(result.power)
    vessel = copy(original_simulation.vessel())
    vessel.set_position(result.point[0], result.point[1])
    vessel.set_course(result.course)
    vessel.set_velocity(result.velocity)
    sat = copy(original_simulation.satellite_at(result.sat_num))
    
    sim = Simulation(area = copy(original_simulation.area()), antenna=antenna, vessel=vessel)
    sim.append_satellite(sat)
    sim.set_current_datetime(result.time)
    sim.update()

    return sim
