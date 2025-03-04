from skyfield.units import Angle, Distance, Velocity

__all__ = ['Angle', 'Distance', 'Velocity', 'Power', 'Frequency']


class Power:
    w: float

    def __init__(self, w: float):
        self.w = w


class Frequency:
    hz: float

    def __init__(self, hz: float):
        self.hz = hz