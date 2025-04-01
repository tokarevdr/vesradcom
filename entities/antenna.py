import abc
from numpy import cos, pi
from scipy.constants import c

from ..units import Angle, Distance, Power, Frequency

__all__ = ['Antenna', 'SmallDipole', 'FiniteLengthDipole']


class Antenna(abc.ABC):
    _frequency: Frequency
    _wavelength: float
    _min_detectable_power: Power
    _max_power: Power
    _power: Power

    def __init__(self, frequency: Frequency, min_detectable_power: Power,
                 max_power: Power, power: Power):
        self._frequency = frequency
        self._wavelength = c / frequency.hz
        self._min_detectable_power = min_detectable_power
        self._max_power = max_power
        self._power = power
    

    @abc.abstractmethod
    def gain(self, theta: float, phi: float) -> float:
        pass
    

    def set_power(self, new_power: Power):
        self._power = new_power


    def power(self) -> Power:
        return self._power
    

    def max_power(self) -> Power:
        return self._max_power
    

    def min_detectable_power(self) -> Power:
        return self._min_detectable_power


    def is_input_detectable(self, p_t: float, g_t: float, azimuth: Angle, elevation: Angle, distance: Distance):
        g_r = self.gain(elevation.radians, azimuth.radians)
        p_r = self.__received_power(p_t, g_t, g_r, distance.m)

        return p_r >= self._min_detectable_power.w
    

    def __received_power(self, p_t: float, g_t: float, g_r: float, d: float):
        '''Friis transmission formula. Euating the power at the terminals of a receive antenna
        as the product of power density of the incident wave and the effective aperture of the receiving antenna
        under idealized conditions given another antenna some distance away transmitting a known amount of power.
        
        Parameters
        ----------
        p_t : float
            The power fed into the transmitting antenna input terminals
        g_t : float
            The transmitting antenna gain
        g_r : float
            The receiving antenna gain
        d : float
            The distance separating the antennas
        
        Returns
        -------
        p_t : float
            The power available at receiving antenna output terminals.
        '''

        return p_t * g_t * g_r * (self._wavelength / (4 * pi * d))**2


class SmallDipole(Antenna):
    def __init__(self, frequency: Frequency, min_detectable_power: Power,
                 max_power: Power, power: Power):
        super().__init__(frequency, min_detectable_power, max_power, power)
    

    def gain(self, theta: float, phi: float):
        return cos(theta)**2
        

class FiniteLengthDipole(Antenna):
    def __init__(self, frequency: Frequency, min_detectable_power: Power,
                 max_power: Power, power: Power):
        super().__init__(frequency, min_detectable_power, max_power, power)
    

    def gain(self, theta: float, phi: float):
        return 0.5 * (cos(theta) * cos(phi) + 1)
        