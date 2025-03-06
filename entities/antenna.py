import abc
from numpy import sin, cos, ones, pi
from scipy.constants import c

from ..units import Angle, Distance, Power, Frequency


__all__ = ['LinearAntenna', 'FiniteLengthDipole']


class LinearAntenna(abc.ABC):
    _frequency: Frequency
    _wavelength: float
    _altitude: Angle
    _azimuth: Angle
    _minimal_detectable_power: Power

    def __init__(self, frequency: Frequency, altitude: Angle, azimuth: Angle, 
                 minimal_detectable_power: Power):
        self._frequency = frequency
        self._altitude = altitude
        self._azimuth = azimuth
        self._wavelength = c / frequency.hz
        self._minimal_detectable_power = minimal_detectable_power
    

    @abc.abstractmethod
    def gain(self, theta: float, phi: float):
        pass


    def altitude(self) -> Angle:
        return self._altitude
    

    def azimuth(self) -> Angle:
        return self._azimuth


    def is_input_detectable(self, p_t: float, g_t: float, azimuth: Angle, elevation: Angle, distance: Distance):
        g_r = self.gain(elevation.radians, azimuth.radians)
        p_r = self.received_power(p_t, g_t, g_r, distance.m)

        print('gain:', g_r)
        print('power:', p_r)
        print('min power:', self._minimal_detectable_power.w)

        return p_r >= self._minimal_detectable_power.w
    

    def received_power(self, p_t: float, g_t: float, g_r: float, d: float):
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


class SmallDipole(LinearAntenna):
    def __init__(self, frequency: Frequency, altitude: Angle, azimuth: Angle, 
                 minimal_detectable_power: Power):
        super().__init__(frequency, altitude, azimuth,
                         minimal_detectable_power)
    

    def gain(self, theta: float, phi: float):
        if hasattr(phi, 'shape'):
            return cos(theta)**2 * ones(phi.shape)
        else:
            return cos(theta)**2
        

class FiniteLengthDipole(LinearAntenna):
    def __init__(self, frequency: Frequency, altitude: Angle, azimuth: Angle, 
                 minimal_detectable_power: Power):
        super().__init__(frequency, altitude, azimuth,
                         minimal_detectable_power)
    

    def gain(self, theta: float, phi: float):
        return 0.5 * (cos(theta) * cos(phi) + 1)
        