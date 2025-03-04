import abc
from numpy import sin, cos, ones, pi
from scipy.constants import c

from ..units import Angle, Distance


__all__ = ['LinearAntenna', 'FiniteLengthDipole']


def received_power(p_t: float, g_t: float, g_r: float, lam: float, d: float):
    '''Friis transmission formula. Euating the power at the terminals of a receive antenna
    as the product of power density of the incident wave and the effective aperture of the receiving antenna
    under idealized conditions given another antenna some distance away transmitting a known amount of power.
    
    Parameters
    ----------
    p_t : float
          The power fed into the transmitting antenna input terminals
    g_t : float
        The transmittin antenna gain
    g_r : float
        The receiving antenna gain
    lam : float
        Wavelength
    d : float
        The distance separating the antennas
    
    Returns
    -------
    p_t : float
        The power available at receiving antenna output terminals.
    '''

    return p_t * g_t * g_r * (lam / (4 * pi * d))**2


class LinearAntenna(abc.ABC):
    def __init__(self, frequency: float, azimuth: Angle, elevation: Angle,
                 minimal_detectable_power: float):
        self.azimuth = Angle(radians=azimuth.radians)
        self.elevation = Angle(radians=elevation.radians)
        self.frequency = frequency
        self._wavelength = c / frequency
        self.minimal_detectable_power = minimal_detectable_power
    

    @abc.abstractmethod
    def gain(self, theta: float, phi: float):
        pass


    def is_input_detectable(self, p_t: float, g_t: float, azimuth: Angle, elevation: Angle, distance: Distance):
        g_r = self.gain(elevation.radians, azimuth.radians)
        p_r = received_power(p_t, g_t, g_r, self._wavelength, distance.m)

        print('gain:', g_r)
        print('power:', p_r)
        print('min power:', self.minimal_detectable_power)

        return p_r >= self.minimal_detectable_power


class SmallDipole(LinearAntenna):
    def __init__(self, frequency: float, azimuth: Angle, elevation: Angle,
                 minimal_detectable_power: float):
        super().__init__(frequency, azimuth, elevation,
                         minimal_detectable_power)
    

    def gain(self, theta: float, phi: float):
        if hasattr(phi, 'shape'):
            return cos(theta)**2 * ones(phi.shape)
        else:
            return cos(theta)**2
        

class FiniteLengthDipole(LinearAntenna):
    def __init__(self, frequency: float, azimuth: Angle, elevation: Angle,
                 minimal_detectable_power: float):
        super().__init__(frequency, azimuth, elevation,
                         minimal_detectable_power)
    

    def gain(self, theta: float, phi: float):
        return 0.5 * (cos(theta) * cos(phi) + 1)
        