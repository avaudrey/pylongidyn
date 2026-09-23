#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
pylongidyn: Longitudinal Vehicle Dynamics Modeling in Python
============================================================

This package provides tools to model and analyze the longitudinal dynamics
of electric road vehicles. It includes representations of ambient conditions,
road profiles, continuous gearboxes and electric motors, enabling simulation
of vehicle performance under realistic driving scenarios.

Main Components
---------------
- rot_rpm_from_speed, rot_speed_from_rpm : Conversion utilities between rpm
and rad/s.
- AmbientAir : Ambient air properties and density calculation, useful for drag
force calculation.
- Road : Road profile handling (slopes, elevations, speed limits).
- LinearContinuousGearbox : Gearbox model with continuously variating gear
ratio.
- SynchronousElectricMotor : Electric motor model with torque/power curves.
- Vehicle : Vehicle body (mass, rotating mass correction, frontal area, drag
coefficient, wheels radius), providing the resistive forces (climbing,
rolling, aerodynamic drag), the motive force needed for a given acceleration
and, conversely, the acceleration reachable with a given motive force.
- VehicleDynamicsModel : Time-stepping simulation (explicit Euler scheme) of
a vehicle following a road profile. The vehicle tries to drive at the speed
limit, using target acceleration and deceleration values, within the torque
and power limits of the electric motor. Braking is shared between the motor
(regenerative braking) and the friction brakes. Results are time series of
speed, acceleration, electrical power, motor operating points and friction
braking power.

Notes
-----
- Units are SI (Pa, °C, rpm, rad/s, m/s, N·m, W).
- Road data must be provided as a CSV file with columns:
  distance_m, altitude_m, speed_limit_kmh
- Stops (rows with a zero speed limit) are currently removed from the road
  profile, and a warning is issued.
- Electrical powers are counted positive when consumed and negative when
  recovered (regenerative braking).
- More information at : https://github.com/avaudrey/pylongidyn

Author: alexandre.vaudrey@pm.me
"""

# ---------------------------------------------------------------------------
#   Copyright (C) 2025 <Alexandre Vaudrey>                                  |
#                                                                           |
#   This program is free software: you can redistribute it and/or modify    |
#   it under the terms of the GNU General Public License as published by    |
#   the Free Software Foundation, either version 3 of the License, or       |
#   (at your option) any later version.                                     |
#                                                                           |
#   This program is distributed in the hope that it will be useful,         |
#   but WITHOUT ANY WARRANTY; without even the implied warranty of          |
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           |
#   GNU General Public License for more details.                            |
#                                                                           |
#   You should have received a copy of the GNU General Public License       |
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.   |
# --------------------------------------------------------------------------|

# NOTE: List of possible/future improvements
# - Taking into account the humidity of air in the calculation of its density.
# - Consider an energy efficiency for the electric motor, that could be a class
# function, given for instance as an empirical function of the torque and speed
# - The vehicle rolling coefficient must be a function of its speed.

import csv
import datetime as dt
import warnings
import numpy as np

# Importation of physical constants required
from scipy.constants import g, R, zero_Celsius
import scipy.optimize as opt

# Molecular weight of air, in g/mol
M_AIR = 29.0


def rot_speed_from_rpm(rot_freq):
    """Conversion of rotational frequency, in rpm, into angular speed,
    in rad/s.
    """
    return rot_freq * np.pi / 30.0


def rot_rpm_from_speed(rot_speed):
    """Conversion of angular speed, in rad/s, into rotational frequency,
    in rpm.
    """
    return rot_speed * 30.0 / np.pi


class AmbientAir:
    """
    Properties of the air into which the vehicle is moving
    ...
    Attributes
    ----------
    pressure : float
        Expressed in pascals (Pa)
    temperature : float
        Expressed in degrees celcius (°C)

    Methods
    -------
    density() : float
        Calculated using the ideal gas law and expressed in kg/m³
    """

    def __init__(self, pressure=1e5, temperature=20.0):
        self.pressure = pressure
        self.temperature = temperature

    def density(self):
        """
        Air density expressed in kg/m³
        """
        # Specific gas constant of air
        gas_constant = R / (M_AIR * 1e-3)
        return self.pressure / (gas_constant * (self.temperature + zero_Celsius))


class Road:
    """
    Road profile the vehicle has to follow.
    ...
    Attributes
    ----------
    max_slope : float
        Maximum value of the road slope used for the regularization of the road
        profile, dimensionless.
    regularization : boolean
        If True, the (positive or negative) slopes of the initial road profile
        are limited by the 'max_slope' value, and the corresponding elevation
        is modified as a consequence.

    Methods
    -------
    climbing_angle(x) : float
        Climbing angle of the road (in °/degrees) at the distance 'x', in m,
        from the starting point.
    elevation(x) : float
        Elevation (in m) of the road at the distance 'x' (in m) from the
        starting point. Each piece of road is supposed to be a linear function.
    export_regularized_profile()
        Road profile with regularized slope is exported to another csv file
        with a slightly modified name but the same structure.
    load_road_data('file_name.csv')
        'file_name.csv' is a csv file composed of three columns, the first
        one with the distance from the starting point (in meters), the seconds
        with the altitude from the sea level (in meters) and the third with the
        speed limit in km/h.
    maximum_authorized_speed() : float
        Maximum value of the authorized speed along the whole path, in m/s.
    maximum_average_speed() : float
        Maximum average speed, in m/s, if the vehicle is at any time at the
        maximum speed authorized.
    maximum_elevation() : float
        Maximum elevation along the whole road, in m.
    minimum_duration() : float
        Minimum duration of the trip, in s, if the vehicle is at any time at the
        maximum speed authorized.
    road_length() : float
        Total length of the road, in m.
    slope(x) : float
        Slope of the road (dimensionless) at the distance 'x' (in m) from the
        starting point.
    speed_limit(x) : float
        Speed limit (in m/s) the vehicle has to respect at the distance 'x'
        (in m) from the starting point.
    speed_statistics : dict
        Reply with a dictionary whom keys are the different levels of speed
        limit along the road and the corresponding values are the sum of all
        road section length corresponding to these speed limits.
    """

    def __init__(self):
        # Name of the file from which solar data are loaded
        self.__data_file_name = "file"
        # Is a data file loaded?
        self.__is_a_data_file_loaded = False
        # List of distance values from the starting point, as np.array
        self.__distance_from_start = np.array([])
        # List of elevation values from the starting point, as np.array
        self.__elevation_from_start = np.array([])
        # List of speed limit values from the starting point, as np.array
        self.__speed_limit = np.array([])
        # Initial value of the road altitude, that must be memorized
        self.__initial_altitude = 0.0
        # Values of the road slope
        self.__road_slope = np.array([])
        # The slope of the road profile is regularize, i.e. limited at the
        # 'max_slope' value if needed
        self.regularization = False
        # Maximum value of the road slope when regularization is applied, here
        # 15% at most
        self.max_slope = 0.15

    # ==== Methods ============================================================

    def load_road_data(self, file_name, delimiter=","):
        """Load the data file containing information about road to follow. For
        the moment, the row separator is a comma ','
        """
        if not isinstance(file_name, str):
            raise TypeError("Path to the data file must be a string!")
        # File opening (first line is avoided). A missing file raises a
        # FileNotFoundError, which is left to the user
        data = np.loadtxt(file_name, skiprows=1, delimiter=delimiter)
        self.__data_file_name = file_name
        self.__is_a_data_file_loaded = True
        # Extraction of the data, in three lists. To avoid any problem of
        # division by zero, any line containing a zero speed is avoided
        mask = data[:, 2] != 0
        # The user is warned, since stops are then ignored
        if not mask.all():
            warnings.warn(
                f"{int((~mask).sum())} line(s) with a zero speed "
                "limit (stops) have been removed from the road "
                "profile."
            )
        # Distance as an array
        self.__distance_from_start = data[:, 0][mask]
        # Transformation of altitude into relative elevation regarding to
        # altitude starting point
        altitude = data[:, 1][mask]
        self.__initial_altitude = altitude[0]
        self.__elevation_from_start = altitude - altitude[0]
        # If needed, the elevation is regularized taking into account the
        # maximum slope authorized
        if self.regularization:
            self.__elevation_from_start = self._slope_regularization()
        # Road slope
        self.__road_slope = np.array(
            [
                float(
                    (
                        self._Road__elevation_from_start[i + 1]
                        - self._Road__elevation_from_start[i]
                    )
                    / (
                        self._Road__distance_from_start[i + 1]
                        - self._Road__distance_from_start[i]
                    )
                )
                for i in range(len(self._Road__distance_from_start) - 1)
            ]
        )
        # Dirty trick in order to give the same length to both arrays
        self.__road_slope = np.append(self.__road_slope, self.__road_slope[-1])
        # Speed limit is converted into m/s
        self.__speed_limit = 1 / 3.6 * data[:, 2][mask]

    def export_regularized_profile(self):
        """Road profile with regularized slope is exported to another csv file
        with a slightly modified name.
        """
        if self.__is_a_data_file_loaded:
            # If the profile regularization is available
            if self.regularization:
                new_file_name = self.__data_file_name[:-4] + "-regularized.csv"
                with open(new_file_name, "w", newline="", encoding="utf-8") as file:
                    write = csv.writer(file)
                    write.writerow(["distance_m", "altitude_m", "speed_limit_kmh"])
                    # The initial altitude is added back, so that the
                    # exported file has the same meaning as the loaded one
                    for distance, altitude, speed in zip(
                        self.__distance_from_start,
                        self.__elevation_from_start + self.__initial_altitude,
                        self.__speed_limit,
                    ):
                        write.writerow(
                            [int(distance), f"{altitude:.2f}", f"{(speed * 3.6):.1f}"]
                        )

    def limit_speed_to_a_maximum_value(self, max_speed):
        """A useful method when you want to reduce the maximum permitted speed
        to a limit value 'max_speed' (in km/h) along the entire route.
        """
        self.__speed_limit = np.clip(
            self.__speed_limit, a_min=None, a_max=max_speed / 3.6
        )

    def road_length(self):
        """Total length of the road, in m."""
        return float(self.__distance_from_start[-1])

    def maximum_elevation(self):
        """Maximum elevation along the road, in m."""
        return float(self.__elevation_from_start.max())

    def maximum_authorized_speed(self):
        """Maximum value of the authorized speed along the whole path."""
        return float(self.__speed_limit.max())

    def maximum_average_speed(self):
        """Maximum average speed if the vehicle is at any time at the
        maximum speed authorized."""
        return float(self.road_length() / self.minimum_duration())

    def minimum_duration(self, format="s"):
        """Minimum duration of the trip if the vehicle is at any time at the
        maximum speed authorized."""
        # Travel time along each piece of the road, in seconds
        if len(self.__distance_from_start) != 1:
            durations = (
                self.__distance_from_start[1:] - self.__distance_from_start[:-1]
            ) / self.__speed_limit[:-1]
        else:
            durations = self.__distance_from_start / self.__speed_limit
        # Total duration in seconds
        duration_s = float(durations.sum())
        if format == "s":
            return duration_s
        # Results possibly printed in "hh:mm:ss" form
        return str(dt.timedelta(seconds=duration_s))

    def _slope_regularization(self):
        """When needed, the elevation of the road can be modified in order to
        maintain, at any distance from the start, the slope at a lower value
        than a maximum one previously defined.
        """
        distance = self.__distance_from_start
        # Previous and new value of the road elevation
        z_old = self.__elevation_from_start
        z_new = [0.0]
        # Calculation of the initial slope
        slope_old = [
            (z_old[i + 1] - z_old[i]) / (distance[i + 1] - distance[i])
            for i in range(len(distance) - 1)
        ]
        # And reconstruction of the elevation considering the limited slope
        # value
        for i in range(len(distance) - 1):
            # We modify the road elevation only if the slope absolute value is
            # larger than the maximum one authorized
            if abs(slope_old[i]) > self.max_slope:
                if slope_old[i] > 0:
                    z_new.append(
                        float(
                            z_new[-1] + (distance[i + 1] - distance[i]) * self.max_slope
                        )
                    )
                else:
                    z_new.append(
                        float(
                            z_new[-1] - (distance[i + 1] - distance[i]) * self.max_slope
                        )
                    )
            else:
                z_new.append(
                    float(z_new[-1] + (distance[i + 1] - distance[i]) * slope_old[i])
                )
        # New elevation values, relative to the starting point like the
        # initial ones
        return np.array(z_new)

    def __where(self, distance):
        """Index 'i' of the road segment [distance_from_start[i],
        distance_from_start[i + 1][ containing 'distance'. The end of the road
        is considered as belonging to the last segment. A binary search is
        used (np.searchsorted), much faster than a scan of the whole array.
        """
        if (distance < 0) or (distance > self.road_length()):
            raise ValueError("Argument must be within the road range distance!")
        i = int(np.searchsorted(self.__distance_from_start, distance, side="right")) - 1
        return min(i, len(self.__distance_from_start) - 2)

    def speed_limit(self, distance):
        """Speed limit (in m/s) the vehicle has to respect at the distance
        'distance' (in m) from the starting point.
        """
        return float(self.__speed_limit[self.__where(distance)])

    def slope(self, distance):
        """Slope of the road (dimensionless) at the distance 'distance' (in m)
        from the starting point.
        """
        # Index of x within the distance list
        i = self.__where(distance)
        slope = (
            self.__elevation_from_start[i + 1] - self.__elevation_from_start[i]
        ) / (self.__distance_from_start[i + 1] - self.__distance_from_start[i])
        return float(slope)

    def climbing_angle(self, distance):
        """Climbing angle of the road (in °/degrees) at the distance 'distance'
        from the starting point.
        """
        return float(np.rad2deg(np.arcsin(self.slope(distance))))

    def elevation(self, distance):
        """Elevation (in m) of the road at the distance 'distance' (in m) from
        the starting point. Each piece of the road is supposed to be a linear
        function."""
        # Index of x within the distance list
        i = self.__where(distance)
        # Percentage of the distance within the concerned piece of road
        elevation = distance - self.__distance_from_start[i]
        return float(elevation * self.slope(distance) + self.__elevation_from_start[i])

    def speed_statistics(self):
        """Reply with a dictionary whom keys are the different levels of
        authorized speed limit along the road and the corresponding values
        are the sum of all road section length corresponding to these speed
        limits.
        """
        speed_statistics = {}
        for i in range(len(self.__speed_limit) - 1):
            # Calculation of each road segment length
            start = self.__distance_from_start[i]
            end = self.__distance_from_start[i + 1]
            segment_lenth = float(end - start)
            # Corresponding speed limit
            speed = self.__speed_limit[i]
            # Add segment length to the corresponding speed limit
            if speed in speed_statistics:
                speed_statistics[float(speed)] += segment_lenth
            else:
                speed_statistics[float(speed)] = segment_lenth
        return speed_statistics


class LinearContinuousGearbox:
    """
    Mechanical transmission of a vehicle that transforms the engine output
    rotational speed into the driven wheels ones. The gear ratio, i.e. the
    ratio of the output to input rotational speed, linearly with the inlet
    speed, from its value at zero speed (minimum_gear_ratio) to its value at
    maximum inlet speed (maximum_gear_ratio). Resulting mechanical torque is
    then computed considering this gear ratio and the energy efficiency of the
    gearbox.
    ...
    Attributes
    ----------
    energy_efficiency : float
        Ratio of the outlet to inlet mechanical (power|torque), dimensionless.
    inlet_frequency : float
        Rotational frequency at the gearbox inlet, in rpm
    inlet_speed : float
        Rotational speed at the gearbox inlet, in rad/s
    inlet_torque : float
        Mechanical torque at the gearbox inlet, in N.m
    maximum_gear_ratio : float
        Maximum value of the output to input speed ratios, corresponding to a
        maximum inlet speed/frequency, dimensionless.
    maximum_inlet_frequency : float
        Maximum value of the inlet rotational frequency, in rpm.
    maximum_inlet_speed : float
        Maximum value of the inlet rotational speed, in rad/s.
    minimum_gear_ratio : float
        Minimum value of the output to input speed ratios, corresponding to a
        minimum inlet speed/frequency, dimensionless.
    outlet_frequency : float
        Rotational frequency at the gearbox outlet, in rpm
    outlet_speed : float
        Rotational frequency at the gearbox outlet, in rad/s
    outlet_torque : float
        Mechanical torque at the gearbox outlet, in N.m

    Methods
    -------
    gear_ratio() : float
        Actual value of the gear ratio.
    inlet_power() : float
        Mechanical power entering the gearbox, in W.
    maximum_outlet_frequency() : float
        Maximum value of the outlet rotational frequency, in rpm.
    maximum_outlet_speed() : float
        Maximum value of the outlet rotational speed, in rad/s.
    outlet_power() : float
        Mechanical power provided at the gearbox outlet, in W.
    """

    def __init__(self):
        # Energy efficiency, at 100% by default
        self.energy_efficiency = 1.0
        # Maximum and minimum values of the gear ratio
        self.minimum_gear_ratio = 1.0
        self.maximum_gear_ratio = 2.0
        # Default max inlet rotational frequency is 3000 rpm, value needed to
        # compute the gear ratio
        self._maximum_inlet_frequency = 3e3
        self._maximum_inlet_speed = rot_speed_from_rpm(3e3)
        # Default inlet rotational frequency is the maximum one, which leads to
        # the maximum transmission ratio
        self._inlet_frequency = 3e3
        self._inlet_speed = rot_speed_from_rpm(3e3)
        # Resulting outlet rotational speed and frequency
        self._outlet_frequency = 3e3 * 2.0
        self._outlet_speed = rot_speed_from_rpm(3e3) * 2.0
        # Inlet torque
        self._inlet_torque = 100.0
        # Outlet torque obtained in considering default values of the energy
        # efficiency and of the gear ratio
        self._outlet_torque = 100.0 * 1.0 / 2.0

    # ==== Minimum and maximum rotation speed and frequency values ============

    def get_maximum_inlet_speed(self):
        """Maximum value of the inlet rotational speed, in rad/s."""
        return self._maximum_inlet_speed

    def set_maximum_inlet_speed(self, speed):
        """Set of the maximum inlet rotational speed, in rad/s."""
        self._maximum_inlet_speed = speed
        self._maximum_inlet_frequency = rot_rpm_from_speed(speed)

    maximum_inlet_speed = property(
        fget=get_maximum_inlet_speed, fset=set_maximum_inlet_speed
    )

    def get_maximum_inlet_frequency(self):
        """Maximum value of the inlet rotational frequency, in rpm."""
        return self._maximum_inlet_frequency

    def set_maximum_inlet_frequency(self, freq):
        """Set of the maximum inlet rotational frequency, in rpm."""
        self._maximum_inlet_frequency = freq
        self._maximum_inlet_speed = rot_speed_from_rpm(freq)

    maximum_inlet_frequency = property(
        fget=get_maximum_inlet_frequency, fset=set_maximum_inlet_frequency
    )

    def maximum_outlet_speed(self):
        """Maximum rotational speed at the gearbox outlet, in rad/s."""
        return self.maximum_inlet_speed * self.maximum_gear_ratio

    def maximum_outlet_frequency(self):
        """Maximum rotational frequency at the gearbox outlet, in rpm."""
        return self.maximum_inlet_frequency * self.maximum_gear_ratio

    def _inlet_speed_from_outlet_speed(self, outlet_speed):
        """Inlet rotational speed, in rad/s, corresponding to a given outlet
        rotational speed, in rad/s. Since the gear ratio r is a linear function
        of the inlet speed w_in, r = r_min + k * w_in with
        k = (r_max - r_min) / w_in_max, the outlet speed w_out = r * w_in is
        a quadratic function of w_in, whose positive root is:
            w_in = 2 * w_out / (r_min + sqrt(r_min**2 + 4 * k * w_out))
        This form remains valid (and accurate) when k = 0.
        """
        # Beyond the maximum outlet speed, the inlet one would exceed its
        # maximum value (and the equation may even have no solution)
        if outlet_speed > self.maximum_outlet_speed() * (1.0 + 1e-9):
            raise ValueError(
                "Outlet speed is higher than the maximum one the gearbox can provide!"
            )
        slope = (
            self.maximum_gear_ratio - self.minimum_gear_ratio
        ) / self.maximum_inlet_speed
        return (
            2.0
            * outlet_speed
            / (
                self.minimum_gear_ratio
                + np.sqrt(self.minimum_gear_ratio**2 + 4.0 * slope * outlet_speed)
            )
        )

    def __get_gear_ratio(self, **kwargs):
        """Value of the actual gear ratio, for any inlet or outlet rotational
        speed or frequency.
        """
        # Allowed input arguments
        allowed_args = {
            "inlet_speed": float,
            "outlet_speed": float,
            "inlet_frequency": float,
            "outlet_frequency": float,
        }

        # Check if input argument is authorized
        # TODO : theses tests may be deleted once the function is called only
        # from the object inside.
        for key, value in kwargs.items():
            if key not in allowed_args:
                raise ValueError(f"Invalid argument '{key}'.")
            if not isinstance(value, allowed_args[key]):
                raise TypeError(
                    f"Argument '{key}' must be of type\
                {allowed_args[key].__name__}, not {type(value).__name__}"
                )
        # And calculation of the resulting gear ratio
        if "inlet_speed" in kwargs:
            variable = kwargs["inlet_speed"] / self.maximum_inlet_speed
        if "inlet_frequency" in kwargs:
            variable = kwargs["inlet_frequency"] / self.maximum_inlet_frequency
        # The gear ratio is a linear function of the INLET speed only, so the
        # inlet speed must first be deduced from the outlet one
        if "outlet_speed" in kwargs:
            variable = (
                self._inlet_speed_from_outlet_speed(kwargs["outlet_speed"])
                / self.maximum_inlet_speed
            )
        if "outlet_frequency" in kwargs:
            variable = (
                self._inlet_speed_from_outlet_speed(
                    rot_speed_from_rpm(kwargs["outlet_frequency"])
                )
                / self.maximum_inlet_speed
            )
        gear_ratio = (
            self.maximum_gear_ratio - self.minimum_gear_ratio
        ) * variable + self.minimum_gear_ratio
        return gear_ratio

    def gear_ratio(self):
        """Actual value of the gear ratio."""
        return self.__get_gear_ratio(inlet_frequency=self.inlet_frequency)

    # ==== Actual values of inlet and outlet speed and frequency ==============

    def get_inlet_frequency(self):
        """Inlet rotational frequency, in rpm."""
        return self._inlet_frequency

    def set_inlet_frequency(self, freq):
        """Set of the inlet rotational frequency, in rpm."""
        # Inlet
        self._inlet_frequency = freq
        self._inlet_speed = rot_speed_from_rpm(freq)
        # Gear ratio
        gear_ratio = self.__get_gear_ratio(inlet_frequency=freq)
        # Outlet
        self._outlet_frequency = gear_ratio * freq
        self._outlet_speed = gear_ratio * rot_speed_from_rpm(freq)

    inlet_frequency = property(fget=get_inlet_frequency, fset=set_inlet_frequency)

    def get_inlet_speed(self):
        """Inlet rotational speed, in rad/s."""
        return self._inlet_speed

    def set_inlet_speed(self, speed):
        """Set of the inlet rotational speed, in rad/s."""
        # Inlet
        self._inlet_speed = speed
        self._inlet_frequency = rot_rpm_from_speed(speed)
        # Gear ratio
        gear_ratio = self.__get_gear_ratio(inlet_speed=speed)
        # Outlet
        self._outlet_frequency = gear_ratio * rot_rpm_from_speed(speed)
        self._outlet_speed = gear_ratio * speed

    inlet_speed = property(fget=get_inlet_speed, fset=set_inlet_speed)

    def get_outlet_frequency(self):
        """Outlet rotational frequency, in rpm."""
        return self._outlet_frequency

    def set_outlet_frequency(self, freq):
        """Set of the outlet rotational frequency, in rpm."""
        # Outlet
        self._outlet_frequency = freq
        self._outlet_speed = rot_speed_from_rpm(freq)
        # Gear ratio
        gear_ratio = self.__get_gear_ratio(outlet_frequency=freq)
        # Inlet
        self._inlet_frequency = freq / gear_ratio
        self._inlet_speed = rot_speed_from_rpm(freq) / gear_ratio

    outlet_frequency = property(fget=get_outlet_frequency, fset=set_outlet_frequency)

    def get_outlet_speed(self):
        """Outlet rotational speed, in rad/s."""
        return self._outlet_speed

    def set_outlet_speed(self, speed):
        """Set of the outlet rotational speed, in rad/s."""
        # Outlet
        self._outlet_speed = speed
        self._outlet_frequency = rot_rpm_from_speed(speed)
        # Gear ratio
        gear_ratio = self.__get_gear_ratio(outlet_speed=speed)
        # Inlet
        self._inlet_frequency = rot_rpm_from_speed(speed) / gear_ratio
        self._inlet_speed = speed / gear_ratio

    outlet_speed = property(fget=get_outlet_speed, fset=set_outlet_speed)

    # ==== Mechanical torques =================================================

    def get_inlet_torque(self):
        """Inlet rotational torque, in N.m."""
        return self._inlet_torque

    def set_inlet_torque(self, torque):
        """Set of the inlet rotational torque, N.m."""
        # Inlet
        self._inlet_torque = torque
        # Outlet. The losses always reduce the power transmitted: when the
        # torque is positive (traction), power flows from inlet to outlet;
        # when it is negative (regenerative braking), from outlet to inlet
        if torque >= 0.0:
            self._outlet_torque = self.energy_efficiency * torque / self.gear_ratio()
        else:
            self._outlet_torque = torque / (self.energy_efficiency * self.gear_ratio())

    inlet_torque = property(fget=get_inlet_torque, fset=set_inlet_torque)

    def get_outlet_torque(self):
        """Outlet rotational torque, in N.m."""
        return self._outlet_torque

    def set_outlet_torque(self, torque):
        """Set of the outlet rotational torque, N.m."""
        # Outlet
        self._outlet_torque = torque
        # Inlet torque, taking into account the gearbox efficiency and the
        # direction of the power flow (see set_inlet_torque)
        if torque >= 0.0:
            self._inlet_torque = torque * self.gear_ratio() / self.energy_efficiency
        else:
            self._inlet_torque = torque * self.gear_ratio() * self.energy_efficiency

    outlet_torque = property(fget=get_outlet_torque, fset=set_outlet_torque)

    # ==== Mechanical power at inlet and outlet ===============================

    def inlet_power(self):
        """Mechanical power entering into the gearbox, in W."""
        return self.inlet_torque * self.inlet_speed

    def outlet_power(self):
        """Mechanical power provided at the gearbox outlet, in W."""
        return self.outlet_torque * self.outlet_speed

    def dissipated_power(self):
        """Mechanical power dissipated into heat inside the gearbox, in W."""
        return self.inlet_power() - self.outlet_power()


class SynchronousElectricMotor:
    """
    Synchronous electric motor used for the vehicle propulsion, which is able
    to provide a constant torque from zero speed to base (speed|frequency). Once
    speed is larger than the based one, maximum torque is decreasing such as the
    resulting maximum power is a constant.

    When either the maximum power or the maximum torque are changed, the base
    (speed|frequency) is modified as a result.
    ...
    Attributes
    ----------
    base_frequency : float
        Rotational frequency, in rpm, beyond which the maximum power is a
        constant and the corresponding torque starts to decrease.
    base_speed : float
        Rotational speed, in rad/s, beyond which the maximum power is a
        constant and the corresponding torque starts to decrease.
    energy_efficiency : float
        Ratio of the mechanical power produced by the electrical one consumed,
        dimensionless.
    maximum_frequency : float
        Maximum value of the motor rotational frequency, in rpm.
    maximum_speed : float
        Maximum value of the motor rotational speed, in rad/s.
    maximum_power : float
        Maximum mechanical power the motor is able to provided once its
        (speed|frequency) is higher than the base one, in W. This power is lower
        thant the electrical one it can consume.
    maximum_torque : float
        Maximum value of the mechanical torque the motor is able to provide
        when its (speed|frequency) is lower than the base one, in N.m.

    Methods
    -------
    maximum_operating_torque : float
        Maximum torque the motor is able to provide at a given value of its
        (speed|frequency), in N.m.
    """

    def __init__(self):
        # Base and maximum (frequency|speed)
        self._base_frequency = 800.0
        self._base_speed = rot_speed_from_rpm(800.0)
        self._maximum_frequency = 2500.0
        self._maximum_speed = rot_speed_from_rpm(2500.0)
        self._maximum_torque = 2500.0
        self._maximum_power = 2500.0 * rot_speed_from_rpm(800.0)
        # Efficiency is a constant value by default
        self.energy_efficiency = 0.9

    def get_base_speed(self):
        """Rotational speed, in rad/s, beyond whch the maximum power is a
        constant and the corresponding torque starts to decrease.
        """
        return self._base_speed

    def set_base_speed(self, speed):
        """Set of the rotational speed, in rad/s, beyond which the maximum
        power is a constant and the corresponding torque starts to decrease.
        """
        # We check whether this speed is lower to the maximum one
        if speed >= self.maximum_speed:
            raise ValueError("Base speed must be lower than the maximum one!")
        self._base_speed = speed
        self._base_frequency = speed * 30.0 / np.pi

    base_speed = property(fget=get_base_speed, fset=set_base_speed)

    def get_base_frequency(self):
        """Rotational frequency, in rpm, beyond which the maximum power is a
        constant and the corresponding torque starts to decrease.
        """
        return self._base_frequency

    def set_base_frequency(self, freq):
        """Set of the rotational frequency, in rpm, beyond which the maximum
        power is a constant and the corresponding torque starts to decrease.
        """
        # We check whether this frequency is lower to the maximum one
        if freq >= self.maximum_frequency:
            raise ValueError("Base frequency must be lower than the maximum one!")
        self._base_frequency = freq
        self._base_speed = freq * np.pi / 30

    base_frequency = property(fget=get_base_frequency, fset=set_base_frequency)

    def get_maximum_speed(self):
        """Maximum value of the motor rotational speed, in rad/s."""
        return self._maximum_speed

    def set_maximum_speed(self, speed):
        """Set of the maximum value of the motor rotational speed, in rad/s."""
        # We check whether this speed is lower to the maximum one
        if speed <= self.base_speed:
            raise ValueError("Maximum speed must be higher than the base one!")
        self._maximum_speed = speed
        self._maximum_frequency = speed * 30.0 / np.pi

    maximum_speed = property(fget=get_maximum_speed, fset=set_maximum_speed)

    def get_maximum_frequency(self):
        """Maximum value of the motor rotational frequency, in rpm."""
        return self._maximum_frequency

    def set_maximum_frequency(self, freq):
        """Set of the maximum value of the motor rotational frequency, in rpm."""
        # We check whether this speed is lower to the maximum one
        if freq <= self.base_frequency:
            raise ValueError("Maximum frequency must be higher than the base one!")
        self._maximum_frequency = freq
        self._maximum_speed = freq * np.pi / 30.0

    maximum_frequency = property(fget=get_maximum_frequency, fset=set_maximum_frequency)

    def get_maximum_torque(self):
        """Maximum value of the mechanical torque the motor is able to provide
        when its (speed|frequency) is lower than the base one, in N.m.
        """
        return self._maximum_torque

    def set_maximum_torque(self, torque):
        """Set of the maximum value of the mechanical torque the motor is able
        to provide when its (speed|frequency) is lower than the base one,
        in N.m.
        """
        self._maximum_torque = torque
        # And the resulting base speed is modified
        self.base_speed = self.maximum_power / torque

    maximum_torque = property(fget=get_maximum_torque, fset=set_maximum_torque)

    def get_maximum_power(self):
        """Maximum mechanical power the motor is able to provided once its
        (speed|frequency) is higher than the base one, in W.
        """
        return self._maximum_power

    def set_maximum_power(self, power):
        """Set of the maximum mechanical power the motor is able to provided
        once its (speed|frequency) is higher than the base one, in W.
        """
        self._maximum_power = power
        # And the resulting base speed is modified
        self.base_speed = power / self.maximum_torque

    maximum_power = property(fget=get_maximum_power, fset=set_maximum_power)

    def maximum_operating_torque(self, **kwargs):
        """Maximum torque the motor is able to provide at a given value of its
        (speed|frequency), in N.m.
        """
        # Entered arguments are either the speed or the frequency
        allowed_args = {"frequency": float, "speed": float}

        # Check if input argument is authorized
        for key, value in kwargs.items():
            if key not in allowed_args:
                raise ValueError(f"Invalid argument '{key}'.")
            if not isinstance(value, allowed_args[key]):
                raise TypeError(
                    f"Argument '{key}' must be of type\
                {allowed_args[key].__name__}, not {type(value).__name__}"
                )
        if "frequency" in kwargs:
            if kwargs["frequency"] <= self.base_frequency:
                return self.maximum_torque
            return self.maximum_power / (kwargs["frequency"] * np.pi / 30)
        if "speed" in kwargs:
            if kwargs["speed"] <= self.base_speed:
                return self.maximum_torque
            return self.maximum_power / kwargs["speed"]


class Vehicle:
    """
    Vehicle whom power consumption is researched for a given road profile.
    ...
    Attributes
    ----------
    drag_coefficient : float
        Drag coefficient of the vehicle, usually noted Cx or Cd, dimensionless
        and considered as a constant within the speed range of the vehicle.
    frontal_area : float
        Expressed in m².
    mass : float
        Gravitational mass of the vehicle in kg.
    rotating_mass_correction_coefficient : float
        Correction coefficient, dimensionless, accounting for the equivalent
        mass increase due to the angular moment of the rotating components.
    wheels_radius : float
        Radius of the vehicle wheels, in m.

    Methods
    -------
    acceleration(motive_force, slope, speed, **kwargs) : float
        Acceleration, in m/s², that could be reached for a given value of a
        motive force, in N. Can also be used with zero motive power to calculate
        negative acceleration.
    climbing_resistance(slope) : float
        Resistance force due to the slope of the road.
    drag_resistance(speed, **density) : float
        Resistance force due to the air drag effect around the vehicle and
        applied in the motion direction, in N. Speed is a mandatory argument
        while air density can be entered as argument if necessary.
    motive_force(acceleration, slope, speed, **kwargs) : float
        Motive force, in N, required to reach a specific acceleration value,
        in m/s². This function can be used also to calculate braking force using
        a negative acceleration as an input.
    power_to_maximum_speed(power, slope, **kwargs) : float
        Calculation of the maximum vehicle speed, in m/s, reachable for a given
        propulsion power, in W. Optional arguments **kwargs can be air_density
        (**density) or rolling coefficient (**rolling_coefficient).
    force_to_maximum_speed(power, slope, **kwargs) : float
        Calculation of the maximum vehicle speed, in m/s, reachable for a given
        propulsion force, in N. Optional arguments **kwargs can be air_density
        (**density) or rolling coefficient (**rolling_coefficient).
    resistance_force(self, slope, speed, **kwargs) : float
        Sum of all the resistance forces, in N, that opposite the vehicle
        motion, as function of the road slope and the vehicle speed at least.
    resistance_power(self, slope, speed, **kwargs) : float
        Total mechanical power, in W, required to maintain the vehicle motion
        against resistance forces.
    rolling_resistance(slope, **rolling_coefficient) : float
        Resistance force due to the rolling of tires on the road and computed
        in taking into account the road slope, in N. A 'rolling_coefficient'
        can be used as argument if necessary.
    normal_weight() : float
        Projection of the vehicle weight along the direction normal to the road,
        in N.
    inertial_mass() : float
        Equivalent mass which is actually multiplied by the acceleration, that
        takes into account the equivalent mass increase due to the angular
        moment of the rotating components.
    weight() : float
        Vehicle weight, in N.
    """

    def __init__(self):
        # Vehicle mass, frontal area and drag coefficient
        self.mass = 1e3
        self.frontal_area = 2.0
        self.drag_coefficient = 0.3
        self.wheels_radius = 0.3
        # Rotating mass correction coefficienta, that takes into account the
        # inertia created by rotating parts in the vehicle
        self.rotating_mass_correction_coefficient = 0.0

    def acceleration(self, motive_force, slope, speed, **kwargs):
        """Acceleration, in m/s², that could be reached for a given value of a
        motive force, in N. If the result is negative, the motive force entered
        as argument is not sufficient to compensate the different resistive
        forces. Can also be used with zero motive power to calculate negative
        acceleration.
        """
        return (
            motive_force - self.resistance_force(slope, speed, **kwargs)
        ) / self.inertial_mass()

    def climbing_resistance(self, slope):
        """Resistance force due to the road slope and applied in the motion
        direction, in newtons (N) : m * g * sin(a) with a the road angle.
        """
        return self.weight() * slope

    def drag_resistance(self, speed, **kwargs):
        """Resistance force due to the air drag effect around the vehicle and
        applied in the motion direction, in newtons (N). Argument 'speed' must
        be entered in m/s. Optional argument 'air_density' can be used to enter
        the specific air density value. A value of 1.2 g/L is used instead.
        """
        if "air_density" in kwargs:
            rho = kwargs.get("air_density")
        else:
            rho = 1.2
        return self.frontal_area * self.drag_coefficient * 0.5 * rho * pow(speed, 2)

    def inertial_mass(self):
        """Equivalent mass which is actually multiplied by the acceleration,
        that takes into account the equivalent mass increase due to the angular
        moment of the rotating components.
        """
        return self.mass * (1 + self.rotating_mass_correction_coefficient)

    def motive_force(self, acceleration, slope, speed, **kwargs):
        """Motive force, in N, required to reach a specific acceleration
        value, in m/s². This function can be used also to calculate braking
        force using a negative acceleration as an input.
        """
        return self.inertial_mass() * acceleration + self.resistance_force(
            slope, speed, **kwargs
        )

    def normal_weight(self, slope):
        """Projection of the vehicle weight along the direction normal to the
        road, in newtons (N) : m * g * cos(a) with a the road angle.
        """
        return float(self.weight() * np.sqrt(1 - pow(slope, 2)))

    def force_to_maximum_speed(self, force, slope, **kwargs):
        """Maximum vehicle speed, in m/s, reachable for a given propulsion
        force, in N. Optional arguments **kwargs can be air air_density
        (**density) or rolling coefficient (**rolling_coefficient).
        """

        # Function to solve to obtain the searched speed, that gives a force
        # as a result
        def f_to_solve(speed):
            return self.resistance_force(slope, speed, **kwargs) - force

        # And solving process starting from an initial zero speed
        sol = opt.root_scalar(f_to_solve, method="secant", x0=0.0)
        return float(sol.root)

    def power_to_maximum_speed(self, power, slope, **kwargs):
        """Maximum vehicle speed, in m/s, reachable for a given propulsion
        power, in W. Optional arguments **kwargs can be air air_density
        (**density) or rolling coefficient (**rolling_coefficient).
        """

        # Function to solve to obtain the searched speed, that gives a power
        # as a result
        def f_to_solve(speed):
            return self.resistance_force(slope, speed, **kwargs) * speed - power

        # And solving process starting from an initial zero speed
        sol = opt.newton(f_to_solve, 0.0)
        return float(sol)

    def resistance_force(self, slope, speed, **kwargs):
        """Sum of all the resistance forces that opposite the vehicle motion,
        in N, as function of the road slope and the vehicle speed at least.
        """
        return (
            self.climbing_resistance(slope)
            + self.rolling_resistance(slope, **kwargs)
            + self.drag_resistance(speed, **kwargs)
        )

    def resistance_power(self, slope, speed, **kwargs):
        """Total mechanical power, in W, required to maintain the vehicle motion
        against resistance forces.
        """
        return self.resistance_force(slope, speed, **kwargs) * speed

    def rolling_resistance(self, slope, **kwargs):
        """Resistance force due to the rolling of tires on the road and
        computed in taking into account the road slope, in newtons (N). Optional
        argument 'rolling_coefficient' can be used to enter a specific value of
        the resistance force to normal force ratio. A value of 0.01 is used
        instead.
        """
        if "rolling_coefficient" in kwargs:
            k = kwargs.get("rolling_coefficient")
        else:
            k = 0.01
        return self.normal_weight(slope) * k

    def weight(self):
        """Vehicle weight, in newtons (N)"""
        return self.mass * g


class VehicleDynamicsModel:
    """
    Simulation of the longitudinal dynamics of a vehicle following a road
    profile. At each time step, the vehicle speed is compared to the local
    speed limit: the vehicle accelerates, decelerates or maintains its speed,
    within the limits of the electric motor maximum torque and power.
    ...
    Attributes
    ----------
    gearbox : LinearContinuousGearbox
    electric_motor : SynchronousElectricMotor
    road : Road
    vehicle : Vehicle
    average_acceleration : float
        Target value of the vehicle positive acceleration, used when a speed
        increase is required, in m/s².
    average_deceleration : float
        Target value (positive) of the vehicle deceleration, used when a speed
        decrease is required, in m/s².
    time_step : float
        Time step used in the calculation, in s. The last step is shortened so
        that the vehicle stops exactly at the end of the road.
    regenerative_braking_fraction : float
        Fraction of the braking torque at the wheels sent to the electric
        motor, before its torque limit is applied, dimensionless. The remaining
        part is dissipated by the friction brakes. Default is 1.
    results : dict
        Time series produced by run_calculation(), as lists:
        'time' (s), 'distance' (m), 'speed' (m/s), 'acceleration' (m/s²),
        'power' (electrical power at the motor terminals, W, positive when
        consumed and negative when recovered), 'motor_torque' (N.m),
        'motor_speed' (rad/s), 'friction_braking_power' (W, positive).
        Each value (but the first one, at t = 0) is applied during the time
        step ending at the corresponding time.

    Methods
    -------
    run_calculation(initial_speed=0., **kwargs)
        Run the simulation. Initial speed in km/h. Optional keyword arguments
        ('rolling_coefficient', 'air_density') are passed to the Vehicle
        methods.
    consumed_energy(regenerative_braking=True) : float
        Electrical energy consumed along the trip, in Wh, with or without the
        energy recovered during braking.
    speeds() : dict
        Minimum, maximum and average speeds along the trip, in m/s.
    trip_duration() : float
        Total duration of the trip, in s.
    """

    def __init__(
        self,
        gearbox: LinearContinuousGearbox,
        electric_motor: SynchronousElectricMotor,
        road: Road,
        vehicle: Vehicle,
        average_acceleration=1.0,
        average_deceleration=1.0,
        time_step=1.0,
        regenerative_braking_fraction=1.0,
    ):
        # Components of the whole model
        self.gearbox = gearbox
        self.motor = electric_motor
        self.road = road
        self.vehicle = vehicle
        # Initial and final speed, in m/s
        self.average_acceleration = average_acceleration
        self.average_deceleration = average_deceleration
        # Time step used along the calculation, in seconds
        self.time_step = time_step
        # Fraction of the braking torque at the wheels that is sent to the
        # electric motor (the remaining part being dissipated by the friction
        # brakes), before the motor torque limit is applied. 1 means that the
        # motor brakes first, as much as it can (upper bound of recovery)
        self.regenerative_braking_fraction = regenerative_braking_fraction
        # Results of the calculation are stored within a dictionary
        self.results = {}

    def consumed_energy(self, regenerative_braking=True):
        """Calculation of the total energy consumed by the vehicle, in Wh, to
        achieve the trip. Whether the option 'regenerative_braking' is set to
        True or False, the negative power values are either included (energy
        recovered) or discarded.
        """
        if not self.results:
            return 0.0
        if regenerative_braking:
            # With regenerative braking, both positive and negative values of
            # the power are taken into account
            power = self.results["power"]
        else:
            # When regenerative braking is not used, only the positive values
            # of power, i.e. the ones actually consumed by vehicle, are taken
            # into account
            power = np.array(self.results["power"]).clip(min=0.0)
        # Claude : each stored power value is applied during the time step that
        # ends at the corresponding time, so a rectangle rule is used (the
        # first value, at t = 0, is only an initialization and is therefore
        # skipped)
        power = np.asarray(power)
        durations = np.diff(self.results["time"])
        return float(np.sum(power[1:] * durations) / 3.6e3)

    def run_calculation(self, initial_speed=0.0, **kwargs):
        """Run the whole vehicle dynamics model. Initial speed must be entered
        here in km/h.
        """
        time, distance, speed, acceleration = 0.0, 0.0, initial_speed / 3.6, 0.0
        self.results["time"] = [time]
        self.results["distance"] = [distance]
        self.results["speed"] = [speed]
        self.results["acceleration"] = [acceleration]
        self.results["power"] = [0.0]
        self.results["motor_torque"] = [0.0]
        self.results["motor_speed"] = [0.0]
        self.results["friction_braking_power"] = [0.0]
        # Calculation is run until the end of the road is achieved
        while distance < self.road.road_length():
            # Local road slope and speed limit
            slope = self.road.slope(distance)
            speed_limit = self.road.speed_limit(distance)
            # Calculation of the gearbox gear ratio from the outlet speed value
            self.gearbox.outlet_speed = speed / self.vehicle.wheels_radius
            # If the current speed is lower than the limit one, we must
            # accelerate
            if speed < speed_limit:
                # Such an acceleration requires a torque, firstly calculated at
                # the wheels, i.e. at the gearbox outlet
                self.gearbox.outlet_torque = (
                    self.vehicle.motive_force(
                        self.average_acceleration, slope, speed, **kwargs
                    )
                    * self.vehicle.wheels_radius
                )
                # The resulting gearbox inlet torque must be compared with the
                # maximum one the electric motor can provide at the given speed
                motor_max_torque = self.motor.maximum_operating_torque(
                    speed=self.gearbox.inlet_speed
                )
                # If the actual torque is lower than the maximum one
                if self.gearbox.inlet_torque < motor_max_torque:
                    # Acceleration is then the nominal one
                    acceleration = self.average_acceleration
                # If not
                else:
                    # The available torque is then the maximum one
                    self.gearbox.inlet_torque = motor_max_torque
                    # So the resulting motive force at the wheels
                    force_wheels = (
                        self.gearbox.outlet_torque / self.vehicle.wheels_radius
                    )
                    # And then the actual acceleration reachable
                    acceleration = self.vehicle.acceleration(
                        force_wheels, slope, speed, **kwargs
                    )
                # We can now calculate the vehicle speed at the next time step
                next_speed = speed + acceleration * self.time_step
                # If the speed limit is actually reached
                if next_speed > speed_limit:
                    # The acceleration needed to reach exactly the speed limit
                    # at the end of the time step is calculated first, from
                    # the current speed
                    acceleration = (speed_limit - speed) / self.time_step
                    # So the corresponding actual gearbox outlet torque,
                    # evaluated at the current speed like the other forces
                    self.gearbox.outlet_torque = (
                        self.vehicle.motive_force(acceleration, slope, speed, **kwargs)
                        * self.vehicle.wheels_radius
                    )
                    # And only then is the speed updated
                    speed = speed_limit
                # If not
                else:
                    speed = next_speed
            # If the current speed is higher than the limit one, we must
            # decrease the vehicle speed
            elif speed > speed_limit:
                # Deceleration required
                acceleration = -self.average_deceleration
                # A braking torque must be applied to follow this speed
                # decrease
                self.gearbox.outlet_torque = (
                    self.vehicle.motive_force(acceleration, slope, speed, **kwargs)
                    * self.vehicle.wheels_radius
                )
                # Speed at the next time step
                next_speed = speed + acceleration * self.time_step
                # As for the speed increase, we can replace the actual speed by
                # the limit one if the latter is reached
                if next_speed < speed_limit:
                    # Same order as above: deceleration first, then torque,
                    # then speed update
                    acceleration = (speed_limit - speed) / self.time_step
                    self.gearbox.outlet_torque = (
                        self.vehicle.motive_force(acceleration, slope, speed, **kwargs)
                        * self.vehicle.wheels_radius
                    )
                    speed = speed_limit
                # If not
                else:
                    speed = next_speed
            # Otherwise, the current speed is equal to the limit one
            else:
                acceleration = 0.0
                # Calculation process is pretty much the same but without any
                # concern of acceleration
                self.gearbox.outlet_torque = (
                    self.vehicle.resistance_force(slope, speed, **kwargs)
                    * self.vehicle.wheels_radius
                )
                # The resulting gearbox inlet torque must be compared with the
                # maximum one the electric motor can provide
                motor_max_torque = self.motor.maximum_operating_torque(
                    speed=self.gearbox.inlet_speed
                )
                if self.gearbox.inlet_torque > motor_max_torque:
                    # Required torque cannot be provided by the electric motor
                    self.gearbox.inlet_torque = motor_max_torque
                    # So the propulsion force at the wheels
                    force_wheels = (
                        self.gearbox.outlet_torque / self.vehicle.wheels_radius
                    )
                    # And the resulting acceleration
                    acceleration = self.vehicle.acceleration(
                        force_wheels, slope, speed, **kwargs
                    )
                    # Actual speed
                    speed += acceleration * self.time_step
            # Braking: whatever the branch above, a negative torque at the
            # wheels is shared between the electric motor (regenerative
            # braking, limited by its maximum torque at the current speed) and
            # the friction brakes, which provide the remaining part. The
            # vehicle dynamics is therefore unchanged.
            friction_braking_power = 0.0
            if self.gearbox.outlet_torque < 0.0:
                total_braking_torque = self.gearbox.outlet_torque
                self.gearbox.outlet_torque = (
                    self.regenerative_braking_fraction * total_braking_torque
                )
                motor_max_torque = self.motor.maximum_operating_torque(
                    speed=self.gearbox.inlet_speed
                )
                if self.gearbox.inlet_torque < -motor_max_torque:
                    self.gearbox.inlet_torque = -motor_max_torque
                # Power dissipated by the friction brakes, positive, in W
                friction_braking_power = float(
                    (self.gearbox.outlet_torque - total_braking_torque)
                    * self.gearbox.outlet_speed
                )
            # Electrical power at the motor terminals: consumed in traction
            # (mechanical power divided by the efficiency), recovered in
            # regenerative braking (mechanical power multiplied by it)
            mechanical_power = float(self.gearbox.inlet_power())
            if mechanical_power >= 0.0:
                electrical_power = (
                    min(self.motor.maximum_power, mechanical_power)
                    / self.motor.energy_efficiency
                )
            else:
                electrical_power = mechanical_power * self.motor.energy_efficiency
            # Speed being now known, distance and time can be calculated. The
            # last time step is shortened so that the vehicle stops exactly at
            # the end of the road (the power is then applied during this
            # shorter duration only)
            step = self.time_step
            if distance + speed * step >= self.road.road_length():
                step = (self.road.road_length() - distance) / speed
                distance = self.road.road_length()
            else:
                distance += speed * step
            time += step
            # And results can now be stored in the dedicated dictionary
            self.results["distance"].append(distance)
            self.results["speed"].append(speed)
            self.results["acceleration"].append(acceleration)
            self.results["time"].append(time)
            # Electrical power at the motor terminals, in W
            self.results["power"].append(electrical_power)
            self.results["motor_torque"].append(self.gearbox.inlet_torque)
            self.results["motor_speed"].append(self.gearbox.inlet_speed)
            self.results["friction_braking_power"].append(friction_braking_power)

    def speeds(self):
        """Minimum, average and maximum speed reached by the vehicle along the
        trip, in m/s. Result is a dictionary.
        """
        if not self.results:
            return {}
        return {
            "min": min(self.results["speed"]),
            "max": max(self.results["speed"]),
            "avg": self.road.road_length() / self.trip_duration(),
        }

    def trip_duration(self):
        """Total duration of the trip, in seconds."""
        if not self.results:
            return 0.0
        return self.results["time"][-1]


if __name__ == "__main__":
    pass
