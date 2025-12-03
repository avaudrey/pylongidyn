# pylongidyn
> A simple, modular vehicle longitudinal dynamics model
## Introduction
**pylongidyn** is a simple and lightweight Python package dedicated to the study of the energy consumption of electric road vehicles.  It is based on a _**longi**tudinal vehicle **dyn**amics_ model, see references [1 page 39, 2 page 15]. Such a model basically consider the concern vehicle as a constant mass object following a specific path. This vehicle is constantly submitted to _climbing_, _rolling_ and _drag resistance_ forces, as well as inertial force when variations of speed are required. When the sum of all these forces is multiplied by the speed at which the vehicle is travelling, it gives us the instantaneous mechanical power needed for maintain or change this motion. 

This package provides simple tools to represent _ambient conditions_, _road profiles_, _continuous gearboxes_ and _electric motors_.

## Installation
The `pylongidyn` module is a single Python file. Simply copy it into your project:
```bash
your_project/
    pylongidyn.py
```
Then import it:
```python
from pylongidyn import *
```

## Main classes
A complete model dynamics calculation requires the creation of different related objects.

### AmbientAir
This class, which represents the vehicle surroundings, is currently used only for determining the **density of air**, required for calculating the [aerodynamic drag](https://en.wikipedia.org/wiki/Drag_(physics)#Aerodynamics) force. 
#### Example
```python
from pylongidyn import AmbientAir

air = AmbientAir(pressure=101325, temperature=25)
print(r"Air density: rho = %.3f kg/m³" % air.density())
```
Pressure must be entered in Pa and temperature in °C. [Ideal gas law](https://en.wikipedia.org/wiki/Ideal_gas_law) is used to model the air behavior.

### Road
This class is dedicated to handle road profiles, slopes, elevation, speed limits and so on.
#### Loading a road profile
Road properties must be loaded for a `csv` file that contains three columns, such as:

|distance (m) | Altitude (m) | Speed limit (km/h)|
|-------------|--------------|-------------------|
|0.|0.|50.|
|200.|4.|50.|
|400.|8.|50.|
|...|...|...|

The **distance** from starting point, expressed in meters, can be discretized regularly or not (so with a possibly non constant step). Here, with a road slope of 2%, the road climb to 4 meters after 200 meters of distance, 8 meters after 400 meters and so on.
#### Example
```python
from pylongidyn import Road

road = Road()
road.load_road_data("myroad.csv")

# Main properties of the road
print("Road length: %.3f km" % (1e-3 * road.road_length()))
print("Max elevation: %.0f m" % road.maximum_elevation())
# Although speed limit is imported in km/h, it is expressed in m/s once manipulated within the Road class
print("Speed at 150m: %.1f km/h" % (3.6 * road.speed_limit(150)))
print("Slope at 150m: %.2f %%" % (1e2 * road.slope(150)))
```
Using a linear interpolation along each piece of road, parameters such as _speed limit_ or _road slope_ can be evaluated at any distance `x` from the starting point, using the methods `speed_limit(x)` and `slope(x)` respectively.
### Vehicle
This class represents the concerned vehicle as an object whom main attributes correspond to practical and typical properties of such a system.
#### Example
```python
from pylongidyn import Vehicle

v = Vehicle()
# Mass is expressed in kg and area in m^2
v.mass = 1400.
v.frontal_area = 2.2
v.drag_coefficient = 0.29
v.wheels_radius = 0.3

# For a 2% road positive slope and a 80 km/h speed:
slope = 0.02     # 5%
speed = 80 / 3.6 # m/s

# Examples of resistive forces
print("Drag force: %.0f N" % v.drag_resistance(speed))
print("Climbing resistance: %.0f N" % v.climbing_resistance(slope))
# The default value of rolling coefficient is 0.01
print("Rolling resistance: %.0f N" % v.rolling_resistance(slope, rolling_coefficient=5e-3))
print("Total resistance: %.0f N" % v.resistance_force(slope, speed))

# Required motive force for acceleration
print("Force for +1 m/s²:", v.motive_force(1, slope, speed))
```
### Electric motor
The model of electric motor used to propel the vehicle is for the moment useful only for the calculation of the maximum mechanical torque it is possible to inject at any moment at the gearbox inlet. It is then defined through only three properties/attributes.
#### Example
```python
from pylongidyn import SynchronousElectricMotor, rot_speed_from_rpm, rot_rpm_from_speed

motor = SynchronousElectricMotor()
# The maximum torque can be provided from zero to base speed and decreases after that
motor.maximum_torque = 250.     # N.m
# The maximum power is reached one the speed reaches the base one
motor.maximum_power = 100e3   	# W
motor.maximum_frequency = 4500. # rpm

print("Maximum speed: %.0f rad/s" % rot_speed_from_rpm(motor.maximum_frequency))
print("Max torque at 1000 rpm: %.0f N.m" %  motor.maximum_operating_torque(frequency=1000.))
```
The rate at which the motor is rotating can be expressed either as a `frequency` in rpm, or as a (rotational) `speed` in rad/s. Conversion from one to the other can be achieve through dedicated functions called `rot_speed_from_rpm` and `rot_rpm_from_speed`, respectively.
### Gearbox
`GearBox` is an object which transforms the mechanical power produced by the electrical engine into the one provided at the vehicle wheels. For the moment, the sole model of gearbox available is a linear one, close to a [CVT](https://en.wikipedia.org/wiki/Continuously_variable_transmission), that linearly transforms the electric motor speed into the vehicle wheels speed. For so, the gearbox [gear ratio](https://en.wikipedia.org/wiki/Gear_train) varies linearly:
1. from the zero motor speed corresponding to the beginning of vehicle motion and to a maximum torque transmitted to the wheels;
2. to a maximum motor speed corresponding to a maximum travel speed for the vehicle.

Once these two values are calculated, any change in either the input of output speeds or torque triggers the calculation of all the resulting other parameters.

For the moment, the _energy efficiency_ is a constant value. 
#### Example
```python
from pylongidyn import LinearContinuousGearbox

gb = LinearContinuousGearbox()

# Maximum torque that must be provided to the wheels
maximum_wheel_torque = 2000. # N.m
# Maximum speed reachable by the vehicle
maximum_vehicle_speed = 130 / 3.6 # m/s
# The gear ratio at minimum speed is imposed by a torque concern
gb.minimum_gear_ratio = motor.maximum_torque / maximum_wheel_torque
# The gear ratio at maximum speed is imposed by the maximum speed of the vehicle
gb.maximum_gear_ratio = (maximum_vehicle_speed / v.wheel_radius) / motor.maximum_speed

gb.inlet_frequency = 1000.      # rpm
print("Current gear ratio: %.3f" % gb.gear_ratio())
print("Outlet speed: %.3f rad/s" % gb.outlet_speed)

gb.inlet_torque = 120          # N.m
print("Outlet torque: %.2f N.m" % gb.outlet_torque)
print("Mechanical power: %.3f W" % gb.outlet_power())
```

### VehicleDynamicsModel
The previously presented objects `Road`, `Vehicle`, `SynchronousElectricMotor` and `LinearContinuousGearbox` are the only ones need to launch a dynamic model calculation. The last parameters required to do so are the values of _nominal acceleration_ and _nominal deceleration_ the vehicle will try to reach anytime its speed has to change. 
#### How does the algorithm works
The whole algorithm basically works as follow : at any time, the vehicle speed is compared to the speed limit and:

1. If the vehicle travels at a speed **lower than the maximum one authorized**, it has to accelerate. The maximum torque avalaible for so is then checked and:

   1.1 If the torque required to reach the _nominal acceleration_ is **lower than the maximum one** providable to the wheels, the vehicle accelerate at this previously defined rate.

   1.2 If not, it accelerates (or even decelerates sometimes, because of the local road slope) according the maximum torque the electric motor is able to provide at this specific moment.

2. If it travels at a **higher speed than the one authorized**, it as to decelerate. For the moment, no **maximum braking torque** is considered and the previously defined _nominal deceleration_ is always reached. 
3. If the vehicle travels at the maximum speed authorized, the same calculation as at step 1.1 is achieved, but **to maintain the speed against the resistive forces** rather than to accelerate.

A very similar model has been published for example in reference [3].
#### Example
```python
from pylongidyn import VehicleDynamicsModel

# Dynamics model
model = VehicleDynamicsModel(
    gearbox=gb,
    electric_motor=motor,
    road=road,
    vehicle=v,
    average_acceleration=1.0, # m/s^2
    average_deceleration=1.0, # m/s^2
    time_step=1.0 # seconds
)

# Run simulation (initial speed = 0 km/h)
model.run_calculation(initial_speed=0.)

print("Trip duration (s):", model.trip_duration())
print("Energy consumed (Wh):", model.consumed_energy())
print("Speed statistics (m/s):", model.speeds())
```
## License
This package is released under the GNU General Public License v3 as described in the source file.
See http://www.gnu.org/licenses/
## References
[1] Jazar, R. N., [_Vehicle Dynamics: Theory and Applications_](https://link.springer.com/book/10.1007/978-3-031-74458-7), Springer, 2008.

[2] Minaker, B. P., [_Fundamentals of Vehicle Dynamics and Modelling_](https://www.wiley.com/en-us/Fundamentals+of+Vehicle+Dynamics+and+Modelling%3A+A+Textbook+for+Engineers+With+Illustrations+and+Examples-p-9781118980095), 2019, Wiley.

[3] Xun, Q., Murgovski, N., & Liu, Y. (2022). _Joint component sizing and energy management for fuel cell hybrid electric trucks_. IEEE Transactions on Vehicular Technology, 71(5), 4863-4878, doi:[10.1109/TVT.2022.3154146](https://doi.org/10.1109/TVT.2022.3154146)
