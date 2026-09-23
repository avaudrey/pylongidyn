# pylongidyn
> A simple, lightweight and modular longitudinal vehicle dynamics model in Python

## Introduction
**pylongidyn** is a small Python module that estimates the **energy consumption** and the **instantaneous power demand** of electric road vehicles. It was first developed to size battery-electric and fuel-cell heavy-duty trucks, but it can be used for any road vehicle.

The module is based on a _**longi**tudinal vehicle **dyn**amics_ model (see references [1, p. 39] and [2, p. 15]). The vehicle is treated as a single mass moving along a road. At every moment it is subject to _climbing_, _rolling_ and _aerodynamic drag_ resistances, plus an inertial force whenever its speed changes. The motive force needed to overcome all these forces, multiplied by the vehicle speed, gives the instantaneous mechanical power required at the wheels. That power is then traced back through the transmission and the electric motor.

Most tools use a prescribed _speed-versus-time_ drive cycle. pylongidyn instead describes the route as **distance, altitude and speed limit**. The simulated vehicle tries to drive at the speed limit, within the limits of its powertrain. The speed profile is therefore a **result** of the simulation, and it depends on the vehicle: an under-powered truck slows down on a steep climb, while a more powerful one keeps its speed.

The whole module is a **single Python file**. It is easy to read, easy to modify, and fast: a 1000 km route with a 1 s time step runs in about one second on a laptop.

## Main features
- Road profiles loaded from a simple CSV file (distance, altitude, speed limit), with optional regularization of unrealistic slopes.
- Vehicle resistive forces: climbing, rolling (with a user-defined rolling coefficient) and aerodynamic drag (with a user-defined air density).
- Electric motor with a constant-torque / constant-power characteristic and a constant energy efficiency.
- Continuously variable transmission (CVT-like) with a constant energy efficiency, applied according to the direction of the power flow.
- Regenerative braking, limited by the motor capability and shared with the friction brakes.
- Time series of speed, acceleration, electrical power, motor operating points and friction braking power, ready to be exported to other tools (for instance, fuel-cell or battery models).

## Installation
pylongidyn is a single Python file. Copy it into your project:
```bash
your_project/
    pylongidyn.py
```
and import it:
```python
import pylongidyn as pld
```
The only requirements are **NumPy** and **SciPy**.

## Physical model
### Resistive and motive forces
For a vehicle of mass $m$, moving at speed $v$ on a road of slope $s$, the resistive forces are

$$F_{res} = \underbrace{m\,g\,s}_{\text{climbing}} + \underbrace{C_{rr}\,m\,g\,\sqrt{1-s^2}}_{\text{rolling}} + \underbrace{\tfrac{1}{2}\,\rho\,C_d\,A\,v^2}_{\text{drag}}$$

where $C_{rr}$ is the rolling coefficient, $\rho$ the air density, $C_d$ the drag coefficient and $A$ the frontal area. The road slope $s = \Delta z / \Delta x$ is used as the sine of the road angle, which is accurate for realistic road slopes.

The motive force needed to reach an acceleration $a$ is

$$F = m\,(1 + k_{rot})\,a + F_{res}$$

where $k_{rot}$ is a correction coefficient that accounts for the inertia of the rotating parts. The torque at the wheels is $T_w = F\,r_w$, with $r_w$ the wheel radius.

### Transmission
The gearbox links the motor (inlet) to the wheels (outlet): $\omega_{out} = r\,\omega_{in}$. Its gear ratio $r$ varies linearly with the motor speed, from its value at zero speed, $r_0$ (`minimum_gear_ratio`), to its value at the maximum motor speed, $r_1$ (`maximum_gear_ratio`):

$$r = r_0 + (r_1 - r_0)\,\frac{\omega_{in}}{\omega_{in,\max}}$$

The gearbox efficiency $\eta_{gb}$ always reduces the transmitted power:
- in **traction** ($T_{out} > 0$), $T_{in} = T_{out}\,r/\eta_{gb}$;
- in **regenerative braking** ($T_{out} < 0$), $T_{in} = T_{out}\,r\,\eta_{gb}$.

### Electric motor
The motor delivers its maximum torque $T_{\max}$ up to its base speed $\omega_b = P_{\max}/T_{\max}$, and a constant maximum power $P_{\max}$ above it. Its electrical power, counted **positive when consumed** and **negative when recovered**, is

$$P_{el} = \begin{cases} P_m / \eta_m & \text{if } P_m \ge 0 \text{ (traction)} \\ P_m\,\eta_m & \text{if } P_m < 0 \text{ (regenerative braking)} \end{cases}$$

where $P_m = T_{in}\,\omega_{in}$ is the mechanical power at the motor shaft and $\eta_m$ its efficiency.

### Simulation algorithm
The simulation uses an explicit time-stepping scheme. At each time step, all forces are evaluated at the current speed, and the vehicle speed is compared to the local speed limit:

1. **Speed below the limit**: the vehicle accelerates at the _target acceleration_ if the motor can provide the required torque. Otherwise, it accelerates (or even decelerates, on a steep climb) with the maximum torque available. If the limit would be exceeded during the step, the acceleration is reduced so that the vehicle reaches the limit exactly.
2. **Speed above the limit**: the vehicle decelerates at the _target deceleration_, without overshooting the limit.
3. **Speed equal to the limit**: the motor provides the torque needed to maintain the speed against the resistive forces. If it cannot, the vehicle slows down.

Whenever the torque required at the wheels is **negative** (braking, or holding the speed downhill), a fraction of it (`regenerative_braking_fraction`) is sent to the motor, within the motor torque limit. The friction brakes provide the rest, so the vehicle motion is not affected. The last time step is shortened so that the vehicle stops exactly at the end of the road.

A similar approach has been published, for instance, in reference [3].

## Main classes
The examples below build, step by step, a model of a 40 t battery-electric truck.

### AmbientAir
Represents the vehicle surroundings. It is currently used to compute the **air density**, from the [ideal gas law](https://en.wikipedia.org/wiki/Ideal_gas_law). Pressure is in Pa and temperature in °C.
```python
import pylongidyn as pld

air = pld.AmbientAir(pressure=101325., temperature=15.)
print("Air density: %.3f kg/m³" % air.density())
```

### Road
Handles road profiles: slopes, elevation, speed limits. Road data are loaded from a CSV file with a header line and three columns:

| distance_m | altitude_m | speed_limit_kmh |
|------------|------------|-----------------|
| 0          | 120        | 90              |
| 500        | 130        | 90              |
| 1000       | 140        | 70              |
| ...        | ...        | ...             |

- The distance step can be irregular.
- Between two points, the altitude varies linearly (constant slope) and the speed limit is constant, equal to the value given at the first point of the segment.
- Altitudes are converted into **elevations** relative to the starting point.
- Rows with a **zero speed limit** (stops) are currently removed, and a warning is issued.

```python
road = pld.Road()
road.load_road_data("myroad.csv")

print("Road length: %.1f km" % (1e-3 * road.road_length()))
print("Maximum elevation: %.0f m" % road.maximum_elevation())
print("Minimum trip duration: %s" % road.minimum_duration(format="hh:mm:ss"))
# Speed limits are given in km/h in the file, but handled in m/s
print("Speed limit at 750 m: %.1f km/h" % (3.6 * road.speed_limit(750.)))
print("Slope at 750 m: %.2f %%" % (1e2 * road.slope(750.)))

# Heavy-duty trucks are limited to 90 km/h, whatever the road
road.limit_speed_to_a_maximum_value(90.)
```
Some elevation profiles, especially those extracted from digital elevation models, contain unrealistic slopes. They can be **regularized**: every slope is limited to `max_slope`, and the elevation is rebuilt accordingly. This option must be set **before** loading the data:
```python
regularized_road = pld.Road()
regularized_road.regularization = True
regularized_road.max_slope = 0.12  # 12 %
regularized_road.load_road_data("myroad.csv")
regularized_road.export_regularized_profile()  # writes myroad-regularized.csv
```

### Vehicle
Represents the vehicle body. Mass is in kg, frontal area in m² and wheel radius in m.
```python
truck = pld.Vehicle()
truck.mass = 40e3
truck.rotating_mass_correction_coefficient = 0.05
truck.frontal_area = 10.
truck.drag_coefficient = 0.6
truck.wheels_radius = 0.5

slope = 0.02      # 2 %
speed = 80 / 3.6  # m/s

print("Drag force: %.0f N" % truck.drag_resistance(speed, air_density=air.density()))
print("Climbing resistance: %.0f N" % truck.climbing_resistance(slope))
# Default rolling coefficient is 0.01
print("Rolling resistance: %.0f N" % truck.rolling_resistance(slope, rolling_coefficient=5e-3))
print("Total resistance: %.0f N" % truck.resistance_force(slope, speed, rolling_coefficient=5e-3))
print("Motive force for +0.5 m/s²: %.0f N" % truck.motive_force(0.5, slope, speed, rolling_coefficient=5e-3))
print("Maximum speed with 300 kW: %.1f km/h" % (3.6 * truck.power_to_maximum_speed(300e3, slope, rolling_coefficient=5e-3)))
```
The optional keyword arguments `rolling_coefficient` and `air_density` are accepted by all force and power methods.

### SynchronousElectricMotor
The motor model gives the maximum torque the motor can deliver at any speed, and converts mechanical power into electrical power with a constant efficiency.
```python
motor = pld.SynchronousElectricMotor()
# Set the maximum frequency first: the base speed, deduced from the maximum
# torque and power, must always remain lower than the maximum speed
motor.maximum_frequency = 4500.  # rpm
motor.maximum_torque = 5e3       # N.m, available from zero to base speed
motor.maximum_power = 400e3      # W, available above base speed
motor.energy_efficiency = 0.9

print("Base frequency: %.0f rpm" % motor.base_frequency)
print("Max torque at 2000 rpm: %.0f N.m" % motor.maximum_operating_torque(frequency=2000.))
```
Rotation rates are given either as a `frequency`, in rpm, or as a (rotational) `speed`, in rad/s. The functions `rot_speed_from_rpm` and `rot_rpm_from_speed` convert from one to the other.

### LinearContinuousGearbox
The gearbox converts the motor rotation into the wheel rotation, with a gear ratio (outlet to inlet speed) that varies linearly with the motor speed. A simple way to size it:
1. at zero speed, the maximum motor torque must provide the maximum torque required at the wheels;
2. at the maximum motor speed, the vehicle must reach its maximum speed.

```python
gearbox = pld.LinearContinuousGearbox()
gearbox.maximum_inlet_frequency = motor.maximum_frequency

maximum_wheel_torque = 3.5e4        # N.m
maximum_vehicle_speed = 100. / 3.6  # m/s
gearbox.minimum_gear_ratio = motor.maximum_torque / maximum_wheel_torque
gearbox.maximum_gear_ratio = (maximum_vehicle_speed / truck.wheels_radius) / motor.maximum_speed
gearbox.energy_efficiency = 0.97

# Setting any inlet or outlet speed updates all the others
gearbox.inlet_frequency = 1000.  # rpm
print("Gear ratio: %.4f" % gearbox.gear_ratio())
print("Vehicle speed: %.1f km/h" % (3.6 * gearbox.outlet_speed * truck.wheels_radius))
# Same for torques
gearbox.inlet_torque = 2000.  # N.m
print("Wheel torque: %.0f N.m" % gearbox.outlet_torque)
print("Power lost in the gearbox: %.2f kW" % (1e-3 * gearbox.dissipated_power()))
```
Here, "minimum" and "maximum" gear ratios mean the ratios at **zero** and **maximum** motor speed. The first one may be larger than the second one.

The vehicle cannot drive faster than `gearbox.maximum_outlet_speed() * wheels_radius`. Asking for a higher speed raises a `ValueError`, so the road speed limits must stay below this value.

### VehicleDynamicsModel
Combines the previous objects and runs the simulation. The last parameters are the _target acceleration_ and _target deceleration_, used whenever the vehicle speed has to change, and the time step.
```python
model = pld.VehicleDynamicsModel(
    gearbox=gearbox,
    electric_motor=motor,
    road=road,
    vehicle=truck,
    average_acceleration=0.5,  # m/s²
    average_deceleration=0.5,  # m/s²
    time_step=1.,              # s
    regenerative_braking_fraction=1.,
)

# Initial speed in km/h; optional arguments are passed to the Vehicle methods
model.run_calculation(initial_speed=0., rolling_coefficient=5e-3, air_density=air.density())

print("Trip duration: %.0f s" % model.trip_duration())
print("Average speed: %.1f km/h" % (3.6 * model.speeds()["avg"]))
print("Energy with regenerative braking: %.2f kWh" % (1e-3 * model.consumed_energy()))
print("Energy without regenerative braking: %.2f kWh" % (1e-3 * model.consumed_energy(regenerative_braking=False)))
```
The results are stored in the `model.results` dictionary, as lists of equal length:

| Key | Unit | Description |
|-----|------|-------------|
| `time` | s | Time since departure |
| `distance` | m | Distance from the starting point |
| `speed` | m/s | Vehicle speed |
| `acceleration` | m/s² | Vehicle acceleration |
| `power` | W | Electrical power at the motor terminals (positive when consumed, negative when recovered) |
| `motor_torque` | N·m | Motor torque |
| `motor_speed` | rad/s | Motor rotational speed |
| `friction_braking_power` | W | Power dissipated by the friction brakes (positive) |

Each value, except the first one at $t = 0$, applies during the time step that **ends** at the corresponding time. For instance, the power profile can be exported for a fuel-cell or battery model:
```python
import numpy as np

np.savetxt("power_profile.csv",
           np.c_[model.results["time"], model.results["power"]],
           delimiter=",", header="time_s,power_W", comments="")
```

## Current limitations
- **No anticipation**: the vehicle starts to brake only once it has entered a section with a lower speed limit.
- **Stops are ignored** (rows with a zero speed limit are removed from the road profile).
- Constant motor and gearbox efficiencies, and a constant rolling coefficient.
- No auxiliary power consumption.
- If the vehicle cannot climb a slope (motor torque too low), the simulation never ends: the road profile and the motor sizing must be checked beforehand.
- The energy source (battery, fuel cell) is not modelled: the output is the electrical power demanded at the motor terminals.

## Roadmap
- Anticipation of speed-limit decreases and stops, with a braking envelope computed along the road.
- Stop durations read from an optional fourth column of the road file.
- Air density taken from an `AmbientAir` object, possibly varying with altitude.
- Speed-dependent target acceleration.
- Auxiliary power consumption.
- Efficiency maps for the motor, and gearboxes with discrete gears.
- Energy source models (battery, PEM fuel-cell system) and energy management strategies.

## Changelog
### September 2026: corrections affecting the results
Several bugs have been fixed. **Results obtained with earlier versions should be recomputed.**
- The optional arguments of `run_calculation` (`rolling_coefficient`, `air_density`) were ignored at constant speed, where the default rolling coefficient (0.01) was used instead. On a flat road at 80 km/h with a rolling coefficient of 0.005, the energy was overestimated by about 50 %.
- `consumed_energy` counted the time step twice. Results were only correct with a 1 s time step. The energy is now integrated with the actual duration of each step.
- Regenerative braking: the motor and gearbox efficiencies were applied in the wrong direction (energy recovered was overestimated), and the recovered power was not limited by the motor capability. Braking is now shared between the motor and the friction brakes.
- The acceleration was set to zero during the time step in which the vehicle reached the speed limit.
- The gear ratio was computed inconsistently when the outlet speed was imposed, so the gearbox did not conserve power (about 1 % error).
- `Road`: the regularization returned absolute altitudes instead of relative elevations, `slope` and `elevation` failed at the end of the road, the segment search was slow on long profiles, a missing file produced an unclear error, and the last time step was not recorded.

With these corrections, the energy computed on a flat road at constant speed matches the analytical value.

## License
This package is released under the GNU General Public License v3. See the `LICENSE` file or http://www.gnu.org/licenses/.

## References
[1] Jazar, R. N., [_Vehicle Dynamics: Theory and Applications_](https://link.springer.com/book/10.1007/978-3-031-74458-7), Springer, 2008.

[2] Minaker, B. P., [_Fundamentals of Vehicle Dynamics and Modelling_](https://www.wiley.com/en-us/Fundamentals+of+Vehicle+Dynamics+and+Modelling%3A+A+Textbook+for+Engineers+With+Illustrations+and+Examples-p-9781118980095), Wiley, 2019.

[3] Xun, Q., Murgovski, N., & Liu, Y. (2022). _Joint component sizing and energy management for fuel cell hybrid electric trucks_. IEEE Transactions on Vehicular Technology, 71(5), 4863-4878, doi:[10.1109/TVT.2022.3154146](https://doi.org/10.1109/TVT.2022.3154146)
