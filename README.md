# pylongidyn

**pylongidyn** is a simple and lightweight Python package dedicated to the study of the energy consumption of electric road vehicles.  It is based on a _longitudinal vehicle dynamics_ model which consider the concern vehicle as constant mass object following a specific path. This vehicle is constantly submitted to _climbing_, _rolling_ and _drag resistance_ forces, as well as inertial force when variations of speed are required. This package provides simple tools to represent _ambient conditions_, _road profiles_, _continuous gearboxes_ and _electric motors_.
## Features
- Ambient air model with density calculation, useful for the calculation of drag effect.
- Road profile handling:
  - Load road data from CSV.
  - Compute slopes, climbing angles, elevations, speed limits.
- Gearbox model with continuous gear ratio and constant efficiency.
- Synchronous electric motor model with maximum torque/power characteristics.

---
## Main classes

### AmbientAir
This class is currently only useful for determining the **density of air**, which is necessary for calculating [aerodynamic drag](https://en.wikipedia.org/wiki/Drag_(physics)#Aerodynamics). For the calculation of the air density, a very simple command as follow can be used.

```python
from pylongidyn import AmbientAir

air = AmbientAir(pressure=101325, temperature=25)
print(r"Air density: rho = %.3f kg/m³" % air.density())
```
Pressure must be entered in Pa and temperature in °C.
> **Future improvements:** taking into account the effect on ambient humidity.
### Road
This class contains all the informations needed by the model to represent the properties of the path the vehicle has to follow. 