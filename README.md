# Simulating the Ring Gaps of Saturn

# Pseudo-code
1. Store the x-y coordinates and masses of Saturn's moons
2. Store x-y coordinates of particles in Saturn's rings
3. Use Euler's method/RK4 to simulate the motion of each particle around Saturn in the presence of Saturn's moons

Preliminary Research:
G (Gravitational Constant) = 6.67E-11 [m^3/(kg*s^2)] 
M (Mass of Saturn) = 5.683E26 [kg]
m (Mass of Mimas) = 3.7493E19 [kg]  *Boost by 10x or so to accelerate effects
d (Distance b/w Saturn and Mimas) = 1.8552E8 [m]
The Cassini Division contains nearly no particles and is the most prominent ‘gap’ in the ring system easily seen from earth. It extends from 117,580 km to 122,170 km from the center of Saturn. (Distances to simulate particles within)
d1 = 1.1758E8 [m]
d2 = 1.2217E8 [m]
