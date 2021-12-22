//  Tyler Kaminski
//  Adham Abdelwahab
//  Vincent Andrade
//  Professor Lukic
//  PEP 351
//
//  I pledge my honor that I have abided by the Stevens Honor System.
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Constants for calculations
float G = 6.67430e-11;
// Mass of the central body (Saturn)
float M = 5.6834e26;
// Mass of Mimas
float m = 3.7493e19;
// Distance from Saturn to Mimas
float d = 1.8552e8;
// time step (1 second)
int t = 1;

// Plot the orbit of Mimas and massless ring particles around Saturn.
// Parameters:
// num_particles: The number of particles to simulate.
// mimas_multiplier: Multiplier of Mimas' mass.
// num_time_steps: The number of time steps to simulate.
void orbit(int num_particles, int mimas_multiplier, int num_time_steps)
{
	// Range of entirety of Saturn's rings.
	// int INNER_LIM = 74500 * pow(10, 3);
	// int OUTER_LIM = 140220 * pow(10, 3);
	// Range around the Cassini Division including B and A rings.
	int INNER_LIM = 92000 * pow(10, 3);
	int OUTER_LIM = 136780 * pow(10, 3);
	// Range of only the Cassini Division.
	// int INNER_LIM = 117580 * pow(10, 3);
	// int OUTER_LIM = 122170 * pow(10, 3);

	float p[num_particles][2];
	float v[num_particles][2];
	for (int i = 0, j = INNER_LIM; j < OUTER_LIM;
	     j += (OUTER_LIM - INNER_LIM) / num_particles, ++i) {
		p[i][0] = j;
		p[i][1] = 0.0;
		v[i][0] = 0.0;
		v[i][1] = sqrt(G * M / p[i][0]);
	}
	// omega used for calculating position of mimas
	float w = sqrt(G * M / pow(d, 3));
	float xm = 0, ym = 0;
	for (int n = 0; n < num_time_steps * t; n += t) {
		xm = d * cos(n * w);
		ym = d * sin(n * w);
		for (int i = 0; i < num_particles; ++i) {
			float rs = sqrt(pow(p[i][0], 2) + pow(p[i][1], 2));
			float rm = sqrt(pow(p[i][0] - xm, 2) +
					pow(p[i][1] - ym, 2));
			float ax = -G * M * p[i][0] / pow(rs, 3) -
				   G * m * mimas_multiplier * (p[i][0] - xm) /
					   pow(rm, 3);
			float ay = -G * M * p[i][1] / pow(rs, 3) -
				   G * m * mimas_multiplier * (p[i][1] - ym) /
					   pow(rm, 3);
			p[i][0] += v[i][0] * t + ax / 2 * pow(t, 2);
			p[i][1] += v[i][1] * t + ay / 2 * pow(t, 2);
			v[i][0] += ax * t;
			v[i][1] += ay * t;
		}
	}
	printf("%f, %f\n", xm, ym);
	for (int i = 0; i < num_particles; ++i)
		printf("%f, %f\n", p[i][0], p[i][1]);
}

int main(int argc, char **argv)
{
	orbit(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
	return 0;
}
