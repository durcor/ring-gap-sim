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
// time step (1 minute)
int t = 60;

void orbit(int num_particles, int mimas_multiplier, int num_time_steps)
{
	// Plot the orbit of Mimas and massless ring particles around Saturn.

	// Parameters:
	// num_particles (int): The number of particles to simulate.
	// mimas_multiplier (int): Multipler of Mimas' mass.
	// num_time_steps (int): The number of time steps to simulate.
	int INNER_LIM = 25000 * pow(10, 3);
	int OUTER_LIM = 75000 * pow(10, 3);
	// Entire range of Saturn's rings
	// int INNER_LIM = 74500 * pow(10, 3);
	// int OUTER_LIM = 140220 * pow(10, 3);
	// Range around the Cassini Division including B and A rings.
	// INNER_LIM = 92000 * 10**3
	// OUTER_LIM = 136780 * 10**3

	double p[num_particles][2];
	double v[num_particles][2];
	for (int i = 0, j = INNER_LIM; j < OUTER_LIM;
	     j += (OUTER_LIM - INNER_LIM) / num_particles, ++i) {
		p[i][0] = (double)j;
		p[i][1] = 0.0;
		v[i][0] = 0.0;
		v[i][1] = sqrt(G * M / p[i][0]);
	}
	// omega used for calculating position of mimas
	double w = sqrt(G * M / pow(d, 3));
	double xm, ym;
	for (int n = 0; n < num_time_steps * t; n += t) {
		xm = d * cos(n * w);
		ym = d * sin(n * w);
		// Plot Saturn and Mimas
		for (int i = 0; i < num_particles; ++i) {
			double rs = sqrt(pow(p[i][0], 2) + pow(p[i][1], 2));
			double rm = sqrt(pow(p[i][0] - xm, 2) +
					 pow(p[i][1] - ym, 2));
			double ax = -G * M * p[i][0] / pow(rs, 3) -
				    G * m * mimas_multiplier * (p[i][0] - xm) /
					    pow(rm, 3);
			double ay = -G * M * p[i][1] / pow(rs, 3) -
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
	if (argc != 4) {
		fprintf(stderr,
			"Usage: %s <num_particles> <mimas_multiplier> <num_time_steps>\n",
			argv[0]);
		return 1;
	}
	orbit(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
	return 0;
}
