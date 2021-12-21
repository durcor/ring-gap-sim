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
#include <vector>

// Constants for calculations
double G = 6.67430e-11;
// Mass of the central body (Saturn)
double M = 5.6834e26;
// Mass of Mimas
double m = 3.7493e19;
// Distance from Saturn to Mimas
double d = 1.8552e8;
// time step (1 minute)
int t = 1;

void orbit(int num_particles, int mimas_multiplier, int num_time_steps)
{
	// Plot the orbit of Mimas and massless ring particles around Saturn.

	// Parameters:
	// num_particles (int): The number of particles to simulate.
	// mimas_multiplier (int): Multipler of Mimas' mass.
	// num_time_steps (int): The number of time steps to simulate.
	// int INNER_LIM = 25000 * pow(10, 3);
	// int OUTER_LIM = 75000 * pow(10, 3);
	// Entire range of Saturn's rings
	// int INNER_LIM = 74500 * pow(10, 3);
	// int OUTER_LIM = 140220 * pow(10, 3);
	// Range around the Cassini Division including B and A rings.
	int INNER_LIM = 92000 * pow(10, 3);
	int OUTER_LIM = 136780 * pow(10, 3);

	std::vector<std::pair<double, double> > p(num_particles);
	std::vector<std::pair<double, double> > v(num_particles);
	for (int i = 0, j = INNER_LIM; j < OUTER_LIM;
	     j += (OUTER_LIM - INNER_LIM) / num_particles, ++i) {
		p[i].first = j;
		p[i].second = 0.0;
		v[i].first = 0.0;
		v[i].second = sqrt(G * M / p[i].first);
	}
	// omega used for calculating position of mimas
	double w = sqrt(G * M / pow(d, 3));
	double xm = 0, ym = 0;
	for (int n = 0; n < num_time_steps * t; n += t) {
		xm = d * cos(n * w);
		ym = d * sin(n * w);
		// Plot Saturn and Mimas
		for (int i = 0; i < num_particles; ++i) {
			double rs =
				sqrt(pow(p[i].first, 2) + pow(p[i].second, 2));
			double rm = sqrt(pow(p[i].first - xm, 2) +
					 pow(p[i].second - ym, 2));
			double ax = -G * M * p[i].first / pow(rs, 3) -
				    G * m * mimas_multiplier *
					    (p[i].first - xm) / pow(rm, 3);
			double ay = -G * M * p[i].second / pow(rs, 3) -
				    G * m * mimas_multiplier *
					    (p[i].second - ym) / pow(rm, 3);
			p[i].first += v[i].first * t + ax / 2 * pow(t, 2);
			p[i].second += v[i].second * t + ay / 2 * pow(t, 2);
			v[i].first += ax * t;
			v[i].second += ay * t;
		}
	}
	printf("%f, %f\n", xm, ym);
	for (int i = 0; i < num_particles; ++i)
		printf("%f, %f\n", p[i].first, p[i].second);
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
