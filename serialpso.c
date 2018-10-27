#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <string.h>
#include <time.h>

double nDimensions, mVelocity, nIterations, seed;
double x_min = -32.768;
double x_max = 32.768;

double ackley(double x[], double nDimensions);
void calculate(double nParticles, double result[]);

int main(int argc, char *argv[]) {
    int i;
    double nParticles;

    //Argument handling START
    for(i=1; i < argc-1; i++) {
        if (strcmp(argv[i], "-D") == 0)
            nDimensions = strtol(argv[i+1],NULL,10);
        else  if (strcmp(argv[i], "-m") == 0)
            nParticles = strtol(argv[i+1],NULL,10);
        else  if (strcmp(argv[i], "-V") == 0)
            mVelocity = strtol(argv[i+1],NULL,10);
        else  if (strcmp(argv[i], "-i") == 0)
            nIterations = strtol(argv[i+1],NULL,10);
        else  if (strcmp(argv[i], "-s") == 0)
            seed = strtol(argv[i+1],NULL,10);
    }
    if (nDimensions == 0)
        nDimensions = 2;
    if (nParticles == 0)
        nParticles = 8;
    if (mVelocity == 0)
        mVelocity = 60;
    if (nIterations == 0)
        nIterations = 1;
    if (seed == 0)
        seed = 1;
    double result[(int)nParticles];
    int j,step;
    double a,b;
    double c1, c2, rho1, rho2, w, fit;
    c1 = c2 = 1.496;
    w = 0.7298;

    //Random number generator initialization
    gsl_rng_env_setup();
    gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, time(0));

    double positions[(int)nParticles][(int)nDimensions];
    double velocities[(int)nParticles][(int)nDimensions];
    double pBestPositions[(int)nParticles][(int)nDimensions];    
    double pBestFitness[(int)nParticles];
    double gBestPosition[(int)nDimensions];
    double gBestFitness = DBL_MAX;

    //particle initialization
    for (i=0; i<nParticles; i++) {
        for (j=0; j<nDimensions; j++) {
            a = x_min + (x_max - x_min) *  gsl_rng_uniform(r);
            b = x_min + (x_max - x_min) *  gsl_rng_uniform(r);
            positions[i][j] = a;
            pBestPositions[i][j] = a;
            velocities[i][j] = (a-b) / 2.;
        }
        pBestFitness[i] = ackley(positions[i],nDimensions);
        if (pBestFitness[i] < gBestFitness) {
            memmove((void *)gBestPosition, (void *)&positions[i], sizeof(double) * nDimensions);
            gBestFitness = pBestFitness[i];
        } 
    }

    //actual calculation
    for (step=0; step<nIterations; step++) {
        for (i=0; i<nParticles; i++) {
            for (j=0; j<nDimensions; j++) {
                // calculate stochastic coefficients
                rho1 = c1 * gsl_rng_uniform(r);
                rho2 = c2 * gsl_rng_uniform(r);
                // update velocity
                velocities[i][j] = w * velocities[i][j] + \
                rho1 * (pBestPositions[i][j] - positions[i][j]) +  \
                rho2 * (gBestPosition[j] - positions[i][j]);
                // update position
                positions[i][j] += velocities[i][j];

                if (positions[i][j] < x_min) {
                    positions[i][j] = x_min;
                    velocities[i][j] = 0;
                } else if (positions[i][j] > x_max) {
                    positions[i][j] = x_max;
                    velocities[i][j] = 0;
                }

            }

            // update particle fitness
            fit = ackley(positions[i], nDimensions);
            // update personal best position?
            if (fit < pBestFitness[i]) {
                pBestFitness[i] = fit;
                // copy contents of positions[i] to pos_b[i]
                memmove((void *)&pBestPositions[i], (void *)&positions[i],
                    sizeof(double) * nDimensions);
            }
            // update gbest??
            if (fit < gBestFitness) {
                // update best fitness
                gBestFitness = fit;
                // copy particle pos to gbest vector
                memmove((void *)gBestPosition, (void *)&positions[i],
                    sizeof(double) * nDimensions);
            }
        }
    }
    printf("Result: %f\n", gBestFitness);
    gsl_rng_free(r);
}

double ackley(double x[], double nDimensions) {
    double c = 2*M_PI;
    double b = 0.2;
    double a = 20;
    double sum1 = 0;
    double sum2 = 0;
    int i;
    for (i=0; i<nDimensions; i++) {
        sum1 = sum1 + gsl_pow_2(x[i]);
        sum2 = sum2 + cos(c*x[i]);
    }
    double term1 = -a * exp(-b*sqrt(sum1/nDimensions));
    double term2 = -exp(sum2/nDimensions);
    return term1 + term2 + a + M_E;
}