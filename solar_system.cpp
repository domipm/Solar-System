#include<cmath>
#include<iostream>
#include<iomanip>

//  Output filenames
#define FICHERO "cond_i.txt" //Initial conditions: Distance to sun (10e6), orbital velocity, mass
#define OUTPUT "out.txt" //Positions (X,Y) for each planet (New row for each iteration)
#define OUTPUT_EM "em.txt" //Total energy and momentum
#define OUTPUT_G "outg.txt" //Geocentric data

//  Simulation parameters
#define MAX_PLANETS 100 //Max. number of planets
#define ITER 10000 //Number of iterations to perform
#define H 0.001 //Time variation delta T

//  Useful constants
#define SUN_MASS 1.989e30 //Solar mass (kg)
#define C 1.496e11 //Astronomical Units (m)
#define G 6.67e-11 //Gravitational constant
#define FACT 5022004.955 //Time conversion factor (1s' = fact s)
#define PI 3.14159265358979323846 //PI

using namespace std;

//  Mass, position, velocity, acceleration, W function arrays (coordinates x,y)
float m[MAX_PLANETS], x[MAX_PLANETS], y[MAX_PLANETS], vx[MAX_PLANETS], vy[MAX_PLANETS], ax[MAX_PLANETS], ay[MAX_PLANETS], wx[MAX_PLANETS], wy[MAX_PLANETS];

//  Arrays to calculate period
float P[MAX_PLANETS] = {0};        // Period of each planet
float yo[MAX_PLANETS] = {0};       // Initial vertical position
bool vuelta[MAX_PLANETS] = {0};    // Whether an orbit has been completed

//  Total energy and momentum of the system
long double E = 0, M = 0;

//  Time
float t = 0;
//  Number of planets
int n = 0;

int main() {


    //  Read initial conditions from file
    FILE *data;
    data = fopen(FICHERO, "r");
    //  Write data into arrays and count number of rows (planets)
    while(fscanf(data, "%f\t%f\t%f", &x[n], &vy[n], &m[n]) != EOF) n++;
    fclose(data);

    //  Rescale

    for (int i = 0; i < n; i++) {
        //Rescale positions
        x[i] *= 1e6*1000/C;
        y[i] *= 1e6*1000/C;
        //Rescale mass
        m[i] *= 1e24/SUN_MASS;
        //Rescale velocity
        vx[i] *= 1000 / C * FACT;
        vy[i] *= 1000 / C * FACT;
    }

    //  Output files

    //  Coordinates
    FILE *out;
    out = fopen(OUTPUT, "w");       // Heliocentric
    FILE *outg;
    outg = fopen(OUTPUT_G, "w");    // Geocentric
    //  Energy and Momentum
    FILE *em;
    em = fopen(OUTPUT_EM, "w"); 

    // Initial loop

    //  Write initial (x,y) values of all planets
    for (int i = 0; i < n; i++) fprintf(out, "%f\t%f\t", x[i], y[i]);                   // Heliocentric coordinates
    for (int i = 0; i < n; i++) fprintf(outg, "%f\t%f\t", (x[i]-x[3]), (y[i]-y[3]));    // Geocentric coordinates

    //  Calculate initial accelerations
    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {

            if (i != j) {

                float R = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) );    // Modulus distance vector R_i,j
                ax[i] += - m[j] * (x[i]-x[j]) / (R*R*R);                                // Initial acceleration (X-axis)

            }

        }

        ay[i] = 0; // Assume initial vertical acceleration null (planets aligned on X-axis)

    }

    //  Initial energy and angular momentum
    long double T = 0, V = 0;
    for (int i = 0; i < n; i++) {

        M += m[i] * (x[i]*vy[i] - y[i]*vx[i]);                                                  // Angular momentum
        E += 0.5 * m[i] * sqrt(vx[i]*vx[i] + vy[i]*vy[i]) * sqrt(vx[i]*vx[i] + vy[i]*vy[i]);    // Kinetic energy

        for (int j = 0; j < n; j++) {

            if (i != j) {

                float R = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) );
                E += - m[i] * m[j] / (2*R);                                                     // Potential energy

            }

        }

    }

    fprintf(em, "%.7f\t%.15Lf\t%.15Lf\n", t, E, M); // Write time, energy and momentum onto file
    
    //  Verlet algorithm

    for (int k = 0; k < ITER; k++) {

        t += H; //Update time
        M = 0;  //Initial energy of each iteration zero
        E = 0;  //Initial momentum of each iteration zero

        //  Begin new line on files for each iteration
        fprintf(out, "\n");
        fprintf(outg, "\n");


        //  Positions

        
        for (int i = 0; i < n; i++) {

            //  Previous vertical position (used to calculate period)
            yo[i] = y[i];

            //  W function (with velocity and acceleration at time T)
            wx[i] = vx[i] + H/2 * ax[i];
            wy[i] = vy[i] + H/2 * ay[i];

            //  New positions
            x[i] += H * vx[i] + H*H/2 * ax[i];
            y[i] += H * vy[i] + H*H/2 * ay[i];

            //  Write new position onto files
            fprintf(out, "%f\t%f\t", x[i], y[i]);
            fprintf(outg, "%f\t%f\t", (x[i]-x[3]), (y[i]-y[3]));

        }

        //  Accelerations

        for (int i = 0; i < n; i++) {

            //  Initialize to zero for each iteration
            ax[i] = 0;
            ay[i] = 0;

            for (int j = 0; j < n; j++) {

                if (i != j) {

                    float R = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) ); // Position vector R_i,j

                    //  Accelerations
                    ax[i] += - m[j] * (x[i] - x[j]) / (R*R*R);
                    ay[i] += - m[j] * (y[i] - y[j]) / (R*R*R);

                }
            }

        }

        //  Velocities

        for (int i = 0; i < n; i++) {

            vx[i] = wx[i] + H/2 * ax[i];
            vy[i] = wy[i] + H/2 * ay[i];

        }

        //  Planet orbital periods

        for (int i = 0; i < n; i++) {
            //  Check if ...
            if (x[i] > 0 && y[i]*yo[i] < 0 && k > 100 && i != 0 && vuelta[i] == 0) {
                vuelta[i] = 1;
                P[i] = k * H * 5022004.955 * 1/3600 * 1/24;
                cout << "PLANET PERIOD " << i << " " << P[i] << " DAYS." << endl;
            }
        }

        //  Energy and momentum

        for (int i = 0; i < n; i++) {

            M += m[i] * (x[i]*vy[i] - y[i]*vx[i]);          // Total angular momentum

            E += 0.5 * m[i] * (vx[i]*vx[i] + vy[i]*vy[i]);  // Kinetic energy

            for (int j = 0; j < n; j++) {

                if (i != j) {

                    float R = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) );
                    E += - m[i] * m[j] / (2*R);             //Potential energy

                }

            }
        } 

        fprintf(em, "%.7f\t%.15Lf\t%.15Lf\n", t, E, M); // Write energy and momentum onto file

    }

    //  Close all files
    fclose(out);
    fclose(outg);
    fclose(em);

    //  Show info on screen
    cout << "SIMULATION TIME: " << setprecision(7) << t * FACT * 1/3600 * 1/24 << " DAYS" << endl; // Real-world time of simulation

    return 0;

}
