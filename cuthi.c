#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct {
    double x;
    double y;
    double z;
} Station;

double cuthi(double x, double y, Station* stations, int num_stations, double s, double alpha);
double crossvalidation(Station* stations, int num_stations, double alpha);
Station* getFunction(int size);
double idw(double x, double y, Station* stations, int num_stations, double alpha);

double cuthi(double x, double y, Station* stations, int num_stations, double s, double alpha) {
    // Check for exact match
    for (int i = 0; i < num_stations; i++) {
        if (stations[i].x == x && stations[i].y == y) {
            return stations[i].z;
        }
    }

    // Precompute distances
    double* r_squared = malloc(num_stations * sizeof(double));
    double** pairwise_distances_squared = malloc(num_stations * sizeof(double*));
    for (int i = 0; i < num_stations; i++) {
        pairwise_distances_squared[i] = malloc(num_stations * sizeof(double));
    }

    for (int i = 0; i < num_stations; i++) {
        double dx = stations[i].x - x;
        double dy = stations[i].y - y;
        r_squared[i] = dx*dx + dy*dy;
    }

    for (int i = 0; i < num_stations; i++) {
        for (int j = 0; j < num_stations; j++) {
            double dx = stations[i].x - stations[j].x;
            double dy = stations[i].y - stations[j].y;
            pairwise_distances_squared[i][j] = dx*dx + dy*dy;
        }
    }

    double wtotal = 0.0;
    double* weights = calloc(num_stations, sizeof(double));

    for (int k = 0; k < num_stations; k++) {
        if (r_squared[k] == 0) {
            free(r_squared);
            for (int i = 0; i < num_stations; i++) free(pairwise_distances_squared[i]);
            free(pairwise_distances_squared);
            free(weights);
            return stations[k].z;
        }
        
        double w = 1.0 / pow(r_squared[k], alpha/2);
        double wcuthi = 1.0;
        
        for (int m = 0; m < num_stations; m++) {
            double adist2 = pairwise_distances_squared[k][m];
            double bdist2 = r_squared[m];
            
            if (adist2 > 0 && bdist2 > 0) {
                double sqrt_ab = sqrt(adist2 * bdist2);
                double cosgama = (adist2 + bdist2 - r_squared[k]) / (2.0 * sqrt_ab);
                double nangle = (cosgama + 1.0) / 2.0;
                
                if (nangle < 0) {
                    if (nangle < -0.01) {
                        printf("nangle ngnalge\n");
                    } else {
                        nangle = 0.0;
                    }
                }
                wcuthi *= pow(nangle, s);
            }
        }
        
        weights[k] = w * wcuthi;
        wtotal += weights[k];
    }

    double res = 0.0;
    for (int k = 0; k < num_stations; k++) {
        res += stations[k].z * weights[k] / wtotal;
    }

    // Cleanup
    free(r_squared);
    for (int i = 0; i < num_stations; i++) free(pairwise_distances_squared[i]);
    free(pairwise_distances_squared);
    free(weights);

    return res;
}

double crossvalidationcuthi(Station* stations, int num_stations, double alpha) {
    double* predictions = malloc(num_stations * sizeof(double));
    double* actuals = malloc(num_stations * sizeof(double));
    
    for (int i = 0; i < num_stations; i++) {
        Station* cv = malloc((num_stations - 1) * sizeof(Station));
        int idx = 0;
        for (int j = 0; j < num_stations; j++) {
            if (j != i) {
                cv[idx++] = stations[j];
            }
        }
        predictions[i] = cuthi(stations[i].x, stations[i].y, cv, num_stations - 1, 0.5, alpha);
        free(cv);
        actuals[i] = stations[i].z;
    }

    // Calculate R^2 score
    double mean = 0.0;
    for (int i = 0; i < num_stations; i++) {
        mean += actuals[i];
    }
    mean /= num_stations;

    double ss_res = 0.0;
    double ss_tot = 0.0;
    for (int i = 0; i < num_stations; i++) {
        ss_res += pow(actuals[i] - predictions[i], 2);
        ss_tot += pow(actuals[i] - mean, 2);
    }

    free(predictions);
    free(actuals);

    return 1.0 - (ss_res / ss_tot);
}


double crossvalidationidw(Station* stations, int num_stations, double alpha) {
    double* predictions = malloc(num_stations * sizeof(double));
    double* actuals = malloc(num_stations * sizeof(double));
    
    for (int i = 0; i < num_stations; i++) {
        Station* cv = malloc((num_stations - 1) * sizeof(Station));
        int idx = 0;
        for (int j = 0; j < num_stations; j++) {
            if (j != i) {
                cv[idx++] = stations[j];
            }
        }
        predictions[i] = idw(stations[i].x, stations[i].y, cv, num_stations - 1, alpha);
        free(cv);
        actuals[i] = stations[i].z;
    }

    // Calculate R^2 score
    double mean = 0.0;
    for (int i = 0; i < num_stations; i++) {
        mean += actuals[i];
    }
    mean /= num_stations;

    double ss_res = 0.0;
    double ss_tot = 0.0;
    for (int i = 0; i < num_stations; i++) {
        ss_res += pow(actuals[i] - predictions[i], 2);
        ss_tot += pow(actuals[i] - mean, 2);
    }

    free(predictions);
    free(actuals);

    return 1.0 - (ss_res / ss_tot);
}


Station* getFunction(int size) {
    Station* stations = malloc(size * sizeof(Station));
    double minlat = -50, maxlat = 150;
    double minlng = -50, maxlng = 150;

    srand(time(NULL));
    for (int i = 0; i < size; i++) {
        double latrand = (double)rand() / RAND_MAX * (maxlat - minlat) + minlat;
        double lngrand = (double)rand() / RAND_MAX * (maxlng - minlng) + minlng;
        stations[i].x = lngrand;
        stations[i].y = latrand;
        stations[i].z = latrand * latrand - lngrand * lngrand;
    }
    return stations;
}

double idw(double x, double y, Station* stations, int num_stations, double alpha) {
    for (int i = 0; i < num_stations; i++) {
        if (stations[i].x == x && stations[i].y == y) {
            return stations[i].z;
        }
    }

    double wtotal = 0.0;
    double* weights = malloc(num_stations * sizeof(double));

    for (int i = 0; i < num_stations; i++) {
        double dx = stations[i].x - x;
        double dy = stations[i].y - y;
        double r2 = dx*dx + dy*dy;
        weights[i] = 1.0 / pow(r2, alpha/2);
        wtotal += weights[i];
    }

    double res = 0.0;
    for (int i = 0; i < num_stations; i++) {
        res += stations[i].z * weights[i] / wtotal;
    }

    free(weights);
    return res;
}

int main() {
    int size = 60;
    Station* stations = getFunction(size);
    
    double result_idw = idw(130, 130, stations, size, 2);
    double result_cuthi = cuthi(0, 0, stations, size, 0.5, 2);

    printf("IDW result is: %f, CUTHI result is: %f, should be 0\n", result_idw, result_cuthi);

    double cv_score = crossvalidationidw(stations, size, 2);
    printf("IDW Cross validation (coefficient of determination) score is: %f (higher is better, 1.0 is perfact)\n", cv_score);
    cv_score = crossvalidationcuthi(stations, size, 2);
    printf("CUTHI Cross validation (coefficient of determination) score is: %f (higher is better, 1.0 is perfact)\n", cv_score);

    free(stations);
    return 0;

}
