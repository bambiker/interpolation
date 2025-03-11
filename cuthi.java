import java.util.Random;

class Station {
    double x;
    double y;
    double z;

    Station(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

public class cuthi {

    public static double cuthi(double x, double y, Station[] stations, double s, double alpha) {
        // Check for exact match
        for (Station station : stations) {
            if (station.x == x && station.y == y) {
                return station.z;
            }
        }

        // Precompute distances
        double[] rSquared = new double[stations.length];
        double[][] pairwiseDistancesSquared = new double[stations.length][stations.length];

        for (int i = 0; i < stations.length; i++) {
            double dx = stations[i].x - x;
            double dy = stations[i].y - y;
            rSquared[i] = dx * dx + dy * dy;
        }

        for (int i = 0; i < stations.length; i++) {
            for (int j = 0; j < stations.length; j++) {
                double dx = stations[i].x - stations[j].x;
                double dy = stations[i].y - stations[j].y;
                pairwiseDistancesSquared[i][j] = dx * dx + dy * dy;
            }
        }

        double wtotal = 0.0;
        double[] weights = new double[stations.length];

        for (int k = 0; k < stations.length; k++) {
            if (rSquared[k] == 0) {
                return stations[k].z;
            }

            double w = 1.0 / Math.pow(rSquared[k], alpha / 2);
            double wcuthi = 1.0;

            for (int m = 0; m < stations.length; m++) {
                double adist2 = pairwiseDistancesSquared[k][m];
                double bdist2 = rSquared[m];

                if (adist2 > 0 && bdist2 > 0) {
                    double sqrtAb = Math.sqrt(adist2 * bdist2);
                    double cosgama = (adist2 + bdist2 - rSquared[k]) / (2.0 * sqrtAb);
                    double nangle = (cosgama + 1.0) / 2.0;

                    if (nangle < 0) {
                        if (nangle < -0.01) {
                            System.out.println("nangle ngnalge");
                        } else {
                            nangle = 0.0;
                        }
                    }
                    wcuthi *= Math.pow(nangle, s);
                }
            }

            weights[k] = w * wcuthi;
            wtotal += weights[k];
        }

        double res = 0.0;
        for (int k = 0; k < stations.length; k++) {
            res += stations[k].z * weights[k] / wtotal;
        }

        return res;
    }

    public static double crossvalidationCuthi(Station[] stations, double alpha) {
        double[] predictions = new double[stations.length];
        double[] actuals = new double[stations.length];

        for (int i = 0; i < stations.length; i++) {
            Station[] cv = new Station[stations.length - 1];
            int idx = 0;
            for (int j = 0; j < stations.length; j++) {
                if (j != i) {
                    cv[idx++] = stations[j];
                }
            }
            predictions[i] = cuthi(stations[i].x, stations[i].y, cv, 0.5, alpha);
            actuals[i] = stations[i].z;
        }

        // Calculate R^2 score
        double mean = 0.0;
        for (double actual : actuals) {
            mean += actual;
        }
        mean /= actuals.length;

        double ssRes = 0.0;
        double ssTot = 0.0;
        for (int i = 0; i < actuals.length; i++) {
            ssRes += Math.pow(actuals[i] - predictions[i], 2);
            ssTot += Math.pow(actuals[i] - mean, 2);
        }

        return 1.0 - (ssRes / ssTot);
    }

    public static double crossvalidationIdw(Station[] stations, double alpha) {
        double[] predictions = new double[stations.length];
        double[] actuals = new double[stations.length];

        for (int i = 0; i < stations.length; i++) {
            Station[] cv = new Station[stations.length - 1];
            int idx = 0;
            for (int j = 0; j < stations.length; j++) {
                if (j != i) {
                    cv[idx++] = stations[j];
                }
            }
            predictions[i] = idw(stations[i].x, stations[i].y, cv, alpha);
            actuals[i] = stations[i].z;
        }

        // Calculate R^2 score
        double mean = 0.0;
        for (double actual : actuals) {
            mean += actual;
        }
        mean /= actuals.length;

        double ssRes = 0.0;
        double ssTot = 0.0;
        for (int i = 0; i < actuals.length; i++) {
            ssRes += Math.pow(actuals[i] - predictions[i], 2);
            ssTot += Math.pow(actuals[i] - mean, 2);
        }

        return 1.0 - (ssRes / ssTot);
    }

    public static Station[] getFunction(int size) {
        Station[] stations = new Station[size];
        double minlat = -50, maxlat = 150;
        double minlng = -50, maxlng = 150;

        Random rand = new Random();
        for (int i = 0; i < size; i++) {
            double latrand = rand.nextDouble() * (maxlat - minlat) + minlat;
            double lngrand = rand.nextDouble() * (maxlng - minlng) + minlng;
            stations[i] = new Station(lngrand, latrand, latrand * latrand - lngrand * lngrand);
        }
        return stations;
    }

    public static double idw(double x, double y, Station[] stations, double alpha) {
        for (Station station : stations) {
            if (station.x == x && station.y == y) {
                return station.z;
            }
        }

        double wtotal = 0.0;
        double[] weights = new double[stations.length];

        for (int i = 0; i < stations.length; i++) {
            double dx = stations[i].x - x;
            double dy = stations[i].y - y;
            double r2 = dx * dx + dy * dy;
            weights[i] = 1.0 / Math.pow(r2, alpha / 2);
            wtotal += weights[i];
        }

        double res = 0.0;
        for (int i = 0; i < stations.length; i++) {
            res += stations[i].z * weights[i] / wtotal;
        }

        return res;
    }

    public static void main(String[] args) {
        int size = 60;
        Station[] stations = getFunction(size);

        double resultIdw = idw(130, 130, stations, 2);
        double resultCuthi = cuthi(0, 0, stations, 0.5, 2);

        System.out.printf("IDW result is: %f, CUTHI result is: %f, should be 0%n", resultIdw, resultCuthi);

        double cvScore = crossvalidationIdw(stations, 2);
        System.out.printf("IDW Cross validation (coefficient of determination) score is: %f (higher is better, 1.0 is perfect)%n", cvScore);
        cvScore = crossvalidationCuthi(stations, 2);
        System.out.printf("CUTHI Cross validation (coefficient of determination) score is: %f (higher is better, 1.0 is perfect)%n", cvScore);
    }
}
