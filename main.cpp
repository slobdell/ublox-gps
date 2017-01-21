#include <iostream>
#include <vector>
#include <sys/time.h>
#include "Ublox.h"

#define NULL_SENTINEL_VALUE -1000
#define EXPECTED_READ_FREQUENCY 4
#define STANDARD_DEVIATION_SAMPLES EXPECTED_READ_FREQUENCY * 10
#define TRY_READ_FREQUENCY 16
#define ONE_SECOND_AS_MICROSECONDS 1000000.0
#define SLEEP_MICROSECONDS ONE_SECOND_AS_MICROSECONDS / TRY_READ_FREQUENCY


// 110m: 0.001
// 11m: 0.0001
// 1.1m: 0.00001
struct timeval tv;
unsigned long getTimestampUs() {
    gettimeofday(&tv,NULL);
    return 1000000 * tv.tv_sec + tv.tv_usec;
}

unsigned long previousTime = -1;

float computeGPSVariance(double *values, int numSamples) {
    float mean = 0.0;
    for(int i=0; i<numSamples;i++) {
        mean += values[i];
    }
    mean /= numSamples;

    float total = 0.0;
    for(int i=0; i<numSamples;i++) {
        float value = values[i];
        total += (value - mean) * (value- mean);
    }
    return total / numSamples;
}

int main() {
	// I used below stuff to test standard deviation. In reality this doesn't work well because
	// the GPS gets a lock on the satellite and is stationary. I think you're better of just knowing
	// that UBlox accuracy is +/- 1 meter, and use that as your standard deviation...square that
	// to get variance
	/*
	 * 
    double *latSamples = new double [STANDARD_DEVIATION_SAMPLES];
    double *lonSamples = new double[STANDARD_DEVIATION_SAMPLES];
    double *altSamples = new double[STANDARD_DEVIATION_SAMPLES];
    double *latVSamples = new double [STANDARD_DEVIATION_SAMPLES];
    double *lonVSamples = new double[STANDARD_DEVIATION_SAMPLES];
    double *altVSamples = new double[STANDARD_DEVIATION_SAMPLES];
	*/


    Ublox gps;
    if(!(gps.testConnection())) {
        printf("GPS CONNECTION UNSTABLE");
    } else {
        printf("Gps connection ok\n");
    }
    // This is a hack...without this the application crashes, I never took the time to figure out
    // root cause

    std::vector<double> pos_data;
    for (int i=0; i<100;i++) {
        gps.decodeSingleMessage(Ublox::NAV_POSLLH, pos_data);
    }

    int i = 0;
    float deltaSeconds = 0;
    double previousLat = NULL_SENTINEL_VALUE;
    double previousLon = NULL_SENTINEL_VALUE;
    double previousAlt = NULL_SENTINEL_VALUE;
    double previousTowData = NULL_SENTINEL_VALUE;
    unsigned long previousTime = NULL_SENTINEL_VALUE;
    while(1) {
        gps.decodeSingleMessage(Ublox::NAV_POSLLH, pos_data);
        double towData = pos_data[0];
        double altitude = pos_data[3] / 1000;
        double latitude = pos_data[2] / 10000000;
        double longitude = pos_data[1] / 10000000;

        // no new value
        if (latitude == previousLat && longitude == previousLon && previousAlt && altitude && previousTowData == towData) {
            usleep(SLEEP_MICROSECONDS);
            continue;
        }
        // This is also a hack. I get 0 values sometimes to indicate junk, so this effectively
        // will not work along the equator
        if (abs(latitude) < 0.1 || abs(longitude) < 0.1) {
            usleep(SLEEP_MICROSECONDS);
            continue;
        }
        if (previousLat == NULL_SENTINEL_VALUE) {
            previousLat = latitude;
            previousLon = longitude;
            previousAlt = altitude;
            previousTowData = towData;
            previousTime = getTimestampUs();
        } else {
            unsigned long currentTime = getTimestampUs();
            deltaSeconds = (currentTime - previousTime) / ONE_SECOND_AS_MICROSECONDS;
            double latitudeVelocity = (latitude - previousLat) / deltaSeconds;
            double longitudeVelocity = (longitude - previousLon) / deltaSeconds;
            double altitudeVelocity = (altitude - previousAlt) / deltaSeconds;
            previousTime = currentTime;
            latVSamples[i] = latitudeVelocity;
            lonVSamples[i] = longitudeVelocity;
            altVSamples[i] = altitudeVelocity;
        }
        printf("Lat: %f, Lon: %f, Alt: %f\n", latitude, longitude, altitude);
		/*
        latSamples[i] = latitude;
        lonSamples[i] = longitude;
        altSamples[i] = altitude;
        i++;
        if(i == STANDARD_DEVIATION_SAMPLES) {
            float latVar = computeGPSVariance(latSamples, i);
            float lonVar = computeGPSVariance(lonSamples, i);
            float altVar = computeGPSVariance(altSamples, i);

            float latVVar = computeGPSVariance(latVSamples, i);
            float lonVVar = computeGPSVariance(lonVSamples, i);
            float altVVar = computeGPSVariance(altVSamples, i);
            printf("lat var: %f, lon var: %f, alt var: %f\n", latVar, lonVar, altVar);
            printf("vlat var: %f, vlon var: %f, valt var: %f\n", latVVar, lonVVar, altVVar);

            i = 0;
        }
		*/
        previousLat = latitude;
        previousLon = longitude;
        previousAlt = altitude;
        previousTowData = towData;

        usleep(SLEEP_MICROSECONDS);
    }
    return 0;
}
