/* 
 * File:   main.cpp
 * Author: Israt
 *
 * Created on October 21, 2016, 7:58 PM
 */

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include "sgp4coord.h"
#include "sgp4unit.h"
#include "sgp4io.h"
#include "sgp4coord.h"

using namespace std;

bool isChanged(double arr1[], double arr2[]) {
    if (arr1[0] == arr2[0] && arr1[1] == arr2[1] && arr1[2] == arr2[2]) {
        return false;
    } else {
        return true;
    }
}

int main(int argc, char* argv[]) {


    //SET UP SOME VARIABLES
    //<editor-fold defaultstate="collapsed" desc="site details and some location variables ">
    double siteLat, siteLon, siteAlt, siteLatRad, siteLonRad;
    //ENTER SITE DETAILS HERE : I put white sands location in decimal which is 
    //White Sands Ground Terminal (WSGT) 32.5007°N 106.6086°W
    siteLat = 32.512027777777774; //+North (Austin ->30.25N)
    siteLon = -105.3886388888888; //+East (Austin ->-97.75) 
    siteAlt = 0.15; //km 
    siteLatRad = siteLat * pi / 180.0;
    siteLonRad = siteLon * pi / 180.0;
    printf("Latitude-- %f, Longitude -- %f \n", siteLatRad, siteLonRad);


    double latlongh[3]; //lat, long in rad, h in km above ellipsoid
    double tdrs8_latlongh[3]; //lat, long in rad, h in km above ellipsoid
    double tdrs9_latlongh[3]; //lat, long in rad, h in km above ellipsoid
    double tdrs10_latlongh[3]; //lat, long in rad, h in km above ellipsoid
    //</editor-fold>

    //ENTER TWO-LINE ELEMENT HERE
    //<editor-fold defaultstate="collapsed" desc="TLE">
    //ISS TLE
    char longstr1[] = "1 25544U 98067A   16308.53695251  .00002415  00000-0  44666-4 0  9996";
    char longstr2[] = "2 25544  51.6423  83.0189 0007076 169.6165 314.0897 15.53516898 26669";

    char tdrs8_longstr1[] = "1 26388U 00034A   16307.79454727 -.00000226  00000-0  00000+0 0  9996";
    char tdrs8_longstr2[] = "2 26388   7.1086  57.8047 0009210 181.8074 178.1014  1.00276995 59941";

    char tdrs9_longstr1[] = "1 27389U 02011A   16308.10657816 -.00000095  00000-0  00000+0 0  9999";
    char tdrs9_longstr2[] = "2 27389   4.7843  82.0463 0019183 220.3119 126.8010  1.00272302 55218";

    char tdrs10_longstr1[] = "1 27566U 02055A   16308.52882761  .00000083  00000-0  00000+0 0  9995";
    char tdrs10_longstr2[] = "2 27566   4.6417  59.3907 0011334 199.5787 160.2494  1.00274596 51007";

    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="R and V vector">
    double ro[3]; // R -position vector for ISS
    double vo[3]; // V-  velocity vector for ISS

    double tdrs8_ro[3]; // R -position vector for tdrs8
    double tdrs8_vo[3]; // V-  velocity vector for tdrs8

    double tdrs9_ro[3]; // R -position vector for tdrs9
    double tdrs9_vo[3]; // V-  velocity vector for tdrs9

    double tdrs10_ro[3]; // R -position vector for tdrs10
    double tdrs10_vo[3]; // V-  velocity vector for tdrs10

    double ground_ro[3]; // R -position vector 
    double prevground_ro[3]; // R -position vector 
    double ground_vo[3]; // V-  velocity vector

    double recef[3]; //R- vector for ISS in Earth Centered Earth Fixed frame or TEME frame
    double prevrecef[3]; //R- vector for ISS in Earth Centered Earth Fixed frame or TEME frame
    double vecef[3]; //V- vector for ISS in Earth Centered Earth Fixed frame or TEME frame

    double tdrs8_recef[3]; //R- vector for TDRS8 in Earth Centered Earth Fixed frame or TEME frame
    double prevtdrs8_recef[3]; //R- vector for TDRS8 in Earth Centered Earth Fixed frame or TEME frame
    double tdrs8_vecef[3]; //V- vector for TDRS8 in Earth Centered Earth Fixed frame or TEME frame

    double tdrs9_recef[3]; //R- vector for TDRS9 in Earth Centered Earth Fixed frame or TEME frame
    double prevtdrs9_recef[3]; //R- vector for TDRS9 in Earth Centered Earth Fixed frame or TEME frame
    double tdrs9_vecef[3]; //V- vector for TDRS9 in Earth Centered Earth Fixed frame or TEME frame

    double tdrs10_recef[3]; //R- vector for TDRS10 in Earth Centered Earth Fixed frame or TEME frame
    double prevtdrs10_recef[3]; //R- vector for TDRS10 in Earth Centered Earth Fixed frame or TEME frame
    double tdrs10_vecef[3]; //V- vector for TDRS10 in Earth Centered Earth Fixed frame or TEME frame

    /*
     * razel           Range, azimuth, and elevation matrix
     * razelrates      Range rate, azimuth rate, and elevation rate matrix
     */

    double razel[3]; // Range, Azimuth and Elevation Vector for ISS
    double razelrates[3]; // Range, Azimuth and Elevation Changing rate Vector for ISS

    double tdrs8_razel[3]; // Range, Azimuth and Elevation Vector for tdrs8
    double tdrs8_razelrates[3]; // Range, Azimuth and Elevation Changing rate Vector for tdrs8

    double tdrs9_razel[3]; // Range, Azimuth and Elevation Vector for tdrs9
    double tdrs9_razelrates[3]; // Range, Azimuth and Elevation Changing rate Vector for tdrs9

    double tdrs10_razel[3]; // Range, Azimuth and Elevation Vector for tdrs10
    double tdrs10_razelrates[3]; // Range, Azimuth and Elevation Changing rate Vector for tdrs10
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Variables used in sgp4 model">
    char typerun, typeinput, opsmode;

    gravconsttype whichconst; // wgs72 or wgs84 World Geodatic System constants

    double sec, secC, juleanDate, juleandateCurrent, tsince, startmfe, stopmfe, deltamin;
    double tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2;


    double tdrs8_juleanDate, tdrs9_juleanDate, tdrs10_juleanDate;
    //</editor-fold>


    //<editor-fold defaultstate="collapsed" desc="TIME related variable initialization">

    //time variables from scenario epoch time
    int year, mon, day, hr, min;
    //current time variables
    int yearC, monC, dayC, hrC, minC;

    typedef char str3[4];
    str3 monstr[13];
    strcpy(monstr[1], "Jan");
    strcpy(monstr[2], "Feb");
    strcpy(monstr[3], "Mar");
    strcpy(monstr[4], "Apr");
    strcpy(monstr[5], "May");
    strcpy(monstr[6], "Jun");
    strcpy(monstr[7], "Jul");
    strcpy(monstr[8], "Aug");
    strcpy(monstr[9], "Sep");
    strcpy(monstr[10], "Oct");
    strcpy(monstr[11], "Nov");
    strcpy(monstr[12], "Dec");

    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Structure Declaration for satellites">
    elsetrec satrec; // ISS related values or structure that holds various data related to ISS
    elsetrec tdrs8_satrec; // TDRS8 related values or structure that holds various data related to TDRS8
    elsetrec tdrs9_satrec; // TDRS9 related values or structure that holds various data related to TDRS9
    elsetrec tdrs10_satrec; // TDRS10 related values or structure that holds various data related to TDRS10
    //</editor-fold>

    float elevation;
    float azimuth; //-180 to 0 to 180


    //SET REAL TIME CLOCK (Set values manually using custom excel function until I find a way to do it automatically)    
    //    set_time(1440763200);

    //SET VARIABLES
    opsmode = 'i';
    typerun = 'c';
    typeinput = 'e';
    whichconst = wgs72;

    // initialize wgs constant based on the which constant value, we are sending reference variables to initialize
    getgravconst(whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2);



    //INITIALIZE SATELLITE TRACKING    

    //printf("Initializing ISS orbit from TLE...\n");
    twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst, startmfe, stopmfe, deltamin, satrec);
    //printf("twoline2rv function complete for ISS...\n");

    //printf("Initializing TDRS8 orbit from it's TLE...\n");
    twoline2rv(tdrs8_longstr1, tdrs8_longstr2, typerun, typeinput, opsmode, whichconst, startmfe, stopmfe, deltamin, tdrs8_satrec);
    //printf("twoline2rv for TDRS8 function complete...\n");

    //printf("Initializing TDRS9 orbit from it's TLE...\n");
    twoline2rv(tdrs9_longstr1, tdrs9_longstr2, typerun, typeinput, opsmode, whichconst, startmfe, stopmfe, deltamin, tdrs9_satrec);
    //printf("twoline2rv for TDRS9 function complete...\n");


    //printf("Initializing TDRS10 orbit from it's TLE...\n");
    twoline2rv(tdrs10_longstr1, tdrs10_longstr2, typerun, typeinput, opsmode, whichconst, startmfe, stopmfe, deltamin, tdrs10_satrec);
    //printf("twoline2rv for TDRS10 function complete...\n");


    //Call propogator to get initial state vector value
    sgp4(whichconst, satrec, 0.0, ro, vo);
    //printf("SGP4 at t = 0 to get initial state vector complete...\n"); 

    //Call propogator to get initial state vector value
    sgp4(whichconst, tdrs8_satrec, 0.0, tdrs8_ro, tdrs8_vo);
    //printf("SGP4 at t = 0 to get initial state vector complete...\n"); 

    //Call propogator to get initial state vector value
    sgp4(whichconst, tdrs9_satrec, 0.0, tdrs9_ro, tdrs9_vo);
    //printf("SGP4 at t = 0 to get initial state vector complete...\n"); 

    //Call propogator to get initial state vector value
    sgp4(whichconst, tdrs10_satrec, 0.0, tdrs10_ro, tdrs10_vo);
    //printf("SGP4 at t = 0 to get initial state vector complete...\n"); 


    juleanDate = satrec.jdsatepoch;
    tdrs8_juleanDate = tdrs8_satrec.jdsatepoch;
    tdrs9_juleanDate = tdrs9_satrec.jdsatepoch;
    tdrs10_juleanDate = tdrs10_satrec.jdsatepoch;
    printf("Julean date %f\n", juleanDate);

    invjday(juleanDate, year, mon, day, hr, min, sec);

    //The epoch defines the time to which all of the time-varying fields in the element set are referenced. 

    printf("Scenario Epoch   %3i %3s%5i%3i:%2i:%12.9f \n", day, monstr[mon], year, hr, min, sec);

    juleandateCurrent = getJulianFromUnix(time(NULL));

    printf("CurrentJulean date %f\n", juleandateCurrent);

    invjday(juleandateCurrent, yearC, monC, dayC, hrC, minC, secC);

    printf("Current Time    %3i %3s%5i%3i:%2i:%12.9f \n", dayC, monstr[monC], yearC, hrC, minC, secC);
    printf("            Time            PosX            PosY            PosZ              Vx              Vy              Vz\n");
    printf("            Time             Lat            Long          Height           Range         Azimuth       Elevation\n");

    //BEGIN SATELLITE TRACKING

    ofstream myfile;
    myfile.open("ISS.txt");
    ofstream geoLoc;
    geoLoc.open("GroundStation.txt");
    ofstream tdrs8;
    tdrs8.open("tdrs8.txt");
    ofstream tdrs9;
    tdrs9.open("tdrs9.txt");
    ofstream tdrs10;
    tdrs10.open("tdrs10.txt");

    site(siteLatRad, siteLonRad, siteAlt, ground_ro, ground_vo);
    /*
     *  // test to see what ijk2ll does
//     */
    //    double test[3];
    //    ijk2ll(ground_ro, test);
    //    printf("ijk2ll-\n");
    //    printf("x- %f, y-%f, z-%f \n", test[0]*180.0/pi, test[1]*180.0/pi, test[2]); //print lat,lon,alt 
    //test ends
    int i = 0;
    std::time_t end = std::time(NULL) + 60.0;
    while (std::time(NULL) <= end) {


        //RUN SGP4 AND COORDINATE TRANSFORMATION COMPUTATIONS
        juleandateCurrent = getJulianFromUnix(time(NULL));

        tsince = (juleandateCurrent - juleanDate) * 24.0 * 60.0;
        sgp4(whichconst, satrec, tsince, ro, vo);
        teme2ecef(ro, vo, juleandateCurrent, recef, vecef);
        ijk2ll(recef, latlongh);
        rv2azel(ro, vo, siteLatRad, siteLonRad, siteAlt, juleandateCurrent, razel, razelrates);


        site(siteLatRad, siteLonRad, siteAlt, ground_ro, ground_vo);
        printf("geo\n");
        printf("x- %f, y-%f, z-%f \n", ground_ro[0], ground_ro[1], ground_ro[2]);

        //CHECK FOR ERRORS
        if (satrec.error > 0) {
            printf("# *** error: t:= %f *** code = %3d\n", satrec.t, satrec.error);
        } else {
            azimuth = razel[1]*180 / pi; // Azimuth * (180/pi)

            elevation = razel[2]*180 / pi; // Elevation  * (180/pi)

            printf("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", satrec.t, recef[0], recef[1], recef[2], vecef[0], vecef[1], vecef[2]);
            printf("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", satrec.t, latlongh[0]*180 / pi, latlongh[1]*180 / pi, latlongh[2], razel[0], razel[1]*180 / pi, razel[2]*180 / pi);

            if (i == 0) {


                for (int j = 0; j < 3; j++) {
                    prevrecef[j] = recef[j];
                    prevground_ro[j] = ground_ro[j];
                }
                myfile << i << " " << recef[0] << " " << recef[1] << " "
                        << recef[2] << " " << vecef[0] << " " << vecef[1] << " "
                        << vecef[2] << " " << latlongh[0] * 180 / pi << " "
                        << latlongh[1] * 180 / pi << " " << latlongh[2] << " "
                        << razel[0] << " " << razel[1] * 180 / pi << " "
                        << razel[2] * 180 / pi << "\n";


                geoLoc << i << " " << ground_ro[0] << " " << ground_ro[1] << " " << ground_ro[2] << "\n";
            } else {
                if (isChanged(recef, prevrecef)) {
                    myfile << i << " " << recef[0] << " " << recef[1] << " "
                            << recef[2] << " " << vecef[0] << " " << vecef[1] << " "
                            << vecef[2] << " " << latlongh[0] * 180 / pi << " "
                            << latlongh[1] * 180 / pi << " " << latlongh[2] << " "
                            << razel[0] << " " << razel[1] * 180 / pi << " "
                            << razel[2] * 180 / pi << "\n";
                    for (int j = 0; j < 3; j++) {
                        prevrecef[j] = recef[j];
                        prevground_ro[j] = ground_ro[j];
                    }
                } else {
                    printf("ISS loc not changed \n");
                }

                if (isChanged(ground_ro, prevground_ro)) {
                    geoLoc << i << " " << ground_ro[0] << " " << ground_ro[1] << " " << ground_ro[2] << "\n";
                }
            }

        }

        //tdrs-8 SGP4 AND COORDINATE TRANSFORMATION COMPUTATIONS

        tsince = (juleandateCurrent - tdrs8_juleanDate) * 24.0 * 60.0;
        sgp4(whichconst, tdrs8_satrec, tsince, tdrs8_ro, tdrs8_vo);
        teme2ecef(tdrs8_ro, tdrs8_vo, juleandateCurrent, tdrs8_recef, tdrs8_vecef);
        ijk2ll(tdrs8_recef, tdrs8_latlongh);
        rv2azel(tdrs8_ro, tdrs8_vo, siteLatRad, siteLonRad, siteAlt, juleandateCurrent, tdrs8_razel, tdrs8_razelrates);
        //CHECK FOR ERRORS
        if (satrec.error > 0) {
            printf("# *** error: t:= %f *** code = %3d\n", tdrs8_satrec.t, tdrs8_satrec.error);
        } else {

            printf("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", tdrs8_satrec.t, tdrs8_recef[0], tdrs8_recef[1], tdrs8_recef[2], tdrs8_vecef[0], tdrs8_vecef[1], tdrs8_vecef[2]);
            printf("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", tdrs8_satrec.t, tdrs8_latlongh[0]*180 / pi, tdrs8_latlongh[1]*180 / pi, tdrs8_latlongh[2], tdrs8_razel[0], tdrs8_razel[1]*180 / pi, tdrs8_razel[2]*180 / pi);

            if (i == 0) {
                for (int j = 0; j < 3; j++) {
                    prevtdrs8_recef[j] = tdrs8_recef[j];
                }
                tdrs8 << i << " " << tdrs8_recef[0] << " " << tdrs8_recef[1] << " "
                        << tdrs8_recef[2] << " " << tdrs8_vecef[0] << " " << tdrs8_vecef[1] << " "
                        << tdrs8_vecef[2] << " " << tdrs8_latlongh[0] * 180 / pi << " "
                        << tdrs8_latlongh[1] * 180 / pi << " " << tdrs8_latlongh[2] << " "
                        << tdrs8_razel[0] << " " << tdrs8_razel[1] * 180 / pi << " "
                        << tdrs8_razel[2] * 180 / pi << "\n";


            } else {
                if (isChanged(tdrs8_recef, prevtdrs8_recef)) {
                    tdrs8 << i << " " << tdrs8_recef[0] << " " << tdrs8_recef[1] << " "
                            << tdrs8_recef[2] << " " << tdrs8_vecef[0] << " " << tdrs8_vecef[1] << " "
                            << tdrs8_vecef[2] << " " << tdrs8_latlongh[0] * 180 / pi << " "
                            << tdrs8_latlongh[1] * 180 / pi << " " << tdrs8_latlongh[2] << " "
                            << tdrs8_razel[0] << " " << tdrs8_razel[1] * 180 / pi << " "
                            << tdrs8_razel[2] * 180 / pi << "\n";
                    for (int j = 0; j < 3; j++) {
                        prevtdrs8_recef[j] = tdrs8_recef[j];
                    }
                } else {
                    printf("tdrs8 loc not changed \n");
                }

            }


        }

        //tdrs-9 SGP4 AND COORDINATE TRANSFORMATION COMPUTATIONS

        tsince = (juleandateCurrent - tdrs9_juleanDate) * 24.0 * 60.0;
        sgp4(whichconst, tdrs9_satrec, tsince, tdrs9_ro, tdrs9_vo);
        teme2ecef(tdrs9_ro, tdrs9_vo, juleandateCurrent, tdrs9_recef, tdrs9_vecef);
        ijk2ll(tdrs9_recef, tdrs9_latlongh);
        rv2azel(tdrs9_ro, tdrs9_vo, siteLatRad, siteLonRad, siteAlt, juleandateCurrent, tdrs9_razel, tdrs9_razelrates);
        //CHECK FOR ERRORS
        if (satrec.error > 0) {
            printf("# *** error: t:= %f *** code = %3d\n", tdrs9_satrec.t, tdrs9_satrec.error);
        } else {

            printf("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", tdrs9_satrec.t, tdrs9_recef[0], tdrs9_recef[1], tdrs9_recef[2], tdrs9_vecef[0], tdrs9_vecef[1], tdrs9_vecef[2]);
            printf("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", tdrs9_satrec.t, tdrs9_latlongh[0]*180 / pi, tdrs9_latlongh[1]*180 / pi, tdrs9_latlongh[2], tdrs9_razel[0], tdrs9_razel[1]*180 / pi, tdrs9_razel[2]*180 / pi);

            if (i == 0) {

                for (int j = 0; j < 3; j++) {
                    prevtdrs9_recef[j] = tdrs9_recef[j];
                }
                tdrs9 << i << " " << tdrs9_recef[0] << " " << tdrs9_recef[1] << " "
                        << tdrs9_recef[2] << " " << tdrs9_vecef[0] << " " << tdrs9_vecef[1] << " "
                        << tdrs9_vecef[2] << " " << tdrs9_latlongh[0] * 180 / pi << " "
                        << tdrs9_latlongh[1] * 180 / pi << " " << tdrs9_latlongh[2] << " "
                        << tdrs9_razel[0] << " " << tdrs9_razel[1] * 180 / pi << " "
                        << tdrs9_razel[2] * 180 / pi << "\n";


            } else {
                if (isChanged(tdrs9_recef, prevtdrs9_recef)) {
                    tdrs9 << i << " " << tdrs9_recef[0] << " " << tdrs9_recef[1] << " "
                            << tdrs9_recef[2] << " " << tdrs9_vecef[0] << " " << tdrs9_vecef[1] << " "
                            << tdrs9_vecef[2] << " " << tdrs9_latlongh[0] * 180 / pi << " "
                            << tdrs9_latlongh[1] * 180 / pi << " " << tdrs9_latlongh[2] << " "
                            << tdrs9_razel[0] << " " << tdrs9_razel[1] * 180 / pi << " "
                            << tdrs9_razel[2] * 180 / pi << "\n";


                    for (int j = 0; j < 3; j++) {
                        prevtdrs9_recef[j] = tdrs9_recef[j];
                    }
                } else {
                    printf("tdrs8 loc not changed \n");
                }

            }



        }

        //tdrs-10 SGP4 AND COORDINATE TRANSFORMATION COMPUTATIONS

        tsince = (juleandateCurrent - tdrs10_juleanDate) * 24.0 * 60.0;
        sgp4(whichconst, tdrs10_satrec, tsince, tdrs10_ro, tdrs10_vo);
        teme2ecef(tdrs10_ro, tdrs10_vo, juleandateCurrent, tdrs10_recef, tdrs10_vecef);
        ijk2ll(tdrs10_recef, tdrs10_latlongh);
        rv2azel(tdrs10_ro, tdrs10_vo, siteLatRad, siteLonRad, siteAlt, juleandateCurrent, tdrs10_razel, tdrs10_razelrates);
        //CHECK FOR ERRORS
        if (satrec.error > 0) {
            printf("# *** error: t:= %f *** code = %3d\n", tdrs10_satrec.t, tdrs10_satrec.error);
        } else {

            printf("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", tdrs10_satrec.t, tdrs10_recef[0], tdrs10_recef[1], tdrs10_recef[2], tdrs10_vecef[0], tdrs10_vecef[1], tdrs10_vecef[2]);
            printf("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", tdrs10_satrec.t, tdrs10_latlongh[0]*180 / pi, tdrs10_latlongh[1]*180 / pi, tdrs10_latlongh[2], tdrs10_razel[0], tdrs10_razel[1]*180 / pi, tdrs10_razel[2]*180 / pi);

            if (i == 0) {
                for (int j = 0; j < 3; j++) {
                    prevtdrs10_recef[j] = tdrs10_recef[j];
                }

                tdrs10 << i << " " << tdrs10_recef[0] << " " << tdrs10_recef[1] << " "
                        << tdrs10_recef[2] << " " << tdrs10_vecef[0] << " " << tdrs10_vecef[1] << " "
                        << tdrs10_vecef[2] << " " << tdrs10_latlongh[0] * 180 / pi << " "
                        << tdrs10_latlongh[1] * 180 / pi << " " << tdrs10_latlongh[2] << " "
                        << tdrs10_razel[0] << " " << tdrs10_razel[1] * 180 / pi << " "
                        << tdrs10_razel[2] * 180 / pi << "\n";


            } else {
                if (isChanged(tdrs10_recef, prevtdrs10_recef)) {
                    tdrs10 << i << " " << tdrs10_recef[0] << " " << tdrs10_recef[1] << " "
                            << tdrs10_recef[2] << " " << tdrs10_vecef[0] << " " << tdrs10_vecef[1] << " "
                            << tdrs10_vecef[2] << " " << tdrs10_latlongh[0] * 180 / pi << " "
                            << tdrs10_latlongh[1] * 180 / pi << " " << tdrs10_latlongh[2] << " "
                            << tdrs10_razel[0] << " " << tdrs10_razel[1] * 180 / pi << " "
                            << tdrs10_razel[2] * 180 / pi << "\n";
                    for (int j = 0; j < 3; j++) {
                        prevtdrs10_recef[j] = tdrs10_recef[j];
                    }
                } else {
                    printf("tdrs8 loc not changed \n");
                }

            }

            

        }
        i++;
        sleep(1);
    } //indefinite loop
    myfile.close();
    geoLoc.close();
    tdrs8.close();
    tdrs9.close();
    tdrs10.close();
    return 0;
}

