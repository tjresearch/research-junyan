#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <io.h>
#include <fstream> 
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;
double juliandate = 2458724.5; //for august 29
double timepast;
//initial speed for saturn v 38,946 km/h
class Spacecraft
{
    public:
     double mass;
     double ax, ay, az;
     double vx, vy, vz;
     double px, py, pz;
    
     Spacecraft() { }
    Spacecraft(double vi, double vj, double vk, double pi, double pj, double pk){
        vx = vi;
        vy = vj;
        vz = vk;
        px = pi;
        py = pj;
        pz = pk;
    
    }
     double getX(){
            
             return px;
    }
    double getY(){
            
             return py;
    }
    double getZ(){
             return pz;
    }
    

     void setX(double x){
             px = x;
    }
    void setY(double y){
            
              py = y;
    }
    void setZ(double z){
           
             pz = z;
    }
    
    double getVx(){
            return vx;
    }
    double getVy(){
            return vy;
    }
    double getVz(){
            return vz;
    }
    
    void setVx(double x){
             vx = x;
    }
    void setVy(double y){
            
              vy = y;
    }
    void setVz(double z){
           
             vz = z;
    }
    };
       double vk2(double a, double r, double v, double t){
           return a * (r + v * t / 2);
       }
       double vk3(double a, double r, double v, double t){
            return a * (r + v * a * r * pow(t, 2) / 4);
       }
       double vk4(double a, double r, double v, double t){
           return a * (r + v * a * (r + v * t / 2) * pow(t, 2) / 2);
       }
       double vf(double a, double r, double v, double t){
           return v + t / 6 * (a * r + 2 * vk2(a,r,v,t) + 2 * vk3(a,r,v,t) + vk4(a,r,v,t));
       }
       double rf(double a, double r, double v, double t){
           return r + t / 6 * (v + v * a * r * t + v * vk2(a,r,v,t) * t + v * vk3(a,r,v,t) * t);
       }
    


class Planet
{
    public:
     double mass;
     double orbitalr;
     double range, rascend, declin, ta, E, ma; //changing variables
     double ascendn, peri, inclin, axis, ec, motion, ima; //constants
     double px, py, pz;
     double vx, vy, vz;
     double d2r = M_PI / 180;
     Planet() { }
    Planet(double m, double o, double i, double r, double ra, double dec, double an, double p, double period, double e, double a){
        mass = m; 
        orbitalr = o;
        inclin = i * d2r;
        range = r;
        rascend = ra * d2r;
        declin = dec * d2r;
        ascendn = an * d2r;
        peri = p * d2r;
        motion = 2 * M_PI /  period ;
        ec = e;
        axis = a;
        double n = acos((1 - r / a ) / e);
        ma = acos((1 - r / a ) / e) - e * sin(acos((1 - r / a ) / e));
        ima = acos((1 - r / a ) / e) - e * sin(acos((1 - r / a ) / e));
        /*
        if(i < 1){
            peri = peri + ascendn;
        }
        */
    }
    double getMeanMotion(){
        return motion;
    }
    double getMass(){
        return mass;
    }
    double getOr(){
        return orbitalr;
    }
    double getMeanAnomaly(){
        return ma;
    }
    double getIMA(){
        return ima;
    }
    double getX(){
             px = range * (cos(ascendn) * cos(ta + peri) - sin(ascendn) * sin(ta + peri) * cos(inclin));
             return px;
    }
    double getY(){
             py = range * (sin(ascendn) * cos(ta + peri) + cos(ascendn) * sin(ta + peri) * cos(inclin));
             return py;
    }
    double getZ(){
             pz = range * (sin(ta + peri) * sin(inclin));
             return pz;
    }
    
    double getVx(){
            return vx;
    }
    double getVy(){
            return vy;
    }
    double getVz(){
            return vz;
    }
    

    void setMA(double d){
        ma = d;
    }
    
    


    double eccentricAnomaly(){
        double pastE = ma;
        double curE = ma;
        do{
            pastE = curE;
            curE = pastE - (pastE - ec * sin(pastE) - ma)/(1 - ec * cos(pastE));
            } while((int)(curE * pow(10, 10)) != (int)(pastE * ( pow(10, 10))));
        E = curE;
        return E;
        }

    double trueAnomaly(){
        double v = 2* atan(sqrt((1 + ec)/(1-ec)) * tan(E/2) );
        if(v < 0)
            v = 2 * M_PI + v;
        ta = v;
        return ta;
                    }

    double ellipticalR(){
        range = abs(axis * (1 - pow(ec, 2))/ (1 + ec * cos(ta)));
        return range;
                    }
    
    double rightAscension(){
    // r = range
    // inclin = inclination
    // v = true Anomaly
    // ascend = Longitude of the ascending node
    // peri = Argument of perigee
        double x = range * (cos(ascendn) * cos(ta + peri) - sin(ascendn) * sin(ta + peri) * cos(inclin));
        double y = range * (sin(ascendn) * cos(ta + peri) + cos(ascendn) * sin(ta + peri) * cos(inclin));
        double ra = atan(y / x);
        if(x < 0){
            ra = M_PI + ra;
                    }
    
        rascend =  fmod(ra, 2 * M_PI);  
        return rascend * 180 / M_PI;
            }

    double Declination(){
        double x = range * (cos(ascendn) * cos(ta + peri) - sin(ascendn) * sin(ta + peri) * cos(inclin));
        double y = range * (sin(ascendn) * cos(ta + peri) + cos(ascendn) * sin(ta + peri) * cos(inclin));
        double z = range * (sin(ta + peri) * sin(inclin));
        double dec = atan(z / sqrt(pow(x,2) + pow(y,2)) );
    
        declin = fmod(dec, 2 * M_PI);  
        return declin * 180 / M_PI;
    }

};
    
void solarsystem(Spacecraft &s, double &a1, double &a2, int &p)
{
    //int year, mon, day, hr, minute;
    double sec;
    
    //planetsAU[3] = 1.66602003028769;
     Planet planets [10];
     planets[0] = Planet(3.302e23, 1, 7.00487, .32971028480559, 130.2016, 6.9345, 48.33167, 29.12703035, 88, 0.20563069,0.38709893); //Mercury
     planets[1] = Planet(4.87e24, 1, 7.0, .32971028480559, 130.2016, 6.9345, 48.33167, 29.12703035, 88, .205,0.38709893); //Venus
     planets[2] = Planet(5.972e24, 1, 0.00005, 1.01014351904246, 335.1159, 0.0009, -11.26064, 114.20783, 365.2, 0.01671022, 1); //Earth
     planets[3] = Planet(0.642e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231); //Mars
      planets[4] = Planet(1898e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231); //Jupiter
       planets[5] = Planet(568e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231); //Saturn
      //planets[3] = Planet(86.8e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231);//Uranus
      //planets[3] = Planet(102e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231); //Neptune
    planets[7] = Planet(0.0146e24, 1, 17.14175, 33.8684306892170, 130.2016, 6.9345, 110.30347, 113.76329, 90560, 0.24880766,39.48168677); //Pluto
    // planets[4] = Planet(3.302e23, 1, 7.0, .32971028480559, 130.2016, 6.9345, 48.33167, 29.12703035, 88, .205,0.38709893); //Jupiter
        double day = 0;
        s = Spacecraft(planets[p].getVx() + 1 * sin(a1) * cos(a2), planets[p].getVy() + 1 * sin(a1) * sin(a2), planets[p].getVz() + 1 * cos(a1), 
        planets[p].getX() + planets[p].getOr() * sin(a1) * cos(a2), planets[p].getY() + planets[p].getOr() * sin(a1) * sin(a2), planets[p].getZ() + 1 * cos(a1));
       // cout << acos((1 - .46669 / .387) / .205) - .205 * sin(acos((1 - .46669 / .387 ) / .205));
    for(double time = juliandate; time <= 88 + juliandate; time += .00069444){
        day += .00069444;
        double spr;
        s.setX(s.getVx() * .00069444 + s.getX());
        s.setY(s.getVy() * .00069444 + s.getY());
        s.setZ(s.getVz() * .00069444 + s.getZ());
        for(int n = 0; n < 1; n++){
            double px = planets[n].getX();
                 double py = planets[n].getY();
                double pz = planets[n].getZ();
            planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));
            planets[n].eccentricAnomaly();
            planets[n].trueAnomaly();
            double range = planets[n].ellipticalR();
            double rascend = planets[n].rightAscension();
            double declin = planets[n].Declination();
            double vx = planets[n].getX() - px;
            double vy = planets[n].getY() - py;
            double vz = planets[n].getZ() - pz;
                px = planets[n].getX();
                 py = planets[n].getY();
                pz = planets[n].getZ();
            cout.precision(10);
           if(day >= 1){
                
                cout << "time: " << time << " ";
                cout << "range: " << range << " ";
                cout << "right Ascension: " << rascend << " ";
                cout << "declination: " << declin << "\n";
                
                day = 0;
                           }    
        }
        
    }
 
}
int main(){
    Spacecraft icarus;
    int i = 0;
    for(double theta = 0; theta < 2 * M_PI; theta+= M_PI / 180){
        for(double phi = 0; phi < M_PI; phi+= M_PI / 180){

            solarsystem(icarus, theta, phi, i);
        }
    }
    return 0;
}