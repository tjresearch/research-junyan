#define _USE_MATH_DEFINES
//#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
//#include <io.h>
#include <fstream> 
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
//initial velocity of spacecraft 0.00645468 AU / day
//earth's orbital radius 6556145.59 meters
//orbital radius is where g field is 9.2672 m/s^2
// 1 AU = 1.496e+11 m
//1 day = 86400 s
using namespace std;
double juliandate = 2458724.5; //for august 29
double timepast;
//initial speed for saturn v 38,946 km/h
class Spacecraft
{
    public:
     double mass;
     double row, theta, phi;
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
        row = sqrt(pow(pi, 2) + pow(pj, 2) + pow(pk, 2));
        theta = atan(pj / pi);
        phi = atan(sqrt(pow(pi, 2) + pow(pj, 2)) / pk );
    }
    double getRow(){
        return row;
    }
    double getTheta(){
        return theta;
    }
    double getPhi(){
        return phi;
    }
    void setSphere(double pi, double pj, double pk){
        row = sqrt(pow(pi, 2) + pow(pj, 2) + pow(pk, 2));
        theta = atan(pj / pi);
        phi = atan(sqrt(pow(pi, 2) + pow(pj, 2)) / pk );
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
    Planet(double m, double o, double i, double r, double ra, double dec, double an, double p, double period, double e, double a, double x, double y, double z, double vxc, double vyc, double vzc){
        mass = m; 
        orbitalr = o;
        inclin = i * d2r;
        range = r;
        //orbitalr = sqrt(6.67 * 10e-11 * m / 9.2672) / 1.496e+11;
        //if(orbitalr < r)
        //orbitalr = 1e-3;
        rascend = ra * d2r;
        declin = dec * d2r;
        ascendn = an * d2r;
        peri = p * d2r;
        motion = 2 * M_PI /  period ;
        ec = e;
        axis = a;
        double n = acos((1 - r / a ) / e);
        px = x;
        py = y;
        pz = z;
        vx = vxc;
        vy = vyc;
        vz = vzc;
       // ma = 0;
       // ma = acos((1 - r / a ) / e) - e * sin(acos((1 - r / a ) / e));
       // ima = acos((1 - r / a ) / e) - e * sin(acos((1 - r / a ) / e));
         ma = atan(((px * vx + py * vy + pz * vz)/ (a * a * motion)) / (1 - sqrt(px * px + py * py + pz * pz) / a)) 
         - ec * sin(atan(((px * vx + py * vy + pz * vz)/ (a * a * motion)) / (1 - sqrt(px * px + py * py + pz * pz) / a)));
       ima = ma;
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
             px = range * sin(M_PI/2 - declin) * cos(rascend);
             return px;
    }
    double getY(){
             py = range * sin(M_PI/2 - declin) * sin(rascend);
             return py;
    }
    double getZ(){
             pz = range * cos(M_PI/2 - declin);
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
    double dist(Spacecraft s, Planet p){
        return sqrt(pow(s.getX() - p.getX(), 2) + pow(s.getY() - p.getY(), 2) + pow(s.getZ() - p.getZ(), 2));
    }
    
int solarsystem(Spacecraft &s, double &a1, double &a2, int &p)
{
    //int year, mon, day, hr, minute;
    double sec;
    ofstream out;
  out.open ("spacecraft.txt");
  
     Planet planets [10];
     planets[0] = Planet(3.302e23, 1e-3, 7.00487, .32971028480559, 130.2016, 6.9345, 48.33167, 29.12703035, 88, 0.20563069,0.38709893,
      -2.112940522890843E-01, 2.499646441178046E-01, 3.980885339059970E-02, -2.717908606532378E-02, -1.704717040208352E-02, 1.100301789608752E-03); //Mercury
    // planets[1] = Planet(4.87e24, 1, 7.0, .32971028480559, 130.2016, 6.9345, 48.33167, 29.12703035, 88, .205,0.38709893); //Venus
   //  planets[2] = Planet(5.972e24, 1, 0.00005, 1.01014351904246, 335.1159, 0.0009, -11.26064, 114.20783, 365.2, 0.01671022, 1); //Earth
   //  planets[3] = Planet(0.642e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231); //Mars
   //   planets[4] = Planet(1898e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231); //Jupiter
    //   planets[5] = Planet(568e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231); //Saturn
    //  planets[6] = Planet(86.8e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231);//Uranus
    //  planets[7] = Planet(102e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231); //Neptune
// planets[8] = Planet(0.0146e24, 1, 17.14175, 33.8684306892170, 130.2016, 6.9345, 110.30347, 113.76329, 90560, 0.24880766,39.48168677); //Pluto
    // planets[4] = Planet(3.302e23, 1, 7.0, .32971028480559, 130.2016, 6.9345, 48.33167, 29.12703035, 88, .205,0.38709893); //Jupiter
        double day = 0;
        s = Spacecraft( 0.00645468 * sin(a2) * cos(a1),  0.00645468 * sin(a2) * sin(a1),  0.00645468 * cos(a2), 
        planets[p].getX() + planets[p].getOr() * sin(a2) * cos(a1), planets[p].getY() + planets[p].getOr() * sin(a2) * sin(a1), planets[p].getZ() + planets[p].getOr() * cos(a2));
       // cout << acos((1 - .46669 / .387) / .205) - .205 * sin(acos((1 - .46669 / .387 ) / .205));
    //for(double time = juliandate; time <= 1000 + juliandate; time += .00069444){
        double time = juliandate;
        while(sqrt(pow(s.getX(), 2) + pow(s.getY(), 2) + pow(s.getZ(), 2)) <= 1){
        day += .00069444;
        double spr;
        /*
        s.setX(rf(1, s.getX(), s.getVx(), .00069444));
       tZ(rf(1, s.getZ(), s.getVz(), .00069444));
        */
        //s.setSp s.setY(rf(1, s.getY(), s.getVy(), .00069444));
        s.setSphere(s.getX(), s.getY(), s.getZ());
        if(sqrt(pow(s.getX(), 2) + pow(s.getY(), 2) + pow(s.getZ(), 2)) < 0.04303120071) //closest we've gotten to the sun
            return 100;
       // cout << a1 << " " << a2 << ": " << s.getRow() << " " << s.getTheta() << " " << s.getPhi() << "\n";
        for(int n = 0; n < 1; n++){
            if(dist(s, planets[n]) < planets[n].getOr() && time != juliandate)
                return n;
            double px = planets[n].getX();
           // if(day >= 1)
            //out << px << "\n";
                 double py = planets[n].getY();
            if(day >= 1)
             out << py << "\n";
                double pz = planets[n].getZ();
            planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));
          // planets[n].setMA(fmod( planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));
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
            if(time == juliandate){
                 s = Spacecraft(s.getVx() + 0.027352689 * sin(a2) * cos(a1),s.getVy() + 0.027352689 * sin(a2) * sin(a1), s.getVz() + 0.027352689 * cos(a2), 
        planets[p].getX() + planets[p].getOr() * sin(a2) * cos(a1), planets[p].getY() + planets[p].getOr() * sin(a2) * sin(a1), planets[p].getZ() + planets[p].getOr() * cos(a2));
            }  
           if(day >= 1){
                
                cout << "time: " << time << " ";
                cout << "range: " << range << " ";
                cout << "right Ascension: " << rascend << " ";
                cout << "declination: " << declin << "\n";
                
               // cout << a1 << " " << a2 << ": " << s.getRow() << " " << s.getTheta() << " " << s.getPhi() << "\n";
              
                //out << s.getX() << "\n";
                // out << s.getY() <<  "\n";
                day = 0;
                           }
        }
       
        s.setX(s.getX() + s.getVx() * .00069444);
        s.setY(s.getY() + s.getVy() * .00069444);
        s.setZ(s.getZ() + s.getVz() * .00069444);
            s.setSphere(s.getX(), s.getY(), s.getZ());
            
        double gravsun = 6.67e-11 * 1.989e30 / pow(s.getRow() * 1.496e+11, 2) ;
        double gravmer = 6.67e-11 * 3.302e23 / pow(dist(s, planets[0]) * 1.496e+11, 2) ;
        double phi = acos((s.getZ() - planets[0].getZ()) / dist(s, planets[0]) );
        double theta = asin((s.getY() - planets[0].getY()) / dist(s, planets[0]));
       s.setVx(vf(gravsun / 1.496e+11 * sin(M_PI - s.getPhi()) * cos(M_PI + s.getTheta()) * 7464960000, s.getX(), s.getVx(), .00069444));
        s.setVy(vf(gravsun / 1.496e+11 * sin(M_PI - s.getPhi()) * sin(M_PI + s.getTheta()) * 7464960000, s.getY(), s.getVy(), .00069444));
        s.setVz(vf(gravsun / 1.496e+11 * cos(M_PI - s.getPhi()) * 7464960000, s.getZ(), s.getVz(), .00069444));
         s.setVx(vf(gravmer / 1.496e+11 * sin( M_PI - phi) * cos(theta + M_PI) * 7464960000, s.getX(), s.getVx(), .00069444));
        s.setVy(vf(gravmer / 1.496e+11 * sin( M_PI - phi) * sin(theta + M_PI) * 7464960000, s.getY(), s.getVy(), .00069444));
        s.setVz(vf(gravmer / 1.496e+11 * cos(M_PI - phi) * 7464960000, s.getZ(), s.getVz(), .00069444));
        time += .00069444;
    }
    out.close();
    return -1;
}
int main(){
  
//int main(int argc, char* argv[]) {}
    //#pragma omp parallel{
    //MPI_Init(NULL, NULL);
    
    // int world_size;
    //MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //int world_rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    /*
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    */
    Spacecraft icarus;
    int i = 0;
    //theta 3.70178 phi 1.53649
    //theta 0 phi 0.05235987756
    //theta  2.19911  phi 1.2217
    // theta 3.787359619
    // 4.136425469 theta
    // theta 5.078903265
    // theta 6.126100816
    //for()
    for(double theta = 3.787359619; theta < 3.787359619 + 1e-10; theta+= M_PI / 180){
        for(double phi = M_PI /2 ; phi < M_PI / 2 + 1e-10; phi+= M_PI / 180){

            //if(solarsystem(icarus, theta, phi, i) != -1)
                cout << theta << " " << phi << ": " << solarsystem(icarus, theta, phi, i) << "\n";
            if(solarsystem(icarus, theta, phi, i) != -1 && solarsystem(icarus, theta, phi, i) != 100){
            phi = M_PI + 1;
            theta = 2 * M_PI + 1;
            }
        }
    }
    //}
    //MPI_Finalize();
        return 0;
}