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
//initial speed for saturn v 38,946 km/h or 40233.5
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
         ma = atan(((px * vx + py * vy + pz * vz)/ (a * a * motion)) / (1 - sqrt(px * px + py * py + pz * pz) / a)) 
         - ec * sin(atan(((px * vx + py * vy + pz * vz)/ (a * a * motion)) / (1 - sqrt(px * px + py * py + pz * pz) / a)));
       ima = ma;
   
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
    
int solarsystem(Spacecraft &s, double &a1, double &a2, int &p, double &t, double &temp)
{
    //int year, mon, day, hr, minute;
    double sec;
    ofstream out;
  out.open ("spacecraft.txt");
  
     Planet planets [10];
     planets[0] = Planet(3.302e23, 1.63083872e-5, 7.00487, .32971028480559, 130.2016, 6.9345, 48.33167, 29.12703035, 88, 0.20563069,0.38709893,
      -2.112940522890843E-01, 2.499646441178046E-01, 3.980885339059970E-02, -2.717908606532378E-02, -1.704717040208352E-02, 1.100301789608752E-03); //Mercury
     planets[1] = Planet(4.87e24, 4.04537843e-5, 3.39471, 0.719257101866, 164.8906, 3.3930, 76.68069, 54.85229, 225, 0.00677323,0.72333199, //Venus
      -6.931931527047550E-01, 1.870889764776485E-01, 4.256943212606793E-02, -5.370335538508818E-03, -1.961855230377539E-02, 4.071529477013288E-05);
    planets[2] = Planet(5.972e24, 4.3825089e-5, 0.00005, 1.01014351904246, 335.1159, 0.0009, -11.26064, 114.20783, 365.2, 0.01671022, 1,
     9.163973382911959E-01, -4.249746745750780E-01, 1.601998994026577E-05, 6.962668804528970E-03, 1.554826356216624E-02, -1.320226017688261E-06); //Earth
     planets[3] = Planet(0.642e24, 2.26574081e-5, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231, //Mars
    -1.538330385243809E+00, 6.376127589280390E-01, 5.110497339903643E-02, -4.834909566264149E-03, -1.173146288124897E-02, -1.271951751534037E-04); 
   planets[4] = Planet(1898e24, 7.81287858264e-4, 1.30530, 5.272299325274, 265.5878, 0.3359, 100.55615, -85.8023, 4328.9, 0.04839266, 5.20336301, //Jupiter
     -4.053769515202177E-01, -5.256590709274248E+00, 3.090371506818809E-02, 7.441268269140511E-03, -2.245620586436468E-04, -1.655281797727711E-04);
      planets[5] = Planet(568e24, 4.2740314218e-4, 2.48446, 1.66602003028769, 288.4594, 0.2227, 113.71504, -41.2831, 10751.8, 0.05415060, 9.53707032,
      3.181044968477705E+00, -9.528596201541061E+00, 3.903112876199565E-02, 4.991280595151790E-03, 1.750654182256989E-03, -2.292205742319197E-04); //Saturn
     planets[6] = Planet(86.8e24, 1.6707924974e-4, 0.76986, 19.83381717157, 33.6702, -0.4996, 74.22988, 96.73436, 30667.5, 0.04716771, 19.19126393,
     1.650566200245873E+01, 1.099604230661651E+01, -1.729240416727421E-01, -2.202500959781696E-03, 3.089681894702222E-03, 3.989664038785872E-05);//Uranus
   planets[7] = Planet(102e24, 1.8111888808e-4, 1.76917, 29.93486047576, 346.9602, -1.0202 , 131.72169, -86.75034, 60225, 0.00858587, 30.06896348,
    2.915843957983342E+01, -6.752549485815638E+00, -5.330168462398258E-01, 6.955819371261789E-04,  3.077650123843038E-03, -7.970622808811658E-05); //Neptune
    double time;
    for(time = juliandate; time <= t; time += .00069444){
        for(int n = 0; n < 8; n++){
          
            double px = planets[n].getX();
           // if(day >= 1)
            //out << px << "\n";
                 double py = planets[n].getY();
           // if(day >= 1)
            // out << py << "\n";
                double pz = planets[n].getZ();
          if(n == 0)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));
          else if(n == 1)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));           
          else if(n == 2)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 185), 2 * M_PI));
            else if(n == 3)     
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 344), 2 * M_PI));
           else if(n == 4)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 2219.526), 2 * M_PI));
            else if(n == 5)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 6111.204), 2 * M_PI));
            else if(n == 6)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 15757.844), 2 * M_PI));
            else if(n == 7)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate - 1695.75), 2 * M_PI));           
            
            planets[n].eccentricAnomaly();
            planets[n].trueAnomaly();
            double range = planets[n].ellipticalR();
            double rascend = planets[n].rightAscension();
            double declin = planets[n].Declination();
            planets[n].vx = (planets[n].getX() - px) / .00069444;
             planets[n].vy = (planets[n].getY() - py) / .00069444;
              planets[n].vz = (planets[n].getZ() - pz) / .00069444;
                px = planets[n].getX();
                 py = planets[n].getY();
                pz = planets[n].getZ();
                cout.precision(10);
        }
    }
        double day = 0;
        s = Spacecraft( 0.00645468 * sin(a2) * cos(a1) + planets[p].getVx(),  0.0064468 * sin(a2) * sin(a1) + planets[p].getVy(),  0.00645468 * cos(a2) + planets[p].getVz(), 
        planets[p].getX() + planets[p].getOr() * sin(a2) * cos(a1), planets[p].getY() + planets[p].getOr() * sin(a2) * sin(a1), planets[p].getZ() + planets[p].getOr() * cos(a2));

       double inittime = time;

         double px = planets[0].getX();
               double  py = planets[0].getY();
              double  pz = planets[0].getZ();
            //  cout << px << " " << py << "\n";
          //    cout << s.getX() << " " << s.getY() << "\n";
        while(sqrt(pow(s.getX(), 2) + pow(s.getY(), 2) + pow(s.getZ(), 2)) <= 450)
        {
            //day += .00069444;
            double spr;
            /*
         s.setX(rf(1, s.getX(), s.getVx(), .00069444));
            tZ(rf(1, s.getZ(), s.getVz(), .00069444));
            */
             s.setX(s.getX() + s.getVx() * .00069444);
                s.setY(s.getY() + s.getVy() * .00069444);
                s.setZ(s.getZ() + s.getVz() * .00069444);
                s.setSphere(s.getX(), s.getY(), s.getZ());
                /*
                cout << planets[0].getX() << " " <<   planets[0].getY() << "\n";
                 cout << s.getX() << " " << s.getY() << "\n";
                 */
            //s.setSp s.setY(rf(1, s.getY(), s.getVy(), .00069444));
            //s.setSphere(s.getX(), s.getY(), s.getZ());
            if(sqrt(pow(s.getX(), 2) + pow(s.getY(), 2) + pow(s.getZ(), 2)) < 0.04303120071) //closest we've gotten to the sun
            return 100;
            // cout << a1 << " " << a2 << ": " << s.getRow() << " " << s.getTheta() << " " << s.getPhi() << "\n";
           // double gravsun = 6.67e-11 * 1.989e30 / pow(s.getRow() * 1.496e+11, 2) ;
            for(int n = 0; n < 8; n++){
                if(dist(s, planets[n]) < planets[n].getOr() && time != inittime){
                temp = time - t;
                return n;
                 }
                double px = planets[n].getX();
                    // if(day >= 1)
                //out << px << "\n";
                 double py = planets[n].getY();
                // if(day >= 1)
                    // out << py << "\n";
                double pz = planets[n].getZ();
                if(n == 0)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));
                else if(n == 1)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));           
                else if(n == 2)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 185), 2 * M_PI));
                else if(n == 3)     
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 344), 2 * M_PI));
                else if(n == 4)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 2219.526), 2 * M_PI));
                else if(n == 5)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 6111.204), 2 * M_PI));
                else if(n == 6)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 15757.844), 2 * M_PI));
                else if(n == 7)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate - 1695.75), 2 * M_PI));           
            
                // planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));           
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

                    if(day >= 1){
                cout << "planet: " << n << " ";
                cout << "time: " << time << " ";
                cout << "range: " << range << " ";
                cout << "right Ascension: " << rascend << " ";
                cout << "declination: " << declin << "\n";
                
               // cout << a1 << " " << a2 << ": " << s.getRow() << " " << s.getTheta() << " " << s.getPhi() << "\n";
              
                //out << s.getX() << "\n";
                // out << s.getY() <<  "\n";
                    }
            }
            if(day >= 1)
                day = 0;
                /*
                s.setX(s.getX() + s.getVx() * .00069444);
                s.setY(s.getY() + s.getVy() * .00069444);
                s.setZ(s.getZ() + s.getVz() * .00069444);
                s.setSphere(s.getX(), s.getY(), s.getZ());
            */
                double gravsun = 6.67e-11 * 1.989e30 / pow(s.getRow() * 1.496e+11, 2) ;
                // double gravmer = 6.67e-11 * 3.302e23 / pow(dist(s, planets[0]) * 1.496e+11, 2) ;
                
             double phi = acos((s.getZ() - planets[0].getZ()) / dist(s, planets[0]) );
                 double theta = asin((s.getY() - planets[0].getY()) / dist(s, planets[0]));
                s.setVx(vf(gravsun / 1.496e+11 * sin(M_PI - s.getPhi()) * cos(M_PI + s.getTheta()) * 7464960000, s.getX(), s.getVx(), .00069444));
                 s.setVy(vf(gravsun / 1.496e+11 * sin(M_PI - s.getPhi()) * sin(M_PI + s.getTheta()) * 7464960000, s.getY(), s.getVy(), .00069444));
                 s.setVz(vf(gravsun / 1.496e+11 * cos(M_PI - s.getPhi()) * 7464960000, s.getZ(), s.getVz(), .00069444));
            
            for(int n = 0; n < 8; n++){
                 phi = acos((s.getZ() - planets[n].getZ()) / dist(s, planets[n]) );
                 theta = asin((s.getY() - planets[n].getY()) / dist(s, planets[n]));
                double gravmer = 6.67e-11 * 3.302e23 / pow(dist(s, planets[n]) * 1.496e+11, 2) ;
              s.setVx(vf(gravmer / 1.496e+11 * sin( M_PI - phi) * cos(theta + M_PI) * 7464960000, s.getX(), s.getVx(), .00069444));
                 s.setVy(vf(gravmer / 1.496e+11 * sin( M_PI - phi) * sin(theta + M_PI) * 7464960000, s.getY(), s.getVy(), .00069444));
                s.setVz(vf(gravmer / 1.496e+11 * cos(M_PI - phi) * 7464960000, s.getZ(), s.getVz(), .00069444));
            }
            
                time += .00069444;
                day += .00069444;
        }
    out.close();
    return -1;
}
int solarsystemp(Spacecraft &s, double &a1, double &a2, int &p, double &t)
{
    //int year, mon, day, hr, minute;
    double sec;
    ofstream out;
  out.open ("spacecraft.txt");
    ofstream outp;
    out.open ("planet.txt");
     Planet planets [10];
     planets[0] = Planet(3.302e23, 1e-3, 7.00487, .32971028480559, 130.2016, 6.9345, 48.33167, 29.12703035, 88, 0.20563069,0.38709893,
      -2.112940522890843E-01, 2.499646441178046E-01, 3.980885339059970E-02, -2.717908606532378E-02, -1.704717040208352E-02, 1.100301789608752E-03); //Mercury
     planets[1] = Planet(4.87e24, .001, 3.39471, 0.719257101866, 164.8906, 3.3930, 76.68069, 54.85229, 225, 0.00677323,0.72333199, //Venus
      -6.931931527047550E-01, 1.870889764776485E-01, 4.256943212606793E-02, -5.370335538508818E-03, -1.961855230377539E-02, 4.071529477013288E-05);
    planets[2] = Planet(5.972e24, 1, 0.00005, 1.01014351904246, 335.1159, 0.0009, -11.26064, 114.20783, 365.2, 0.01671022, 1,
     9.163973382911959E-01, -4.249746745750780E-01, 1.601998994026577E-05, 6.962668804528970E-03, 1.554826356216624E-02, -1.320226017688261E-06); //Earth
     planets[3] = Planet(0.642e24, 1, 1.85061, 1.66602003028769, 130.2016, 6.9345, 49.57854, 286.4623, 687.0, 0.09341233, 1.52366231, //Mars
    -1.538330385243809E+00, 6.376127589280390E-01, 5.110497339903643E-02, -4.834909566264149E-03, -1.173146288124897E-02, -1.271951751534037E-04); 
   planets[4] = Planet(1898e24, 1, 1.30530, 5.272299325274, 265.5878, 0.3359, 100.55615, -85.8023, 4328.9, 0.04839266, 5.20336301, //Jupiter
     -4.053769515202177E-01, -5.256590709274248E+00, 3.090371506818809E-02, 7.441268269140511E-03, -2.245620586436468E-04, -1.655281797727711E-04);
      planets[5] = Planet(568e24, 1, 2.48446, 1.66602003028769, 288.4594, 0.2227, 113.71504, -41.2831, 10751.8, 0.05415060, 9.53707032,
      3.181044968477705E+00, -9.528596201541061E+00, 3.903112876199565E-02, 4.991280595151790E-03, 1.750654182256989E-03, -2.292205742319197E-04); //Saturn
     planets[6] = Planet(86.8e24, 1, 0.76986, 19.83381717157, 33.6702, -0.4996, 74.22988, 96.73436, 30667.5, 0.04716771, 19.19126393,
     1.650566200245873E+01, 1.099604230661651E+01, -1.729240416727421E-01, -2.202500959781696E-03, 3.089681894702222E-03, 3.989664038785872E-05);//Uranus
   planets[7] = Planet(102e24, 1, 1.76917, 29.93486047576, 346.9602, -1.0202 , 131.72169, -86.75034, 60225, 0.00858587, 30.06896348,
    2.915843957983342E+01, -6.752549485815638E+00, -5.330168462398258E-01, 6.955819371261789E-04,  3.077650123843038E-03, -7.970622808811658E-05); //Neptune
    double time;
    for(time = juliandate; time <= t; time += .00069444){
        for(int n = 0; n < 8; n++){
          
            double px = planets[n].getX();
           // if(day >= 1)
            //out << px << "\n";
                 double py = planets[n].getY();
           // if(day >= 1)
            // out << py << "\n";
                double pz = planets[n].getZ();
          if(n == 0)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));
          else if(n == 1)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));           
          else if(n == 2)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 185), 2 * M_PI));
            else if(n == 3)     
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 344), 2 * M_PI));
           else if(n == 4)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 2219.526), 2 * M_PI));
            else if(n == 5)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 6111.204), 2 * M_PI));
            else if(n == 6)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 15757.844), 2 * M_PI));
            else if(n == 7)
                planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate - 1695.75), 2 * M_PI));           
            
            planets[n].eccentricAnomaly();
            planets[n].trueAnomaly();
            double range = planets[n].ellipticalR();
            double rascend = planets[n].rightAscension();
            double declin = planets[n].Declination();
            planets[n].vx = (planets[n].getX() - px) / .00069444;
             planets[n].vy = (planets[n].getY() - py) / .00069444;
              planets[n].vz = (planets[n].getZ() - pz) / .00069444;
                px = planets[n].getX();
                 py = planets[n].getY();
                pz = planets[n].getZ();
                cout.precision(10);
        }
    }
        double day = 0;
        s = Spacecraft( 0.00645468 * sin(a2) * cos(a1) + planets[p].getVx(),  0.0064468 * sin(a2) * sin(a1) + planets[p].getVy(),  0.00645468 * cos(a2) + planets[p].getVz(), 
        planets[p].getX() + planets[p].getOr() * sin(a2) * cos(a1), planets[p].getY() + planets[p].getOr() * sin(a2) * sin(a1), planets[p].getZ() + planets[p].getOr() * cos(a2));

       double inittime = time;

         double px = planets[0].getX();
               double  py = planets[0].getY();
              double  pz = planets[0].getZ();
            //  cout << px << " " << py << "\n";
          //    cout << s.getX() << " " << s.getY() << "\n";
        while(sqrt(pow(s.getX(), 2) + pow(s.getY(), 2) + pow(s.getZ(), 2)) <= 450)
        {
            //day += .00069444;
            double spr;
            /*
         s.setX(rf(1, s.getX(), s.getVx(), .00069444));
            tZ(rf(1, s.getZ(), s.getVz(), .00069444));
            */
             s.setX(s.getX() + s.getVx() * .00069444);
                s.setY(s.getY() + s.getVy() * .00069444);
                s.setZ(s.getZ() + s.getVz() * .00069444);
                s.setSphere(s.getX(), s.getY(), s.getZ());
                /*
                cout << planets[0].getX() << " " <<   planets[0].getY() << "\n";
                 cout << s.getX() << " " << s.getY() << "\n";
                 */
            //s.setSp s.setY(rf(1, s.getY(), s.getVy(), .00069444));
            //s.setSphere(s.getX(), s.getY(), s.getZ());
            if(sqrt(pow(s.getX(), 2) + pow(s.getY(), 2) + pow(s.getZ(), 2)) < 0.04303120071) //closest we've gotten to the sun
            return 100;
            // cout << a1 << " " << a2 << ": " << s.getRow() << " " << s.getTheta() << " " << s.getPhi() << "\n";
           // double gravsun = 6.67e-11 * 1.989e30 / pow(s.getRow() * 1.496e+11, 2) ;
            for(int n = 0; n < 8; n++){
                if(dist(s, planets[n]) < planets[n].getOr() && time != inittime){
                    
                return n;
                 }
                double px = planets[n].getX();
                    // if(day >= 1)
                //out << px << "\n";
                 double py = planets[n].getY();
                // if(day >= 1)
                    // out << py << "\n";
                double pz = planets[n].getZ();
                if(n == 0)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));
                else if(n == 1)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));           
                else if(n == 2)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 185), 2 * M_PI));
                else if(n == 3)     
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 344), 2 * M_PI));
                else if(n == 4)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 2219.526), 2 * M_PI));
                else if(n == 5)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 6111.204), 2 * M_PI));
                else if(n == 6)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate + 15757.844), 2 * M_PI));
                else if(n == 7)
                    planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate - 1695.75), 2 * M_PI));           
            
                // planets[n].setMA(fmod(planets[n].getIMA() + planets[n].getMeanMotion() * (time - juliandate), 2 * M_PI));           
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

                    if(n == p && day >= 1){
                        /*
                cout << "planet: " << n << " ";
                cout << "time: " << time << " ";
                cout << "range: " << range << " ";
                cout << "right Ascension: " << rascend << " ";
                cout << "declination: " << declin << "\n";
                */
               // cout << a1 << " " << a2 << ": " << s.getRow() << " " << s.getTheta() << " " << s.getPhi() << "\n";
              
                out << s.getX() << " ";
                 out << s.getY() <<  "\n";
                 outp << planets[n].getX() << " ";
                 outp << planets[n].getY() << "\n";
                    }
            }
            if(day >= 1)
                day = 0;
                /*
                s.setX(s.getX() + s.getVx() * .00069444);
                s.setY(s.getY() + s.getVy() * .00069444);
                s.setZ(s.getZ() + s.getVz() * .00069444);
                s.setSphere(s.getX(), s.getY(), s.getZ());
            */
                double gravsun = 6.67e-11 * 1.989e30 / pow(s.getRow() * 1.496e+11, 2) ;
                // double gravmer = 6.67e-11 * 3.302e23 / pow(dist(s, planets[0]) * 1.496e+11, 2) ;
                
             double phi = acos((s.getZ() - planets[0].getZ()) / dist(s, planets[0]) );
                 double theta = asin((s.getY() - planets[0].getY()) / dist(s, planets[0]));
                s.setVx(vf(gravsun / 1.496e+11 * sin(M_PI - s.getPhi()) * cos(M_PI + s.getTheta()) * 7464960000, s.getX(), s.getVx(), .00069444));
                 s.setVy(vf(gravsun / 1.496e+11 * sin(M_PI - s.getPhi()) * sin(M_PI + s.getTheta()) * 7464960000, s.getY(), s.getVy(), .00069444));
                 s.setVz(vf(gravsun / 1.496e+11 * cos(M_PI - s.getPhi()) * 7464960000, s.getZ(), s.getVz(), .00069444));
            
            for(int n = 0; n < 8; n++){
                 phi = acos((s.getZ() - planets[n].getZ()) / dist(s, planets[n]) );
                 theta = asin((s.getY() - planets[n].getY()) / dist(s, planets[n]));
                double gravmer = 6.67e-11 * 3.302e23 / pow(dist(s, planets[n]) * 1.496e+11, 2) ;
              s.setVx(vf(gravmer / 1.496e+11 * sin( M_PI - phi) * cos(theta + M_PI) * 7464960000, s.getX(), s.getVx(), .00069444));
                 s.setVy(vf(gravmer / 1.496e+11 * sin( M_PI - phi) * sin(theta + M_PI) * 7464960000, s.getY(), s.getVy(), .00069444));
                s.setVz(vf(gravmer / 1.496e+11 * cos(M_PI - phi) * 7464960000, s.getZ(), s.getVz(), .00069444));
            }
            
                time += .00069444;
                day += .00069444;
        }
    out.close();
    return -1;
}
int main(){
  double time = 2458730.5;
  
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
     cout << "Where are you starting in our solar system? The planets are numbered (0 as Mercury, 1 as Venus...). " << "\n";
     cin >> i;
     int d = 0;
     cout << "What is your planetary destination?" << "\n";
     cin >> d;
     cout << "what is the earliest julian date you can start your journey there? (No earlier than 8/29/2019 which is 2458724.5)" << "\n";
     cin >> juliandate;
     cout << "what is the latest julian date you can start your journey there?" << "\n";
     cin >> time;
    //theta 3.70178 phi 1.53649
    //theta 0 phi 0.05235987756
    //theta  2.19911  phi 1.2217
    // theta 3.787359619
    // 4.136425469 theta
    // theta 5.078903265
    // theta 6.126100816
    //for()
    /*
    double temp;
    double theta = 3.787359619;
    double phi = M_PI /2;
    solarsystem(icarus, theta, phi , i, time, temp);
    */
    
    double mintt;
    double minta;
    double mintp;
    double mint = 3460; //time it took new horizons to reach pluto in days
    for(double t = juliandate; t < time; t++){
    for(double theta = 3.787359619; theta < 3.787359619 + 1e-10; theta+= M_PI / 180){
        for(double phi = M_PI /2 ; phi < M_PI / 2 + 1e-10; phi+= M_PI / 180){
            double temp;
            //if(solarsystem(icarus, theta, phi, i) != -1)
                cout << time << " " << theta << " " << phi << ": " << solarsystem(icarus, theta, phi, i, t, temp) << "\n";
            //if(solarsystem(icarus, theta, phi, i) != i && solarsystem(icarus, theta, phi, i) != 100){
            if(solarsystem(icarus, theta, phi, i, t, temp) != d){
            if(temp < mint){
                mint = temp;
                mintt = t;
                minta = theta;
                mintp = phi;
            }
            phi = M_PI + 1;
            theta = 2 * M_PI + 1;
            }
        }
    }
    }
    //}
    //MPI_Finalize();
    cout << "here are your best conditions to travel: " << "\n";
    cout << mint << " " << mintt << " " << minta << " " << " " << mintp;
    solarsystemp(icarus, minta, mintp, i, mintt);
        return 0;
        
}