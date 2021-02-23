#include "Simulation.h"

void Run(int pot, double m, double w, int nP, double dr, int nS, double dt, int nPrint, int nAnalysis){
    std::stringstream s;
    if (pot==0) s << "Free";
    else if (pot==1) s << "HOsc_" << w;
    else s << "SN";
    s << "_" << m;
    string name = s.str();
    
    std::stringstream nErr;
    nErr << "Analysis_" << m << "_Error";
    
    pair< Sim, vector< complex<double> > > init = initialise(nP, nS, dr, dt, m, pot, w, name, nPrint, nAnalysis);
    
    Sim sim = init.first;
    vector< complex<double> > U = init.second;
    
    setupPrintS(name,nP,dr);
    setupPrintS(nErr.str(),nP,dr);
    
    printSim(U,name,0);
    printAnalysis(U,sim,0,m,true);
    
    for (int i=1; i<=nS; i++){
        U = Update(U,sim);
        if ((i%nPrint)==0){
            printSim(U,name,i*dt);  
            cout << "Done step " << i << endl;
        }
        if ((i%nAnalysis)==0)
            printAnalysis(U,sim,i*dt,m,false);
    }    
}

pair< Sim, vector< complex<double> > > initialise(int nP, int nS, double dr, double dt, double m, int Pot, double w, string name, int nPrint, int nAnalysis){
    
    Sim sim;    
    sim.nPnts = nP;
    sim.nStps = nS;
    sim.dr = dr;
    sim.dt = dt;
    sim.M = m;
    sim.pot = Pot;
    sim.W = w;
    
    sim.nPrint = nPrint;
    sim.nAnalysis = nAnalysis;
    
    complex<double> expon = exp(complex<double>(0,-0.5*dt/(m*dr*dr)));
    complex<double> a = 0.5*(1.0+expon);
    complex<double> b = 0.5*(1.0-expon);
    expon = exp(complex<double>(0,-0.25*dt/(m*dr*dr)));
    
    sim.aa = a*a;
    sim.ab = a*b;
    sim.bb = b*b;
    sim.ac = a*expon;
    sim.bc = b*expon;
    
    vector< complex<double> > U;
    U.push_back(complex<double>(0,0));
    for (int i=1; i<nP-1; i++){double dist = i*dr;
        U.push_back(dist*exp(-0.5*dist*dist)*pow(M_PI,-0.75));
    } 
    
    U.push_back(complex<double>(0,0));
    
    printSimDetails(sim,name);
   
   return pair< Sim, vector< complex<double> > > (sim, U);    
}

vector< complex<double> > appM1M2(vector< complex<double> > U, Sim sim){
    vector< complex<double> > UNew;
    
    UNew.push_back(U[0]);
    UNew.push_back(sim.ac*U[1]+sim.ab*U[2]+sim.bb*U[3]);
    UNew.push_back(sim.bc*U[1]+sim.aa*U[2]+sim.ab*U[3]);
    
    for (int i=3; i<=sim.nPnts-5; i+=2){
        UNew.push_back(sim.ab*U[i-1]+sim.aa*U[i]+sim.ab*U[i+1]+sim.bb*U[i+2]);
        UNew.push_back(sim.bb*U[i-1]+sim.ab*U[i]+sim.aa*U[i+1]+sim.ab*U[i+2]);
    }
    
    int i = sim.nPnts-3;
    UNew.push_back(sim.ab*U[i-1]+sim.aa*U[i]+sim.bc*U[i+1]);
    UNew.push_back(sim.bb*U[i-1]+sim.ab*U[i]+sim.ac*U[i+1]);
    UNew.push_back(U[i+2]);
    
    return UNew;
}

vector< complex<double> > appM2M1(vector< complex<double> > U, Sim sim){
    vector< complex<double> > UNew;
    
    UNew.push_back(U[0]);
    UNew.push_back(sim.ac*U[1]+sim.bc*U[2]);
    
    for (int i=2; i<=sim.nPnts-4; i+=2){
        UNew.push_back(sim.ab*U[i-1]+sim.aa*U[i]+sim.ab*U[i+1]+sim.bb*U[i+2]);
        UNew.push_back(sim.bb*U[i-1]+sim.ab*U[i]+sim.aa*U[i+1]+sim.ab*U[i+2]);
    }
    
    int i = sim.nPnts-2;
    UNew.push_back(sim.bc*U[i-1]+sim.ac*U[i]);
    UNew.push_back(U[i+1]);
    
    return UNew;
}

vector<double> V(vector< complex<double> > U, Sim sim){
    vector<double> v;
    
    if (sim.pot==0){ //Free particle
        return vector<double>(sim.nPnts,0);
    } 
    
    else if (sim.pot==1){ //harmonic oscillator
        v.push_back(0);
        for (int i=1; i<sim.nPnts-1; i++){
            v.push_back(0.5*sim.M*sim.W*sim.W*i*i*sim.dr*sim.dr);
        }
        v.push_back(0);
    }
    
    else if (sim.pot==2){ //SN 
         v = fastSNV(U,sim);        
    }    
    
    return v;
}

vector<double> fastSNV(vector< complex<double> > U, Sim sim){
    vector<double> v(sim.nPnts,0);
    
    for (int i=1; i<sim.nPnts-1; i++){
        double ri = i*sim.dr;
        
        double nsUi = normSq(U[i]);
        v[i] += nsUi/ri;
        
        for (int j=i+1; j<sim.nPnts-1; j++){
            double rj = j*sim.dr;
            
            v[i] += normSq(U[j])/rj;
            v[j] += nsUi/rj;
        }
        
        v[i] *= -4*M_PI*sim.M*sim.M*sim.dr;
    }
    
    return v;
}

vector< complex<double> > appExpVU(vector< complex<double> > U, Sim sim){      
    if (sim.pot==0) return vector< complex<double> > (sim.nPnts, complex<double>(1,0));
    
    vector<double> v = V(U,sim);
    
    vector< complex<double> > UNew;
    
    UNew.push_back(complex<double>(0,0));
    for (int i=1; i<sim.nPnts-1; i++){
        UNew.push_back(U[i]*exp(complex<double>(0,-1*sim.dt*v[i])));
    }
    UNew.push_back(complex<double>(0,0));
    
    return UNew;
}

vector< complex<double> > Update(vector< complex<double> > U, Sim sim){
    vector< complex<double> > UNew = appM2M1(U,sim);
    UNew = appExpVU(UNew,sim);
    UNew = appM1M2(UNew,sim);
    return UNew;
}

vector< complex<double> > exact(Sim sim, double t){
    vector< complex<double> > ex;
    
    if (sim.pot==0 || sim.pot==2){
        double ratio = t/sim.M;
        for (int i=0; i<sim.nPnts; i++){
            double x2 = i*i*sim.dr*sim.dr;
            ex.push_back(i*sim.dr*pow(M_PI,-0.75)*pow(complex<double>(1,ratio),-1.5)*exp(complex<double>(-1,ratio)*x2*0.5/(1+ratio*ratio)));
        }
    }
    
    else if (sim.pot==1){
        double MW = sim.M*sim.W;
        double cWt = cos(sim.W*t);
        double sWt = sin(sim.W*t);
        complex<double> frac = MW/complex<double>(MW*cWt,sWt);
        for (int i=0; i<sim.nPnts; i++){
            double x = i*sim.dr;
            ex.push_back(x*pow(M_PI,-0.75)*pow(frac,1.5)*exp(-0.5*frac*(complex<double>(cWt,MW*sWt)*x*x)));
        }
    }    
    
    return ex;
}

void printSimDetails(Sim sim, string name){
    string fname = "Data/SimDetails_"+name+".txt";
    ofstream output(fname.c_str());
    
    if (output.is_open()){
        if (sim.pot==0) output << "3D Free Particle" << endl;
        else if (sim.pot==1) output << "3D Harmonic Potential" << endl;
        else output << "3D SN" << endl;
        
        output  << "nPnts: " << sim.nPnts << endl 
                << "nStps: " << sim.nStps << endl
                << "dr: " << sim.dr << endl
                << "dt: " << sim.dt << endl
                << "M: " << sim.M << endl;
        if (sim.pot==1)
            output << "W: " << sim.W << endl;
        output << "nPrint: " << sim.nPrint << endl
               << "nAnalysis: " << sim.nAnalysis << endl;
        output.close();
    }
}

void setupPrintS(string name, int nP, double dx){
    string r = "Data/"+name+"_re.txt";
    string i = "Data/"+name+"_im.txt";
    string n = "Data/"+name+"_n2.txt";
        
    ofstream re (r.c_str());
    ofstream im (i.c_str());
    ofstream n2 (n.c_str());

    if (re.is_open() && im.is_open() && n2.is_open()){
        
        re << "R"; im << "R";  n2 << "R";
        
        for (int i=0; i<nP; i++){
            double x = i*dx;
            re << ", " << x;
            im << ", " << x;
            n2 << ", " << x;
        }
        
        re.close(); im.close(); n2.close();
    }    
}

void printSim(vector< complex<double> > U, string name, double t){
    string r = "Data/"+name+"_re.txt";
    string i = "Data/"+name+"_im.txt";
    string n = "Data/"+name+"_n2.txt";
        
    ofstream re (r.c_str(), fstream::out | fstream::app);
    ofstream im (i.c_str(), fstream::out | fstream::app);
    ofstream n2 (n.c_str(), fstream::out | fstream::app);

    if (re.is_open() && im.is_open() && n2.is_open()){

        re << endl << "t = " << t << " ";
        im << endl << "t = " << t << " ";
        n2 << endl << "t = " << t << " ";

        for (int j=0; j<U.size(); j++){
            re << ", " << U[j].real();
            im << ", " << U[j].imag();
            n2 << ", " << normSq(U[j]);
        }

        re.close();
        im.close();
        n2.close();
    }
}