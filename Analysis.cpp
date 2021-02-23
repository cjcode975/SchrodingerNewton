#include "Analysis.h"

pair< vector<double>, vector< complex<double> > > analysis(vector< complex<double> > U, Sim sim, double t){
        vector< complex<double> > ex = exact(sim, t);
        vector<double> v = V(U,sim);
        
        vector<bool> found (5,false);
        double maxProb = 0;
        
        double frac = 1/(2*sim.M*sim.dr*sim.dr);
        
        pair< vector<double>, vector< complex<double> > > an = make_pair(vector<double>(12,0),vector< complex<double> >(sim.nPnts,0));
        //Energy, Total prob , E(r), Var(r), Max(P(r)), d_TV, d_HELL, r_50, r_60, r_70, r_80, r_90
        
        for (int i=1; i<sim.nPnts-1; i++){
            an.second[i] = ex[i]-U[i];
            
            double prob = 4*M_PI*normSq(U[i]);
            double prob_exact = 4*M_PI*normSq(ex[i]);            
            
            an.first[0] += (bar(U[i])*((2*frac+v[i])*U[i]-frac*U[i-1]-frac*U[i+1])).real()*sim.dr; //calculate the energy. 
            
            an.first[1] += 0.5*prob*sim.dr; //add half to the probability integral. This gives the probability of r<i
            if (an.first[1] >= 0.5 && !found[0]){ an.first[7] = i*sim.dr; found[0]=true; }
            else if (an.first[1] >= 0.6 && !found[1]){ an.first[8] = i*sim.dr; found[1]=true; }
            else if (an.first[1] >= 0.7 && !found[2]){ an.first[9] = i*sim.dr; found[2]=true; }
            else if (an.first[1] >= 0.8 && !found[3]){ an.first[10] = i*sim.dr; found[3]=true; }
            else if (an.first[1] >= 0.9 && !found[4]){ an.first[11] = i*sim.dr; found[4]=true; }
            an.first[1] += 0.5*prob*sim.dr; // add the rest
            
            an.first[2] += i*sim.dr*prob*sim.dr; // expected position. 
            
            if (prob>maxProb){
                maxProb = prob;
                an.first[4] = i*sim.dr;                
            }
            
            an.first[5] += abs(prob-prob_exact)*sim.dr; //d_TV.
            an.first[6] += (prob_exact-2*sqrt(prob*prob_exact)+prob)*sim.dr; // d_HELL.           
        }
         
        an.first[5] /= 2;
        an.first[6] = sqrt(0.5*an.first[6]);
        
        for (int i=1; i<sim.nPnts-1; i++){
                an.first[3] += (an.first[2]-i*sim.dr)*(an.first[2]-i*sim.dr)*4*M_PI*normSq(U[i])*sim.dr;
        }
        
        return an;
}

void setupPrintA(int nAn, int nS, double Dt){
    ofstream en("Data/Analysis_Energy.txt");
    ofstream un ("Data/Analysis_Unitary.txt");
    ofstream ep ("Data/Analysis_Expected_Position.txt");
    ofstream vp ("Data/Analysis_Variance_Position.txt");
    ofstream mp ("Data/Analysis_Loc_Max_Prob.txt");
    ofstream dt ("Data/Analysis_dTV.txt");
    ofstream dh ("Data/Analysis_dHELL.txt");
    ofstream r5 ("Data/Analysis_r_50.txt");
    ofstream r6 ("Data/Analysis_r_60.txt");
    ofstream r7 ("Data/Analysis_r_70.txt");
    ofstream r8 ("Data/Analysis_r_80.txt");
    ofstream r9 ("Data/Analysis_r_90.txt");
    
    if (en.is_open() && un.is_open() && ep.is_open() && vp.is_open() && mp.is_open() 
            && dt.is_open() && dh.is_open()
            && r5.is_open() && r6.is_open() && r7.is_open() && r8.is_open() && r9.is_open()){
        
        en << "Time"; un << "Time"; ep << "Time"; vp << "Time"; mp << "Time"; 
        dt << "Time"; dh << "Time";
        r5 << "Time"; r6 << "Time"; r7 << "Time"; r8 << "Time"; r9 << "Time";
        
        for (int i=0; i<=nS; i+=nAn){
            double t = i*Dt;
            en << ", " << t; un << ", " << t; ep << ", " << t; vp << ", " << t; mp << ", " << t; 
            dt << ", " << t; dh << ", " << t;
            r5 << ", " << t; r6 << ", " << t; r7 << ", " << t; r8 << ", " << t; r9 << ", " << t;
        }
        
        en.close(); un.close(); ep.close(); vp.close(); mp.close();
        dt.close(); dh.close();
        r5.close(); r6.close(); r7.close(); r8.close(); r9.close();
    }
    
}

void printAnalysis(vector< complex<double> > U, Sim sim, double t, double m, bool m_first){
    pair< vector<double>, vector< complex<double> > > an = analysis(U,sim, t);
    
    ofstream en ("Data/Analysis_Energy.txt", fstream::out | fstream::app);
    ofstream un ("Data/Analysis_Unitary.txt", fstream::out | fstream::app);
    ofstream ep ("Data/Analysis_Expected_Position.txt", fstream::out | fstream::app);
    ofstream vp ("Data/Analysis_Variance_Position.txt", fstream::out | fstream::app);
    ofstream mp ("Data/Analysis_Loc_Max_Prob.txt", fstream::out | fstream::app);
    ofstream dt ("Data/Analysis_dTV.txt", fstream::out | fstream::app);
    ofstream dh ("Data/Analysis_dHELL.txt", fstream::out | fstream::app);
    ofstream r5("Data/Analysis_r_50.txt", fstream::out | fstream::app);
    ofstream r6 ("Data/Analysis_r_60.txt", fstream::out | fstream::app);
    ofstream r7 ("Data/Analysis_r_70.txt", fstream::out | fstream::app);
    ofstream r8("Data/Analysis_r_80.txt", fstream::out | fstream::app);
    ofstream r9 ("Data/Analysis_r_90.txt", fstream::out | fstream::app);
    
    if (en.is_open() && un.is_open() && ep.is_open() && vp.is_open() && mp.is_open() 
            && dt.is_open() && dh.is_open()
            && r5.is_open() && r6.is_open() && r7.is_open() && r8.is_open() && r9.is_open()){
        
        if (m_first){
            en << endl << "Mass " << m; un << endl << "Mass " << m; 
            ep << endl << "Mass " << m; vp << endl << "Mass " << m;
            mp << endl << "Mass " << m; 
            dt << endl << "Mass " << m; dh << endl << "Mass " << m;
            r5 << endl << "Mass " << m; r6 << endl << "Mass " << m; 
            r7 << endl << "Mass " << m; r8 << endl << "Mass " << m; 
            r9 << endl << "Mass " << m;            
        }
        
        en << ", " << an.first[0]; un << ", " << an.first[1]; ep << ", " << an.first[2];
        vp << ", " << an.first[3]; mp << ", " << an.first[4]; 
        dt << ", " << an.first[5]; dh << ", " << an.first[6];
        r5 << ", " << an.first[7]; r6 << ", " << an.first[8]; r7 << ", " << an.first[9]; 
        r8 << ", " << an.first[10]; r9 << ", " << an.first[11];
        
        en.close(); un.close(); ep.close(); vp.close(); mp.close();
        dt.close(); dh.close();
        r5.close(); r6.close(); r7.close(); r8.close(); r9.close();        
    }  
    
    std::stringstream s;
    s << "Analysis_" << m << "_Error";
    printSim(an.second,s.str(),t);
}