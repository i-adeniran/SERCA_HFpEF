/*
 * This program implements an electromechanical cardiac myocyte by coupling 
 * the Bartolucci et al. 2020 Human ventricular cardiac cell model, which is an extension of 
 * the O'Hara-Virag-Varro-Rudy 2011 Human ventricular cardiac cell model, the 
 * extended Rice myofilament mechanics model (exended with ATP energetics by Tran et al. 2010) 
 * and the Cortassa 2006 mitochondrial bioenergetics model.
 *
 * O'Hara-Virag-Varro-Rudy Ventricular Cell Model
 * ==============================================
 * O'Hara T, Virág L, Varró A, Rudy Y., “Simulation of the undiseased human cardiac
 * ventricular action potential: model formulation and experimental validation”.
 * PLoS Comput Biol. 2011 May;7(5):e1002061. doi: 10.1371/journal.pcbi.1002061.
 *
 * Rice Myofilament Mechanics Model
 * ================================
 * J. Rice, F. Wang, D. Bers, and P. Tombe, “Approximate model of cooperative 
 * activation and crossbridge cycling in cardiac muscle using ordinary differential 
 * equations”. Biophysical Journal, vol. 95, pp. 2368– 2390, 2008.
 *
 * Tran et al. Extended Model with ATP Energetics
 * ==============================================
 * Tran K, Smith NP, Loiselle DS, Crampin EJ., "A metabolite-sensitive, thermodynamically 
 * constrained model of cardiac cross-bridge cycling: implications for force development 
 * during ischemia". Biophys J. 2010 Jan 20;98(2):267-76. doi: 10.1016/j.bpj.2009.10.011.
 * 
 * Cortassa Mitochondrial bioenergetics Model
 * ==========================================
 * Cortassa S, Aon MA, O'Rourke B, Jacques R, Tseng HJ, Marbán E, Winslow RL. A computational 
 * model integrating electrophysiology, contraction, and mitochondrial bioenergetics in the 
 * ventricular myocyte. Biophys J. 2006 Aug 15;91(4):1564-89. doi: 10.1529/biophysj.105.076174. 
 *
 * Compile with: g++ -std=c++23 -O3 -m64 -o ordcell myofilament_ode.cpp mitochondria_ode.cpp main_ordcell.cpp
 */

#include <numbers>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <myofilament.hpp>
#include <mitochondria.hpp>

enum class CellType { EPI, MCELL, ENDO };

class CardioMyocyte {
public:
  CardioMyocyte( CellType cell_type = CellType::ENDO, bool isHFpEF = false ) : lsoda_integrator_electrophysiology( this ) {
    is_HFpEF = isHFpEF;

    if (is_HFpEF) {
      std::cout << " -> Creating HFpEF Cardiomyocyte ...  = " << std::endl;  
      
      INAL_scaling_factor   = 2.34;
      INA_scaling_factor    = 0.6;
      VOLUME_scaling_factor = 1.25;
      vcell                 = 1000.0*std::numbers::pi*rad*rad*L * VOLUME_scaling_factor;
      Ageo                  = 2.0*std::numbers::pi*rad*rad+2.0*std::numbers::pi*rad*L;
      Acap                  = 2.0*Ageo;
      vmyo                  = 0.68*vcell;
      vmito                 = 0.26*vcell;
      vnsr                  = 0.0552*vcell;
      vjsr                  = 0.0048*vcell;
      vsr                   = 0.95*(vnsr+vjsr); //0.06*vcell
      vss                   = 0.02*vcell;
    }
    else
      std::cout << " -> Creating Healthy Cardiomyocyte ...  = " << std::endl;  

    std::cout << "      VOLUME_scaling_factor is " << VOLUME_scaling_factor << std::endl;

    std::cout << "Pacing Frequency = " << 1000.0/CL << "Hz" << std::endl;

    celltype = cell_type;

    y_electrophysiology[1]  = 1.0;        //jnca
    y_electrophysiology[2]  = 0.0;        //nca
    y_electrophysiology[3]  = 0.0;        //Ok       
    y_electrophysiology[4]  = 0.0;        //I2k      
    y_electrophysiology[5]  = 0.0;        //I1k      
    y_electrophysiology[6]  = 1.0;        //Ck       
    y_electrophysiology[7]  = 0.0;        //Okp      
    y_electrophysiology[8]  = 0.0;        //I2kp     
    y_electrophysiology[9]  = 0.0;        //I1kp     
    y_electrophysiology[10] = 1.0;        //Ckp        
    y_electrophysiology[11] = 0.0;        //I2Cak    
    y_electrophysiology[12] = 0.0;        //I1Cak    
    y_electrophysiology[13] = 0.0;        //CCak      
    y_electrophysiology[14] = 0.0;        //I2Cakp    
    y_electrophysiology[15] = 0.0;        //I1Cakp     
    y_electrophysiology[16] = 0.0;        //CCakp     
    y_electrophysiology[17] = 0.999637;   //IC1_dutta
    y_electrophysiology[18] = 6.83208e-5; //IC2_dutta
    y_electrophysiology[19] = 1.80145e-8; //C1_dutta 
    y_electrophysiology[20] = 8.26619e-5; //C2_dutta 
    y_electrophysiology[21] = 0.00015551; //O_dutta  
    y_electrophysiology[22] = 5.67623e-5; //IO_dutta 
    y_electrophysiology[23] = 0.0;        //IOb_dutta
    y_electrophysiology[24] = 0.0;        //Ob_dutta 
    y_electrophysiology[25] = 0.0;        //Cb_dutta 
    y_electrophysiology[26] = 0.0;        //D_dutta     

    rtol_electrophysiology[0]  = 0.0;
    atol_electrophysiology[0]  = 0.0;
    for ( std::size_t i = 1; i <= neq_electrophysiology; ++i ) {
      rtol_electrophysiology[i] = 1.0E-10;//1.0E-10;
      atol_electrophysiology[i] = 1.0E-10;
    }
  }


  void run( std::ofstream& out ) {
    auto NOUT         = static_cast<std::size_t>(ft/dt);

    Myofilament myofilament( is_HFpEF );
    Mitochondria mitochondria( vcell );
    
    auto count        = 1;
    auto storage_time = 1050.0; //CL + 50.0; 

    auto t_start      = std::chrono::high_resolution_clock::now();
    // for (int step = 0; step < NOUT; ++step) {
    while (t <= ft) {
      t += dt;

      //=====================
      // MITOCHONDRIA: Part 1
      //=====================
      mitochondria.calculateConstantModelParameters( vmyo, vsr, vnsr, vjsr, vss, cao, ko, nao );
      mitochondria.sharedVariableUpdateGPC_1( v );		

      // Reversal Potentials
      reversalPotentials();

      mitochondria.sharedVariableUpdateGPC_2( v, EK );

      
      
      lsoda_integrator_electrophysiology.lsoda(&CardioMyocyte::fexElectrophysiology, neq_electrophysiology, y_electrophysiology, &t_begin_electrophysiology, 
                                     tout_electrophysiology, itol_electrophysiology, rtol_electrophysiology, atol_electrophysiology, 
                                     itask_electrophysiology, &istate_electrophysiology, iopt_electrophysiology, jt_electrophysiology, 
                                     iwork1_electrophysiology, iwork2_electrophysiology, iwork5_electrophysiology, iwork6_electrophysiology, 
                                     iwork7_electrophysiology, iwork8_electrophysiology, iwork9_electrophysiology, rwork1_electrophysiology, 
                                     rwork5_electrophysiology, rwork6_electrophysiology, rwork7_electrophysiology, nullptr);
      if (istate_electrophysiology <= 0) {
        printf("Error istate = %d\n", istate_electrophysiology);
        std::exit(0);
      }
      tout_electrophysiology += dt;

      RGC();
      stimulus();
      vo = v;
      voltage();
      dVdt_APD();
      FBC( myofilament.TropTot_dt, mitochondria.VnaCa, mitochondria.Vuni );

      //=====================
      // MITOCHONDRIA: Part 2
      //=====================
			/*After ATP mod**/
			mitochondria.calculateV_AM( myofilament.fappT, myofilament.Pxb );
			/**after Mitochondria mod**/
			mitochondria.calculateATPm();
			mitochondria.calculateDmuH();
			mitochondria.calculateNAD();

			mitochondria.sharedVariableUpdateMitochondria();

			mitochondria.calculateVCS();
			mitochondria.calculateVACO();
			mitochondria.calculateVIDH();
			mitochondria.calculateVKGDH();
			mitochondria.calculateVSL();
			mitochondria.calculateVSDH();
			mitochondria.calculateVFH();
			mitochondria.calculateVMDH();
			mitochondria.calculateVAAT();
			mitochondria.calculateVNO_VHNe();
			mitochondria.calculateVFO_VHFe();
			mitochondria.calculateVATPase_Vhu();
			mitochondria.calculateVANT_Vhleak();

      computeMyofilamentForce( myofilament, mitochondria );

      //=====================
      // MITOCHONDRIA: Part 3
      //=====================
			mitochondria.calculateVuni( cai );		    // Moved up
			mitochondria.calculateVnaCa( cai, nai );	// from below (in the Mitene folder)

      computeMitochondria( mitochondria );

      if ( count %  10000000 == 0 ) std:: cout << t*100.0/ft << "% complete" << std::endl; // Display progress
      
      if ( t >= (ft - storage_time) ) {
        out	<< t 	                                      << std::setw(15)  //0
            << v				                                << std::setw(15)  //1
            << cai				                              << std::setw(15)  //2
            << casr                                     << std::setw(15)  //3
            << myofilament.SL                           << std::setw(15)  //4
            << myofilament.F_passive                    << std::setw(15)  //5
            << myofilament.Fnorm                        << std::setw(15)  //6
            << myofilament.F_total                      << std::setw(15)  //7
            << mitochondria.Dpsi      		              << std::setw(15)  //8
            << mitochondria.CrPi_cyto/mitochondria.ATPi	<< std::setw(15)  //9
            << mitochondria.CrPi_mito/mitochondria.ATPi	<< std::setw(15)  //10
            << mitochondria.VNO				                  << std::setw(15)  //11
            << mitochondria.V_AM         	              << std::setw(15)  //12   
            << mitochondria.Vuni                        << std::setw(15)  //13
            << Jup         	                            << std::setw(15)  //14
            << ICaL        	                            << std::setw(15)  //15
            << Jrel        	                            << std::setw(15)  //16
            << INaCa       	                            << std::setw(15)  //17
            << CaMKa       	                            << std::setw(15)  //18
            << mitochondria.Cam                         << std::setw(15)  //19 
            << mitochondria.VANT                        << std::setw(15)  //20 
            << mitochondria.VATPase                     << std::setw(15)  //21 
            << mitochondria.NADH                        << std::setw(15)  //22
            << mitochondria.NAD                         << std::setw(15)  //23
            << mitochondria.NAD/mitochondria.NADH       << std::setw(15)  //24 
            << mitochondria.ADPm                        << std::setw(15)  //25
            << mitochondria.ATPm                        << std::setw(15)  //26
            << mitochondria.ADP                         << std::setw(15)  //27
            << mitochondria.ATPi                        << std::setw(15)  //28
            << mitochondria.ATPi/mitochondria.ADP       << std::setw(15)  //29
            << mitochondria.VCK_mito                    << std::setw(15)  //30
            << mitochondria.VCK_cyto                    << std::setw(15)  //31
            << std::endl;
      }

      ++count;
    }

    out.close();
    myofilament.lsoda_integrator.n_lsoda_terminate();
    mitochondria.lsoda_integrator.n_lsoda_terminate();
    lsoda_integrator_electrophysiology.n_lsoda_terminate();

    auto t_stop       = std::chrono::high_resolution_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::minutes>(t_stop - t_start);
    std::cout << "Total computing time = " << time_elapsed.count() << " minutes." << std::endl;

    std::cout << "APD: " << APD << std::endl;
  }


  void voltage() {
    v += -dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+IClCa+IClbk+Ist);
  }

  
  void stimulus() {
    if ( (t > (start+n*CL)) and (t < (start+duration+n*CL-dt)) ) {
      if ( Ist == 0 )
        vrest=v;
  
      Ist=amp;
    }    
    else if ( t > (start+duration+n*CL-dt) ) {   
      Ist = 0.0;
      n   = n+1;
    }
  }


  void dVdt_APD() {
    vdot_old  = vdot;
    vdot      = (v-vo)/dt;
    if ( (APD_flag == 0) and (v > -40.0) and (vdot < vdot_old) ) {
      vdot_max    = vdot_old;
      t_vdot_max  = t-dt;
      APD_flag    = 1;
    }
    if ( (APD_flag == 1) and (v < 0.9*vrest) ) {
      APD=t-t_vdot_max;
      APD_flag=0;
      }
  }

  void computeMyofilamentForce( Myofilament &myofilament, Mitochondria &mitochondria ) {
    std::array<double, 3> myofilament_data = { cai, mitochondria.ADP, mitochondria.ATPi }; 
    myofilament.lsoda_integrator.lsoda(&Myofilament::fexMyofilament, myofilament.neq_myofilament, myofilament.y_myofilament, &myofilament.t_begin_myofilament, 
                            myofilament.tout_myofilament, myofilament.itol_myofilament, myofilament.rtol_myofilament, myofilament.atol_myofilament, 
                            myofilament.itask_myofilament, &myofilament.istate_myofilament, myofilament.iopt_myofilament, myofilament.jt_myofilament, 
                            myofilament.iwork1_myofilament, myofilament.iwork2_myofilament, myofilament.iwork5_myofilament, myofilament.iwork6_myofilament, 
                            myofilament.iwork7_myofilament, myofilament.iwork8_myofilament, myofilament.iwork9_myofilament, myofilament.rwork1_myofilament, 
                            myofilament.rwork5_myofilament, myofilament.rwork6_myofilament, myofilament.rwork7_myofilament, myofilament_data.data());
    if (myofilament.istate_myofilament <= 0) {
        printf("Error istate = %d\n", myofilament.istate_myofilament);
        std::exit(0);
    }
    myofilament.tout_myofilament += dt;
  }

  void computeMitochondria( Mitochondria &mitochondria ) {
    std::array<double, 3> mitochondrial_data = { Jup, INaK, IpCa }; 
    mitochondria.lsoda_integrator.lsoda(&Mitochondria::fexMitochondria, mitochondria.neq_mitochondria, mitochondria.y_mitochondria, &mitochondria.t_begin_mitochondria, 
                            mitochondria.tout_mitochondria, mitochondria.itol_mitochondria, mitochondria.rtol_mitochondria, mitochondria.atol_mitochondria, 
                            mitochondria.itask_mitochondria, &mitochondria.istate_mitochondria, mitochondria.iopt_mitochondria, mitochondria.jt_mitochondria, 
                            mitochondria.iwork1_mitochondria, mitochondria.iwork2_mitochondria, mitochondria.iwork5_mitochondria, mitochondria.iwork6_mitochondria, 
                            mitochondria.iwork7_mitochondria, mitochondria.iwork8_mitochondria, mitochondria.iwork9_mitochondria, mitochondria.rwork1_mitochondria, 
                            mitochondria.rwork5_mitochondria, mitochondria.rwork6_mitochondria, mitochondria.rwork7_mitochondria, mitochondrial_data.data());
    if (mitochondria.istate_mitochondria <= 0) {
        printf("Error istate = %d\n", mitochondria.istate_mitochondria);
        std::exit(0);
    }
    mitochondria.tout_mitochondria += dt;
  }


  void fexElectrophysiology( double t, double *y, double *ydot, void *data ) {
    jnca      = y[0];
    nca       = y[1];
    Ok        = y[2]; 
    I2k       = y[3]; 
    I1k       = y[4]; 
    Ck        = y[5]; 
    Okp       = y[6]; 
    I2kp      = y[7]; 
    I1kp      = y[8]; 
    Ckp       = y[9];   
    I2Cak     = y[10]; 
    I1Cak     = y[11]; 
    CCak      = y[12]; 
    I2Cakp    = y[13]; 
    I1Cakp    = y[14];  
    CCakp     = y[15]; 
    IC1_dutta = y[16];
    IC2_dutta = y[17];
    C1_dutta  = y[18];
    C2_dutta  = y[19];
    O_dutta   = y[20];
    IO_dutta  = y[21];
    IOb_dutta = y[22];
    Ob_dutta  = y[23];
    Cb_dutta  = y[24];
    D_dutta   = y[25];
    

    // ICaL, ICaNa, ICaK current %%%
    // ICaL parameters
    double bGCaL  = 1.0;
    double cPCa   = 0.9;
    double kCDI   = 9.0;
    double fICaLp = (1.0/(1.0+KmCaMK/CaMKa))*(1.0-undo_p);
    // Up/Down rates
    double r_down = (1.0e-1)*(1.0-undo_CDI);
    double r_up   = (r_down*nca/(1.0-nca))*(1.0-undo_CDI);

    ///===
    double jncass = 1.0/(1.0+exp((v+19.58+25.0)/3.696));
    double tjnca  = 1.0;
    
    ydot[0]       = ((jncass-jnca)/tjnca); // jnca
    double Kmn    = 0.05;
    double k2n    = 1000.0;
    double km2n   = jnca*150.0;
    double anca   = (1.0-nca)/std::pow((1.0+Kmn/cass), 4.0); 

    ydot[1]       = (anca*k2n-nca*km2n); // nca
    double dss    = 1.0/(1.0+exp((-(v+3.940))/4.230));
    double td     = 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
    
    alpha         = dss/td;
    beta          = (1.0-dss) / td;
    ///===

    double vffrt  = v*F*F/(R*T);
    double vfrt   = v*F/(R*T);

    // Recovery (jca)
    double jcass_new = 1.0/(1.0+exp((v+19.58)/3.696));
    double jcass_VD  = jcass_new;
    double jcass_CD  = jcass_new;
    double jcass_VDp = jcass_new;
    double jcass_CDp = jcass_new;
    double tjca_new  = 35.0 + 350.0*exp(-std::pow((v-(-20.0)), 2.0)/std::pow(2.0*10.0, 2.0));
    double tjca_VD   = tjca_new;
    double tjca_VDp  = tjca_new;
    double tjca_CD   = tjca_new;
    double tjca_CDp  = tjca_new;
    // psi and omega rates
    double psi_VD=jcass_VD/tjca_VD;
    double psi_VDp=jcass_VDp/tjca_VDp;
    double psi_CD=jcass_CD/tjca_CD;
    double psi_CDp=jcass_CDp/tjca_CDp;
    double omega_VD =(1.0-jcass_VD)/tjca_VD;
    double omega_VDp=(1.0-jcass_VDp)/tjca_VDp;
    double omega_CD =(1.0-jcass_CD)/tjca_CD;
    double omega_CDp=(1.0-jcass_CDp)/tjca_CDp;

    double f1ss_0    = 0.8 / (1.0+exp((v+19.58)/3.696)) + 0.2;
    double tf1_0     = 1.0*(70.0 + 1.2/ (0.0045*exp((v+20.0)/(-50.0))+0.0045*exp((v+30.0)/10.0)));
    double ktaup     = 2.5;
    double gamma_VD  = (1.0-f1ss_0)/ tf1_0;
    double delta_VD  = f1ss_0  / tf1_0;
    double gamma_VDp = gamma_VD/ktaup;
    double delta_VDp = delta_VD/ktaup;
    double gamma_CD  = gamma_VD*kCDI;
    double delta_CD  = delta_VD*kCDI;
    double gamma_CDp = gamma_VDp*kCDI;
    double delta_CDp = delta_VDp*kCDI;
    double tf1_VD    = 1.0/(gamma_VD+delta_VD);
    double tf1_CD    = 1.0/(gamma_CD+delta_CD);
    double f1ss_VD   = gamma_VD / (gamma_VD+delta_VD);
    double f1ss_CD   = gamma_CD / (gamma_CD+delta_CD);

    double tf2_new = 1.0*(100.0 + 0.0/ (0.0035*exp((v+5.0)/(-84.0))+0.0035*exp((v+5.0)/4.0)));
    double tf2_VD  = tf2_new;
    double tf2_CD  = tf2_VD/kCDI;
    double tf2_VDp = tf2_new*ktaup;
    double tf2_CDp = tf2_VD/kCDI*ktaup;
    // Reversibility
    double theta_VD  = alpha*gamma_VD*psi_VD/tf2_VD/(alpha*gamma_VD*psi_VD+beta*delta_VD*omega_VD);
    double theta_CD  = alpha*gamma_CD*psi_CD/tf2_CD/(alpha*gamma_CD*psi_CD+beta*delta_CD*omega_CD);
    double theta_VDp = alpha*gamma_VDp*psi_VDp/tf2_VDp/(alpha*gamma_VDp*psi_VDp+beta*delta_VDp*omega_VDp);
    double theta_CDp = alpha*gamma_CDp*psi_CDp/tf2_CDp/(alpha*gamma_CDp*psi_CDp+beta*delta_CDp*omega_CDp);
    double eta_VD    = 1.0/tf2_VD - theta_VD;
    double eta_VDp   = 1.0/tf2_VDp - theta_VDp;
    double eta_CD    = 1.0/tf2_CD - theta_CD;
    double eta_CDp   = 1.0/tf2_CDp - theta_CDp;

    tf2_VD           = 1.0/(eta_VD+theta_VD);
    tf2_CD           = 1.0/(eta_CD+theta_CD);
    double f2ss_VD   = eta_VD/(eta_VD+theta_VD);
    double f2ss_CD   = eta_CD/(eta_CD+theta_CD);

    // Driving Forces
    double PhiCaL  = 4.0*vffrt*(1.2*cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
    double PhiCaNa = 1.0*vffrt*(0.75*nass*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
    double PhiCaK  = 1.0*vffrt*(0.75*kss*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
    double PCa     = 0.0001*cPCa;
    if ( celltype == CellType::EPI )
      PCa = PCa*1.4;
    else if ( celltype == CellType::MCELL )
      PCa = PCa*2.0;
    double PCap    = 1.1*PCa;
    double PCaNa   = 0.00125*PCa;
    double PCaK    = 3.574e-4*PCa;
    double PCaNap  = 0.00125*PCap;
    double PCaKp   = 3.574e-4*PCap;

    // Markov Model: VDI states
    double OCak   = 1.0-CCak-I1Cak-I2Cak-Ck-I1k-I2k-Ok;
    double OCakp  = 1.0-CCakp-I1Cakp-I2Cakp-Ckp-I1kp-I2kp-Okp;
    ydot[2]       = (alpha*Ck       + delta_VD*I1k    - (beta+gamma_VD)*Ok       - r_up*Ok   + r_down*OCak);          // 2: dOk     
    ydot[3]       = (eta_VD*I1k     + omega_VD*Ck     - (theta_VD+psi_VD)*I2k    - r_up*I2k  + r_down*I2Cak);         // 3: dI2k    
    ydot[4]       = (theta_VD*I2k   + gamma_VD*Ok     - (eta_VD+delta_VD)*I1k    - r_up*I1k  + r_down*I1Cak);        // 4: dI1k    
    ydot[5]       = (beta*Ok        + psi_VD*I2k      - (omega_VD+alpha)*Ck      - r_up*Ck   + r_down*CCak);          // 5: dCk     
    ydot[6]       = (alpha*Ckp      + delta_VDp*I1kp  - (beta+gamma_VDp)*Okp     - r_up*Okp  + r_down*OCakp);         // 6: dOkp    
    ydot[7]       = (eta_VDp*I1kp   + omega_VDp*Ckp   - (theta_VDp+psi_VDp)*I2kp - r_up*I2kp + r_down*I2Cakp);        // 7: dI2kp   
    ydot[8]       = (theta_VDp*I2kp + gamma_VDp*Okp   - (eta_VDp+delta_VDp)*I1kp - r_up*I1kp + r_down*I1Cakp);        // 8: dI1kp     
    ydot[9]       = (beta*Okp       + psi_VDp*I2kp    - (omega_VDp+alpha)*Ckp    - r_up*Ckp  + r_down*CCakp);         // 9: dCkp    
    // Markov Model: CDI states
    ydot[10]      = (eta_CD*I1Cak     + omega_CD*CCak   - (theta_CD+psi_CD)*I2Cak + r_up*I2k - r_down*I2Cak);        // 10: I2Cak   
    ydot[11]      = (theta_CD*I2Cak   + gamma_CD*OCak   - (eta_CD+delta_CD)*I1Cak + r_up*I1k - r_down*I1Cak);        // 11: I1Cak   
    ydot[12]      = (beta*OCak     + psi_CD*I2Cak    - (omega_CD+alpha)*CCak + r_up*Ck - r_down*CCak);               // 12: CCak   
    ydot[13]      = (eta_CDp*I1Cakp   + omega_CDp*CCakp - (theta_CDp+psi_CDp)*I2Cakp + r_up*I2kp - r_down*I2Cakp);   // 13: I2Cakp 
    ydot[14]      = (theta_CDp*I2Cakp + gamma_CDp*OCakp - (eta_CDp+delta_CDp)*I1Cakp + r_up*I1kp - r_down*I1Cakp);   // 14: I1Cakp 
    ydot[15]      = (beta*OCakp   + psi_CDp*I2Cakp  - (omega_CDp+alpha)*CCakp + r_up*Ckp - r_down*CCakp);            // 15: CCakp  

    // ICaL ICaNa ICaK currents
    double ICaL_VD   = PCa    * PhiCaL  * Ok;
    double ICaL_VDp  = PCap   * PhiCaL  * Okp;
    double ICaL_CD   = PCa    * PhiCaL  * OCak;
    double ICaL_CDp  = PCap   * PhiCaL  * OCakp;
    double ICaNa_VD  = PCaNa  * PhiCaNa * Ok;
    double ICaNa_VDp = PCaNap * PhiCaNa * Okp;
    double ICaNa_CD  = PCaNa  * PhiCaNa * OCak;
    double ICaNa_CDp = PCaNap * PhiCaNa * OCakp;
    double ICaK_VD   = PCaK   * PhiCaK  * Ok;
    double ICaK_VDp  = PCaKp  * PhiCaK  * Okp;
    double ICaK_CD   = PCaK   * PhiCaK  * OCak;
    double ICaK_CDp  = PCaKp  * PhiCaK  * OCakp;
    // ICaL VD vs CD & ICaL p vs np
    double ICaLnp  = ICaL_VD  + ICaL_CD;
    double ICaLp   = ICaL_VDp + ICaL_CDp;
    double ICaLVD  = ICaL_VD  * (1.0-fICaLp)  + ICaL_VDp*fICaLp;
    double ICaLCD  = ICaL_CD  * (1.0-fICaLp)  + ICaL_CDp*fICaLp;

    double ICaNanp = ICaNa_VD  + ICaNa_CD;
    double ICaNap  = ICaNa_VDp + ICaNa_CDp;
    double ICaKnp  = ICaK_VD  + ICaK_CD;
    double ICaKp   = ICaK_VDp + ICaK_CDp;

    ICaL  = (ICaLp*fICaLp + ICaLnp  *(1.0-fICaLp))*bGCaL;
    ICaNa = (ICaNap*fICaLp + ICaNanp*(1.0-fICaLp))*bGCaL;
    ICaK  = (ICaKp*fICaLp + ICaKnp  *(1.0-fICaLp))*bGCaL;

    // ICaL conductance*
    gICaL = ICaL/PhiCaL;

    // IKr current from DUTTA - Markovian formulation
    // IKr parameters
    double bGKr        = 1.2;
    double k1_dutta    = 1.0;
    double k2_dutta    = 1.0;
    double k3_dutta    = 1.0;
    double k4_dutta    = 1.0;
    double k11_dutta   = 1.0;
    double k21_dutta   = 1.0;
    double k31_dutta   = 1.0;
    double k41_dutta   = 1.0;
    double k51_dutta   = 1.0;
    double k61_dutta   = 1.0;
    double k52_dutta   = 1.0;
    double k62_dutta   = 1.0;
    double k53_dutta   = 1.0;
    double k63_dutta   = 1.0;
    double kD_dutta    = 1.0;
    
    double A1_dutta    = 0.0264;      // per_millisecond (in IKr)
    double A11_dutta   = 0.0007868;   // per_millisecond (in IKr)
    double A2_dutta    = 4.986e-6;    // per_millisecond (in IKr)
    double A21_dutta   = 5.455e-6;    // per_millisecond (in IKr)
    double A3_dutta    = 0.001214;    // per_millisecond (in IKr)
    double A31_dutta   = 0.005509;    // per_millisecond (in IKr)
    double A4_dutta    = 1.854e-5;    // per_millisecond (in IKr)
    double A41_dutta   = 0.001416;    // per_millisecond (in IKr)
    double A51_dutta   = 0.4492;      // per_millisecond (in IKr)
    double A52_dutta   = 0.3181;      // per_millisecond (in IKr)
    double A53_dutta   = 0.149;       // per_millisecond (in IKr)
    double A61_dutta   = 0.01241;     // per_millisecond (in IKr)
    double A62_dutta   = 0.3226;      // per_millisecond (in IKr)
    double A63_dutta   = 0.008978;    // per_millisecond (in IKr)
    double B1_dutta    = 4.631e-5;    // per_millivolt (in IKr)
    double B11_dutta   = 1.535e-8;    // per_millivolt (in IKr)
    double B2_dutta    = -0.004226;   // per_millivolt (in IKr)
    double B21_dutta   = -0.1688;     // per_millivolt (in IKr)
    double B3_dutta    = 0.008516;    // per_millivolt (in IKr)
    double B31_dutta   = 7.771e-9;    // per_millivolt (in IKr)
    double B4_dutta    = -0.04641;    // per_millivolt (in IKr)
    double B41_dutta   = -0.02877;    // per_millivolt (in IKr)
    double B51_dutta   = 0.008595;    // per_millivolt (in IKr)
    double B52_dutta   = 3.613e-8;    // per_millivolt (in IKr)
    double B53_dutta   = 0.004668;    // per_millivolt (in IKr)
    double B61_dutta   = 0.1725;      // per_millivolt (in IKr)
    double B62_dutta   = -0.0006575;  // per_millivolt (in IKr)
    double B63_dutta   = -0.02215;    // per_millivolt (in IKr)
    // GKr_b       = 0.046585   % milliS_per_microF (in IKr)
    double Kmax_IKr      = 0.0;       // dimensionless (in IKr)
    double Kt_dutta      = 3.5e-5;    // per_millisecond (in IKr)
    double Ku_dutta      = 0.0;       // per_millisecond (in IKr)
    double Temp_dutta    = 37.0;      // celsius (in IKr)
    double halfmax_dutta = 1.0;       // dimensionless (in IKr)
    double n_dutta       = 1.0;       // dimensionless (in IKr)
    double q1_dutta      = 4.843;     // dimensionless (in IKr)
    double q11_dutta     = 4.942;     // dimensionless (in IKr)
    double q2_dutta      = 4.23;      // dimensionless (in IKr)
    double q21_dutta     = 4.156;     // dimensionless (in IKr)
    double q3_dutta      = 4.962;     // dimensionless (in IKr)
    double q31_dutta     = 4.22;      // dimensionless (in IKr)
    double q4_dutta      = 3.769;     // dimensionless (in IKr)
    double q41_dutta     = 1.459;     // dimensionless (in IKr)
    double q51_dutta     = 5.0;       // dimensionless (in IKr)
    double q52_dutta     = 4.663;     // dimensionless (in IKr)
    double q53_dutta     = 2.412;     // dimensionless (in IKr)
    double q61_dutta     = 5.568;     // dimensionless (in IKr)
    double q62_dutta     = 5.0;       // dimensionless (in IKr)
    double q63_dutta     = 5.682;     // dimensionless (in IKr)
    double vhalf_dutta   = 1.0;       // millivolt (in IKr)

    double rate1_dutta   = A1_dutta*exp(B1_dutta*v)*exp((Temp_dutta-20.0)*log(q1_dutta)/10.0)*k1_dutta;
    double rate2_dutta   = A2_dutta*exp(B2_dutta*v)*exp((Temp_dutta-20.0)*log(q2_dutta)/10.0)*k2_dutta;
    double rate3_dutta   = A3_dutta*exp(B3_dutta*v)*exp((Temp_dutta-20.0)*log(q3_dutta)/10.0)*k3_dutta;
    double rate4_dutta   = A4_dutta*exp(B4_dutta*v)*exp((Temp_dutta-20.0)*log(q4_dutta)/10.0)*k4_dutta;
    double rate11_dutta  = A11_dutta*exp(B11_dutta*v)*exp((Temp_dutta-20.0)*log(q11_dutta)/10.0)*k11_dutta;
    double rate21_dutta  = A21_dutta*exp(B21_dutta*v)*exp((Temp_dutta-20.0)*log(q21_dutta)/10.0)*k21_dutta;
    double rate31_dutta  = (A31_dutta*exp(B31_dutta*v)*exp((Temp_dutta-20.0)*log(q31_dutta)/10.0))*k31_dutta;
    double rate41_dutta  = (A41_dutta*exp(B41_dutta*v)*exp((Temp_dutta-20.0)*log(q41_dutta)/10.0))*k41_dutta;
    double rate51_dutta  = A51_dutta*exp(B51_dutta*v)*exp((Temp_dutta-20.0)*log(q51_dutta)/10.0)*k51_dutta;
    double rate52_dutta  = A52_dutta*exp(B52_dutta*v)*exp((Temp_dutta-20.0)*log(q52_dutta)/10.0)*k52_dutta;
    double rate53_dutta  = A53_dutta*exp(B53_dutta*v)*exp((Temp_dutta-20.0)*log(q53_dutta)/10.0)*k53_dutta;
    double rate61_dutta  = A61_dutta*exp(B61_dutta*v)*exp((Temp_dutta-20.0)*log(q61_dutta)/10.0)*k61_dutta;
    double rate62_dutta  = A62_dutta*exp(B62_dutta*v)*exp((Temp_dutta-20.0)*log(q62_dutta)/10.0)*k62_dutta;
    double rate63_dutta  = A63_dutta*exp(B63_dutta*v)*exp((Temp_dutta-20.0)*log(q63_dutta)/10.0)*k63_dutta;
    double rateD_dutta   = Ku_dutta*exp(n_dutta*log(D_dutta))/(exp(n_dutta*log(D_dutta))+halfmax_dutta)*kD_dutta;
    double vrect_dutta   = 1.0/(1.0+exp(-(v-vhalf_dutta)/6.789));

    ydot[16]             = (-(rate11_dutta*IC1_dutta-rate21_dutta*IC2_dutta)+rate51_dutta*C1_dutta-rate61_dutta*IC1_dutta);                                                                       // dIC1_dutta 
    ydot[17]             = (rate11_dutta*IC1_dutta-rate21_dutta*IC2_dutta-(rate3_dutta*IC2_dutta-rate4_dutta*IO_dutta)+rate52_dutta*C2_dutta-rate62_dutta*IC2_dutta);                             // dIC2_dutta 
    ydot[18]             = (-(rate1_dutta*C1_dutta-rate2_dutta*C2_dutta)-(rate51_dutta*C1_dutta-rate61_dutta*IC1_dutta));                                                                         // dC1_dutta  
    ydot[19]             = (rate1_dutta*C1_dutta-rate2_dutta*C2_dutta-(rate31_dutta*C2_dutta-rate41_dutta*O_dutta)-(rate52_dutta*C2_dutta-rate62_dutta*IC2_dutta));                               // dC2_dutta  
    ydot[20]             = (rate31_dutta*C2_dutta-rate41_dutta*O_dutta-(rate53_dutta*O_dutta-rate63_dutta*IO_dutta)-(Kmax_IKr*rateD_dutta*O_dutta-Ku_dutta*Ob_dutta));                            // dO_dutta   
    ydot[21]             = (rate3_dutta*IC2_dutta-rate4_dutta*IO_dutta+rate53_dutta*O_dutta-rate63_dutta*IO_dutta-(Kmax_IKr*rateD_dutta*IO_dutta-Ku_dutta*rate53_dutta/rate63_dutta*IOb_dutta));  // dIO_dutta  
    ydot[22]             = (Kmax_IKr*rateD_dutta*IO_dutta-Ku_dutta*rate53_dutta/rate63_dutta*IOb_dutta+Kt_dutta*vrect_dutta*Cb_dutta-Kt_dutta*IOb_dutta);                                         // dIOb_dutta 
    ydot[23]             = (Kmax_IKr*rateD_dutta*O_dutta-Ku_dutta*Ob_dutta+Kt_dutta*vrect_dutta*Cb_dutta-Kt_dutta*Ob_dutta);                                                                      // dOb_dutta  
    ydot[24]             = (-(Kt_dutta*vrect_dutta*Cb_dutta-Kt_dutta*Ob_dutta)-(Kt_dutta*vrect_dutta*Cb_dutta-Kt_dutta*IOb_dutta));                                                               // dCb_dutta  
    ydot[25]             = 0.0;                                                                                                                                                                   // dD_dutta   

    double GKr = 0.046*bGKr;
    if (celltype == CellType::EPI)
      GKr = GKr*1.1; // Modified respect ORd (*1.3)
    else if (celltype == CellType::MCELL)
      GKr = GKr*0.8;
    IKr = GKr*sqrt(ko/5.4)*O_dutta*(v-EK);
  }


  void reversalPotentials() {
    ENa               = (R*T/F)*log(nao/nai);
    EK                = (R*T/F)*log(ko/ki);
    const double PKNa = 0.01833;
    EKs               = (R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai));
    EKshift           = 8.0; // Shift for compesation in Ito,IKs,IK1 due to LJP
  }


  void RGC () {
    CaMKb = CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
    CaMKa = ( CaMKb+CaMKt ) * CaMKa_scaling_factor;
    double vffrt = v*F*F/(R*T);
    double vfrt = v*F/(R*T);
    
    // INa
    double bGNa = 0.27;
    double thf_shift = 0.075;
    double mss = 1.0/(1.0+exp((-(v+39.57))/9.871));
    double tm = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
    m = mss-(mss-m)*exp(-dt/tm);
    double hss = 1.0/(1.0+exp((v+78.5)/6.22));
    double thf = 1.0/(3.6860e-6*exp(-(v+3.8875)/7.8579)+16.0*exp((v-0.4963)/9.1843))+thf_shift;
    double ths = 1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
    double Ahf = 0.99;
    double Ahs = 1.0-Ahf;
    hf = hss-(hss-hf)*exp(-dt/thf);
    hs = hss-(hss-hs)*exp(-dt/ths);
    double h = Ahf*hf+Ahs*hs;
    double jss = hss;
    double tj = (4.8590+1.0/(0.8628*exp(-(v+116.7258)/7.6005)+1.1096*exp((v+6.2719)/9.0358)));
    j = jss-(jss-j)*exp(-dt/tj);
    double hssp = 1.0/(1.0+exp((v+84.7)/6.22));
    double thsp = 3.0*ths;
    hsp = hssp-(hssp-hsp)*exp(-dt/thsp);
    double hp = Ahf*hf+Ahs*hsp;
    double tjp = 1.46*tj;
    jp = jss-(jss-jp)*exp(-dt/tjp);
    double GNa = 75.0 * bGNa;
    double fINap = (1.0/(1.0+KmCaMK/CaMKa));
    INa = GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);
    if (is_HFpEF)
      INa *= INA_scaling_factor;
    
    // INaL
    double bGnal = 2.8;
    double mLss = 1.0/(1.0+exp((-(v+42.85))/5.264));
    double tmL = tm;
    mL = mLss-(mLss-mL)*exp(-dt/tmL);
    double hLss=1.0/(1.0+exp((v+87.61)/7.488));
    double thL=200.0;
    hL = hLss-(hLss-hL)*exp(-dt/thL);
    double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
    double thLp=3.0*thL;
    hLp = hLssp-(hLssp-hLp)*exp(-dt/thLp);
    double GNaL=0.0075 * bGnal;
    if ( celltype == CellType::EPI ) GNaL = GNaL * 0.7;
    double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
    INaL = GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
    if (is_HFpEF) 
      INaL *= INAL_scaling_factor;
    
    // Ito
    double ass=1.0/(1.0+exp((-(v+EKshift-14.34))/14.82));
    double ta=1.0515 / (1.0/(1.2089*(1.0+exp(-(v+EKshift-18.4099)/29.3814)))+3.5/(1.0+exp((v+EKshift+100.0)/29.3814)));
    a = ass-(ass-a)*exp(-dt/ta);
    double iss=1.0/(1.0+exp((v+EKshift+43.94)/5.711));
    double delta_epi = 0.0;
    if ( celltype == CellType::EPI )
      delta_epi = 1.0-(0.95/(1.0+exp((v+EKshift+70.0)/5.0)));
    else
      delta_epi = 1.0;
    double tiF=4.562+1.0/(0.3933*exp((-(v+EKshift+100.0))/100.0)+0.08004*exp((v+EKshift+50.0)/16.59));
    double tiS=23.62+1.0/(0.001416*exp((-(v+EKshift+96.52))/59.05)+1.780e-8*exp((v+EKshift+114.1)/8.079));
    tiF = tiF * delta_epi;
    tiS = tiS * delta_epi;
    double AiF=1.0/(1.0+exp((v-213.6+EKshift)/151.2));
    double AiS=1.0-AiF;
    iF = iss-(iss-iF)*exp(-dt/tiF);
    iS = iss-(iss-iS)*exp(-dt/tiS);
    double i=AiF*iF+AiS*iS;
    double assp=1.0/(1.0+exp((-(v+EKshift-24.34))/14.82));
    ap = assp-(assp-ap)*exp(-dt/ta);
    double dti_develop=1.354+1.0e-4/(exp((v+EKshift-167.4)/15.89)+exp(-(v+EKshift-12.23)/0.2154));
    double dti_recover=1.0-0.5/(1.0+exp((v+EKshift+70.0)/20.0));
    double tiFp=dti_develop*dti_recover*tiF;
    double tiSp=dti_develop*dti_recover*tiS;
    iFp = iss-(iss-iFp)*exp(-dt/tiFp);
    iSp = iss-(iss-iSp)*exp(-dt/tiSp);
    double ip=AiF*iFp+AiS*iSp;
    double Gto=0.02;
    if ( (celltype == CellType::EPI) or (celltype == CellType::MCELL) ) Gto = Gto * 4.0;
    double fItop=(1.0/(1.0+KmCaMK/CaMKa));
    Ito = Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);
    
    // IKs
    double bGKs = 2.0;
    double xs1ss=1.0/(1.0+exp((-(v+11.60 + EKshift))/8.932));
    double txs1=817.3+1.0/(2.326e-4*exp((v+48.28+ EKshift)/17.80)+0.001292*exp((-(v+210.0+ EKshift))/230.0));
    xs1 = xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
    double xs2ss=xs1ss;
    double txs2=1.0/(0.01*exp((v-50.0+ EKshift)/20.0)+0.0193*exp((-(v+66.54+ EKshift))/31.0));
    xs2 = xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
    double KsCa=1.0+0.6/(1.0+std::pow((3.8e-5/cai), 1.4));
    double GKs = 0.0034 * bGKs;
    if ( celltype == CellType::EPI ) GKs = GKs * 1.4;
    IKs = GKs*KsCa*xs1*xs2*(v-EKs);
    
    // IK1
    double kslope_rk1 = 1.09;
    double bGK1 = 0.71;
    double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59+EKshift)/(1.5692*ko+3.8115)));
    double txk1=122.2/(exp((-(v+EKshift+127.2))/20.36)+exp((v+EKshift+236.8)/69.33));
    xk1 = xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
    double rk1=1.0/(1.0+exp((v+105.8-2.6*ko+EKshift)/(kslope_rk1*9.493)));
    double GK1 = 0.1908 * bGK1;
    if ( celltype == CellType::EPI )
      GK1 = GK1 * 1.2;
    else if ( celltype == CellType::MCELL )
      GK1 = GK1 * 1.3;
    
    IK1 = GK1*sqrt(ko)*rk1*xk1*(v-EK);
    
    // INaCa Current
    double bGncx = 2.4;
    double kna1=15.0;
    double kna2=5.0;
    double kna3=88.12;
    double kasymm=12.5;
    double wna=6.0e4;
    double wca=6.0e4;
    double wnaca=5.0e3;
    double kcaon=1.5e6;
    double kcaoff=5.0e3;
    double qna=0.5224;
    double qca=0.1670;
    double hca=exp((qca*v*F)/(R*T));
    double hna=exp((qna*v*F)/(R*T));
    double h1=1.0+nai/kna3*(1.0+hna);
    double h2=(nai*hna)/(kna3*h1);
    double h3=1.0/h1;
    double h4=1.0+nai/kna1*(1.0+nai/kna2);
    double h5=nai*nai/(h4*kna1*kna2);
    double h6=1.0/h4;
    double h7=1.0+nao/kna3*(1.0+1.0/hna);
    double h8=nao/(kna3*hna*h7);
    double h9=1.0/h7;
    double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
    double h11=nao*nao/(h10*kna1*kna2);
    double h12=1.0/h10;
    double k1=h12*cao*kcaon;
    double k2=kcaoff;
    double k3p=h9*wca;
    double k3pp=h8*wnaca;
    double k3=k3p+k3pp;
    double k4p=h3*wca/hca;
    double k4pp=h2*wnaca;
    double k4=k4p+k4pp;
    double k5=kcaoff;
    double k6=h6*cai*kcaon;
    double k7=h5*h2*wna;
    double k8=h8*h11*wna;
    double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
    double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
    double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
    double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
    double E1=x1/(x1+x2+x3+x4);
    double E2=x2/(x1+x2+x3+x4);
    double E3=x3/(x1+x2+x3+x4);
    double E4=x4/(x1+x2+x3+x4);
    double KmCaAct=150.0e-6;
    double allo=1.0/(1.0+std::pow((KmCaAct/cai), 2.0));
    double zna=1.0;
    double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    double JncxCa=E2*k2-E1*k1;
    double Gncx=0.0008 * bGncx;
    if ( celltype == CellType::EPI )
      Gncx = Gncx * 1.2;
    else if ( celltype == CellType::MCELL )
      Gncx = Gncx * 1.4;
    
    INaCa_i = 0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);
    
    h1 = 1.0+nass/kna3*(1.0+hna);
    h2 = (nass*hna)/(kna3*h1);
    h3 = 1.0/h1;
    h4 = 1.0+nass/kna1*(1.0+nass/kna2);
    h5 = nass*nass/(h4*kna1*kna2);
    h6 = 1.0/h4;
    h7 = 1.0+nao/kna3*(1.0+1.0/hna);
    h8 = nao/(kna3*hna*h7);
    h9 = 1.0/h7;
    h10 = kasymm+1.0+nao/kna1*(1.0+nao/kna2);
    h11 = nao*nao/(h10*kna1*kna2);
    h12 = 1.0/h10;
    k1   = h12*cao*kcaon;
    k2   = kcaoff;
    k3p  = h9*wca;
    k3pp = h8*wnaca;
    k3   = k3p+k3pp;
    k4p  = h3*wca/hca;
    k4pp = h2*wnaca;
    k4 = k4p+k4pp;
    k5 = kcaoff;
    k6 = h6*cass*kcaon;
    k7 = h5*h2*wna;
    k8 = h8*h11*wna;
    x1 = k2*k4*(k7+k6)+k5*k7*(k2+k3);
    x2 = k1*k7*(k4+k5)+k4*k6*(k1+k8);
    x3 = k1*k3*(k7+k6)+k8*k6*(k2+k3);
    x4 = k2*k8*(k4+k5)+k3*k5*(k1+k8);
    E1 = x1/(x1+x2+x3+x4);
    E2 = x2/(x1+x2+x3+x4);
    E3 = x3/(x1+x2+x3+x4);
    E4 = x4/(x1+x2+x3+x4);
    KmCaAct = 150.0e-6;
    allo = 1.0/(1.0+std::pow((KmCaAct/cass), 2.0));
    JncxNa = 3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    JncxCa = E2*k2-E1*k1;
    INaCa_ss = 0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);
    
    INaCa = INaCa_i+INaCa_ss;
    
    // INaK
    double bGnak = 2.0;
    double k1p=949.5;
    double k1m=182.4;
    double k2p=687.2;
    double k2m=39.4;
    k3p = 1899.0;
    double k3m=79300.0;
    k4p = 639.0;
    double k4m=40.0;
    double Knai0=9.073;
    double Knao0=27.78;
    double delta= -0.1550;
    double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
    double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
    double Kki=0.5;
    double Kko=0.3582;
    double MgADP=0.05;
    double MgATP=9.8;
    double Kmgatp=1.698e-7;
    double H=1.0e-7;
    double eP=4.2;
    double Khp=1.698e-7;
    double Knap=224.0;
    double Kxkur=292.0;
    double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
    double a1=(k1p*std::pow((nai/Knai), 3.0))/(std::pow((1.0+nai/Knai), 3.0)+std::pow((1.0+ki/Kki), 2.0)-1.0);
    double b1=k1m*MgADP;
    double a2=k2p;
    double b2=(k2m*std::pow((nao/Knao), 3.0))/(std::pow((1.0+nao/Knao), 3.0)+std::pow((1.0+ko/Kko), 2.0)-1.0);
    double a3=(k3p*std::pow((ko/Kko), 2.0))/(std::pow((1.0+nao/Knao), 3.0)+std::pow((1.0+ko/Kko), 2.0)-1.0);
    double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
    double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
    double b4=(k4m*std::pow((ki/Kki), 2.0))/(std::pow((1.0+nai/Knai), 3.0)+std::pow((1.0+ki/Kki), 2.0)-1.0);
    x1 = a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
    x2 = b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
    x3 = a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
    x4 = b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
    E1 = x1/(x1+x2+x3+x4);
    E2 = x2/(x1+x2+x3+x4);
    E3 = x3/(x1+x2+x3+x4);
    E4 = x4/(x1+x2+x3+x4);
    double zk=1.0;
    double JnakNa=3.0*(E1*a3-E2*b3);
    double JnakK=2.0*(E4*b1-E3*a1);
    double Pnak=30.0 * bGnak;
    if ( celltype == CellType::EPI )
      Pnak = Pnak * 0.9;
    else if ( celltype == CellType::MCELL )
      Pnak = Pnak * 0.7;
    
    INaK = Pnak*(zna*JnakNa+zk*JnakK);
    
    // ICl
    double IClCa_si = 0.0;
    double IClb_si = 0.0;
    double Cli = 15.0;   // Intracellular Cl  [mM]
    double Clo = 150.0;  // Extracellular Cl  [mM]
    double ecl = (R*T/F)*log(Cli/Clo); // [mV]
    double GClCa =0.5* 0.109625;   // [mS/uF]
    double GClB = 1.0*9.0e-3;        // [mS/uF]
    double KdClCa = 100.0e-3;    // [mM]
    IClCa = IClCa_si*GClCa/(1.0+KdClCa/cass)*(v-ecl);
    IClbk = IClb_si*GClB*(v-ecl);

    // IKb
    double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
    double GKb=0.003;
    if ( celltype == CellType::EPI ) GKb = GKb * 0.6;
    IKb = GKb*xkb*(v-EK);
    
    // INab
    double PNab=3.75e-10;
    INab = PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);
    
    // ICab
    double bPcab = 4.0;
    double PCab=2.5e-8*bPcab;
    ICab = (1.0-undo_ICab)*PCab*4.0*vffrt*(1.2*cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
    
    // IpCa
    double GpCa=0.0005;
    IpCa = GpCa*cai/(0.0005+cai);
  }


  void FBC( double TropTot_dt, double VNaCa, double Vuni ) {
    double CaMKb = CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
    CaMKa = ( CaMKb+CaMKt ) * CaMKa_scaling_factor;
    CaMKt = CaMKt + (dt*(aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt));

    // Diffusion Fluxes
    double bJdiff= 1.7;
    JdiffNa = (nass-nai)/2.0;
    JdiffK  = (kss-ki)/2.0;
    Jdiff   = (cass-cai)*bJdiff/0.2;

    //==================
    // Bartolucci et al.
    //==================
    // === RyRs CICR from SR ===
    // Jrel parameters
    double g_irel_max = coefficients[0]; // 20*10^-3; % millimolar_per_second (in calcium_dynamics)
    double RyRa1       = 0.05;;  // uM
    double RyRa2       = 0.03;;  // uM
    double RyRohalf    = 0.12-(RyRa1-RyRa2/2.0); // uM
    double RyRchalf    = 0.10-(RyRa1-RyRa2/2.0); // uM

    double fJrelp      = (1.0/(1.0+KmCaMK/CaMKa));
    double RyRSRCass   = (1.0 - 1.0/(1.0 + exp((casr-0.3)/0.1)));
    double RyRainfss   = RyRa1-RyRa2/(1.0 + exp((1000.0*cass-(0.043))/0.0082));
    double RyRtauadapt = 1000.0; // ms
    RyRa = RyRainfss - (RyRainfss - RyRa) * exp(-dt/RyRtauadapt);

    double RyRoinfss   = (1.0 - 1.0/(1.0 +  exp((1000.0*cass-(RyRa+ RyRohalf))/0.003)));
    double RyRtauact   = 18.75*std::pow(10.0, (-1.0))/1.875;       // ms
    RyRo = RyRoinfss - (RyRoinfss - RyRo) * exp(-dt/RyRtauact);

    double RyRcinfss   = (1.0/(1.0 + exp((1000.0*cass-(RyRa+RyRchalf))/0.001)));
    double RyRtauinact = 2.0*87.5*1.0/10.0;    // ms
    RyRc = RyRcinfss - (RyRcinfss - RyRc) * exp(-dt/RyRtauinact);

    if (celltype == CellType::MCELL ) {
      double g_irel_max_M = g_irel_max*1.7;
      Jrelnp = g_irel_max_M*RyRSRCass*RyRo*RyRc*(casr-cass);
    }
    else
      Jrelnp = g_irel_max*RyRSRCass*RyRo*RyRc*(casr-cass);

    double RyRtauinactp = RyRtauinact*1.25;
    RyRcp = RyRcinfss - (RyRcinfss - RyRcp) * exp(-dt/RyRtauinactp);

    double g_irel_max_p = g_irel_max*1.25;
    if (celltype == CellType::MCELL ) g_irel_max_p = g_irel_max*1.25*1.7;
    
    Jrelp = g_irel_max_p*RyRSRCass*RyRo*RyRcp*(casr-cass);
    Jrel = ((1.0-fJrelp)*Jrelnp+fJrelp*Jrelp);

    // == Ca2+ Uptake Flux ==
    // SERCA Jup parameter 
    double cJup = coefficients[1]; // 3.13
    double Jupnp = 0.004375*cai/(cai+0.00092);
    double Jupp  = 2.75*0.004375*cai/(cai+0.00092-0.00017);
    if (celltype == CellType::EPI ) {
      Jupnp = Jupnp*1.3;
      Jupp  = Jupp*1.3;
    }
    
    double fJupp = (1.0/(1.0+KmCaMK/CaMKa));
    Jleak = 0.0123*casr/15.0;
    Jup   = cJup*((1.0-fJupp)*Jupnp+fJupp*Jupp);
    if (is_HFpEF)
      Jup *= SERCA_scaling_factor;

    double Vmax_SRCaP  = 1.0*5.3114e-3; // [mM/msec] (286 umol/L cytosol/sec)
    double Kmf         = 0.246e-3;          // [mM] default
    double Kmr         = 1.7;               // [mM]L cytosol
    double hillSRCaP   = 1.787;       // [mM]
    Jup2 = Vmax_SRCaP*(std::pow((cai/Kmf), hillSRCaP)-std::pow((casr/Kmr), hillSRCaP))/(1.0+std::pow((cai/Kmf), hillSRCaP)+std::pow((casr/Kmr), hillSRCaP));

    // Tranlocation Flux (useful if considering JSR+NSR)
    Jtr = 0.0; //(cansr-cajsr)/100.0

    nai = nai + (dt*(-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo));
    nass = nass + (dt*(-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa));

    ki = ki + (dt*(-(Ito+IKr+IKs+IK1+IKb+Ist-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo));
    kss = kss + (dt*(-(ICaK)*Acap/(F*vss)-JdiffK));

    // Calcium Buffer Constants
    cmdnmax = 0.05;
    if (celltype == CellType::EPI ) cmdnmax = cmdnmax*1.2; // Modified  respect ORd (*1.3)
    
    double Bcai=1.0/(1.0+cmdnmax*kmcmdn/std::pow((kmcmdn+cai), 2.0));
    
    // Feedback from mitochondria
    double mito_cai_contribution = (VNaCa - Vuni)*vmito/vmyo;
    cai = cai + (dt*(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vsr/vmyo+Jleak*vsr/vmyo+Jdiff*vss/vmyo-(2.0*TropTot_dt/1000.0) + mito_cai_contribution)));

    double Bcass = 1.0/(1.0+BSRmax*KmBSR/std::pow((KmBSR+cass), 2.0)+BSLmax*KmBSL/std::pow((KmBSL+cass), 2.0));
    cass = cass + (dt*(Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vsr/vss-Jdiff)));

    double Bcasr = 1.0/(1.0+csqnmax*kmcsqn/std::pow((kmcsqn+casr), 2.0));
    casr = casr + (dt*(Bcasr*(Jup-Jleak-Jrel)));
    //=====END Bartolucci et al.============================
  }

// private:
  bool is_HFpEF        = false;
  CellType celltype    = CellType::ENDO;

  double num_beats     = 100.0;
  double CL            = 1000.0;//pacing cycle length 
  double ft            = num_beats*CL;//final time
  
  double amp           = -53.0; //-80//stimulus amplitude in uA/uF
  double start         = 0.0;   //start time of the stimulus, relative to each beat		
  double duration      = 0.5;   //duration of the stimulus in ms

  double maxEGTA       = 0.0;
  double undo_p        = 0.0; // Phosphorylation on/off
  double undo_CDI      = 0.0; // CDI on/off
  double undo_ICab     = 0.0; // ICab on/off

  // Introduce APD, timing, and counting parameters
  double dt            = 0.005;
  double APD_flag      = 0;
  double APD           = 0.0;
  double vo            = v;
  double t             = 0.0;
  double vdot          = 0.0;
  double n             = 0.0;
  double t_vdot_max    = 0.0;
  double vrest         = 0.0;
  double vdot_old      = 0.0;
  double vdot_max      = 0.0;
  
  // Initial values for state variables, there are 41 of them
  double v             = -87.0;
  double nai           = 7.0;
  double nass          = nai;
  double ki            = 145.0;
  double kss           = ki;
  double cai           = 1.0e-4;
  double cass          = cai;
  double cansr         = 1.2;
  double cajsr         = cansr;
  double casr          = 1.2;
  double m             = 0.0;
  double hf            = 1.0;
  double hs            = 1.0;
  double j             = 1.0;
  double hsp           = 1.0;
  double jp            = 1.0;
  double mL            = 0.0;
  double hL            = 1.0;
  double hLp           = 1.0;
  double a             = 0.0;
  double iF            = 1.0;
  double iS            = 1.0;
  double ap            = 0.0;
  double iFp           = 1.0;
  double iSp           = 1.0;
  double d             = 0.0;
  double ff            = 1.0;
  double fs            = 1.0;
  double fcaf          = 1.0;
  double fcas          = 1.0;
  double jca           = 1.0;
  double nca           = 0.0;
  double ffp           = 1.0;
  double fcafp         = 1.0;
  double xrf           = 0.0;
  double xrs           = 0.0;
  double xs1           = 0.0;
  double xs2           = 0.0;
  double xk1           = 1.0;
  double Jrelnp        = 0.0;
  double Jrelp         = 0.0;
  double CaMKt         = 0.0;
  double jnca          = 1.0;
  double alpha         = 0.0;
  double beta          = 0.0;
  double RyRa          = 0.03;
  double RyRo          = 0.0;
  double RyRc          = 1.0;
  double RyRcp         = 1.0;
  double Ta            = 0.0; // Active Force

  // Constants
  double nao                   = 144.0;//extracellular sodium in mM
  double cao                   = 1.8;//extracellular calcium in mM
  double ko                    = 5.4;//extracellular potassium in mM

  // Buffer paramaters
  double BSRmax                = 0.047;
  double KmBSR                 = 0.00087;
  double BSLmax                = 1.124;
  double KmBSL                 = 0.0087;
  double cmdnmax               = 0.05;
  double kmcmdn                = 0.00238;
  double trpnmax               = 0.07;
  double kmtrpn                = 0.0005;
  double csqnmax               = 10.0 * std::pow(10.0,(-1.0));
  double kmcsqn                = 0.8;

  // CaMK paramaters
  double aCaMK                 = 0.05;
  double bCaMK                 = 0.00068;
  double CaMKo                 = 0.05;
  double KmCaM                 = 0.0015;
  double KmCaMK                = 0.15;

  // Physical constants
  double R                     = 8314.0;
  double T                     = 310.0;
  double F                     = 96485.0;

  // Cell geometry
  double L                     = 0.01;
  double rad                   = 0.0011;
  double vcell                 = 1000.0*std::numbers::pi*rad*rad*L;
  double Ageo                  = 2.0*std::numbers::pi*rad*rad+2.0*std::numbers::pi*rad*L;
  double Acap                  = 2.0*Ageo;
  double vmyo                  = 0.68*vcell;
  double vmito                 = 0.26*vcell;
  double vnsr                  = 0.0552*vcell;
  double vjsr                  = 0.0048*vcell;
  double vsr                   = 0.95*(vnsr+vjsr); //0.06*vcell
  double vss                   = 0.02*vcell;

  // Introduce variables for reversal potentials, currents, fluxes, and CaMK
  double EKshift       = 8.0; // Shift for compensation in Ito,IKs,IK1 due to LJP
  double zca           = 2.0;
  double ENa           = 0.0;
  double EK            = 0.0;
  double EKs           = 0.0;
  double INa           = 0.0;
  double INaL          = 0.0;
  double Ito           = 0.0;
  double ICaL          = 0.0;
  double ICaNa         = 0.0;
  double ICaK          = 0.0;
  double IKr           = 0.0;
  double IKs           = 0.0;
  double IK1           = 0.0;
  double INaCa_i       = 0.0;
  double INaCa_ss      = 0.0;
  double INaCa         = 0.0;
  double INaK          = 0.0;
  double IKb           = 0.0;
  double INab          = 0.0;
  double IpCa          = 0.0;
  double ICab          = 0.0;
  double Ist           = 0.0;
  double Jrel          = 0.0;
  double Jup           = 0.0;
  double Jup2          = 0.0;
  double Jtr           = 0.0;
  double Jdiff         = 0.0;
  double JdiffNa       = 0.0;
  double JdiffK        = 0.0;
  double Jleak         = 0.0;
  double CaMKa         = 0.0;
  double CaMKb         = 0.0;
  double gICaL         = 0.0;
  double IClCa         = 0.0;
  double IClbk         = 0.0;

  double CCak          = 0.0;
  double I1Cak         = 0.0;
  double I2Cak         = 0.0;
  double Ck            = 1.0;
  double I1k           = 0.0;
  double I2k           = 0.0;
  double CCakp         = 0.0;
  double I1Cakp        = 0.0;
  double I2Cakp        = 0.0;
  double Ckp           = 1.0;
  double I1kp          = 0.0;
  double I2kp          = 0.0;
  double Ok            = 0.0;
  double Okp           = 0.0;
  double IC1_dutta     = 0.999637;
  double IC2_dutta     = 6.83208e-5;
  double C1_dutta      = 1.80145e-8;
  double C2_dutta      = 8.26619e-5;
  double O_dutta       = 0.00015551;
  double IO_dutta      = 5.67623e-5;
  double IOb_dutta     = 0.0;
  double Ob_dutta      = 0.0;
  double Cb_dutta      = 0.0;
  double D_dutta       = 0.0;
  double XS            = 0.0;
  double XW            = 0.0;
  double Ca_TRPN       = 0.0;
  double TmBlocked     = 1.0;
  double ZETAS         = 0.0;
  double ZETAW         = 0.0;
  double dCa_TRPN      = 0.0;

  //=============================================================
  // Initial values for ICaL Markov State Variables (ODEs: LSODA)
  //=============================================================
  static const int neq_electrophysiology = 26;
  int iwork1_electrophysiology           = 0;
  int iwork2_electrophysiology           = 0;
  int iwork5_electrophysiology           = 0;
  int iwork6_electrophysiology           = 0;
  int iwork7_electrophysiology           = 0;
  int iwork8_electrophysiology           = 0;
  int iwork9_electrophysiology           = 0;
  double rwork1_electrophysiology        = 0.0;
  double rwork5_electrophysiology        = 0.0;
  double rwork6_electrophysiology        = 0.0;
  double rwork7_electrophysiology        = 0.0;
  double y_electrophysiology[neq_electrophysiology+1]; 

  double t_begin_electrophysiology       = 0.0;
  double dt_electrophysiology            = dt;
  double tout_electrophysiology          = t_begin_electrophysiology + dt;
  int itol_electrophysiology             = 2;
  double rtol_electrophysiology[neq_electrophysiology+1];
  double atol_electrophysiology[neq_electrophysiology+1];
  int itask_electrophysiology            = 1;
  int istate_electrophysiology           = 1;
  int iopt_electrophysiology             = 0;
  int jt_electrophysiology               = 2;


  const double coefficients[2] = {0.022, 3.0};

  //======================
  // HFpEF Scaling Factors
  //======================
  double VOLUME_scaling_factor;
  double INAL_scaling_factor            = 2.34;
  double INA_scaling_factor             = 0.6;
  double SERCA_scaling_factor           = 1.0;
  double CaMKa_scaling_factor           = 1.0;
  
  LSODA<CardioMyocyte> lsoda_integrator_electrophysiology;
};
