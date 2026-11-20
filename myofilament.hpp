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

#include <iostream>
#include <array>
#include <algorithm>
#include <lsoda.hpp>

#define sign(a) ((a) < (0.0) ? (-1.0) : (1.0))
#define heav(a) ((a) < (0.0) ? (0.0)  : (1.0))

class Myofilament {
public:
  Myofilament( bool isHFpEF = false) : lsoda_integrator( this ), is_HFpEF(isHFpEF) { 
    if (isHFpEF)
      std::cout << " -> -> Constructing HFpEF Myofilament!!!" << std::endl;
    else
      std::cout << " -> -> Constructing Healthy Myofilament!!!" << std::endl;

    // Intiial values for myofilament ODEs
    y_myofilament[0]    = 0.0;        //Dummy
    y_myofilament[0+1]  = Nxb;        // Non-permissive fraction in overlap region
    y_myofilament[1+1]  = XBpreR;     // Fraction of XBs in pre-rotated state (Stiffness generation)
    y_myofilament[2+1]  = XBpostR;    // Fraction of XBs in post-rotated state (force generation)
    y_myofilament[3+1]  = x_XBpreR;   // Strain of XBs in the pre-rotated state
    y_myofilament[4+1]  = x_XBpostR;  // Strain of XBs in the post-rotated state
    y_myofilament[5+1]  = TropCaL;    // Ca bound to the low affinity troponin site
    y_myofilament[6+1]  = TropCaH;    // Ca bound to the high affinity troponin site
    y_myofilament[7+1]  = SL;         // Sarcomere length (um)
    y_myofilament[8+1]  = intf;       // Inegrated force

    rtol_myofilament[0] = 0.0;
    atol_myofilament[0] = 0.0;
    for ( std::size_t i = 1; i <= neq_myofilament; ++i ) {
      rtol_myofilament[i] = 1.0E-7;
      atol_myofilament[i] = 1.0E-6;
    }

    perm50        = 0.5;
    n_perm        = 15.0;
    if ( is_HFpEF ) {
      PCon_titin    = 0.002 * 2.0;
      gapp        = 70.0e-3 * 0.7;
    }
  }

  auto calculatePassiveForces( double sarcomere_length ) -> double {
    // Passive forces: Trabeculae: Titin and collagen, Single cells: Titin
    PF_titin = sign(sarcomere_length-SL_rest)*PCon_titin*std::exp(std::abs(sarcomere_length-SL_rest)*PExp_titin-1.0);
    
    if ( sarcomere_length >= SL_collagen )
      PF_collagen = PCon_collagen*std::exp(PExp_collagen*(sarcomere_length-SL_collagen)-1.0);
    else
      PF_collagen = 0;
    
    return PF_titin + PF_collagen; //Trabeculae
  }

  void fexMyofilament( double t, double *y, double *ydot, void *data ) {
    const double* pdata = static_cast<const double*>(data);

    double cai  = pdata[0];
    MgADPm      = pdata[1];
    MgATPm      = pdata[2];

    Nxb     	  = y[0];
    XBpreR  	  = y[1];
    XBpostR 	  = y[2];
    x_XBpreR	  = y[3];
    x_XBpostR	  = y[4];
    TropCaL 	  = y[5];
    TropCaH 	  = y[6];
    SL      	  = y[7];
    intf    	  = y[8];

    // Convert pH to mM
    Hm          = pow(10.0,(-pHm))*1.0e3; // H concentration in mM
    
    // Ca binding to troponin C
    konT        = kon*1.0*pow(Q_kon,((TmpC-37.0)/10.0));
    konT_app    = konT*(pow(kdHCa,m_atp) + pow(Href,m_atp))/(pow(kdHCa,m_atp) + pow(Hm,m_atp));
    
    koffLT      = koffL*1.0*1.0*pow(Q_koff,((TmpC-37.0)/10.0)); // Intact as opposed to skinned
    koffHT      = koffH*1.0*1.0*pow(Q_koff,((TmpC-37.0)/10.0));
    
    ydot[5]/*d_TropCaL*/ = konT_app*cai*1000.0*(1.0-y[5]) - koffLT*y[5];
    ydot[6]/*d_TropCaH*/ = konT_app*cai*1000.0*(1.0-y[6]) - koffHT*y[6];
    
    // *****************************************************
    // Thin filament activation rates
    
    // Sarcomere geometry
    sovr_ze  = std::min(L_thick/2.0, SL/2.0);
    sovr_cle = std::max( SL/2.0 - (SL-L_thin),L_hbare/2.0 );
    L_sovr   = sovr_ze - sovr_cle; // Length of single overlap region
    
    // Overlap fraction for thick filament
    sov_thick   = L_sovr*2.0/(L_thick - L_hbare);
    // Overlap fraction for thin filament
    sov_thin    = L_sovr/L_thin;
    
    TropReg     = (1.0-sov_thin)*TropCaL + sov_thin*TropCaH;
    permtot     = pow( (1.0/(1.0+pow((perm50/TropReg),n_perm))),0.5 );
    permtot_p_n = std::min(100.0,(1.0/permtot));
    
    // Rate constants governing the transition btw Permissive and Non
    kn_pT       = kn_p*permtot*pow(Q_kn_p,((TmpC-37.0)/10.0));
    kp_nT       = kp_n*permtot_p_n*pow(Q_kp_n,((TmpC-37.0)/10.0));
    
    
    // *****************************************************
    // Cross-bridge cycling rates
    // Pxb to XBpreR
    fappT       = fapp*Xsp*pow(Q_fapp,((TmpC-37.0)/10.0));
    
    // Pi-dependent transition XBpreR to Pxb
    gappslmd    = 1.0 + (1.0-sov_thick)*gslmod;
    gappT       = gapp*gappslmd*Xsp*pow(Q_gapp,((TmpC-37.0)/10.0));
    gappT_true  = gappT/Pi_ref;  // True first order rate constant
    
    // XBpreR to XBpostR
    hfmd        = std::exp(-sign(x_XBpreR)*hfmdc*pow((x_XBpreR/x0),2.0));
    hfT         = hf*hfmd*Xsp*pow(Q_hf,((TmpC-37.0)/10.0));
    
    // H-dependent transition XBpostR to XBpreR
    hbT         = hb*Xsp*pow(Q_hb,((TmpC-37.0)/10.0));
    hbT_true    = hbT/Href;
    hbT_app     = ((kdADP+MgADP_ref)/MgADP_ref)*(MgADPm/(kdADP+MgADPm))*hbT_true;
    
    // MgATP-dependent transition from XBpostR to Pxb
    if (x_XBpostR < x0) 
      gxbmd = exp(sigma_p*pow(((x0-x_XBpostR)/x0),2));
    else
      gxbmd = exp(sigma_n*pow( ((x0-x_XBpostR)/x0),2 ));
    
    gxbT        = gxb*std::max(gxbmd,1.0)*Xsp*std::pow(Q_gxb,((TmpC-37)/10));
    gxbT_true   = gxbT/MgATP_ref;
    gxbT_app    = ((kdADP+MgADP_ref)/(kdADP+MgADPm))*gxbT_true;
    
    // Pxb to XBpostR - Introduced for thermodynamic efficiency
    G0          = -29600.0 + log(10.0)*R_myo*(TmpC+273.0)*(-log10(1.0e-7));
    K_MgATP     = exp(-G0/(R_myo*(TmpC+273.0)))*1.0e6;  // 1e6 to convert from M2 to mM2
    fxbT        = (kdADP*fappT*hfT*gxbT_true)/(gappT_true*hbT_true*K_MgATP);
    fxb         = (kdADP*fapp*hf*(gxb/MgATP_ref))/((gapp/Pi_ref)*(hb/Href)*K_MgATP); // Used for calculating max
        
    ap1         = fappT;
    ap2         = hfT;
    ap3         = MgATPm*gxbT_app;
    am1         = Pi*gappT_true;
    am2         = Hm*hbT_app;
    am3         = fxbT;
    
    // *****************************************************
    ////  Cross_bridge transitions
    
    // Update all RUs
    Pxb                   = 1.0 - y[0] - y[1] - y[2];
    ydot[0]/*d_Nxb*/      = -kn_pT*y[0] + kp_nT*Pxb;
    ydot[1]/*d_XBpreR*/   = ap1*Pxb - am1*y[1] - ap2*y[1] + am2*y[2];
    ydot[2]/*d_XBpostR*/  = ap2*y[1] + am3*Pxb - am2*y[2] - ap3*y[2];
    
    ATPase      = ap3*XBpostR - am3*Pxb;
    // *****************************************************
    // Mean strain of strongly-bound states
    
    // Get velocity
    dSL         = 0.0; //Isometric
    
    // Steady State occupancies of the bound states - Duty fractions
    Sum_duty    = am3*am2 + ap3*ap1 + am2*ap1 + ap1*ap2 + am3*am1 + am3*ap2 + ap2*ap3 + am3*am1 + ap3*am1;
    
    dutyPreR    = (am3*am2 + ap3*ap1 + am2*ap1)/Sum_duty;
    dutyPostR   = (ap1*ap2 + am3*am1 + am3*ap2)/Sum_duty;
    
    // Rate of change of the average distortions
    ydot[3]/*d_x_XBpreR*/  = dSL/2.0 + (Psi/dutyPreR)*(-ap1*y[3]) + am2*(y[4]-x0-y[3]);
    ydot[4]/*d_x_XBpostR*/ = dSL/2.0 + (Psi/dutyPostR)*(ap2*(y[3] + x0 - y[4]));
    
    // Compute the active force
    F_active    = kxb*sov_thick*(x_XBpreR*XBpreR + x_XBpostR*XBpostR);                       // /*d_x_XBpostR*/
    
    // Maximal state occupancies under optimal conditions
    max_XBpreR  = (fapp*(hb+gxb)+fxb*fapp)/(gxb*hf + fapp*hf + gapp*hb + gapp*gxb + fapp*hb + fapp*gxb + fxb*(hb+gapp+hb));
    max_XBpostR = (fapp*hf+fxb*(gapp+hb))/(gxb*hf + fapp*hf + gapp*hb + gapp*gxb + fapp*hb + fapp*gxb + fxb*(hb+gapp+hb));
    
    // Factor to normalise force
    FnormD      = kxb*x0*max_XBpostR;
    Fnorm       = F_active/FnormD;
    
    // Get passive forces
    F_preload   = sign(SLset-SL_rest)*PCon_titin*(exp(PExp_titin*abs(SLset-SL_rest))-1.0);//CalculatePassiveForces(1.95)
    F_passive   = calculatePassiveForces(SL); //CalculatePassiveForces(SL)
    F_afterload = 0.0;// Either constant or due to series elastic element
    if (SEon)
      F_afterload = KSE*(SLset-SL);
    
    //  Compute the integral of forces
    ydot[8]/*dintf*/ = (-F_passive+F_preload-Fnorm+F_afterload);     // Total force
    ydot[7]/*dSL*/   = contrflag*((y[8]+(SLset-y[7])*viscosity)/mass)*heav(y[7]-SL_min)*heav(SL_max-y[7]);

    if ( is_HFpEF )
      F_total  = Fnorm + F_passive;
    else
      F_total          = Fnorm + F_passive;
    
    // Ca buffering by low-affinity troponin C (LTRPNCa)
    // Compute fraction of strongly bound for compute effective binding of Ca
    FrSBXB       = (y[2]+y[1])/(max_XBpostR + max_XBpreR);
    FrSBXB_dt    = (ydot[2]+ydot[1])/(max_XBpostR + max_XBpreR);
    
    // Compute derivative of z-line end of single overlap region
    sovr_ze_dt   = -ydot[7]/2*heav(L_thick-y[7]);
    
    // Compute derivative of centerline of end of single overlap region
    sovr_cle_dt  = -ydot[7]/2*heav((2.0*L_thin-y[7])-L_hbare);
    
    // Compute the derivative of the length of the single overlap
    len_sovr_dt  = sovr_ze_dt-sovr_cle_dt;
    
    // Compute the derivative of the overlap fraction for thin filament
    SOVFThin_dt  = len_sovr_dt/L_thin;
    
    // Compute the derivative of the overlap fraction for thick filament
    SOVFThick_dt = 2.0*len_sovr_dt/(L_thick-L_hbare);
    
    TropTot      = Trop_conc*((1-sov_thin)*y[5] + sov_thin*(FrSBXB*y[6]+(1-FrSBXB)*y[5]));
    TropTot_dt   = Trop_conc*(-SOVFThin_dt*y[5]+(1-sov_thin)*ydot[5] + SOVFThin_dt*(FrSBXB*y[6]+(1-FrSBXB)*y[5]) + sov_thin*(FrSBXB_dt*y[6]+FrSBXB*ydot[6]-FrSBXB_dt*y[5] + (1-FrSBXB)*ydot[5]));
  }
  
// private:
  // Set temperature and SL
  double TmpC          = 37.0;

  // Metabolite reference concentrations
  double Pi_ref        = 2.0;   //mM
  double Pi            = 2.0;   //mM
  double MgATP_ref     = 5.0;   //mM
  double MgADP_ref     = 0.036; //mM

  // MgADP dissociation constant
  double kdADP         = 0.004; // mM

  // Competitive H binding parameters
  double kdHCa         = 1.0e-5; // mM
  double m_atp         = 1.0;
  double pH_ref        = 7.15;
  double Href          = pow( 10.0,(-pH_ref) )*1.0e3;

  // Gas constant
  double R_myo         = 8.314; //J/K/mol

  // Sarcomere geometry (um)
  double SL_max      = 2.4;   // (um) minimum sarcomere length
  double SL_min      = 1.4;   // (um) maximum sarcomere length
  double L_thick     = 1.65;  // (um) thick filament length
  double L_hbare     = 0.1;   // (um) length of bare portion of thick filament
  double L_thin      = 1.2;   // (um) thin filament length

  // Temperature dependence
  double Q_kon       = 1.5;
  double Q_koff      = 1.3;
  double Q_kn_p      = 1.6;
  double Q_kp_n      = 1.6;
  double Q_fapp      = 6.25;
  double Q_gapp      = 2.5;
  double Q_hf        = 6.25;
  double Q_hb        = 6.25;
  double Q_gxb       = 6.25;

  // Species parameter
  double Xsp         = 0.2;         // 1 for rat, 0.2 for Guinea Pig, Human?
      
  // Ca binding to troponin and thin filament
  double kon         = 50.0e-3;     // (uM-1 ms-1)
  double koffL       = 250.0e-3;    // (ms-1)
  double koffH       = 25.0e-3;     // (ms-1)
  double perm50      = 0.5;         //  perm variable that controls n to p transition
  double n_perm      = 15.0;        //  in Hill-like fashion
    
  double kn_p        = 500.0e-3;    // (ms-1)
  double kp_n        = 50.0e-3;     // (ms-1)

  // Thin filament regulation and XB cycling
  double fapp        = 500.0e-3;    // (ms-1) XB on rate
  double gapp        = 70.0e-3;     // (ms-1) XB off rate
  double gslmod      = 6.0;         // Controls SL effect on gapp
  double hf          = 2000.0e-3;   // (ms-1) Rate between pre-force and force states
  double hfmdc       = 5.0;         //
  double hb          = 400.0e-3;    // (ms-1) Rate between pre-force and force states
  double gxb         = 70.0e-3;     // (ms-1) ATP consuming transition rate
  double sigma_p     = 8.0;         // Distortion dependence of STP using transition gxb
  double sigma_n     = 1.0;         //

  // Mean strain of strongly-bond states
  double x0          = 0.007;       // (um) Strain induced by head rotation
  double Psi         = 2.0;         // Scaling factor balancing SL motion and XB cycling

  // Normalised active and passive force
  double SL_rest       = 1.9+0.3;   // (um) Rest SL length for 0 passive force
  double PCon_titin    = 0.002;     // (Normalised Force) Passive force due to titin
  double PExp_titin    = 10.0;      // (um-1) This applies to trabeculae and single cells only
  double SL_collagen   = 2.25;      // (uM) Resting length for collagen
  double PCon_collagen = 0.02;      // (normalised force) Passive force due to collagen
  double PExp_collagen = 70.0;      // (um-1) This applies to trabeculae and single cells only

  // Calculation of compdoublee muscle response
  double mass          = 0.00025;   // ([norm Force ms^2]/um) Muscle mass
  double viscosity     = 0.003e3;   // ([normalised Force ms]/um) Muscle viscosity

  double KSE           = 1.0;       // 50 // (normalised Force/um) Series elastic element
  double kxb           = 120.0;     // (mNewtons per mm^2) Maximal force
  double Trop_conc     = 70.0;      // (uM) Troponin concentration

  // Set metabolite concentrations, pH levels and SL
  double MgATPm        = 5.0;       // mM
  double MgADPm        = 36.0e-3;   // mM
  double pHm           = 7.15;

  // Time-varying parameters
  double contrflag     = 1.0;       // 0 for isometric contraction 1 for active contraction
  double SLset         = 2.2+0.01;  // Initial SL or SL for isometric contractions
  double SEon          = true;      // 0 w/o series elastic element, 1 w/element

  //===========================================
  // Initial values for myofilament ODEs: LSODA
  //===========================================
  static const int neq_myofilament  = 9;
  int iwork1_myofilament            = 0;
  int iwork2_myofilament            = 0;
  int iwork5_myofilament            = 0;
  int iwork6_myofilament            = 0;
  int iwork7_myofilament            = 0;
  int iwork8_myofilament            = 0;
  int iwork9_myofilament            = 0;
  double rwork1_myofilament         = 0.0;
  double rwork5_myofilament         = 0.0;
  double rwork6_myofilament         = 0.0;
  double rwork7_myofilament         = 0.0;
  double y_myofilament[neq_myofilament+1]; 

  double Nxb      = 0.9997;    
  double XBpreR   = 0.0001;    
  double XBpostR  = 0.0001;    
  double x_XBpreR = 0.0;         
  double x_XBpostR= 0.007;     
  double TropCaL  = 0.01447254;
  double TropCaH  = 0.1276;    
  double SL       = 2.2;       
  double intf     = 0.0;       

  double t_begin_myofilament = 0.0;
  double dt_myofilament      = 0.005;
  double tout_myofilament    = t_begin_myofilament + dt_myofilament;
  int itol_myofilament       = 2;
  double rtol_myofilament[neq_myofilament+1];
  double atol_myofilament[neq_myofilament+1];
  int itask_myofilament      = 1;
  int istate_myofilament     = 1;
  int iopt_myofilament       = 0;
  int jt_myofilament         = 2;

  // Set metabolite concentrations, pH levels and SL
  double Hm              = 0.0;  
  double konT            = 0.0;
  double konT_app        = 0.0;  
  double koffLT          = 0.0;
  double koffHT          = 0.0;  
  double sovr_ze         = 0.0;  
  double sovr_cle        = 0.0;  
  double L_sovr          = 0.0;
  double sov_thick       = 0.0;  
  double sov_thin        = 0.0;  
  double TropReg         = 0.0;
  double permtot         = 0.0;
  double permtot_p_n     = 0.0;
  double kn_pT           = 0.0;
  double kp_nT           = 0.0;
  double fappT           = 0.0;
  double fxb             = 0.0;
  double fxbT            = 0.0;
  double gappslmd        = 0.0;
  double gappT_true      = 0.0;
  double gappT           = 0.0;
  double hfmd            = 0.0;
  double hfT             = 0.0;
  double hbT             = 0.0;
  double hbT_true        = 0.0;
  double gxbmd           = 0.0;
  double gxbT            = 0.0;
  double gxbT_true       = 0.0;
  double gxbT_app        = 0.0;
  double G0              = 0.0;
  double K_MgATP         = 0.0;
  double ap1             = 0.0;
  double ap2             = 0.0;
  double ap3             = 0.0;
  double am1             = 0.0;
  double am2             = 0.0;
  double am3             = 0.0;
  double hbT_app         = 0.0;
  double Pxb             = 0.0; 
  double ATPase          = 0.0; 
  double dSL             = 0.0; 
  double Sum_duty        = 0.0;     
  double dutyPreR        = 0.0;   
  double dutyPostR       = 0.0;     
  double F_active        = 0.0;   
  double max_XBpreR      = 0.0;       
  double max_XBpostR     = 0.0;     
  double FnormD          = 0.0;   
  double Fnorm           = 0.0;   
  double F_preload       = 0.0;   
  double F_passive       = 0.0;   
  double F_total         = 0.0;   
  double F_afterload     = 0.0;     
  double PF_titin        = 0.0;     
  double PF_collagen     = 0.0;      
  double FrSBXB          = 0.0;   
  double FrSBXB_dt       = 0.0;   
  double sovr_ze_dt      = 0.0;     
  double sovr_cle_dt     = 0.0;     
  double len_sovr_dt     = 0.0;       
  double SOVFThin_dt     = 0.0;     
  double SOVFThick_dt    = 0.0;     
  double TropTot         = 0.0;   
  double TropTot_dt      = 0.0; 
  double mKTR            = 0.0; 

  bool is_HFpEF          = false;
  LSODA<Myofilament> lsoda_integrator;
};
