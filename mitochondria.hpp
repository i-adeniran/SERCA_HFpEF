/*
 * This program implements an electromechanical cardiac myocyte by coupling 
 * the O'Hara-Virag-Varro-Rudy 2011 Human ventricular cardiac cell model and the 
 * extended Rice myofilament mechanics model (exended with ATP energetics by Tran et al. 2010).
 *
 * O'Hara-Virag-Varro-Rudy Ventricular Cell Model
 * ==============================================
 * O'Hara T, Virág L, Varró A, Rudy Y., “Simulation of the undiseased human cardiac
 * ventricular action potential: model formulation and experimental validation”.
 * PLoS Comput Biol. 2011 May;7(5):e1002061. doi: 10.1371/journal.pcbi.1002061.
 *
 * Rice Myofilament Mechanics Model
 * ================================
 * J. Rice, F. Wang, D. Bers, and P. Tombe, “Aproximate model of cooperative 
 * activation and crossbridge cycling in cardiac muscle using ordinary differential 
 * equations”. Biophysical Journal, vol. 95, pp. 2368– 2390, 2008.
 *
 * Tran et al. Extended Model with ATP Energetics
 * ==============================================
 * Tran K, Smith NP, Loiselle DS, Crampin EJ., "A metabolite-sensitive, thermodynamically 
 * constrained model of cardiac cross-bridge cycling: implications for force development 
 * during ischemia". Biophys J. 2010 Jan 20;98(2):267-76. doi: 10.1016/j.bpj.2009.10.011.
 *
 * Compile with: g++ -std=c++23 -O3 -m64 -o ordcell myofilament_ode.cpp mitochondria_ode.cpp main_ordcell.cpp
 */

#include <numbers>
#include <lsoda.hpp>

class Mitochondria {
public:
  Mitochondria( double cell_volume ) : lsoda_integrator( this ) { 
    std::cout << " -> -> Constructing Mitochondria!!!" << std::endl; 

    mito_L                        = 0.01;
    mito_rad                      = 0.0011;
    mito_vcell                    = cell_volume;
    mito_Ageo                     = 2.0*std::numbers::pi*mito_rad*mito_rad+2.0*std::numbers::pi*mito_rad*mito_L;
    Acap                          = 2.0*mito_Ageo;
    vmito                         = 0.26*mito_vcell;
    Vmyo                          = 0.68*mito_vcell;
    
    // Intiial values for mitochondria ODEs
    y_mitochondria[0]                         = 0.0;
    y_mitochondria[1]                         = ATPi_cyto;
    y_mitochondria[2]                         = ATPi;      
    y_mitochondria[3]                         = CrPi_cyto; 
    y_mitochondria[4]                         = CrPi_mito; 
    y_mitochondria[5]                         = Cam;       
    y_mitochondria[6]                         = ADPm;      
    y_mitochondria[7]                         = Dpsi;     
    y_mitochondria[8]                         = NADH;      
    y_mitochondria[9]                         = ISOC;      
    y_mitochondria[10]                        = AKG;       
    y_mitochondria[11]                        = SCoA;      
    y_mitochondria[12]                        = Succ;      
    y_mitochondria[13]                        = FUM;       
    y_mitochondria[14]                        = MAL;       
    y_mitochondria[15]                        = Oaa;       
    y_mitochondria[16]                        = ASP; 
    
    rtol_mitochondria[0] = 0.0;
    atol_mitochondria[0] = 0.0;
    for ( std::size_t i = 1; i <= neq_mitochondria; ++i ) {
      rtol_mitochondria[i] = 1.0E-7;//1.0E-10;
      atol_mitochondria[i] = 1.0E-6;
    }
  }

  void fexMitochondria( double t, double *y, double *ydot, void *data ) {
    const double* pdata = static_cast<const double*>(data);

    double Jup  = pdata[0];
    double INaK = pdata[1];
    double IpCa = pdata[2];

    ATPi_cyto       = y[0];
    ATPi            = y[1];
    CrPi_cyto       = y[2];
    CrPi_mito       = y[3];
    Cam             = y[4];
    ADPm            = y[5];
    Dpsi            = y[6];
    NADH            = y[7];
    ISOC            = y[8];
    AKG             = y[9]; 
    SCoA            = y[10];
    Succ            = y[11];  
    FUM             = y[12];
    MAL             = y[13];
    Oaa             = y[14];  
    ASP             = y[15];

    //==============
    // CalculateCK()
    //==============
    Vt_CRP2  = kt_2 * ( y[3] - y[2] );
    VCK_cyto = kf_2 * ( ( CRT_cyto - y[2] ) * y[0] -  y[2] * ( 8.0 - y[0] ) * inv_keq);
    VCK_mito = kf_3 * ( ( CRT_mito - y[3] ) * y[1] - y[3] * ADP * inv_keq);
    
    ydot[3]/*dCrPi_mito*/	= VCK_mito - Vt_CRP2;
    ydot[2]/*dCrPi_cyto*/	= Vt_CRP2  + VCK_cyto;
    ydot[1]/*dATPi*/		  = 0.615 * VANT - V_AM - 0.5*Jup - ( 6.371e-5 * ( INaK + IpCa ) ) - VCK_mito; // Defined normally in CalculateFATPi
    ydot[0]/*dATPi_cyto*/	= -VCK_cyto - VATPase_cyto;
    //==============================================================>>> END: CalculateCK();
      
    //=================
    // CalculateFATPi()
    //=================
    //ydot[1]/*dATPi*/ = 0.615 * VANT - V_AM - 0.5* Jup - ( 6.371e-5 * ( INaK + IpCa ) );
      
    //=========================
    // CalculateF_mitene()
    //=========================
    ydot[4]/*dCam*/  = fm * (Vuni - VnaCa);
    ydot[5]/*dADPm*/ = VANT - VATPase - VSL;
    ydot[6]/*dDpsi*/ = -(-VHNe - VHFe + Vhu + VANT + Vhleak + two_b * VnaCa + 2.0 * Vuni ) * inv_Cmito;
    ydot[7]/*dNADH*/ = -VNO + VIDH + VKGDH + VMDH;
    ydot[8]/*dISOC*/ = VACO - VIDH;
    ydot[9]/*dAKG*/  = VIDH + VAAT - VKGDH;
    ydot[10]/*dSCoA*/= VKGDH - VSL;
    ydot[11]/*dSucc*/= VSL - VSDH;
    ydot[12]/*dFUM*/ = VSDH - VFH;
    ydot[13]/*dMAL*/ = VFH - VMDH;
    ydot[14]/*dOaa*/ = VMDH - VCS - VAAT;
    ydot[15]/*dASP*/ = VAAT - kcnsASP * y[15];
  }

  void calculateConstantModelParameters( double vmyo, double vsr, double vnsr, double vjsr, double vss, double cao, double Ko, double nao ) {
    Vtotal                        = (vmyo+vsr+vnsr+vjsr+vss)/0.74;
    Vmito                         = vmito;                              // ***Vtotal * 0.36	// assuming Vmito is 36% of the cellular volume

    // Optimization Precalcs
    Vmyo_1000_F_Acap_C_m          = Vmyo * ( Faraday * 1000.0) / (mAcap * C_m);
    Co                            = Vtotal * ( 2.0 * cao + Ko + nao) / Vmyo;
    RTF_05                        = RT_over_F * 0.5;
    RTF_log_Nao                   = RT_over_F * log(nao);
    RTF_log_Ko                    = RT_over_F * log(Ko);
    RTF_log_Nao_Ko                = RT_over_F * log(Ko + 0.01833 * nao);
    RTF_05_log_Cao                = RTF_05 * log(cao);
    FaradayE3                     = (1000.0*Faraday);
    F_over_RT                     = 1.0 / RT_over_F;
    KmNaiP1p5                     = sqrt(KmNai*KmNai*KmNai); // KmNai^1.5
    inv_KiADP_CaP                 = 1.0 / KiADP_CaP;
    KmATP_AM_Ki_AM                = KmATP_AM / Ki_AM;

    // Mitochondria Constants : ATPm
    DmuH_Constant                 = -2.303 * RT_over_F * DpH;
    VCS_C1                        = (KCS * EtCS * AcCoA) / (KmAcCoA + AcCoA);
    one_inv_KACOeq                = 1.0 + 1.0 / KACOeq;
    VIDH_Constant                 = 1.0 + mH / kh_1 + kh_2 / mH;
    kIDH_EtID                     = kIDH * EtID;
    inv_KADP                      = 1.0 / KADP;
    inv_KaCa                      = 1.0 / KaCa;
    inv_KidhNADH                  = 1.0 / KidhNADH;
    KmKGNAD_KmIDNAD               = KmKGNAD / KmIDNAD;
    Mg_Kmg_1                      = Mg / Kmg + 1.0;
    Mg_Kmg_1_Kca                  = Mg_Kmg_1 / Kca;
    kKGDH_EtKG                    = kKGDH * EtKG;
    CoA_KSLeq                     = CoA / KSLeq;
    kSDH_EtSDH                    = kSDH * EtSDH;
    KmSucc_KiFUM                  = KmSucc / KiFUM;
    inv_KiOxaa                    = 1.0 / KiOxaa;
    kfFH_KFHeq                    = kfFH / KFHeq;
    kMDH_Fh_EtMD                  = pow( ( 1.0 / ( 1.0 + Kh3 / mH + Kh3 * Kh4 / pow(mH,2) ) ) ,2) * ( 1.0 / ( 1.0 + mH / Kh1 + pow(mH,2) / ( Kh1 * Kh2 ) ) + Koff) * kMDH * EtMD;;
    Kmal_Kioaa                    = Kmal / Kioaa;
    VAAT_Constant                 = kfAAT * GLU * kcnsASP * KAATeq / kfAAT;
    kcnsASP_KAATeq_kfAAT          = kcnsASP * KAATeq / kfAAT;
    KfAAT_GLU                     = kfAAT * GLU;
    KfAAT_KAATeq                  = kfAAT/KAATeq;
    kres_sq_KmIDNAD               = kres * kres / KmIDNAD;
    exp_6_FRT_Dpsio               = exp( 6.0 * Dpsio * F_over_RT);
    FRT_6_g                       = 6.0 * mmg * F_over_RT;
    ra_rc1_exp_6_FRT_Dpsio        = ra + rc1 * exp_6_FRT_Dpsio;
    r1_exp_6_FRT_Dpsio            = r1 * exp_6_FRT_Dpsio;
    rhoREN_ra_rc1_exp_6_FRT_Dpsio = 0.5 * rhoREN * ra_rc1_exp_6_FRT_Dpsio;
    rhoREN_rc2                    = 0.5 * rhoREN * rc2;
    rhoREN_ra                     = 0.5 * rhoREN * ra;
    rhoRen_6_ra                   = 6.0 * rhoREN * ra;
    rhoRen_6_ra_rb                = 6.0 * rhoREN * ( ra + rb );
    AREF                          = RT_over_F * log( kresf * sqrt( FADH2 / FAD ) );
    exp_AREF_FRT                  = exp( AREF * F_over_RT );
    ra_rc2_exp_AREF_FRT           = 0.5 * (ra + rc2 * exp_AREF_FRT);
    VFO_C1                        = ( ra + rc1 * exp_6_FRT_Dpsio ) * exp_AREF_FRT * 0.5;
    ra_exp_AREF_FRT               = 4.0 * ra * exp_AREF_FRT;
    ra_rb                         = 4.0 * (ra + rb);
    VFO_VHFe_C1                   = ( 1.0 + r1 * exp_AREF_FRT) * exp_6_FRT_Dpsio;
    r2_r3_exp_AREF_FRT            = r2 + r3 * exp_AREF_FRT;
    exp_3_FRT_Dpsio               = exp( 3.0 * Dpsio * F_over_RT);
    FRT_3                         = 3.0 * F_over_RT;
    kf1_Pi                        = kf1 / Pi;
    VATPase_C1                    = ( 100.0 * pa + pc1 * exp_3_FRT_Dpsio);
    pa_pb_3                       = 3.0 * (pa + pb);
    pa_300                        = 300.0 * pa;
    p1_exp_3_FRT_Dpsio            = p1 * exp_3_FRT_Dpsio;
    hm_F_over_RT                  = hm * F_over_RT;
    VmDT_75                       = 0.75 * VmDT;
    VmDT_20                       = 20.0 * VmDT;
    tenDiv9                       = 10.0 / 9.0;
    // Pause doing Mitochondria : VHLeak

    KmATP_SR_Ki_SR                = KmATP_SR / Ki_SR;
    inv_Ki_prime_SR               = 1.0 / Ki_prime_SR;

    // Start Mitochondria Again : Vuni
    inv_ktrans                    = 1.0 / ktrans;
    inv_kact                      = 1.0 / kact;
    Vmuni_ktrans                  = Vmuni / ktrans;
    FRT2                          = 2.0 * F_over_RT;
    b_05                          = b * 0.5;
    // Pause Mitochondria Again : Vanca

    inv_7p5                       = 1.0 / 7.5;
    inv_6                         = 1.0 / 6.0;
    inv_9p5                       = 1.0 / 9.5;

    // Mitochondria Last Time : F_Mitene
    inv_Cmito                     = 1.0 / Cmito;
    two_b                         = 2.0 * b;

    // CK Mod
    inv_keq                       = 1.0 / keq;
  }

  // V_AM is the rate of ATP hydrolysis by myofibrils (AM ATPase)
  void calculateV_AM( double fappT, double Pxb ) {
    // Parameters: V_AM_scaler, V_AM_max, f_01, f_12, f_23, KmATP_AM, Ki_AM
    V_AM = V_AM_max * (fappT * Pxb) / (1.0 + inv_ATPi * ( KmATP_AM + KmATP_AM_Ki_AM * ADP ));
  }

  void calculateATPm() {
    ATPm = mCm - (ADPm);
  }	
	
	void calculateDmuH() {
		DmuH = DmuH_Constant + (Dpsi);
	}

	void calculateNAD() {
		NAD = CPN - (NADH);
	}

	void calculateVCS() {
		//Parameter: KmOaa, KCS, EtCS, AcCoA, KmAcCoA
		VCS = VCS_C1 * (Oaa) / ((Oaa) +  KmOaa);
	}

	void calculateVACO() {
		VACO = kfACO * ( CIK - (AKG) - (SCoA) - (Succ) - (FUM) - (MAL) - (Oaa) - (ISOC)* one_inv_KACOeq );
	}

	void calculateVIDH() {
		double Fa = 1.0 / (( 1.0 + (ADPm) * inv_KADP) * (1.0 + (Cam) * inv_KaCa));
		double Fi = 1.0 + (NADH) * inv_KidhNADH;
		VIDH      = kIDH_EtID / (VIDH_Constant + KmIDNAD_NAD * Fi + std::pow( Kmiso/(ISOC), nID) * Fa * (1.0 + KmIDNAD_NAD  * Fi));
	}

	void calculateVKGDH() {
		double mma = ( (Mg_Kmg_1 + Mg_Kmg_1_Kca*(Cam)) );
		VKGDH      = kKGDH_EtKG * mma / ( mma + pow(KmKG/(AKG),nKG) + KmKGNAD_KmIDNAD * KmIDNAD_NAD );
	}

	//Mitochondia;Also CIT; Intemediate Variable;  Dependent on calculateATPm
	void calculateVSL() {
		VSL = kfSL * ( (SCoA) * (ADPm) - CoA_KSLeq * (Succ) * ATPm);
	}

	//Mitochondia; Intemediate Variable;  Optimize Things that use this;1/Succ
	void calculateVSDH() {
		//Parameters: kSDH, EtSDH, KmSucc, KiFUM, KiOxaa
		VSDH = kSDH_EtSDH * (Succ) / ((Succ) + (KmSucc + KmSucc_KiFUM*(FUM)) * (1.0 + inv_KiOxaa*(Oaa)) );
	}

	//Mitochondia; Intemediate Variable;  Optimize Things that use this;Dependent on NAD(KmIDNAD_NAD)
	void calculateVFH() {
		//Parameters: kfFH, KFHeq
		VFH = kfFH * (FUM) - kfFH_KFHeq * (MAL);
	}

	//Mitochondia; Intemediate Variable;  Optimize Things that use this;
	void calculateVMDH() {
		//Parameters: Kh1, Kh2, Kh3, Kh4, Koff, mH, Kioaa, KmmNAD, Kmal, EtMD, kMDH
		VMDH = kMDH_Fh_EtMD * (MAL) * NAD / ( ( (MAL) + Kmal + (Oaa) * Kmal_Kioaa ) * ( KmmNAD + NAD ) );
	}

	//Mitochondia; Intemediate Variable;  Optimize Things that use this;
	void calculateVAAT() {
		//Parameters: kfAAT, GLU, kcnsASP, KAATeq, 
		VAAT = VAAT_Constant * (Oaa) / ( kcnsASP_KAATeq_kfAAT + (AKG) );
	}

	//Mitochondia; Matlab model Calculation of VAAT. Uses ASP
	void calculateVAAT_ASP() {
		//Parameters: kfAAT, GLU, KAATeq
		VAAT = KfAAT_GLU * (Oaa) - KfAAT_KAATeq * (ASP) * (AKG);
	}

	//Mitochondia; Intemediate Variable;  Optimize Things that use this;RT_over_F; F_over_RT;KmIDNAD_NAD;DmuH;exp_6_FRT_Dpsio;exp_FRT_6_g_DmuH;
	void calculateVNO_VHNe() {
		//Parameters: kres, rhoREN, Dpsio, mmg, ra, rc1, r1, rc2, rb
		double AREN         =  std::sqrt( (NADH) * kres_sq_KmIDNAD * KmIDNAD_NAD );
		double denominator  = 1.0 / ( (exp_6_FRT_Dpsio + r1_exp_6_FRT_Dpsio *  AREN ) + ( r2 + r3 *  AREN ) * exp_FRT_6_g_DmuH );
	
		VNO                 = ( (rhoREN_ra_rc1_exp_6_FRT_Dpsio + rhoREN_rc2 * exp_FRT_6_g_DmuH) *  AREN  - rhoREN_ra * exp_FRT_6_g_DmuH ) * denominator;
		VHNe                = (rhoRen_6_ra * AREN  - rhoRen_6_ra_rb * exp_FRT_6_g_DmuH ) *	denominator;
	}

	//Mitochondia; Intemediate Variable;  Optimize Things that use this;RT_over_F;exp_6_FRT_Dpsio;exp_FRT_6_g_DmuH;
	void calculateVFO_VHFe() {
		//Parameters: kresf, FADH2, FAD, r2, r3, rhoREF
		double denominator = rhoREF / (VFO_VHFe_C1 + r2_r3_exp_AREF_FRT * exp_FRT_6_g_DmuH);
		//Unused
		//VFO  = ( VFO_C1 - ra_rc2_exp_AREF_FRT * exp_FRT_6_g_DmuH) * denominator;
		VHFe = ( ra_exp_AREF_FRT - ra_rb*exp_FRT_6_g_DmuH ) * denominator;
	}

	//Mitochondia; Intemediate Variable;  Optimize Things that use this; RT_over_F; exp_3_FRT_Dpsio;F_over_RT; 1/ADPm; ATPm
	void calculateVATPase_Vhu() {
		//This has a numeric instability due to round off (I guess) dependent on the order of evaluation of variables in VATPase
		// DO NOT combine the exp_3FRT_DmuH statements, breaks the model
		//Parameters: kf1, Pi, pa, pc1, pc2, p1, p2, p3 ,rhoF1
		//Bad naming on some of the parameters
		double exp_3FRT_DmuH  = std::exp( FRT_3 * DmuH );
		double AF1            =  kf1_Pi * ATPm / (ADPm);
		double denominator    = -rhoF1 / ( exp_3_FRT_Dpsio + p1_exp_3_FRT_Dpsio * AF1 + ( p2 + p3 * AF1 ) * exp_3FRT_DmuH);
	//	VATPase = ( VATPase_C1 * AF1 - (pa + pc2 * AF1) * exp_3FRT_DmuH ) * denominator;
		VATPase               = ( (VATPase_C1 + pc2 * exp_3FRT_DmuH) * AF1 - pa * exp_3FRT_DmuH ) * denominator;
		Vhu                   = ( pa_300 + pa_300 * AF1  - pa_pb_3 * exp_3FRT_DmuH ) * denominator; // Rate of H uptake by mitochondrial via F1FoATPase
	}

	//Mitochondia;Different from paper; Intemediate Variable;  Optimize Things that use this; ATPm;ADP
	void calculateVANT_Vhleak() {
		//Parameters: gh, VmDT; hm
		double ATPi_ADP   = (ATPi) / ADP;
		double ADPm_ATPm  = (ADPm) / ATPm;
		VANT              = ( VmDT_75 - VmDT_20 * ATPi_ADP * ADPm_ATPm * exp(-F_over_RT * (Dpsi)) ) /	( ( 1.0 + tenDiv9 * ATPi_ADP * exp( -hm_F_over_RT * (Dpsi) ) ) * (1.0 + 18 * ADPm_ATPm ) );
		Vhleak            = gh * DmuH;
	}

	//Mitochondia;Vuni; F_over_RT; FRT2_Dpsi; Treat Cai as parameter in Mitochondia only model
	void calculateVuni( double cai) {
		//Parameters: Vmuni; ktrans; mmL; kact; na
		double Cai_ktrans_plus1     = 1.0 + (cai) * inv_ktrans;
		double Cai_ktrans_plus1_p3  = Cai_ktrans_plus1 * Cai_ktrans_plus1 * Cai_ktrans_plus1;
		Vuni                        =  ( Vmuni_ktrans * (cai) * FRT2_Dpsi * Cai_ktrans_plus1_p3  / ( ( Cai_ktrans_plus1_p3 * Cai_ktrans_plus1 + mmL / std::pow( 1.0 + (cai) * inv_kact, na) ) * ( 1.0 - std::exp(-FRT2_Dpsi) ) ) );
	}

	///Mitochondia; F_over_RT; FRT2_Dpsi; 1/Nai; 1/Cam; 1/Cai ; Cai Nai as params in Mito only
	void calculateVnaCa( double cai, double nai ) {
		//Parameters: b; VmNC; Kna; mmn; Knca;
		VnaCa = VmNC * std::exp( b_05 * FRT2_Dpsi ) * (Cam) / ( (cai) * std::pow( ( 1.0 + Kna / (nai) ) , mmn) * ( 1.0 + Knca / (Cam) ) );
	}

	//C_K Mod, Comes before calculateFATPi. Combination of getFCrPi(),getCK()
	void calculateCK() {
		//Parameters: kt_2, kf_2, kf_3, keq, CRT_cyto, CRT_cyto, CRT_mito, VATPase_cyto
		//Precomputes: inv_keq = 1.0 / keq;
		double Vt_CRP2  = kt_2 * ( (CrPi_mito) - (CrPi_cyto) );
		double VCK_cyto = kf_2 * ( ( CRT_cyto - (CrPi_cyto) ) * (ATPi_cyto) -  (CrPi_cyto) * ( 8.0 - (ATPi_cyto) ) * inv_keq);
		double VCK_mito = kf_3 * ( ( CRT_mito - (CrPi_mito) ) * (ATPi) - (CrPi_mito) * ADP * inv_keq);

    //		(*dCrPi_mito) = VCK_mito - Vt_CRP2;
    //		(*dCrPi_cyto) = Vt_CRP2 + VCK_cyto;
    //		(*dATPi) = 0.615 * VANT - V_AM - 0.5*Jup - ( 6.371e-5 * ( INaK + IpCa ) ) - VCK_mito; // defined normally in calculateFATPi
    //		(*dATPi_cyto) = - VCK_cyto - VATPase_cyto;
	}

	void sharedVariableUpdateGPC_1( double v ) {
		//Before IpCa
		ADP             = 8.0 - (ATPi);
		inv_ATPi        = 1.0 / (ATPi);

		//Done after V -> Algebraic Potential, before ICA
		VF_over_RT      = (v) * F_over_RT;
		exp_VF_over_RT  = exp( VF_over_RT );
		VFsq_over_RT    = FaradayE3 * VF_over_RT;	//VF_over_RT is unitless mV/mV
		exp2VFRT        = exp_VF_over_RT * exp_VF_over_RT;

		V_30            = (v) + 30.0;
	}
	
	//Common variable computations for the GPC Model Part 2
	//Done after calculateReversalPotentials and before calculateIk1
	void sharedVariableUpdateGPC_2( double v, double EK ) {
		//Done After calculateReversalPotentials, before Ik1
		V_E_K = (v) - EK;
	}

	//Common variable computations for the Mitochondria Model
	void sharedVariableUpdateMitochondria() {
		//Must be done after calculateNAD is Called
		KmIDNAD_NAD       = KmIDNAD / NAD;
		//After calculateDmuH
		exp_FRT_6_g_DmuH  = exp(FRT_6_g * DmuH);
		//Anytime
		FRT2_Dpsi         = FRT2 * ( (Dpsi) - 91.0 );
	}

// private:
  // Variables from cardiomyocyte
  double mito_L                        = 0.01;
  double mito_rad                      = 0.0011;
  double mito_vcell                    = 1000.0*3.14*mito_rad*mito_rad*mito_L;
  double mito_Ageo                     = 2.0*std::numbers::pi*mito_rad*mito_rad+2.0*std::numbers::pi*mito_rad*mito_L;
  double Acap                          = 2.0*mito_Ageo;

  // Physical constants
  double R                             = 8314.0;
  double T                             = 310.0;
  double F                             = 96485.0;

  //==========================
  // Metabolism and Energetics
  //==========================
  double vmito                         = 0.26*mito_vcell;

  double V_AM_scaler                   = 1.0; // variable?
  double Cmito                         = 1.812e-3;
  double C_m                           = 1.0;
  double mAcap                         = Acap;// 1.534E-4
  double mCm                           = 1.5; // Should this be Cmito or cell?
  double DpH                           = -0.6;
  double RT_over_F                     = R*T/F;
  double F_over_RT                     = F/(R*T);
  double CPN                           = 10.0; //mM
  double KCS                           = 0.05; //ms-1
  double EtCS                          = 0.4; // mM
  double AcCoA                         = 1.0; // mM
  double KmAcCoA                       = 1.26e-2; // mM
  double KmOaa                         = 6.4e-4; //mM
  double kfACO                         = 1.25e-2; // ms-1
  double CIK                           = 1.0; //mM
  double KACOeq                        = 2.22;
  double KADP                          = 0.62;
  double KaCa                          = 0.0005;
  double KidhNADH                      = 0.19;
  double kIDH                          = 0.05;
  double EtID                          = 0.109;
  double mH                            = 2.5e-5;
  double kh_1                          = 8.1e-5;
  double kh_2                          = 5.98e-5;
  double KmIDNAD                       = 0.923;
  double Kmiso                         = 1.52;
  double nID                           = 2.0;
  double Kmg                           = 0.0308;
  double Mg                            = 0.4;
  double Kca                           = 0.00127;
  double kKGDH                         = 7.5e-2;
  double EtKG                          = 0.5;
  double KmKG                          = 1.94;
  double nKG                           = 1.2;
  double KmKGNAD                       = 38.7;
  double kfSL                          = 0.005;
  double kSDH                          = 0.005;
  double EtSDH                         = 0.5;
  double KmSucc                        = 0.03;
  double KmSucc_KiFUM                  = 0.03/1.3;
  double KiOxaa                        = 0.15;
  double kfFH                          = 0.00332;
  double KFHeq                         = 1.0;
  double Kh1                           = 1.131e-5;
  double Kh2                           = 26.7;
  double Kh3                           = 6.68e-9;
  double Kh4                           = 5.62e-6;
  double Koff                          = 0.0399;
  double Kioaa                         = 0.0031;
  double KmmNAD                        = 0.2244;
  double Kmal                          = 1.493;
  double EtMD                          = 0.154;
  double kMDH                          = 0.111;
  double kfAAT                         = 6.44e-4;
  double GLU                           = 10.0;
  double kcnsASP                       = 1.5e-6;
  double KAATeq                        = 6.6;
  double kres                          = 1.35e18;
  double rhoREN                        = 1.0e-1;
  double Dpsio                         = 50.0;
  double mmg                           = 0.85;
  double ra                            = 6.394e-13;
  double rb                            = 1.762e-16;
  double rc1                           = 2.656e-22;
  double rc2                           = 8.632e-30;
  double r1                            = 2.077e-18;
  double r2                            = 1.728e-9;
  double r3                            = 1.059e-26;
  double kresf                         = 5.765e13;
  double FADH2                         = 1.24;
  double FAD                           = 0.01;
  double rhoREF                        = 3.75e-4;
  double kf1                           = 1.71e6;
  double Pi                            = 2.0;
  double p1                            = 1.346e-8;
  double p2                            = 7.739e-7;
  double p3                            = 6.65e-15;
  double pa                            = 1.656e-8;
  double pb                            = 3.373e-10;
  double pc1                           = 9.651e-17;
  double pc2                           = 4.585e-17;
  double rhoF1                         = 0.05;
  double gh                            = 2.0e-7;
  double VmDT                          = 0.015;
  double hm                            = 0.5;
  double ktrans                        = 0.019;
  double kact                          = 0.00038;
  double Vmuni                         = 0.0275;
  double b                             = 0.5;
  double na                            = 2.8;
  double VmNC                          = 0.625e-4;
  double Kna                           = 9.4;
  double mmn                           = 3.0;
  double Knca                          = 3.75e-4;
  double kt_2                          = 2.0e-3;
  double kf_2                          = 1.4e-4;
  double kf_3                          = 1.33e-6;
  double keq                           = 0.0095;
  double CRT_cyto                      = 25.0;
  double CRT_mito                      = 25.0;
  double VATPase_cyto                  = 1.0e-5;
  double fm                            = 0.0003;
  double KmATP_SR                      = 0.01;
  double Ki_SR                         = 0.14;
  double Ki_prime_SR                   = 5.1;
  double Km1ATP_CaP                    = 0.012;
  double Km2ATP_CaP                    = 0.23;
  double KiADP_CaP                     = 1.0;
  double KmpCa                         = 0.0005;
  double Faraday                       = F/1000.0;
  double CoA                           = 0.02;
  double aL                            = 2.0;
  double bL                            = 2.0;
  double Kfb                           = 0.00024;
  double Ki_AM                         = 0.26;
  double KiFUM                         = 1.3;
  double KmATP_AM                      = 0.03;
  double KmNai                         = 10.0;
  double Krb                           = 1.64269;
  double KSLeq                         = 3.115;
  double mmL                           = 110.0;
  double V_AM_max                      = 0.00048;
  double Vmyo                          = 0.68*mito_vcell;
  double FaradayE3                     = 1000.0*Faraday;
  double inv_keq                       = 0.0; //CK Mod
  double Vmito                         = 0.0;
  double Vtotal                        = 0.0;
  double one_inv_KACOeq                = 0.0; //Mitochondria Model
  double inv_11                        = 0.0;
  double inv_11p1                      = 0.0;
  double inv_6                         = 0.0;
  double inv_6p8                       = 0.0;
  double inv_7p5                       = 0.0;
  double inv_9p5                       = 0.0;
  double inv_ATPi                      = 0.0;
  double inv_Cmito                     = 0.0; //Mitochondria Model
  double inv_KADP                      = 0.0;  //Mitochondria Model
  double inv_KaCa                      = 0.0;  //Mitochondria Model
  double inv_kact                      = 0.0;  //Mitochondria Model
  double inv_Ki_prime_SR               = 0.0;
  double inv_Ki1AD_NaK                 = 0.0;
  double inv_KiADP_CaP                 = 0.0;
  double inv_KidhNADH                  = 0.0; //Mitochondria Model
  double inv_KiOxaa                    = 0.0; //Mitochondria Model
  double KmNaiP1p5                     = 0.0;
  double inv_ktrans                    = 0.0;  //Mitochondria Model
  double inv_LTRPNtot_Ktrop_half       = 0.0;   //Force Model
  double tenDiv9                       = 0.0; //Mitochondria Model
  double twoThirds                     = 0.0;  //Force Model
  double two_b                         = 0.0;  //Mitochondria Model
  double FRT2                          = 0.0; //Mitochondria Model
  double ADP                           = 0.0;
  double AREF                          = 0.0; //Mitochondria Model
  double b_05                          = 0.0; //Mitochondria Model
  double Co                            = 0.0;
  double CoA_KSLeq                     = 0.0;  //Mitochondria Model
  double DmuH_Constant                 = 0.0; //Mitochondria Model
  double exp_AREF_FRT                  = 0.0; //Mitochondria Model
  double exp_3_FRT_Dpsio               = 0.0; //Mitochondria Model
  double exp_6_FRT_Dpsio               = 0.0; //Mitochondria Model
  double f_23_g_12_mod                 = 0.0; //Force Model
  double FRT_3                         = 0.0; //Mitochondria Model
  double FRT_6_g                       = 0.0; //Mitochondria Model
  double hm_F_over_RT                  = 0.0; //Mitochondria Model
  double kcnsASP_KAATeq_kfAAT          = 0.0;   //Mitochondria Model
  double kf1_Pi                        = 0.0; //Mitochondria Model
  double KfAAT_GLU                     = 0.0;  //Mitochondria Model
  double KfAAT_KAATeq                  = 0.0; //Mitochondria Model
  double kfFH_KFHeq                    = 0.0;  //Mitochondria Model
  double kIDH_EtID                     = 0.0;  //Mitochondria Model
  double kKGDH_EtKG                    = 0.0;  //Mitochondria Model
  double Kmal_Kioaa                    = 0.0;  //Mitochondria Model
  double KmATP_AM_Ki_AM                = 0.0;
  double KmATP_SR_Ki_SR                = 0.0;
  double kMDH_Fh_EtMD                  = 0.0; //Mitochondria Model
  double KmKGNAD_KmIDNAD               = 0.0; //Mitochondria Model
  double KmnsCa_p3                     = 0.0;
  double kres_sq_KmIDNAD               = 0.0; //Mitochondria Model
  double kSDH_EtSDH                    = 0.0;  //Mitochondria Model
  double Mg_Kmg_1                      = 0.0;  //Mitochondria Model
  double Mg_Kmg_1_Kca                  = 0.0; //Mitochondria Model
  double p1_exp_3_FRT_Dpsio            = 0.0;    //Mitochondria Model
  double pa_pb_3                       = 0.0; //Mitochondria Model
  double pa_300                        = 0.0; //Mitochondria Model
  double r1_exp_6_FRT_Dpsio            = 0.0;    //Mitochondria Model
  double r2_r3_exp_AREF_FRT            = 0.0;    //Mitochondria Model
  double ra_rc1_exp_6_FRT_Dpsio        = 0.0;   //Mitochondria Model
  double ra_exp_AREF_FRT               = 0.0; //Mitochondria Model
  double ra_rb                         = 0.0; //Mitochondria Model
  double ra_rc2_exp_AREF_FRT           = 0.0;    //Mitochondria Model
  double rhoRen_6_ra                   = 0.0;  //Mitochondria Model
  double rhoRen_6_ra_rb                = 0.0; //Mitochondria Model
  double rhoREN_ra                     = 0.0;  //Mitochondria Model
  double rhoREN_ra_rc1_exp_6_FRT_Dpsio = 0.0; //Mitochondria Model
  double rhoREN_rc2                    = 0.0;  //Mitochondria Model
  double RTF_05                        = 0.0;
  double RTF_05_log_Cao                = 0.0;
  double RTF_log_Ko                    = 0.0;
  double RTF_log_Nao                   = 0.0;
  double RTF_log_Nao_Ko                = 0.0;
  double V_30                          = 0.0;
  double V_AM_scaler_max_1_f_01_12_23  = 0.0;
  double VAAT_Constant                 = 0.0; //Mitochondria Model
  double VATPase_C1                    = 0.0;  //Mitochondria Model
  double VCS_C1                        = 0.0; //Mitochondria Model
  double VFO_C1                        = 0.0; //Mitochondria Model
  double VFO_VHFe_C1                   = 0.0;  //Mitochondria Model
  double VIDH_Constant                 = 0.0; //Mitochondria Model
  double VmDT_75                       = 0.0; //Mitochondria Model
  double VmDT_20                       = 0.0; //Mitochondria Model
  double Vmuni_ktrans                  = 0.0; //Mitochondria Model
  double Vmyo_1000_F_Acap_C_m          = 0.0;
  double Vmyo_VNSR                     = 0.0;
  double Vmyo_VSS                      = 0.0;
  double KmIDNAD_NAD                   = 0.0;   //Mitochondria Model
  double exp_FRT_6_g_DmuH              = 0.0; //Mitochondria Model
  double FRT2_Dpsi                     = 0.0;   //Mitochondria Model
  double kltrpn_plus_Cai               = 0.0;  //Force Model
  double V_E_K                         = 0.0;
  double VFsq_over_RT                  = 0.0;
  double exp2VFRT                      = 0.0;
  double VF_over_RT                    = 0.0;
  double exp_VF_over_RT                = 0.0;
  double V_AM                          = 0.0;
  double ATPm                          = 0.0;  //Mitochondria Model
  double DmuH                          = 0.0;  //Mitochondria Model
  double NAD                           = 0.0;   //Mitochondria Model
  double VCS                           = 0.0;   //Mitochondria Model
  double VACO                          = 0.0;  //Mitochondria Model
  double VIDH                          = 0.0;  //Mitochondria Model
  double VKGDH                         = 0.0;  //Mitochondria Model
  double VSL                           = 0.0;   //Mitochondria Model
  double VSDH                          = 0.0;  //Mitochondria Model
  double VFH                           = 0.0;   //Mitochondria Model
  double VMDH                          = 0.0;  //Mitochondria Model
  double VAAT                          = 0.0;  //Mitochondria Model
  double VNO                           = 0.0;   //Mitochondria Model
  double VHNe                          = 0.0;  //Mitochondria Model
  double VHFe                          = 0.0;  //Mitochondria Model
  double VATPase                       = 0.0;  //Mitochondria Model
  double Vhu                           = 0.0;   //Mitochondria Model
  double VANT                          = 0.0;  //Mitochondria Model
  double Vhleak                        = 0.0;  //Mitochondria Model
  double FN_Ca                         = 0.0;  //Force Model
  double Jtrpn                         = 0.0;
  double beta_SS                       = 0.0;
  double beta_JSR                      = 0.0;
  double beta_i                        = 0.0;
  double Vuni                          = 0.0;  //Mitochondria Model
  double VnaCa                         = 0.0;  //Mitochondria Model
  double Vt_CRP2                       = 0.0;
  double VCK_cyto                      = 0.0;
  double VCK_mito                      = 0.0;

  //==============================
  // Mitochondrial State Variables
  //==============================
  // From CellML
  double ATPi_cyto                     = 7.97824;
  double ATPi                          = 7.98691;
  double CrPi_cyto                     = 19.4302;
  double CrPi_mito                     = 19.4385;
  double Cam                           = 0.000137836;
  double ADPm                          = 0.0258454;
  double Dpsi                          = 164.073;
  double NADH                          = 7.76311;
  double ISOC                          = 0.496682;
  double AKG                           = 0.000610554;
  double SCoA                          = 0.247364;
  double Succ                          = 0.000394005;
  double FUM                           = 0.0191924;
  double MAL                           = 0.00925885;
  double Oaa                           = 1.08962e-07;
  //Mito Matlab version, get a good value from Sonia later
  double ASP                           = 0.0502614;


  //===========================================
  // Initial values for mitochondria ODEs: LSODA
  //===========================================
  static const int neq_mitochondria            = 16;
  int iwork1_mitochondria                      = 0;
  int iwork2_mitochondria                      = 0;
  int iwork5_mitochondria                      = 0;
  int iwork6_mitochondria                      = 0;
  int iwork7_mitochondria                      = 0;
  int iwork8_mitochondria                      = 0;
  int iwork9_mitochondria                      = 0;
  double rwork1_mitochondria                   = 0.0;
  double rwork5_mitochondria                   = 0.0;
  double rwork6_mitochondria                   = 0.0;
  double rwork7_mitochondria                   = 0.0;
  double y_mitochondria[neq_mitochondria+1];

  double t_begin_mitochondria                  = 0.0;
  double dt_mitochondria                       = 0.005;
  double tout_mitochondria                     = t_begin_mitochondria + dt_mitochondria;
  int itol_mitochondria                        = 2;
  double rtol_mitochondria[neq_mitochondria+1];
  double atol_mitochondria[neq_mitochondria+1];
  double itask_mitochondria                    = 1;
  int istate_mitochondria                      = 1;
  double iopt_mitochondria                     = 0;
  double jt_mitochondria                       = 2;

  LSODA<Mitochondria> lsoda_integrator;
};
