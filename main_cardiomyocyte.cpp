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
#include <fstream>
#include <sstream>
#include "cardiomyocyte.hpp"

auto main() -> int {
  auto celltype = CellType::ENDO;
  auto isHFpEF  = true;

  double inal_therapy = 0.21;

  std::ostringstream file_name;

  std::ofstream output_filename;
  if ( celltype == CellType::EPI ) {
    std::cout << "Running an Epicardial cell ..." << std::endl;
    if (isHFpEF) 
      output_filename.open("outputs/epi_hfpef.txt");
    else
      output_filename.open("outputs/epi.txt");
  }
  else if ( celltype == CellType::MCELL ) {
    std::cout << "Running an Mid-myocardial cell ..." << std::endl;
    if (isHFpEF) 
      output_filename.open("outputs/mcell_hfpef.txt");
    else
      output_filename.open("outputs/mcell.txt");
  }
  else { //( celltype == LVENDO )
    std::cout << "Running an Endocardial cell ..." << std::endl;
    file_name << "outputs/endo";  
    if (isHFpEF) {
      // file_name << "_hfpef_pacing" << (pacing_rate*100.0) << "percent.txt";
      file_name << "_hfpef_inal_therapy" << (inal_therapy*100.0) << "percent.txt";
      output_filename.open(file_name.str());
    }
    else {
      file_name << "_inal_therapy" << inal_therapy << "percent.txt";
      output_filename.open(file_name.str());
    }
  }

  // Run cell
  CardioMyocyte ordcell( celltype, isHFpEF );
  ordcell.run( output_filename );
}
