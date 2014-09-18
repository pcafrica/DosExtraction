#include "dosModel.h"

using namespace std::chrono;
using namespace constants;
using namespace utility;

DosModel::DosModel()
	: initialized_(false) {}

DosModel::DosModel(const ParamList & params)
	: initialized_(true), params_(params) {}

void DosModel::simulate(const GetPot & config, const std::string & input_experim, const std::string & output_directory,
                        const std::string & output_plot_subdir, const std::string & output_filename) const
{
	if ( !initialized_ ) {
		throw std::logic_error("ERROR: parameters list of DosModel has not been properly initialized.");
	}
	
	// Define output filenames.
	const std::string output_fitting_filename = output_filename + "_fitting.txt";
	const std::string output_CV_filename = output_filename + "_CV.csv";
	
	// Open output files.
	std::ofstream output_fitting;
	std::ofstream output_CV;
	
	output_fitting.open(output_directory + output_fitting_filename, std::ios_base::out);
	
	output_CV.open(output_directory + output_CV_filename, std::ios_base::out);
	output_CV.setf(std::ios_base::scientific);
	output_CV.precision(std::numeric_limits<Real>::digits10);
	
	if ( !output_fitting.is_open() || !output_CV.is_open() ) {
		throw std::ofstream::failure("ERROR: output files cannot be opened or directory does not exist.");
	}
	
	output_fitting << "Running on thread: " << omp_get_thread_num() << "." << std::endl;
	
	// Timing.
	high_resolution_clock::time_point initTime = high_resolution_clock::now();
	
	print_block( ( "Simulation No. " + std::to_string(params_.simulationNo_) + " started.").c_str(), output_fitting );
	
	VectorXr V = VectorXr::LinSpaced(params_.nSteps_, params_.V_min_, params_.V_max_);
	
	Index semicNodesNo = floor( 0.6 * params_.nNodes_ );
	Index   insNodesNo = params_.nNodes_ - semicNodesNo;
	
	// Mesh creation.
	output_fitting << "Creating mesh...";
	VectorXr x = VectorXr::Zero( params_.nNodes_ );	// The mesh.
	
	{
		VectorXr temp1 = VectorXr::LinSpaced(  semicNodesNo, -params_.t_semic_,              0);
		VectorXr temp2 = VectorXr::LinSpaced(insNodesNo + 1,                 0, params_.t_ins_);
		
		x << temp1, temp2.segment(1, temp2.size() - 1);
	}
	
	VectorXr xm = 0.5 * ( x.segment(1, x.size() - 1) + x.segment(0, x.size() - 1) );
	print_done(output_fitting);
	
	// System assembly.
	output_fitting << "Assembling system matrices...";
	VectorXr eps = EPS0 * params_.eps_semic_ * VectorXr::Ones( xm.size() );
	
	for ( int i = 0; i < eps.size(); ++i ) {
		if ( xm(i) > 0.0 ) {
			eps(i) = EPS0 * params_.eps_ins_;
		}
	}
	
	Bim1D bimSolver(x);
	bimSolver.assembleStiff( eps, VectorXr::Ones( params_.nNodes_ ) );
	
	{
		VectorXr temp = VectorXr::Zero( xm.size() );
		
		for ( int i = 0; i < temp.size(); ++i ) {
			if ( xm(i) < 0.0 ) {
				temp(i) = 1.0;
			}
		}
		
		bimSolver.assembleMass( temp, VectorXr::Ones( params_.nNodes_ ) );
	}
	
	SparseXr A = bimSolver.Stiff();	// Stiffness matrix.
	SparseXr M = bimSolver.Mass ();	// Mass matrix.
	
	print_done(output_fitting);
	
	// Computing nodes and weights of quadrature rule.
	output_fitting << "Computing nodes and weights of quadrature rule";
	GaussHermiteRule rule( config("GaussHermite/nNodes", 101) );
	output_fitting << " using " << rule.nNodes() << " nodes...";
	
	try {
		rule.apply(config);
	} catch ( const std::runtime_error & runtimeError ) {
		throw runtimeError;
	}
	
	print_done(output_fitting);
	
	// Variables initialization.
	output_fitting << "Initializing variables...";
	GaussianCharge charge_fun(params_, rule);
	
	MatrixXr Phi  = MatrixXr::Zero(     x.size(), V.size() );
	MatrixXr Dens = MatrixXr::Zero( semicNodesNo, V.size() );
	
	VectorXr cTot     = VectorXr::Zero( V.size() );
	VectorXr charge_n = VectorXr::Zero( V.size() );
	
	print_done(output_fitting);
	
	output_fitting << "Running Newton solver for non-linear Poisson equation..." << std::endl;
	output_fitting << "\tMax No. of iterations set: " << config("NLP/maxIterationsNo", 1000) << std::endl;
	output_fitting << "\tTolerance set: " << config("NLP/tolerance", 1.0e-4) << std::endl;
	
	// Start simulation.
	for ( Index i = 0; i < V.size(); ++i ) {
		// Print current iteration.
		if ( i == 0 || (i + 1) % 10 == 0 || i == V.size() - 1 ) {
			output_fitting << std::endl << "\titeration: " << (i + 1) << "/" << params_.nSteps_;
		}
		
		VectorXr phiOld = VectorXr::Zero( x.size() );
		
		if ( i == 0 ) {
			phiOld = -VectorXr::LinSpaced(phiOld.size(), params_.Wf_ / Q - params_.Ea_ / Q, params_.Wf_ / Q - params_.Ea_ / Q - V(i));
		} else {
			phiOld = Phi.col(i - 1) + VectorXr::LinSpaced(phiOld.size(), 0, V(i) - V(i - 1));
		}
		
		NonLinearPoisson1D nlpSolver( bimSolver, config("NLP/maxIterationsNo", 1), config("NLP/tolerance", 1.0e-4) );
		nlpSolver.apply(x, phiOld, charge_fun);
		
		Phi.col(i) = nlpSolver.phi();
		
		VectorXr charge = charge_fun.charge( Phi.col(i).segment(0, semicNodesNo) );
		Dens.col(i) = - charge / Q;
		
		cTot(i) = nlpSolver.cTot();
		
		charge_n(i) = numerics::trapz((VectorXr) x.segment(0, semicNodesNo), charge);
	}
	
	print_done(output_fitting);
	
	// Timing.
	high_resolution_clock::time_point finalTime = high_resolution_clock::now();
	output_fitting << "Simulation took " << duration_cast<seconds>(finalTime - initTime).count()
	               << " seconds." << std::endl;
	               
	// Post-processing and creation of output files.
	
	post_process(config, input_experim, output_fitting, output_CV, (VectorXr) x.segment(0, semicNodesNo),
	             (VectorXr) Dens.col(Dens.cols() - 1), V, cTot);
	             
	output_fitting.close();
	output_CV     .close();
	
	// Create output Gnuplot files.
	try {
		save_plot(output_directory, output_plot_subdir, output_CV_filename, output_filename);
	} catch ( const std::exception & genericException ) {
		throw genericException;
	}
	
	return;
}

void DosModel::post_process(const GetPot & config, const std::string & input_experim, std::ostream & output_fitting,
                            std::ostream & output_CV, const VectorXr & x_semic, const VectorXr & dens,
                            const VectorXr & V_simulated, const VectorXr & C_simulated) const
{
	assert( x_semic    .size() == dens       .size() );
	assert( V_simulated.size() == C_simulated.size() );
	
	CsvParser parser_experim(input_experim, config("hasHeaders", true));
	
	Real A_semic = parser_experim.importCell(1, 3);	// Semiconductor area [m^2].
	Real C_sb    = parser_experim.importCell(1, 4);	// Sbord capacitance [F].
	
	VectorXr V_experim = parser_experim.importCol(1);
	VectorXr C_experim = parser_experim.importCol(2);
	
	assert( V_experim.size() == C_experim.size() );
	
	// Sorting "V_experim" and "C_experim". The order is established by "V_experim".
	{
		VectorXpair<Real> order = numerics::sort_pair( V_experim );
		VectorXr V_sort = VectorXr::Zero( V_experim.size() );
		VectorXr C_sort = VectorXr::Zero( C_experim.size() );
		
		for ( int i = 0; i < order.size(); ++i ) {
			V_sort(i) = V_experim( order(i).second );
			C_sort(i) = C_experim( order(i).second );
		}
		
		V_experim = V_sort;
		C_experim = C_sort;
	}
	
	VectorXr dC_dV_experim   = numerics::deriv(C_experim,                            V_experim  );
	VectorXr dC_dV_simulated = numerics::deriv(C_simulated.array() * A_semic + C_sb, V_simulated);
	
	Real charge_center_of_mass = numerics::trapz( x_semic.cwiseProduct(dens) ) / numerics::trapz( dens );
	Real cAccStar = C_simulated.maxCoeff();	// Simulated.
	
	Index j_e = 0;
	dC_dV_experim  .maxCoeff(&j_e);
	
	Index j_s = 0;
	dC_dV_simulated.maxCoeff(&j_s);
	
	Real V_shift = V_simulated(j_s) - V_experim(j_e);
	
	VectorXr     C_interp = numerics::interp1(V_experim,     C_experim, V_simulated.array() - V_shift);
	VectorXr dC_dV_interp = numerics::interp1(V_experim, dC_dV_experim, V_simulated.array() - V_shift);
	
	Real error_L2 = std::sqrt( numerics::error_L2(C_interp, C_simulated.array() * A_semic + C_sb, V_simulated, V_shift) );
	Real error_H1 = std::sqrt( error_L2 * error_L2 +
	                           numerics::error_L2(dC_dV_interp, dC_dV_simulated, V_simulated, V_shift)
	                         );
	                         
	// Print to output.
	output_fitting << std::endl;
	output_fitting << "V_shift = " << V_shift << std::endl;
	output_fitting << "Charge center of mass = " << charge_center_of_mass << std::endl;
	output_fitting << "C_acc* = " << cAccStar << std::endl;
	output_fitting << std::endl;
	output_fitting << "Distance between experimental and simulated voltage values:" << std::endl;
	output_fitting << "\t L2-distance = " << error_L2 << std::endl;
	output_fitting << "\t H1-distance = " << error_H1 << std::endl;
	
	output_CV << "V_experim, C_experim, dC/dV_experim, V_simulated, C_simulated, dC/dV_simulated" << std::endl;
	
	for ( int i = 0; i < std::max( V_simulated.size(), V_experim.size() ); ++i ) {
		if ( i < V_experim.size() ) {
			output_CV << V_experim(i) << ", " << C_experim(i) << ", " << dC_dV_experim(i) << ", ";
		} else {
			output_CV << ",,, ";
		}
		
		if ( i < V_simulated.size() ) {
			output_CV << V_simulated(i) - V_shift << ", " << C_simulated(i) * A_semic + C_sb << ", " << dC_dV_simulated(i);
		} else {
			output_CV << ",,";
		}
		
		output_CV << std::endl;
	}
	
	return;
}

void DosModel::gnuplot_commands(const std::string & output_CV_filename, std::ostream & os) const
{
	os << "set datafile separator \",\";\n";
	os << "set format y \"%.2te%+03T\";\n";
	os << "\n";
	os << "stats \"" + output_CV_filename + "\" using 1 name \"V\" nooutput;\n";
	os << "\n";
	os << "set multiplot;\n";
	os << "set size 1, 0.5;\n";
	os << "\n";
	os << "set origin 0.0, 0.5;\n";
	os << "\tset xlabel \"V_gate [V]\";\n";
	os << "\tset ylabel \"dC/dV [F/V]\";\n";
	os << "\tset title \"dC/dV [F/V], N0 = " << params_.N0_ << "\";\n";
	os << "\tplot [V_min:V_max] \"" + output_CV_filename + "\" using 1:3 title \"Experimental\" with lines lw 2, \\\n";
	os << "\t                   \"" + output_CV_filename + "\" using 4:6 title \"Simulated\"    with lines lw 2;\n";
	os << "\n";
	os << "set origin 0.0, 0.0;\n";
	os << "\tset xlabel \"V_gate [V]\";\n";
	os << "\tset ylabel \"C [F]\";\n";
	os << "\tset title \"C [F]\";\n";
	os << "\tplot [V_min:V_max] \"" + output_CV_filename + "\" using 1:2 title \"Experimental\" with lines lw 2, \\\n";
	os << "\t                   \"" + output_CV_filename + "\" using 4:5 title \"Simulated\"    with lines lw 2;\n";
	os << "\n";
	os << "unset multiplot;\n";
	
	return;
}

void DosModel::save_plot(const std::string & output_directory, const std::string & output_plot_subdir,
                         const std::string & output_CV_filename, const std::string & output_filename) const
{
	// Save script for later reuse.
	const std::string output_plot_filename = output_plot_subdir + output_filename + "_plot.gp";
	
	std::ofstream output_plot;
	output_plot.open(output_directory + output_plot_filename, std::ios_base::out);
	
	if ( !output_plot.is_open() ) {
		throw std::runtime_error("ERROR: Gnuplot output file cannot be opened or directory does not exist.");
	}
	
	gnuplot_commands("../" + output_CV_filename, output_plot);
	output_plot << "\n";
	output_plot << "pause mouse;\n";
	output_plot.close();
	
	// Create .pdf plot file.
	Gnuplot output_pdf;
	output_pdf << "set terminal pdf;\n";
	output_pdf << "set output \"" + output_directory + output_filename + "_plot.pdf\";\n";
	output_pdf << "\n";
	gnuplot_commands(output_directory + output_CV_filename, output_pdf);
	output_pdf << "\n";
	output_pdf << "set output;\n";
	
	return;
}
