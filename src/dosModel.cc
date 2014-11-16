/* C++11 */

/**
 * @file   dosModel.cc
 * @author Pasquale Claudio Africa <pasquale.africa@gmail.com>
 * @date   2014
 *
 * This file is part of the "DosExtraction" project.
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 */

#include "dosModel.h"

using namespace gnuplotio;

using namespace std::chrono;
using namespace constants;
using namespace utility;

DosModel::DosModel()
    : initialized_(false), params_(), V_shift_(0.0), error_L2_(0.0), error_H1_(0.0), error_Peak_(0.0),
      C_acc_experim_(0.0), C_acc_simulated_(0.0), C_dep_experim_(0.0) {}

DosModel::DosModel(const ParamList & params)
    : initialized_(true), params_(params), V_shift_(0.0), error_L2_(0.0), error_H1_(0.0), error_Peak_(0.0),
      C_acc_experim_(0.0), C_acc_simulated_(0.0), C_dep_experim_(0.0) {}

void DosModel::simulate(const GetPot & config, const std::string & input_experim, const std::string & output_directory,
                        const std::string & output_plot_subdir, const std::string & output_filename)
{
    if ( !initialized_ )
    {
        throw std::logic_error("ERROR: list of parameters in DosModel has not been properly initialized.");
    }
    
    // Define output filenames.
    const std::string output_info_filename = output_filename + "_info.txt";
    const std::string output_CV_filename = output_filename + "_CV.csv";
    
    // Open output files.
    std::ofstream output_info;
    std::ofstream output_CV;
    
    output_info.open(output_directory + output_info_filename, std::ios_base::out);
    
    output_CV.open(output_directory + output_CV_filename, std::ios_base::out);
    output_CV.setf(std::ios_base::scientific);
    output_CV.precision(std::numeric_limits<Real>::digits10);
    
    if ( !output_info.is_open() || !output_CV.is_open() )
    {
        throw std::ofstream::failure("ERROR: output files cannot be opened or directory does not exist.");
    }
    
    output_info << "Running on thread: " << omp_get_thread_num() << "." << std::endl;
    
    // Timing.
    high_resolution_clock::time_point initTime = high_resolution_clock::now();
    
    print_block( ( "Simulation No. " + std::to_string(params_.simulationNo_) + " started.").c_str(), output_info );
    
    VectorXr V = VectorXr::LinSpaced(params_.nSteps_, params_.V_min_, params_.V_max_);
    
    Index semicNodesNo = floor( 0.6 * params_.nNodes_ );
    Index   insNodesNo = params_.nNodes_ - semicNodesNo;
    
    // Mesh creation.
    output_info << "Creating mesh...";
    VectorXr x = VectorXr::Zero( params_.nNodes_ );    // The mesh.
    
    {
        VectorXr temp1 = VectorXr::LinSpaced(  semicNodesNo, -params_.t_semic_,              0);
        VectorXr temp2 = VectorXr::LinSpaced(insNodesNo + 1,                 0, params_.t_ins_);
        
        x << temp1, temp2.segment(1, temp2.size() - 1);
    }
    
    VectorXr xm = 0.5 * ( x.segment(1, x.size() - 1) + x.segment(0, x.size() - 1) );
    print_done(output_info);
    
    // System assembly.
    output_info << "Assembling system matrices...";
    VectorXr eps = params_.eps_semic_ * VectorXr::Ones( xm.size() );
    
    for ( Index i = 0; i < eps.size(); ++i )
    {
        if ( xm(i) > 0.0 )
        {
            eps(i) = params_.eps_ins_;
        }
    }
    
    Bim1D bimSolver(x);
    
    bimSolver.assembleStiff( eps, VectorXr::Ones( params_.nNodes_ ) );
    
    {
        VectorXr temp = VectorXr::Zero( xm.size() );
        
        for ( Index i = 0; i < temp.size(); ++i )
        {
            if ( xm(i) < 0.0 )
            {
                temp(i) = 1.0;
            }
        }
        
        bimSolver.assembleMass( temp, VectorXr::Ones( params_.nNodes_ ) );
    }
    
    SparseXr A = bimSolver.Stiff();    // Stiffness matrix.
    SparseXr M = bimSolver.Mass ();    // Mass matrix.
    
    print_done(output_info);
    
    // Computing nodes and weights of quadrature.
    output_info << "Computing nodes and weights of quadrature";
    
    QuadratureRule * quadRule = nullptr;
    
    {
        QuadratureRuleFactory * quadRuleFactory;
        
        Index rule = config("QuadratureRule/rule", 1);
        
        switch ( rule )
        {
            case 1:
                output_info << " (Gauss-Hermite rule)";
                quadRuleFactory = new GaussHermiteRuleFactory;
                break;
                
            case 0:
                output_info << " (Gauss-Laguerre rule)";
                quadRuleFactory = new GaussLaguerreRuleFactory;
                break;
                
            default:
                throw std::runtime_error("ERROR: wrong variable \"rule\" set in the configuration file (only 1 or 0 allowed).");
                break;
        }
        
        quadRule = quadRuleFactory->BuildRule( config("QuadratureRule/nNodes", 101) );
        
        delete quadRuleFactory;
    }
    
    output_info << " using " << quadRule->nNodes() << " nodes...";
    
    try
    {
        quadRule->apply(config);
    }
    catch ( const std::exception & genericException )
    {
        throw;
    }
    
    print_done(output_info);
    
    // Constitutive relation.
    output_info << "Initializing constitutive relation for the Density of States";
    
    Charge * charge_fun = nullptr;
    
    {
        ChargeFactory * chargeFactory;
        
        Index constitutive_relation = config("DOS", 1);
        
        switch ( constitutive_relation )
        {
            case 1:
                output_info << " (Gaussian)";
                chargeFactory = new GaussianChargeFactory;
                break;
                
            case 0:
                output_info << " (Exponential)";
                chargeFactory = new ExponentialChargeFactory;
                break;
                
            default:
                throw std::runtime_error("ERROR: wrong variable \"DOS\" set in the configuration file (only 1 or 0 allowed).");
                break;
        }
        
        charge_fun = chargeFactory->BuildCharge(params_, *quadRule);
        
        delete chargeFactory;
    }
    
    output_info << "...";
    print_done(output_info);
    
    // Initialize Newton solver for the non-linear Poisson equation.
    Index maxIterationsNo = config("NLP/maxIterationsNo", 100);
    Real  tolerance       = config("NLP/tolerance", 1.0e-4);
    
    NonLinearPoisson1D nlpSolver(bimSolver, maxIterationsNo, tolerance);
    
    // Variables initialization.
    output_info << "Initializing variables...";
    
    MatrixXr Phi  = MatrixXr::Zero(     x.size(), V.size() );
    MatrixXr Dens = MatrixXr::Zero( semicNodesNo, V.size() );
    
    VectorXr cTot     = VectorXr::Zero( V.size() );
    VectorXr charge_n = VectorXr::Zero( V.size() );
    
    print_done(output_info);
    
    output_info << "Running Newton solver for non-linear Poisson equation..." << std::endl;
    output_info << "\tMax No. of iterations set: " << maxIterationsNo << std::endl;
    output_info << "\tTolerance set: " << tolerance << std::endl;
    
    // Start simulation.
    for ( Index i = 0; i < V.size(); ++i )
    {
        // Print current step number.
        if ( i == 0 || (i + 1) % 10 == 0 || i == V.size() - 1 )
        {
            output_info << std::endl << "\tstep: " << (i + 1) << "/" << params_.nSteps_;
        }
        
        VectorXr phiOld = VectorXr::Zero( x.size() );
        
        if ( i == 0 )
        {
            phiOld = -VectorXr::LinSpaced(phiOld.size(), params_.Wf_ / Q - params_.Ea_ / Q, params_.Wf_ / Q - params_.Ea_ / Q - V(i));
        }
        else
        {
            phiOld = Phi.col(i - 1) + VectorXr::LinSpaced(phiOld.size(), 0, V(i) - V(i - 1));
        }
        
        nlpSolver.apply(phiOld, *charge_fun);
        
        Phi.col(i) = nlpSolver.phi();
        
        VectorXr charge = charge_fun->charge( Phi.col(i).segment(0, semicNodesNo) );
        Dens.col(i) = - charge / Q;
        
        cTot(i) = nlpSolver.cTot();
        
        charge_n(i) = numerics::trapz((VectorXr) x.segment(0, semicNodesNo), charge);
    }
    
    print_done(output_info);
    
    // Timing.
    high_resolution_clock::time_point finalTime = high_resolution_clock::now();
    output_info << "Simulation took " << duration_cast<seconds>(finalTime - initTime).count()
                << " seconds." << std::endl;
                
    // Free up memory to avoid leaks.
    delete quadRule;
    quadRule = nullptr;
    
    delete charge_fun;
    charge_fun = nullptr;
    
    // Post-processing and creation of output files.
    post_process(config, input_experim, output_info, output_CV, params_.A_semic_, params_.C_sb_,
                 (VectorXr) x.segment(0, semicNodesNo), (VectorXr) Dens.col(Dens.cols() - 1), V, cTot);
                 
    output_info << std::endl;
    output_info << "C_sb = " << params_.C_sb_ << std::endl;
    output_info << "t_semic = " << params_.t_semic_ << std::endl;
    
    output_info.close();
    output_CV.close();
    
    // Create output Gnuplot files.
    try
    {
        save_plot(output_directory, output_plot_subdir, output_CV_filename, output_filename);
    }
    catch ( const std::exception & genericException )
    {
        throw;
    }
    
    return;
}

void DosModel::post_process(const GetPot & config, const std::string & input_experim, std::ostream & output_info,
                            std::ostream & output_CV, const Real & A_semic, const Real & C_sb,
                            const VectorXr & x_semic, const VectorXr & dens, const VectorXr & V_simulated,
                            const VectorXr & C_simulated)
{
    assert( x_semic    .size() == dens       .size() );
    assert( V_simulated.size() == C_simulated.size() );
    
    CsvParser parser_experim(input_experim, config("skipHeaders", true));
    
    VectorXr V_experim = parser_experim.importCol(1);
    VectorXr C_experim = parser_experim.importCol(2);
    
    assert( V_experim.size() == C_experim.size() );
    
    // Sorting "V_experim" and "C_experim". The order is established by "V_experim".
    {
        VectorXpair<Real> sort = numerics::sort_pair( V_experim );
        VectorXr V_sort = VectorXr::Zero( V_experim.size() );
        VectorXr C_sort = VectorXr::Zero( C_experim.size() );
        
        for ( Index i = 0; i < sort.size(); ++i )
        {
            V_sort(i) = V_experim( sort(i).second );
            C_sort(i) = C_experim( sort(i).second );
        }
        
        V_experim = V_sort;
        C_experim = C_sort;
    }
    
    VectorXr dC_dV_experim   = numerics::deriv(C_experim,                            V_experim  );
    VectorXr dC_dV_simulated = numerics::deriv(C_simulated.array() * A_semic + C_sb, V_simulated);
    
    Real center_of_charge = numerics::trapz( x_semic.cwiseProduct(dens) ) / numerics::trapz( dens );    // Center of charge.
    Real cAccStar = C_simulated.maxCoeff();    // Simulated.
    
    // Compute V_shift.
    {
        Index j_e = 0;
        dC_dV_experim  .maxCoeff(&j_e);
        
        Index j_s = 0;
        dC_dV_simulated.maxCoeff(&j_s);
        
        V_shift_ = V_simulated(j_s) - V_experim(j_e);
    }
    
    VectorXr     C_interp = numerics::interp1(V_experim,     C_experim, V_simulated.array() - V_shift_);
    VectorXr dC_dV_interp = numerics::interp1(V_experim, dC_dV_experim, V_simulated.array() - V_shift_);
    
    // Save for automatic fitting.
    {
        C_acc_experim_ = C_experim(C_experim.size() - 1);
        
        // Find the value in (V_simulated - V_shift) nearest to V_experim(end).
        Index i = 0;
        (V_simulated.array() - V_shift_ - V_experim(V_experim.size() - 1)).abs().minCoeff(&i);
        
        C_acc_simulated_ = C_simulated(i) * A_semic + C_sb;
        
        C_dep_experim_ = C_experim(0);
    }
    
    // Compute errors.
    error_L2_ = std::sqrt( numerics::error_L2(C_interp, C_simulated.array() * A_semic + C_sb, V_simulated.array() - V_shift_) );
    error_H1_ = std::sqrt( error_L2_ * error_L2_ +
                           numerics::error_L2(dC_dV_interp, dC_dV_simulated, V_simulated.array() - V_shift_)
                         );
    error_Peak_ = std::abs( numerics::nonNaN(dC_dV_interp).maxCoeff() - numerics::nonNaN(dC_dV_simulated).maxCoeff() );
    
    // Print to output.
    output_info << std::endl;
    output_info << "V_shift = " << V_shift_ << std::endl;
    output_info << "Center of charge = " << center_of_charge << std::endl;
    output_info << "C_acc* = " << cAccStar << std::endl;
    output_info << std::endl;
    output_info << "Distance between experimental and simulated capacitance values:" << std::endl;
    output_info << "\tL2-distance = " << error_L2_ << std::endl;
    output_info << "\tH1-distance = " << error_H1_ << std::endl;
    output_info << "\tPeak-distance (on dC/dV) = " << error_Peak_ << std::endl;
    
    output_CV << "V_experim, C_experim, dC/dV_experim, V_simulated, C_simulated, dC/dV_simulated" << std::endl;
    
    for ( Index i = 0; i < std::max( V_simulated.size(), V_experim.size() ); ++i )
    {
        if ( i < V_experim.size() )
        {
            output_CV << V_experim(i) << ", " << C_experim(i) << ", " << dC_dV_experim(i) << ", ";
        }
        else
        {
            output_CV << ",,, ";
        }
        
        if ( i < V_simulated.size() )
        {
            output_CV << V_simulated(i) - V_shift_ << ", " << C_simulated(i) * A_semic + C_sb << ", " << dC_dV_simulated(i);
        }
        else
        {
            output_CV << ",,";
        }
        
        output_CV << std::endl;
    }
    
    return;
}

void DosModel::save_plot(const std::string & output_directory, const std::string & output_plot_subdir,
                         const std::string & csv_filename, const std::string & output_filename) const
{
    // Save script for later reuse.
    const std::string output_plot_filename = output_plot_subdir + output_filename + "_plot.gp";
    
    std::ofstream output_plot;
    output_plot.open(output_directory + output_plot_filename, std::ios_base::out);
    
    if ( !output_plot.is_open() )
    {
        throw std::runtime_error("ERROR: Gnuplot output file cannot be opened or directory does not exist.");
    }
    
    gnuplot_commands("../" + csv_filename, output_plot);
    
    output_plot << std::endl;
    output_plot << "pause mouse;" << std::endl;
    output_plot.close();
    
    // Create .png plot file.
    Gnuplot output_png;
    output_png << "set terminal pngcairo enhanced size 891, 614;" << std::endl;
    output_png << "set output \"" + output_directory + output_filename + "_plot.png\";" << std::endl;
    output_png << std::endl;
    
    gnuplot_commands(output_directory + csv_filename, output_png);
    
    output_png << std::endl;
    output_png << "set output;" << std::endl;
    
    return;
}

void DosModel::gnuplot_commands(const std::string & csv_filename, std::ostream & os) const
{
    os << "set datafile separator \",\";" << std::endl;
    os << "set format y \"%.2te%+03T\";" << std::endl;
    os << std::endl;
    os << "set key right center;" << std::endl;
    os << std::endl;
    os << "stats \"" + csv_filename + "\" using 1 name \"V\" nooutput;" << std::endl;
    os << std::endl;
    os << "set multiplot layout 2, 1 title \"";
    
    os.setf(std::ios_base::scientific);
    os.precision(4);
    os << "N0=" << params_.N0_ << ", σ=" << params_.sigma_ / KB_T << ", T=" << params_.T_ <<  ",  Phi_B=" << (params_.Wf_ - params_.Ea_) / Q;
    os << "\\nN0_2=" << params_.N0_2_ << ", σ_2=" << params_.sigma_2_ / KB_T << ", shift_2=" << params_.shift_2_;
    os << "\\nN0_3=" << params_.N0_3_ << ", σ_3=" << params_.sigma_3_ / KB_T << ", shift_3=" << params_.shift_3_;
    os << "\\nN0_4=" << params_.N0_4_ << ", σ_4=" << params_.sigma_4_ / KB_T << ", shift_4=" << params_.shift_4_;
    os << "\\nN0_e=" << params_.N0_exp_ << ", λ_e=" << params_.lambda_exp_ / KB_T;
    os << "\\nV_{shift}=" << V_shift_ << ", nNodes=" << params_.nNodes_ << ", nSteps=" << params_.nSteps_;
    os << "\" font \", 10\";" << std::endl;
    
    os << "\tset xlabel \"V_{gate} - V_{shift} [V]\" offset 0, 0.75;" << std::endl;
    os << std::endl;
    os << "\tset ylabel \"dC/dV [F/V]\";" << std::endl;
    os << "\tplot [V_min:V_max] \"" + csv_filename + "\" using 1:3 title \"Experimental\" with lines lw 2, \\" << std::endl;
    os << "\t                   \"" + csv_filename + "\" using 4:6 title \"Simulated\"    with lines lw 2;" << std::endl;
    os << std::endl;
    os << "\tset ylabel \"C [F]\";" << std::endl;
    os << "\tplot [V_min:V_max] \"" + csv_filename + "\" using 1:2 title \"Experimental\" with lines lw 2, \\" << std::endl;
    os << "\t                   \"" + csv_filename + "\" using 4:5 title \"Simulated\"    with lines lw 2;" << std::endl;
    os << std::endl;
    os << "unset multiplot;" << std::endl;
    
    return;
}
