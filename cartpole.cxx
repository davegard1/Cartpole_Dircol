//////////////////////////////////////////////////////////////////////////
////////////////        low_thrust.cxx               /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Low thrust orbit transfer problem         ////////////////
//////// Last modified: 16 February 2009                  ////////////////
//////// Reference:     Betts  (2001)             	  ////////////////
//////// (See PSOPT handbook for full reference)           ///////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
#include <iostream>
#include "structs.h"

using namespace PSOPT;



//////////////////////////////////////////////////////////////////////////
///////////////////  Define auxiliary functions                 //////////
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{

    return 0.;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{

    if (iphase == 1) {
	   adouble u = controls[0];
	   return (u*u);
   }
   else {
	   return (0);
   }
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   // Local integers
   int i, j;

	const double m1_ {1.0};
	const double m2_ {0.3};
	const double l_ {0.5};
	const double g_ {-9.81};


	// Input states
    adouble* u = controls;

	adouble x1 = states[0]; // cart position
	adouble x2 = states[1]; // cart speed
	adouble x3 = states[2]; // pendulum angle
	adouble x4 = states[3]; // pendulum angle rate of change

	// Intermediate states
	adouble denom;

	denom = (1. - cos(x3) * cos(x3)) ;
	
    adouble u1 = u[0];

	// Output states
	adouble y1, y2, y3, y4;

	y1 = x2;

	y2 = (l_*m2_ * sin(x3) * x4*x4 + u1 + m2_*g_* cos(x3)*sin(x3) )/ (m1_ + m2_ * denom);

	y3 = x4;

	y4 = -(l_*m2_*cos(x3)*sin(x3) * x4*x4 + u1*cos(x3) + (m1_+m2_)*g_*sin(x3) ) / (l_*m1_ + l_*m2_ * denom);





   derivatives[ 0 ] = y1;
   derivatives[ 1 ] = y2;
   derivatives[ 2 ] = y3;
   derivatives[ 3 ] = y4;


//    path[ 0 ] = powu[0];

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{


   int offset;

   adouble x0   = initial_states[ 0 ];
   adouble xd0  = initial_states[ 1 ];
   adouble th0  = initial_states[ 2 ];
   adouble thd0 = initial_states[ 3 ];


   adouble xf   = final_states[ 0 ];
   adouble xdf  = final_states[ 1 ];
   adouble thf  = final_states[ 2 ];
   adouble thdf = final_states[ 3 ];



   if (iphase==1) {
   	e[ 0 ]  = x0;
   	e[ 1 ]  = xd0;
   	e[ 2 ]  = th0;
   	e[ 3 ]  = thd0;
   }

   if (1 == 1) offset = 4;
   else offset = 0;

   if (iphase == 1 ) {
   	e[ offset + 0 ]  = xf;
   	e[ offset + 1 ]  = xdf;
   	// e[ offset + 2 ]  = asin ( sin ( (thf - M_PI) / 2.0 ) ) / M_PI;
   	e[ offset + 2 ]  = thf;
   	e[ offset + 3 ]  = thdf;
   }

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
 //  auto_link_multiple(linkages, xad, 1);
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    MSdata msdata;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        				= "Cartpole problem";
    problem.outfilename                 = "cartpole.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			          = 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   				= 4;
    problem.phases(1).ncontrols 				= 1;
    problem.phases(1).nparameters               = 0;
    problem.phases(1).nevents   	        	= 8;
    problem.phases(1).npath     				= 0;
    problem.phases(1).nodes                     << 10;

    psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    double x0 = 0.;
    double xd0 = 0.0;
    double th0 = 0.0;
    double thd0 = 0.0;
    double tf = 4.0;



    // problem.phases(1).bounds.lower.parameters <<  tauL;

    // problem.phases(1).bounds.upper.parameters <<  tauU;

    problem.phases(1).bounds.lower.states << -100., -10., -M_PI*2, -M_PI*2;

    problem.phases(1).bounds.upper.states << 100., 10., M_PI*2, M_PI*2;

    problem.phases(1).bounds.lower.controls << -20.0;

    problem.phases(1).bounds.upper.controls << 20.0;

    problem.phases(1).bounds.lower.events << x0, xd0, th0, thd0, -10., 0.0, M_PI, 0.0;

    problem.phases(1).bounds.upper.events << x0, xd0, th0, thd0, 10., 0.0, M_PI, 0.0;

    double EQ_TOL = 0.001;


    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.0;
    problem.phases(1).bounds.upper.EndTime      = tf * 2.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 		    = &integrand_cost;
    problem.endpoint_cost 			= &endpoint_cost;
    problem.dae             		= &dae;
    problem.events 					= &events;
    problem.linkages				= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			= 20;
    int ncontrols          = problem.phases(1).ncontrols;
    int nstates            = problem.phases(1).nstates;

    MatrixXd u_guess    =  zeros(ncontrols,nnodes);
    MatrixXd x_guess    =  zeros(nstates,nnodes);
    MatrixXd time_guess =  linspace(0.0,5.0,nnodes);

    MatrixXd param_guess = zeros(1,1);


    // x_guess.block(3,0,1,nnodes) = Eigen::MatrixXd::LinSpaced(nnodes,th0,0.0);
    x_guess.block(3,0,1,nnodes) =  linspace(th0, M_PI, nnodes);
    x_guess.block(1,0,1,nnodes) =  linspace(x0, tf, nnodes);

    // std::cout << u_guess.rows() << std::endl;
    // std::cout << u_guess.cols() << std::endl;
    // std::cout << u_guess.block(3,0,1,nnodes) << std::endl;

    // u_guess.block(3,0,1,nnodes) = Eigen::MatrixXd::Ones(1,nnodes);
    // x_guess    = load_data("../../../examples/low_thrust/X0.dat",nstates  , nnodes );
    // time_guess = load_data("../../../examples/low_thrust/T0.dat",1        , nnodes );


    auto_phase_guess(problem, u_guess, x_guess, param_guess, time_guess);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.ipopt_linear_solver         = "ma97";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    // algorithm.defect_scaling              = "jacobian-based";
    // algorithm.jac_sparsity_ratio          =  0.11; // 0.05;
    algorithm.collocation_method          = "Hermite-Simpson";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_increment_factor     = 0.2;
    algorithm.mr_max_iterations           = 8;





////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd x, u, t;

    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);


    // t = t/3600.0;


    // MatrixXd tau = solution.get_parameters_in_phase(1);

    // Print(tau,"tau");

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////


    OutputObs output;
    output.fname = "../cartpole_dircol.dat";
    
    output.f_output.open(output.fname);

    double n = x.cols();
    for (int i = 0; i < n; ++i)
    {
        output.f_output << t(i) << ", " << x.col(i).format(output.fmt) << std::endl;
    }

    output.f_output.close();


    OutputObs output_u;
    output_u.fname = "../cartpole_dircol_u.dat";
    
    output_u.f_output.open(output_u.fname);

    for (int i = 0; i < n; ++i)
    {
        output_u.f_output << t(i) << ", " << u.col(i).format(output_u.fmt) << std::endl;
    }

    output_u.f_output.close();
////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    // MatrixXd x1 = x.row(0); 
    // MatrixXd x2 = x.row(1); 
    // MatrixXd x3 = x.row(2); 
    // MatrixXd x4 = x.row(3); 

    // // MatrixXd u1 = u.row(0); 


    // plot(t,x1,problem.name+": states", "time (s)", "x (m)","x (m)");

    // plot(t,x2,problem.name+": states", "time (s)", "xdot (m/s)","xdot");

    // plot(t,x3,problem.name+": states", "time (s)", "theta (deg)","theta");

    // plot(t,x4,problem.name+": states", "time (s)", "theta dot (deg/s)","theta dot");



    // plot(t,u,problem.name+": controls","time (s)", "u", "u");




}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
