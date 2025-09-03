#ifndef NDEBUG
#define NDEBUG
#endif
// #undef NDEBUG

#include <assert.h>
#include <cmath>
#include <random>
#include <iostream>

#include "mex.h"
#include "mex_eig_dense.h"

// simulation type
enum simtype {Euler, Heun};

inline constexpr simtype SIMTYPE = Euler;

#include "parameters.h"
#include "reactions.h"

// number of threads
// inline constexpr size_t ntheads = 8;

// matlab fun: [c,t] = reaction_sim(p,c0);
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    // set the number of threads
    // omp_set_num_threads(ntheads);
    // Eigen::setNbThreads(ntheads);

    // Eigen::initParallel();

    // load parameters
    load_parameters(prhs);

    // simulation time step
    dt = p.t_sim/double(p.n_step);

    // initial concentrations
    Eigen::Map<Eigen::VectorXd> c0 = mat2eig_vector(prhs[1]);
    // followed species
    Eigen::Matrix<size_t, Eigen::Dynamic, 1> spec_out = mat2eig_size_t_vector(prhs[2]);

    // output: c(t)
    Eigen::Map<Eigen::MatrixXd> c_out = eig2mat_matrix(plhs[0],p.n_data+1,spec_out.rows());
    // output: t
    Eigen::Map<Eigen::VectorXd> t_out = eig2mat_vector(plhs[1],p.n_data+1);

    double c[NSPEC];//,dc[NSPEC];
    // n_out_pos: position in the output
    // n_out: how often we would like to write data to the output
    size_t n_out_pos,n_out = p.n_step/p.n_data;

    //#pragma omp parallel private(c,t,n_out_pos) // dc
    {

        // random distribution
        // random starting point
        unsigned seed = time(NULL);
        std::mt19937 generator(seed);
        //         std::normal_distribution<double> distribution(0,sqrt(dt/p.omega));
        std::normal_distribution<double> distribution(0.0, p.sigma_inf * std::sqrt(1.0 - std::exp(-2.0 * dt / p.tau)));
        //         std::normal_distribution<double> distribution(0.0, 1.0);

        // initial concentration
        for (size_t i=0;i<NSPEC;++i) {
            c[i] = c0(i);
        }
        t = 0.0;
        n_out_pos=0;

        // output concentrations for the initial state
        t_out(n_out_pos) = t;
        for (size_t nspec=0;nspec<(size_t)spec_out.rows();++nspec){
            c_out(n_out_pos,nspec) = c[spec_out(nspec)];
        }

        // simulation steps
        for(size_t step=1;step<=p.n_step;++step) {
            t += dt;

            // reactions => calculate dc
            reactions(c,distribution,generator);

            // output concentrations
            if ((step % n_out) == 0) {
                ++n_out_pos;
                t_out(n_out_pos) = t;
                for (size_t nspec=0;nspec<(size_t)spec_out.rows();++nspec){
                    c_out(n_out_pos,nspec) = c[spec_out(nspec)];
                }
            }
        }
    }


    return;
}