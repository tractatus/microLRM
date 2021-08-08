/******************************************************************************
* Coherent Point Drift
* Copyright (C) 2014 Pete Gadomski <pete.gadomski@gmail.com>
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
******************************************************************************/

#include <cpdArma/nonrigid_lowrank.hpp>

#include "affinity_eigenvectors.hpp"
#include "debug.hpp"
#include "sigma2.hpp"
#include "spdiag_locations.hpp"


namespace cpdArma
{


NonrigidLowrank::NonrigidLowrank(float tol, int max_it, float outliers,
                                 bool use_fgt,
                                 float epsilon, float beta, float lambda, arma::uword numeig)
    : Nonrigid(tol, max_it, outliers, use_fgt, epsilon, beta, lambda)
    , m_numeig(numeig)
{}



Registration::ResultPtr NonrigidLowrank::execute(const arma::mat& X,
        const arma::mat& Y, double sigma2) const
{
    const arma::uword N = X.n_rows;
    const arma::uword M = Y.n_rows;
    const arma::uword D = Y.n_cols;

    DEBUG("parameters, tol: " << get_tol() <<
          ", max_it: " << get_max_it() <<
          ", outliers: " << get_outliers() <<
          ", epsilon: " << get_epsilon() <<
          ", beta: " << get_beta() <<
          ", lambda: " << get_lambda());

    const arma::uword numeig = (get_numeig() == 0) ?
                               arma::uword(std::sqrt(Y.n_rows)) :
                               get_numeig();
    DEBUG("number of eigenvectors: " << numeig);

    const double sigma2_init = sigma2;

    arma::mat T = Y;
    arma::mat W = arma::zeros<arma::mat>(M, D);

    int iter = 0;
    double ntol = get_tol() + 10;
    double L = 0;

    arma::mat Q, S;
    find_affinity_eigenvectors(Y, get_beta(), numeig, get_epsilon(), Q, S);

    arma::sp_mat invS(spdiag_locations(numeig), 1 / S.diag());

    double L_old, Np;
    arma::vec P1(M), Pt1(M);
    arma::mat PX(M, D);

    while (iter < get_max_it() &&
           ntol > get_tol() &&
           sigma2 > 10 * std::numeric_limits<double>::epsilon())
    {
        //std::cout << "nonrigid_lowrank iteration: dL= " << ntol << ", iter= " << iter << ", sigma2= " << sigma2 <<  std::endl;
        L_old = L;
        arma::mat QtW = Q.t() * W;

        L = find_P(X, T, sigma2, P1, Pt1, PX);

        L = L + get_lambda() / 2 * arma::trace(QtW.t() * S * QtW);
        ntol = std::abs((L - L_old) / L);

        DEBUG("nonrigid_lowrank iteration: dL= " << ntol <<
              ", iter= " << iter << ", sigma2= " << sigma2);

        arma::sp_mat dP(spdiag_locations(M), P1);
        arma::mat dPQ = dP * Q;
        arma::mat F = PX - dP * Y;

        W = 1 / (get_lambda() * sigma2) *
            (F - dPQ * (arma::solve(get_lambda() * sigma2 * invS + Q.t() * dPQ,
                                    Q.t() * F)));

        T = Y + (Q * (S * (Q.t() * W)));

        Np = arma::sum(P1);
        sigma2 = std::abs(
                     arma::accu(arma::pow(X, 2) % arma::repmat(Pt1, 1, D)) +
                     arma::accu(arma::pow(T, 2) % arma::repmat(P1, 1, D)) -
                     2 * arma::trace(PX.t() * T)) /
                 (Np * D);

        ++iter;
    }

    ResultPtr result(new Result());
    result->Y = T;
    return result;
}


}
