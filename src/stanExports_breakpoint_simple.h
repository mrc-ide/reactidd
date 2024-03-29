// Generated by rstantools.  Do not edit by hand.

/*
    reactidd is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    reactidd is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with reactidd.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_breakpoint_simple_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_breakpoint_simple");
    reader.add_event(81, 79, "end", "model_breakpoint_simple");
    return reader;
}
#include <stan_meta_header.hpp>
class model_breakpoint_simple
  : public stan::model::model_base_crtp<model_breakpoint_simple> {
private:
        int num_data;
        int num_all;
        row_vector_d Y;
        row_vector_d N;
        std::vector<double> X;
        double max_X;
        std::vector<double> X_all;
        std::vector<int> indexes;
        std::vector<double> restriction_dates;
        double tau1;
        double tau2;
        double tau3;
        double tau4;
        double tau5;
public:
    model_breakpoint_simple(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_breakpoint_simple(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_breakpoint_simple_namespace::model_breakpoint_simple";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "num_data", "int", context__.to_vec());
            num_data = int(0);
            vals_i__ = context__.vals_i("num_data");
            pos__ = 0;
            num_data = vals_i__[pos__++];
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "num_all", "int", context__.to_vec());
            num_all = int(0);
            vals_i__ = context__.vals_i("num_all");
            pos__ = 0;
            num_all = vals_i__[pos__++];
            current_statement_begin__ = 6;
            validate_non_negative_index("Y", "num_data", num_data);
            context__.validate_dims("data initialization", "Y", "row_vector_d", context__.to_vec(num_data));
            Y = Eigen::Matrix<double, 1, Eigen::Dynamic>(num_data);
            vals_r__ = context__.vals_r("Y");
            pos__ = 0;
            size_t Y_j_1_max__ = num_data;
            for (size_t j_1__ = 0; j_1__ < Y_j_1_max__; ++j_1__) {
                Y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("N", "num_data", num_data);
            context__.validate_dims("data initialization", "N", "row_vector_d", context__.to_vec(num_data));
            N = Eigen::Matrix<double, 1, Eigen::Dynamic>(num_data);
            vals_r__ = context__.vals_r("N");
            pos__ = 0;
            size_t N_j_1_max__ = num_data;
            for (size_t j_1__ = 0; j_1__ < N_j_1_max__; ++j_1__) {
                N(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("X", "num_data", num_data);
            context__.validate_dims("data initialization", "X", "double", context__.to_vec(num_data));
            X = std::vector<double>(num_data, double(0));
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_k_0_max__ = num_data;
            for (size_t k_0__ = 0; k_0__ < X_k_0_max__; ++k_0__) {
                X[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 9;
            context__.validate_dims("data initialization", "max_X", "double", context__.to_vec());
            max_X = double(0);
            vals_r__ = context__.vals_r("max_X");
            pos__ = 0;
            max_X = vals_r__[pos__++];
            current_statement_begin__ = 10;
            validate_non_negative_index("X_all", "num_all", num_all);
            context__.validate_dims("data initialization", "X_all", "double", context__.to_vec(num_all));
            X_all = std::vector<double>(num_all, double(0));
            vals_r__ = context__.vals_r("X_all");
            pos__ = 0;
            size_t X_all_k_0_max__ = num_all;
            for (size_t k_0__ = 0; k_0__ < X_all_k_0_max__; ++k_0__) {
                X_all[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("indexes", "num_data", num_data);
            context__.validate_dims("data initialization", "indexes", "int", context__.to_vec(num_data));
            indexes = std::vector<int>(num_data, int(0));
            vals_i__ = context__.vals_i("indexes");
            pos__ = 0;
            size_t indexes_k_0_max__ = num_data;
            for (size_t k_0__ = 0; k_0__ < indexes_k_0_max__; ++k_0__) {
                indexes[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 12;
            validate_non_negative_index("restriction_dates", "5", 5);
            context__.validate_dims("data initialization", "restriction_dates", "double", context__.to_vec(5));
            restriction_dates = std::vector<double>(5, double(0));
            vals_r__ = context__.vals_r("restriction_dates");
            pos__ = 0;
            size_t restriction_dates_k_0_max__ = 5;
            for (size_t k_0__ = 0; k_0__ < restriction_dates_k_0_max__; ++k_0__) {
                restriction_dates[k_0__] = vals_r__[pos__++];
            }
            // initialize transformed data variables
            current_statement_begin__ = 16;
            tau1 = double(0);
            stan::math::fill(tau1, DUMMY_VAR__);
            current_statement_begin__ = 17;
            tau2 = double(0);
            stan::math::fill(tau2, DUMMY_VAR__);
            current_statement_begin__ = 18;
            tau3 = double(0);
            stan::math::fill(tau3, DUMMY_VAR__);
            current_statement_begin__ = 19;
            tau4 = double(0);
            stan::math::fill(tau4, DUMMY_VAR__);
            current_statement_begin__ = 20;
            tau5 = double(0);
            stan::math::fill(tau5, DUMMY_VAR__);
            // execute transformed data statements
            current_statement_begin__ = 22;
            stan::math::assign(tau1, get_base1(restriction_dates, 1, "restriction_dates", 1));
            current_statement_begin__ = 23;
            stan::math::assign(tau2, get_base1(restriction_dates, 2, "restriction_dates", 1));
            current_statement_begin__ = 24;
            stan::math::assign(tau3, get_base1(restriction_dates, 3, "restriction_dates", 1));
            current_statement_begin__ = 25;
            stan::math::assign(tau4, get_base1(restriction_dates, 4, "restriction_dates", 1));
            current_statement_begin__ = 26;
            stan::math::assign(tau5, get_base1(restriction_dates, 5, "restriction_dates", 1));
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 31;
            num_params_r__ += 1;
            current_statement_begin__ = 32;
            validate_non_negative_index("beta", "6", 6);
            num_params_r__ += 6;
            current_statement_begin__ = 33;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_breakpoint_simple() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 31;
        if (!(context__.contains_r("prev0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable prev0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("prev0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "prev0", "double", context__.to_vec());
        double prev0(0);
        prev0 = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, prev0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable prev0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 32;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "6", 6);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(6));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(6);
        size_t beta_j_1_max__ = 6;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lub_unconstrain(-(0.5), 0.5, beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 33;
        if (!(context__.contains_r("delay")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable delay missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("delay");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "delay", "double", context__.to_vec());
        double delay(0);
        delay = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0, 14, delay);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable delay: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 31;
            local_scalar_t__ prev0;
            (void) prev0;  // dummy to suppress unused var warning
            if (jacobian__)
                prev0 = in__.scalar_lb_constrain(0, lp__);
            else
                prev0 = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 32;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_lub_constrain(-(0.5), 0.5, 6, lp__);
            else
                beta = in__.vector_lub_constrain(-(0.5), 0.5, 6);
            current_statement_begin__ = 33;
            local_scalar_t__ delay;
            (void) delay;  // dummy to suppress unused var warning
            if (jacobian__)
                delay = in__.scalar_lub_constrain(0, 14, lp__);
            else
                delay = in__.scalar_lub_constrain(0, 14);
            // transformed parameters
            current_statement_begin__ = 38;
            validate_non_negative_index("beta_t", "num_all", num_all);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta_t(num_all);
            stan::math::initialize(beta_t, DUMMY_VAR__);
            stan::math::fill(beta_t, DUMMY_VAR__);
            current_statement_begin__ = 39;
            validate_non_negative_index("prev_t", "num_all", num_all);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> prev_t(num_all);
            stan::math::initialize(prev_t, DUMMY_VAR__);
            stan::math::fill(prev_t, DUMMY_VAR__);
            current_statement_begin__ = 40;
            validate_non_negative_index("Y_hat", "num_data", num_data);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> Y_hat(num_data);
            stan::math::initialize(Y_hat, DUMMY_VAR__);
            stan::math::fill(Y_hat, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 42;
            for (int i = 1; i <= num_all; ++i) {
                current_statement_begin__ = 43;
                if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau1 + delay)))) {
                    current_statement_begin__ = 44;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 1, "beta", 1)), 
                                "assigning variable beta_t");
                } else if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau2 + delay)))) {
                    current_statement_begin__ = 47;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 2, "beta", 1)), 
                                "assigning variable beta_t");
                } else if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau3 + delay)))) {
                    current_statement_begin__ = 50;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 3, "beta", 1)), 
                                "assigning variable beta_t");
                } else if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau4 + delay)))) {
                    current_statement_begin__ = 53;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 4, "beta", 1)), 
                                "assigning variable beta_t");
                } else if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau5 + delay)))) {
                    current_statement_begin__ = 56;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 5, "beta", 1)), 
                                "assigning variable beta_t");
                } else {
                    current_statement_begin__ = 59;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 6, "beta", 1)), 
                                "assigning variable beta_t");
                }
            }
            current_statement_begin__ = 64;
            stan::model::assign(prev_t, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        prev0, 
                        "assigning variable prev_t");
            current_statement_begin__ = 66;
            for (int i = 2; i <= num_all; ++i) {
                current_statement_begin__ = 67;
                stan::model::assign(prev_t, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (get_base1(prev_t, (i - 1), "prev_t", 1) * get_base1(beta_t, (i - 1), "beta_t", 1)), 
                            "assigning variable prev_t");
            }
            current_statement_begin__ = 70;
            stan::math::assign(Y_hat, stan::model::rvalue(prev_t, stan::model::cons_list(stan::model::index_multi(indexes), stan::model::nil_index_list()), "prev_t"));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 38;
            size_t beta_t_j_1_max__ = num_all;
            for (size_t j_1__ = 0; j_1__ < beta_t_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(beta_t(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: beta_t" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable beta_t: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 39;
            size_t prev_t_j_1_max__ = num_all;
            for (size_t j_1__ = 0; j_1__ < prev_t_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(prev_t(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: prev_t" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable prev_t: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            current_statement_begin__ = 40;
            size_t Y_hat_j_1_max__ = num_data;
            for (size_t j_1__ = 0; j_1__ < Y_hat_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(Y_hat(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: Y_hat" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable Y_hat: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 78;
            lp_accum__.add((multiply(Y, stan::math::log(Y_hat)) + multiply(subtract(N, Y), stan::math::log(subtract(1, Y_hat)))));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("prev0");
        names__.push_back("beta");
        names__.push_back("delay");
        names__.push_back("beta_t");
        names__.push_back("prev_t");
        names__.push_back("Y_hat");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(6);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(num_all);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(num_all);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(num_data);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_breakpoint_simple_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double prev0 = in__.scalar_lb_constrain(0);
        vars__.push_back(prev0);
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_lub_constrain(-(0.5), 0.5, 6);
        size_t beta_j_1_max__ = 6;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        double delay = in__.scalar_lub_constrain(0, 14);
        vars__.push_back(delay);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 38;
            validate_non_negative_index("beta_t", "num_all", num_all);
            Eigen::Matrix<double, Eigen::Dynamic, 1> beta_t(num_all);
            stan::math::initialize(beta_t, DUMMY_VAR__);
            stan::math::fill(beta_t, DUMMY_VAR__);
            current_statement_begin__ = 39;
            validate_non_negative_index("prev_t", "num_all", num_all);
            Eigen::Matrix<double, Eigen::Dynamic, 1> prev_t(num_all);
            stan::math::initialize(prev_t, DUMMY_VAR__);
            stan::math::fill(prev_t, DUMMY_VAR__);
            current_statement_begin__ = 40;
            validate_non_negative_index("Y_hat", "num_data", num_data);
            Eigen::Matrix<double, Eigen::Dynamic, 1> Y_hat(num_data);
            stan::math::initialize(Y_hat, DUMMY_VAR__);
            stan::math::fill(Y_hat, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 42;
            for (int i = 1; i <= num_all; ++i) {
                current_statement_begin__ = 43;
                if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau1 + delay)))) {
                    current_statement_begin__ = 44;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 1, "beta", 1)), 
                                "assigning variable beta_t");
                } else if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau2 + delay)))) {
                    current_statement_begin__ = 47;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 2, "beta", 1)), 
                                "assigning variable beta_t");
                } else if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau3 + delay)))) {
                    current_statement_begin__ = 50;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 3, "beta", 1)), 
                                "assigning variable beta_t");
                } else if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau4 + delay)))) {
                    current_statement_begin__ = 53;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 4, "beta", 1)), 
                                "assigning variable beta_t");
                } else if (as_bool(logical_lt(get_base1(X_all, i, "X_all", 1), (tau5 + delay)))) {
                    current_statement_begin__ = 56;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 5, "beta", 1)), 
                                "assigning variable beta_t");
                } else {
                    current_statement_begin__ = 59;
                    stan::model::assign(beta_t, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                stan::math::exp(get_base1(beta, 6, "beta", 1)), 
                                "assigning variable beta_t");
                }
            }
            current_statement_begin__ = 64;
            stan::model::assign(prev_t, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        prev0, 
                        "assigning variable prev_t");
            current_statement_begin__ = 66;
            for (int i = 2; i <= num_all; ++i) {
                current_statement_begin__ = 67;
                stan::model::assign(prev_t, 
                            stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                            (get_base1(prev_t, (i - 1), "prev_t", 1) * get_base1(beta_t, (i - 1), "beta_t", 1)), 
                            "assigning variable prev_t");
            }
            current_statement_begin__ = 70;
            stan::math::assign(Y_hat, stan::model::rvalue(prev_t, stan::model::cons_list(stan::model::index_multi(indexes), stan::model::nil_index_list()), "prev_t"));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t beta_t_j_1_max__ = num_all;
                for (size_t j_1__ = 0; j_1__ < beta_t_j_1_max__; ++j_1__) {
                    vars__.push_back(beta_t(j_1__));
                }
                size_t prev_t_j_1_max__ = num_all;
                for (size_t j_1__ = 0; j_1__ < prev_t_j_1_max__; ++j_1__) {
                    vars__.push_back(prev_t(j_1__));
                }
                size_t Y_hat_j_1_max__ = num_data;
                for (size_t j_1__ = 0; j_1__ < Y_hat_j_1_max__; ++j_1__) {
                    vars__.push_back(Y_hat(j_1__));
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_breakpoint_simple";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "prev0";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = 6;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "delay";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t beta_t_j_1_max__ = num_all;
            for (size_t j_1__ = 0; j_1__ < beta_t_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta_t" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t prev_t_j_1_max__ = num_all;
            for (size_t j_1__ = 0; j_1__ < prev_t_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "prev_t" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t Y_hat_j_1_max__ = num_data;
            for (size_t j_1__ = 0; j_1__ < Y_hat_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Y_hat" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "prev0";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = 6;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "delay";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t beta_t_j_1_max__ = num_all;
            for (size_t j_1__ = 0; j_1__ < beta_t_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "beta_t" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t prev_t_j_1_max__ = num_all;
            for (size_t j_1__ = 0; j_1__ < prev_t_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "prev_t" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t Y_hat_j_1_max__ = num_data;
            for (size_t j_1__ = 0; j_1__ < Y_hat_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Y_hat" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_breakpoint_simple_namespace::model_breakpoint_simple stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
