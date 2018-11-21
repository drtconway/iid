#include <iostream>
using std::cout; using std::endl;
using std::left; using std::fixed; using std::right; using std::scientific;
#include <iomanip>
using std::setw;
using std::setprecision;
#include <limits>

#include <algorithm>
#include <string>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>
using boost::multiprecision::cpp_dec_float_50;

#include <boost/random.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/gamma.hpp>

void gen_x_signed_exp(const double r, boost::random::mt19937& rng, int N, std::vector<cpp_dec_float_50>& vec)
{
    boost::random::uniform_real_distribution<> runif;
    boost::random::exponential_distribution<cpp_dec_float_50> rexp(r);

    for (int i = 0; i < N; ++i)
    {
        cpp_dec_float_50 x = rexp(rng);
        double u = runif(rng);
        if (u < 0.5)
        {
            x = -x;
        }
        vec.push_back(x);
    }
    std::sort(vec.begin(), vec.end());
}

void gen_x_unsigned_exp(const double r, boost::random::mt19937& rng, int N, std::vector<cpp_dec_float_50>& vec, bool srtd = true)
{
    boost::random::uniform_real_distribution<> runif;
    boost::random::exponential_distribution<cpp_dec_float_50> rexp(r);

    for (int i = 0; i < N; ++i)
    {
        cpp_dec_float_50 x = rexp(rng);
        vec.push_back(x);
    }
    if (srtd)
    {
        std::sort(vec.begin(), vec.end());
    }
}

void gen_x_unsigned_unif(boost::random::mt19937& rng, int N, std::vector<cpp_dec_float_50>& vec, bool srtd = true)
{
    boost::random::uniform_real_distribution<> runif;

    for (int i = 0; i < N; ++i)
    {
        cpp_dec_float_50 x = runif(rng);
        vec.push_back(x);
    }
    if (srtd)
    {
        std::sort(vec.begin(), vec.end());
    }
}

void gen_k_unsigned_exp(const double r, boost::random::mt19937& rng, int N, std::vector<unsigned long>& vec, bool srtd = true)
{
    boost::random::uniform_real_distribution<> runif;
    boost::random::exponential_distribution<double> rexp(r);

    for (int i = 0; i < N; ++i)
    {
        double x = rexp(rng);
        vec.push_back(static_cast<unsigned long>(x));
    }
    if (srtd)
    {
        std::sort(vec.begin(), vec.end());
    }
}

int main(int argc, const char* argv[])
{
    std::string fun(argv[1]);

    int N = 200;
    if (argc >= 4)
    {
        N = boost::lexical_cast<unsigned long>(argv[2]);
    }

    unsigned long S = 17;
    if (argc >= 5)
    {
        S = boost::lexical_cast<unsigned long>(argv[3]);
        std::cerr << S << endl;
    }
    boost::random::mt19937 rng(S);

    std::cout.precision(std::numeric_limits<cpp_dec_float_50>::digits10);

    if (fun == "beta")
    {
        std::vector<cpp_dec_float_50> A;
        gen_x_unsigned_exp(0.05, rng, N, A, false);
        std::vector<cpp_dec_float_50> B;
        gen_x_unsigned_exp(0.05, rng, N, B, false);

        std::vector<cpp_dec_float_50> Y;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 y = boost::math::beta(A[i], B[i]);
            Y.push_back(y);
        }

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << A[i] << ", " << B[i] << ", " << Y[i] << ", " << log(Y[i]) << "]" << endl;
        }
        return 0;
    }

    if (fun == "betaInt")
    {
        std::vector<unsigned long> A;
        gen_k_unsigned_exp(0.01, rng, N, A, false);
        std::vector<unsigned long> B;
        gen_k_unsigned_exp(0.01, rng, N, B, false);

        std::vector<cpp_dec_float_50> Y;
        for (int i = 0; i < N; ++i)
        {
	    A[i] += 1;
	    B[i] += 1;
            cpp_dec_float_50 y = boost::math::beta(A[i], B[i]);
            Y.push_back(y);
        }

        cout << "seed: " << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << A[i] << ", " << B[i] << ", " << Y[i] << ", " << log(Y[i]) << "]" << endl;
        }
        return 0;
    }

    if (fun == "ibetaInt")
    {
        std::vector<unsigned long> A;
        gen_k_unsigned_exp(0.01, rng, N, A, false);
        std::vector<unsigned long> B;
        gen_k_unsigned_exp(0.01, rng, N, B, false);

        std::vector<cpp_dec_float_50> X;
        gen_x_unsigned_unif(rng, N, X, false);

        std::vector<cpp_dec_float_50> Y;
        std::vector<cpp_dec_float_50> Z;
        boost::random::uniform_real_distribution<> runif;
        for (int i = 0; i < N; ++i)
        {
            ++A[i];
            ++B[i];
            cpp_dec_float_50 y = boost::math::ibeta(A[i], B[i], X[i]);
            cpp_dec_float_50 z = boost::math::ibetac(A[i], B[i], X[i]);
            Y.push_back(y);
            Z.push_back(z);
        }

        cout << "seed: " << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << A[i] << ", " << B[i] << ", " << X[i] << ", " << Y[i] << ", " << log(Y[i]) << ", "
                                                          << Z[i] << ", " << log(Z[i]) << "]" << endl;
        }
        return 0;
    }

    if (fun == "ibeta")
    {
        std::vector<cpp_dec_float_50> A;
        gen_x_unsigned_exp(0.05, rng, N, A, false);
        std::vector<cpp_dec_float_50> B;
        gen_x_unsigned_exp(0.05, rng, N, B, false);

        std::vector<cpp_dec_float_50> X;
        gen_x_unsigned_unif(rng, N, X, false);

        std::vector<cpp_dec_float_50> Y;
        std::vector<cpp_dec_float_50> Z;
        boost::random::uniform_real_distribution<> runif;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 y = boost::math::ibeta(A[i], B[i], X[i]);
            cpp_dec_float_50 z = boost::math::ibetac(A[i], B[i], X[i]);
            Y.push_back(y);
            Z.push_back(z);
        }

        cout << "seed: " << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << A[i] << ", " << B[i] << ", " << X[i] << ", " << Y[i] << ", " << log(Y[i]) << ", "
                                                          << Z[i] << ", " << log(Z[i]) << "]" << endl;
        }
        return 0;
    }

    if (fun == "choose")
    {
        std::vector<unsigned long> J;
        std::vector<unsigned long> K;
        gen_k_unsigned_exp(0.01, rng, N, J);
        gen_k_unsigned_exp(0.01, rng, N, K, false);

        std::vector<cpp_dec_float_50> C;
        std::vector<cpp_dec_float_50> D;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 c = boost::math::binomial_coefficient<cpp_dec_float_50>(J[i] + K[i], K[i]);
            cpp_dec_float_50 d = log(c);
            C.push_back(c);
            D.push_back(d);
        }

        cout << "seed: " << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << (J[i] + K[i]) << ", " << K[i] << ", " << C[i] << ", " << D[i] << "]" << endl;
        }
        return 0;
    }

    if (fun == "erf")
    {
        std::vector<cpp_dec_float_50> X;
        gen_x_signed_exp(0.2, rng, N, X);

        std::vector<cpp_dec_float_50> Y;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 y = erf(X[i]);
            Y.push_back(y);
        }
        cout << "seed: " << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << X[i] << ", " << Y[i] << "]" << endl;
        }
        return 0;
    }

    if (fun == "erfc")
    {
        std::vector<cpp_dec_float_50> X;
        gen_x_signed_exp(0.2, rng, N, X);

        std::vector<cpp_dec_float_50> Y;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 y = erfc(X[i]);
            Y.push_back(y);
        }

        cout << "seed: " << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << X[i] << ", " << Y[i] << ", " << log(Y[i]) << "]" << endl;
        }
        return 0;
    }

    if (fun == "gamma")
    {
        std::vector<cpp_dec_float_50> X;
        gen_x_unsigned_exp(0.05, rng, N, X);

        std::vector<cpp_dec_float_50> Y;
        std::vector<cpp_dec_float_50> Z;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 y = tgamma(X[i]);
            cpp_dec_float_50 z = lgamma(X[i]);
            Y.push_back(y);
            Z.push_back(z);
        }

        cout << "seed: " << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << X[i] << ", " << Y[i] << ", " << Z[i] << "]" << endl;
        }
        cout << "]" << endl;
        return 0;
    }

    if (fun == "igamma")
    {
        std::vector<cpp_dec_float_50> A;
        gen_x_unsigned_exp(0.1, rng, N, A);
        std::vector<cpp_dec_float_50> X;
        gen_x_unsigned_exp(0.1, rng, N, X, false);

        std::vector<cpp_dec_float_50> Y;
        std::vector<cpp_dec_float_50> Z;
        std::vector<cpp_dec_float_50> U;
        std::vector<cpp_dec_float_50> V;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 y = boost::math::gamma_p(A[i], X[i]);
            cpp_dec_float_50 z = boost::math::gamma_q(A[i], X[i]);
            Y.push_back(y);
            Z.push_back(z);
        }

        cout << "seed: " << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cout << "- [" << A[i] << ", " << X[i] << ", " << Y[i] << ", " << log(Y[i])
                                                          << ", " << Z[i] << ", " << log(Z[i]) << "]" << endl;
        }
        return 0;
    }

   return 1;
}
