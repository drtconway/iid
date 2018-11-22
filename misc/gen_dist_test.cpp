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

#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/students_t.hpp>

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
    boost::random::uniform_real_distribution<> runif;

    std::cout.precision(std::numeric_limits<cpp_dec_float_50>::digits10);

    if (fun == "beta")
    {
        boost::random::exponential_distribution<double> rexp(0.01);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            unsigned long a = 1 + static_cast<unsigned long>(rexp(rng));
            unsigned long b = 1 + static_cast<unsigned long>(rexp(rng));
            cpp_dec_float_50 x = runif(rng);
            boost::math::beta_distribution<cpp_dec_float_50> dst(a, b);
            cout << "- [" << a << ", " << b << ", " << x
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, x)
                               << ", " << cdf(dst, x)
                               << ", " << cdf(complement(dst, x))
                               << "]" << endl;
        }
        return 0;
    }

    if (fun == "binom")
    {
        boost::random::exponential_distribution<double> rexp(0.005);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 p = runif(rng);
            unsigned long n = 1 + static_cast<unsigned long>(rexp(rng));
            unsigned long k = static_cast<unsigned long>((n+1)*runif(rng));
            boost::math::binomial_distribution<cpp_dec_float_50> dst(n, p);
            cout << "- [" << n << ", " << p << ", " << k
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, k)
                               << ", " << cdf(dst, k)
                               << ", " << cdf(complement(dst, k))
                               << "]" << endl;
        }
        return 0;
    }

    if (fun == "chisq")
    {
        boost::random::exponential_distribution<double> rexp1(0.05);
        boost::random::exponential_distribution<double> rexp2(0.01);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            unsigned long n = 1 + static_cast<unsigned long>(rexp1(rng));
            cpp_dec_float_50 x = rexp2(rng);
            boost::math::chi_squared_distribution<cpp_dec_float_50> dst(n);
            cout << "- [" << n << ", " << x
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, x)
                               << ", " << cdf(dst, x)
                               << ", " << cdf(complement(dst, x))
                               << "]" << endl;
        }
        return 0;
    }

    if (fun == "gamma")
    {
        boost::random::exponential_distribution<double> rexp1(0.05);
        boost::random::exponential_distribution<double> rexp2(0.01);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 a = 1 + static_cast<unsigned long>(rexp1(rng));
            cpp_dec_float_50 x = rexp2(rng);
            boost::math::gamma_distribution<cpp_dec_float_50> dst(a);
            cout << "- [" << a << ", " << x
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, x)
                               << ", " << cdf(dst, x)
                               << ", " << cdf(complement(dst, x))
                               << "]" << endl;
        }
        return 0;
    }

    if (fun == "geom")
    {
        boost::random::exponential_distribution<double> rexp(0.02);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 p = runif(rng);
            unsigned long k = static_cast<unsigned long>(rexp(rng));
            boost::math::geometric_distribution<cpp_dec_float_50> dst(p);
            cout << "- [" << p << ", " << k
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, k)
                               << ", " << cdf(dst, k)
                               << ", " << cdf(complement(dst, k))
                               << "]" << endl;
        }
        return 0;
    }

    if (fun == "hyper")
    {
        boost::random::exponential_distribution<double> rexp(0.02);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            unsigned long M = 1 + static_cast<unsigned long>(rexp(rng));
            unsigned long K = 1 + static_cast<unsigned long>(rexp(rng));
            unsigned long N = M + K;
            unsigned long n = static_cast<unsigned long>((N+1)*runif(rng));
            unsigned long k = static_cast<unsigned long>((n+1)*runif(rng));
            while (k > K || n - k > M) {
                k = static_cast<unsigned long>((n+1)*runif(rng));
            }
            boost::math::hypergeometric_distribution<cpp_dec_float_50> dst(K, n, N);
            cout << "- [" << N << ", " << K << ", " << n << ", " << k
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, k)
                               << ", " << cdf(dst, k)
                               << ", " << cdf(complement(dst, k))
                               << "]" << endl;
        }
        return 0;
    }

    if (fun == "norm")
    {
        boost::random::exponential_distribution<double> rexp(0.02);
        boost::random::exponential_distribution<double> rexp2(0.005);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 mu = rexp(rng);
            if (runif(rng) < 0.5) {
                mu = -mu;
            }
            cpp_dec_float_50 sig = rexp(rng);
            cpp_dec_float_50 x = rexp2(rng);
            if (runif(rng) < 0.5) {
                x = -x;
            }
            boost::math::normal_distribution<cpp_dec_float_50> dst(mu, sig);
            cout << "- [" << mu << ", " << sig << ", " << x
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, x)
                               << ", " << cdf(dst, x)
                               << ", " << cdf(complement(dst, x))
                               << "]" << endl;
        }
        return 0;
    }

    if (fun == "pois")
    {
        boost::random::exponential_distribution<double> rexp(0.02);
        boost::random::exponential_distribution<double> rexp2(0.005);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 lam = rexp(rng);
            unsigned long k = static_cast<unsigned int>(rexp(rng));
            if (runif(rng) < 0.5)
            {
                k = static_cast<unsigned int>(rexp2(rng));
            }
            boost::math::poisson_distribution<cpp_dec_float_50> dst(lam);
            cout << "- [" << lam << ", " << k
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, k)
                               << ", " << cdf(dst, k)
                               << ", " << cdf(complement(dst, k))
                               << "]" << endl;
        }
        return 0;
    }

    if (fun == "stud")
    {
        boost::random::exponential_distribution<double> rexp(0.075);
        boost::random::exponential_distribution<double> rexp2(0.005);

        cout << "seed:" << S << endl;
        cout << "data:" << endl;
        for (int i = 0; i < N; ++i)
        {
            cpp_dec_float_50 nu = 2 + rexp(rng);
            cpp_dec_float_50 x = static_cast<unsigned int>(rexp2(rng));
            if (runif(rng) < 0.5)
            {
                x = -x;
            }
            boost::math::students_t_distribution<cpp_dec_float_50> dst(nu);
            cout << "- [" << nu << ", " << x
                               << ", " << mean(dst)
                               << ", " << variance(dst)
                               << ", " << pdf(dst, x)
                               << ", " << cdf(dst, x)
                               << ", " << cdf(complement(dst, x))
                               << "]" << endl;
        }
        return 0;
    }

   return 1;
}
