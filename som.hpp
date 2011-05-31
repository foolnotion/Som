#ifndef SOM_HPP
#define SOM_HPP

#include <vector>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

#define NDEBUG
#define BOOST_DISABLE_ASSERTS

#include <boost/function.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

using namespace boost::numeric; // for convenience: boost::numeric::ublas:: becomes ublas::

/* Using a global generator (we don't want to create a new one at every call */
namespace // anonymous namespace - no one shall get their filthy hands on the twister!
{
boost::mt19937 twister;
}

namespace som
{
namespace random
{
double double_range(double start, double end)
{
    boost::uniform_real<> dist(start, end);
    return dist(twister);
}
} // namespace random

namespace util
{
namespace
{
/*! \brief Split a space-separated string into tokens.

    Streams in C++ have the special ability of reading until a whitespace. The stringstream
    will use the output operator (>>) and put a string into \variable buf everytime a
    whitespace is met; buf is then used to push_back() into the vector.
    The vector tokens will contain all the words in str.

    @param str The string to be split */
std::vector<std::string>
split_by_spaces(std::string& str)
{
    std::string buf;
    std::stringstream ss(str);

    std::vector<std::string> tokens;

    while (ss >> buf)
    {
        tokens.push_back(buf);
    }
    return tokens;
}
} // anonymous namespace

typedef std::vector<std::vector<double> > samples;

std::vector<std::vector<double> >
load_samples_from_file(const char* filename)
{
    std::vector<std::vector<double> > samples;
    std::string line;
    std::ifstream istr(filename);
    while (istr.good())
    {
        std::getline(istr,line);
        std::vector<std::string> tokens = split_by_spaces(line);

        if (tokens.size() > 0)
        {
            std::vector<double> v;
            for (std::vector<std::string>::size_type i = 0; i != tokens.size(); ++i)
            {
                std::istringstream iss(tokens[i]);
                double val;
                iss >> val;
                v.push_back(val);
            }
            samples.push_back(v);
        }
    }
    std::cout << "Samples size: " << samples.size() << std::endl;
    return samples;
}
} // namespace util

class point3
{
    friend std::ostream& operator<< (std::ostream& os, const som::point3& p)
    {
        os << "(" << p.x << "," << p.y << "," << p.z << ")";
        return os;
    }
public:
    point3() : x(0), y(0), z(0) {}
    point3(unsigned x, unsigned y, unsigned z) : x(x), y(y), z(z) {}
    unsigned x;
    unsigned y;
    unsigned z;
};

template <std::size_t N, std::size_t S>
class map
{
    friend std::ostream& operator<< (std::ostream& os, const som::map<N,S>& map)
    {
        for (unsigned i = 0; i != N; ++i)
        {
            for (unsigned j = 0; j != N; ++j)
            {
                for (unsigned k = 0; j != N; ++j)
                {
                    os << map(i,j,k);
                }
            }
            os << std::endl;
        }
        return os;
    }

    typedef typename boost::multi_array<ublas::bounded_vector<double,S>, 3>::index index;

public:
    /*! \brief Map constructor

        For now the map uses the euclidean distance (boost::numeric::ublas::norm_2)

        \todo Enable usage of custom metrics

        \tparam N Size (number of elements) for one dimension
        \tparam S Size of the sample vector
     */
    map<N,S>()
    {
        grid3_.resize(boost::extents[N][N][N]);

        std::cout << "Map dimensions: " << N << "x" << N << "x" << N << " (" << N*N*N << " elements)" << std::endl;

        for (index x = 0; x != N; ++x)
            for (index y = 0; y != N; ++y)
                for (index z = 0; z != N; ++z)
                {
                    for (unsigned i = 0; i != S; ++i)
                    {
                        grid3_[x][y][z](i) = som::random::double_range(0.f,1.f);
                    }
                }
    }

    /*! \brief Get best matching unit

        Returns the position (in 3d coordinates) of the closest vector (in terms of distance to the sample vector,
        depending on what metric the map uses).

        By default the euclidean distance is used:
        \f[ d(p,q)=\sqrt{(p_1-q_1)^2+(p_2-q_2)^2+...+(p_i-q_i)^2+...+(p_n-q_n)^2}
        \f]

        \param sample The sample vector
        \returns A som::point3 representing the coordinates of the closest node
     */
    som::point3
    best_matching_unit(const ublas::bounded_vector<double,S>& sample)
    {
        som::point3 p;
        double min_distance = 2.0;
        for (index x = 0; x != N; ++x)
            for (index y = 0; y != N; ++y)
                for (index z = 0; z != N; ++z)
                {
                    double distance = ublas::norm_2(grid3_[x][y][z]-sample);
                    if (min_distance > distance)
                    {
                        min_distance = distance;
                        p.x = x;
                        p.y = y;
                        p.z = z;
                    }
                }
        return p;
    }

    /*! \brief Load input samples from a som::util::sample to the map's internal sample vector

        This function basically transforms each \c std::vector (corresponding to an input sample read from file) into an \c ublas::bounded_vector. The dimensions must match.

        \todo Work out a way to make sure the dimensions will match.

        \param samples The input sample vector
     */
    void load_samples(const som::util::samples& samples)
    {
        samples_.resize(samples.size());
        for (typename std::vector<ublas::bounded_vector<double,S> >::size_type i = 0; i != samples_.size(); ++i)
            for (unsigned j = 0; j != samples[i].size(); ++j)
            {
                samples_[i](j) = samples[i][j];
            }
    }

    /*! \brief Convenience operator for accessing the map's elements

        \returns The element at position \f$(i,j,k)\f$
    */
    const ublas::bounded_vector<double,S>& operator()(unsigned i, unsigned j, unsigned k) const { return grid3_[i][j][k]; }
    std::vector<ublas::bounded_vector<double,S> > samples_; //!< The input samples
private:
    boost::multi_array<ublas::bounded_vector<double,S>, 3> grid3_; //!< 3d grid of nodes
//    boost::function<double (ublas::vector_expression<double>& e)> norm_;
};

} // namespace som

#endif // SOM_HPP
