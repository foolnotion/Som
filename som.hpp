#ifndef SOM_HPP
#define SOM_HPP

#include <ostream>
#include <vector>

#define NDEBUG
#define BOOST_DISABLE_ASSERTS
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
using namespace boost::numeric; // for convenience: boost::numeric::ublas::vector becomes ublas::vector

namespace som
{
namespace random
{
    inline double double_range(double start, double end);
} // namespace random

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

typedef boost::function <double (const ublas::vector<double>& v1, const ublas::vector<double>& v2)> metric;

class map
{
    friend std::ostream& operator<< (std::ostream&, const som::map&);

public:
    /* constructor */
    explicit map(unsigned map_size, unsigned sample_size);

    /* getters and setters */
    unsigned size() const { return map_size_; } //!< Returns the map size.
    /*! \brief Convenience operator for accessing the map's nodes

        \returns The som::node_ptr at position (i,j,k)
    */
    const ublas::vector<double>& operator()(unsigned i, unsigned j, unsigned k) const { return grid3_[i][j][k]; }

    /* application logic */
    som::point3 best_matching_unit(const ublas::vector<double>& sample);
    void train();
    void load_samples(const std::vector<ublas::vector<double> >&);

private:
    unsigned map_size_; //!< Size of the self-organizing map
    unsigned sample_size_; //!< Size of the sample vectors
    unsigned n_samples_; //!< Number of input samples

    boost::array<ublas::vector<double>, 0> input_samples_; //!< The input samples
    boost::multi_array<ublas::vector<double>, 3> grid3_; //!< 3d grid of nodes
};

} // namespace som

#endif // SOM_HPP
