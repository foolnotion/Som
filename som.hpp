#ifndef SOM_HPP
#define SOM_HPP

#include <ostream>
#include <vector>

#define NDEBUG
#define BOOST_DISABLE_ASSERTS // production mode
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/future.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/ref.hpp>

using namespace boost::numeric; // for convenience: boost::numeric::ublas::vector becomes ublas::vector

namespace som
{
namespace random
{
double double_range(double start, double end);
}

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

/*! \brief Wrapper around the euclidean norm from boost::numeric::ublas

    Takes two weight vectors and computes the euclidean distance between them.

    \param v1 The first vector
    \param v2 The second vector
    \returns The euclidean distance between the two vectors
  */
template <class T>
struct norm_2
{
    double operator() (const ublas::vector<T>& v1, const ublas::vector<T>& v2)
    {
        return ublas::norm_2(v1-v2);
    }
};

class node
{
    friend std::ostream& operator<< (std::ostream&, const som::node&);

public:
    node(const unsigned input_size);
    ublas::vector<double>& get_weights() { return weights_; }

private:
    ublas::vector<double> weights_;
};

typedef boost::shared_ptr<node> node_ptr;
typedef boost::function <double (const ublas::vector<double>& v1, const ublas::vector<double>& v2)> metric;

typedef std::vector<std::vector<som::node_ptr> > grid2;

class map
{
    friend std::ostream& operator<< (std::ostream&, const som::map&);

public:
    /* constructor */
    explicit map(unsigned map_size, unsigned sample_size, const som::metric& metric);

    /* getters and setters */
    unsigned size() const { return map_size_; } //!< Returns the map size.
    const som::metric& get_distance() const { return metric_; }
    /*! \brief Convenience operator for accessing the map's nodes

        \returns The som::node_ptr at position (i,j,k)
    */
    som::node_ptr operator()(unsigned i, unsigned j, unsigned k) const { return grid3_[i][j][k]; }

    /* application logic */
    som::point3 best_matching_unit(ublas::vector<double>& sample);
    void train();
    void load_samples(const std::vector<ublas::vector<double> >&);

private:
    unsigned map_size_; //!< Size of the self-organizing map
    unsigned sample_size_; //!< Size of the sample vectors
    unsigned n_samples_; //!< Number of input samples

    boost::array<ublas::vector<double>, 0> input_samples_; //!< The input samples
    boost::multi_array<som::node_ptr, 3> grid3_; //!< 3d grid of nodes
    som::metric metric_; //!< Distance function (boost::function)
};

struct thread_result
{
    double min_distance;
    som::point3 position;
};

} // namespace som

#endif // SOM_HPP
