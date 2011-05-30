#ifndef SOM_HPP
#define SOM_HPP

#include <ostream>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/threadpool.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/ref.hpp>

using namespace boost::numeric; // for convenience: boost::numeric::ublas::vector becomes ublas::vector

namespace som
{
namespace random
{
double double_range(double start, double end);
}

struct position
{
    friend std::ostream& operator<< (std::ostream& os, const som::position& p)
    {
        os << "(" << p.x << "," << p.y << ")";
        return os;
    }

    unsigned x;
    unsigned y;
};

struct abstract_distance
{
    virtual double
    operator() (const ublas::vector<double>& v1, const ublas::vector<double>& v2) =0;
};

struct euclidean_distance : abstract_distance
{
    double operator() (const ublas::vector<double>& v1, const ublas::vector<double>& v2)
    {
        return ublas::norm_2(v1-v2);
    }
};

struct euclidean_distance1
{
    double operator() (const ublas::vector<double>& v1, const ublas::vector<double>& v2)
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
    boost::shared_ptr<abstract_distance> distance_;
};

typedef boost::shared_ptr<node> node_ptr;
typedef boost::function <double (const ublas::vector<double>& v1,
                                 const ublas::vector<double>& v2)> distance_function;

class map
{
    friend std::ostream& operator<< (std::ostream&, const som::map&);

public:
    map(unsigned map_size, unsigned sample_size, const boost::shared_ptr<som::abstract_distance>& distance);
    unsigned get_size() { return map_size_; } //!< Returns the map size.
    const distance_function& get_distance() const { return distance_function_; }
    const boost::threadpool::pool& threadpool() const { return tp_; }
    som::node_ptr operator() (unsigned i, unsigned j);
    som::position get_bmu(ublas::vector<double>& sample);
    som::position get_bmu1(ublas::vector<double>& sample);
    som::position get_bmu2(ublas::vector<double>& sample);
    som::position get_bmu3(ublas::vector<double>& sample);
    void train();
    void load_samples(const std::vector<ublas::vector<double> >&);

private:
    unsigned map_size_; //!< Size of the self-organizing map
    unsigned sample_size_; //!< Size of the sample vectors
    unsigned n_samples_; //!< Number of input samples
    std::vector<ublas::vector<double> > input_samples_; //!< The input samples
    ublas::matrix<som::node_ptr> nodes_; //!< The node (som::node_ptr) matrix.
    boost::shared_ptr<som::abstract_distance> distance_; //!< The distance functor
    distance_function distance_function_; //!< Distance function (boost::function)
    /* threadpool and workers */
    boost::threadpool::pool tp_;
    boost::mutex mutex_;
};

struct thread_result
{
    double min_distance;
    som::position position;
};

} // namespace som

#endif // SOM_HPP
