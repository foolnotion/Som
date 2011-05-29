#ifndef SOM_HPP
#define SOM_HPP

#include <ostream>
#include <vector>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/math/special_functions/hypot.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

namespace som
{
/* This function takes a closed range */
double random(double start, double end);

struct position
{
    unsigned x;
    unsigned y;
};

struct abstract_distance
{
    virtual double
    operator() (const boost::numeric::ublas::vector<double>& v1, const boost::numeric::ublas::vector<double>& v2) =0;
};

struct euclidean_distance : abstract_distance
{
    double operator() (const boost::numeric::ublas::vector<double>& v1, const boost::numeric::ublas::vector<double>& v2)
    {
        return boost::numeric::ublas::norm_2(v1-v2);
    }
};

class node
{
    friend std::ostream& operator<< (std::ostream&, const som::node&);

public:
    node(const unsigned input_size);
    boost::numeric::ublas::vector<double>& get_weights() { return weights_; }

private:
    boost::numeric::ublas::vector<double> weights_;
    boost::shared_ptr<abstract_distance> distance_;
};

typedef boost::shared_ptr<node> node_ptr;

class map
{
    friend std::ostream& operator<< (std::ostream&, const som::map&);

public:
    map(unsigned map_size, unsigned sample_size, const boost::shared_ptr<som::abstract_distance>& distance);
    unsigned get_size() { return map_size_; } //!< Returns the map size.
    som::node_ptr operator() (unsigned i, unsigned j);
    som::position get_bmu(boost::numeric::ublas::vector<double>& sample);
    void train();

private:
    std::vector<boost::numeric::ublas::vector<double> > input_samples_; //!< The input samples
    boost::shared_ptr<som::abstract_distance> distance_; //!< The distance functor
    boost::numeric::ublas::matrix<som::node_ptr> nodes_; //!< The node (som::node_ptr) matrix.
    unsigned map_size_; //!< Size of the self-organizing map
    unsigned sample_size_; //!< Size of the sample vectors
    unsigned n_samples_; //!< Number of input samples
};

} // namespace som

#endif // SOM_HPP
