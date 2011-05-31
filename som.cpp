#include "som.hpp"

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/detail/iterator.hpp>


/* Using a global generator (we don't want to create a new one at every call */
boost::mt19937 twister;

std::ostream&
som::operator<< (std::ostream& os, const som::node &n)
{
    os << n.weights_;
    std::cout << std::endl;
    return os;
}

std::ostream&
som::operator<< (std::ostream& os, const som::map& map)
{    
    const unsigned map_size = map.size();
    for (unsigned i = 0; i != map_size; ++i)
    {
        for (unsigned j = 0; j != map_size; ++j)
        {
            for (unsigned k = 0; j != map_size; ++j)
            {
                os << *map(i,j,k);
            }
        }
        os << std::endl;
    }
    return os;
}

double som::random::double_range(double start, double end)
{
    boost::uniform_real<> dist(start, end);
    return dist(twister);
}

/*! \brief Node constructor

  Initializes a map node with random weights
  \param sample_size Size of the sample vector, which is also the size of the node's internal weights vector
  */
som::node::node(const unsigned sample_size)
{
    weights_.resize(sample_size);
    for (unsigned i = 0; i != sample_size; ++i)
    {
        weights_.insert_element(i, som::random::double_range(0.f,1.f));
    }
}

/*! \brief Map constructor
    \param map_size Matrix size (the map has all the nodes organized as a matrix)
    \param sample_size Size of the sample vector
    \param distance A reference to the distance functor used to determine the best mathing unit and its neighbours
    \sa abstract_distance, euclidean_distance
 */
som::map::map(unsigned map_size, unsigned sample_size, const som::metric& metric)
    : map_size_(map_size),
      sample_size_(sample_size),
      metric_(metric)
{
    typedef boost::multi_array<som::node_ptr, 3>::index index;

    metric_ = som::norm_2<double>(); // TODO: initialize this properly

    grid3_.resize(boost::extents[map_size_][map_size_][map_size_]);

    for (index i = 0; i != map_size_; ++i)
    {
        for (index j = 0; j != map_size_; ++j)
        {
            for (index k = 0; k != map_size; ++k)
            {
                grid3_[i][j][k] = boost::make_shared<som::node>(sample_size_);
            }
        }
    }
}

/*! \brief Get best matching unit

    Returns the position of the closest node to the sample vector (in terms of distance, depending on what metric the map uses).
    By default the euclidean distance is used.

    \param sample The sample vector
    \returns A som::position representing the matrix indexes of the closest node
    \sa abstract_distance, euclidean_distance
  */
som::point3
som::map::best_matching_unit(ublas::vector<double>& sample)
{
    typedef boost::multi_array<som::node_ptr, 3>::index index;
    som::point3 p;
    double min_distance = 2.0;
    for (index i = 0; i != map_size_; ++i)
    {
        for (index j = 0; j != map_size_; ++j)
        {
            for (index k = 0; k != map_size_; ++k)
            {
                double distance = metric_(grid3_[i][j][k]->get_weights(), sample);
                if (min_distance > distance)
                {
                    min_distance = distance;
                    p.x = i;
                    p.y = j;
                    p.z = k;
                }
            }
        }
    }
    return p;
}
