#include "som.hpp"

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
    for (unsigned i = 0; i != map.nodes_.size1(); ++i)
    {
        for (unsigned j = 0; j != map.nodes_.size2(); ++j)
        {
            os << *map.nodes_(i,j);
        }
        os << std::endl;
    }
    return os;
}

double som::random(double start, double end)
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
        weights_.insert_element(i, som::random(0.f,1.f));
    }
}

/*! \brief Map constructor
    \param map_size Matrix size (the map has all the nodes organized as a matrix)
    \param sample_size Size of the sample vector
    \param distance A reference to the distance functor used to determine the best mathing unit and its neighbours
    \sa abstract_distance, euclidean_distance
 */
som::map::map(u_int size, u_int input_size, const boost::shared_ptr<som::abstract_distance>& distance)
    : distance_(distance)
{
    sample_size_ = input_size;
    map_size_ = size;
    nodes_ = boost::numeric::ublas::matrix<som::node_ptr>(size, size);

    for (unsigned i = 0; i != map_size_; ++i)
    {
        for (unsigned j = 0; j != map_size_; ++j)
        {
            nodes_(i,j) = (boost::make_shared<som::node>(sample_size_));
        }
    }
}

/*! \brief Convenience operator for accessing the map's nodes

    \param i The line of the nodes matrix
    \param j The column of the nodes matrix
    \returns The som::node_ptr at position (i,j)
*/
som::node_ptr
som::map::operator() (unsigned i, unsigned j)
{
    return nodes_(i,j);
}

/*! \brief Get best matching unit

    Returns the position of the closest node to the sample vector (in terms of distance, depending on what metric the map uses). By default the euclidean distance is used.

    \param sample The sample vector
    \returns A som::position representing the matrix indexes of the closest node
    \sa abstract_distance, euclidean_distance
  */
som::position
som::map::get_bmu(boost::numeric::ublas::vector<double> &sample)
{
    som::position p; // use this trick to avoid assignment of smart pointers (expensive copy)
    double min_distance = (*distance_)(nodes_(0,0)->get_weights(), sample);
    for (unsigned i = 0; i != map_size_; ++i)
    {
        for (unsigned j = 0; j != map_size_; ++j)
        {
            double distance = (*distance_)(nodes_(i,j)->get_weights(), sample);
            if (min_distance > distance)
            {
                min_distance = distance;
                p.x = i;
                p.y = j;
            }
        }
    }
    return p;
}

