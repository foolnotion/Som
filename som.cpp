#include "som.hpp"

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/detail/iterator.hpp>

namespace
{
/* Using a global generator (we don't want to create a new one at every call
   Making it anonymous so no one would touch it with their filthy hands */
boost::mt19937 twister;
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
                os << map(i,j,k);
            }
        }
        os << std::endl;
    }
    return os;
}

double
som::random::double_range(double start, double end)
{
    boost::uniform_real<> dist(start, end);
    return dist(twister);
}

/*! \brief Map constructor
    \param map_size Matrix size (the map has all the nodes organized as a matrix)
    \param sample_size Size of the sample vector
    \param distance A reference to the distance functor used to determine the best mathing unit and its neighbours
    \sa abstract_distance, euclidean_distance
 */
som::map::map(unsigned map_size, unsigned sample_size) : map_size_(map_size), sample_size_(sample_size)
{
    typedef boost::multi_array<ublas::vector<double>, 3>::index index;

    grid3_.resize(boost::extents[map_size_][map_size_][map_size_]);

    for (index x = 0; x != map_size_; ++x)
        for (index y = 0; y != map_size_; ++y)
            for (index z = 0; z != map_size_; ++z)
            {
                grid3_[x][y][z].resize(sample_size);
                for (unsigned i = 0; i != sample_size; ++i)
                {
                    grid3_[x][y][z](i) = som::random::double_range(0.f, 1.f);
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
som::map::best_matching_unit(const ublas::vector<double>& sample)
{
    typedef boost::multi_array<ublas::vector<double>, 3>::index index;
    som::point3 p;
    double min_distance = 2.0;
    for (index i = 0; i != map_size_; ++i)
    {
        for (index j = 0; j != map_size_; ++j)
        {
            for (index k = 0; k != map_size_; ++k)
            {
                double distance = ublas::norm_2(grid3_[i][j][k] - sample);
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
