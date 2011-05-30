#include "som.hpp"

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/detail/iterator.hpp>

/* Using a global generator (we don't want to create a new one at every call */
boost::mt19937 twister;

/* Global threadpool (too lazy to pass around a global context, and singletons are evil!) */
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
som::map::map(unsigned map_size, unsigned sample_size, const boost::shared_ptr<som::abstract_distance>& distance)
    : map_size_(map_size),
      sample_size_(sample_size),
      distance_(distance)
{
    tp_.size_controller ().resize (boost::thread::hardware_concurrency ());
    distance_function_ = som::euclidean_distance1();
    nodes_ = ublas::matrix<som::node_ptr>(map_size_, map_size_);

    for (unsigned i = 0; i != map_size_; ++i)
    {
        for (unsigned j = 0; j != map_size_; ++j)
        {
            nodes_(i,j) = (boost::make_shared<som::node>(sample_size_));
        }
    }
}

void
som::map::load_samples(const std::vector<ublas::vector<double> >& samples)
{
    input_samples_ = samples;
    std::cout << "Loaded " << input_samples_.size() << " input samples." << std::endl;
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

    Returns the position of the closest node to the sample vector (in terms of distance, depending on what metric the map uses).
    By default the euclidean distance is used.

    \param sample The sample vector
    \returns A som::position representing the matrix indexes of the closest node
    \sa abstract_distance, euclidean_distance
  */
som::position
som::map::get_bmu(ublas::vector<double> &sample)
{
    som::position p; // use this trick to avoid assignment of smart pointers (expensive copy)
    double min_distance = 2.0;
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

som::position
som::map::get_bmu1(ublas::vector<double> &sample)
{
    som::position p; // use this trick to avoid assignment of smart pointers (expensive copy)
    double min_distance = 2.0;
    for (unsigned i = 0; i != map_size_; ++i)
    {
        for (unsigned j = 0; j != map_size_; ++j)
        {
            double distance = distance_function_(nodes_(i,j)->get_weights(), sample);
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

som::thread_result
calculate_distance(som::map& map, ublas::vector<double>& sample, unsigned offset, unsigned length)
{
    som::distance_function distance_function = map.get_distance();
    som::thread_result r;
    r.min_distance = 2.0;
    for (unsigned i = offset; i != offset+length; ++i)
    {
        for (unsigned j = 0; j != map.get_size(); ++j)
        {
            double distance = distance_function(map(i,j)->get_weights(), sample);
            if (r.min_distance > distance)
            {
                r.min_distance = distance;
                r.position.x = i;
                r.position.y = j;
            }
        }
    }
    std::cout << r.position << std::endl;
    return r;
}

som::position
som::map::get_bmu2(ublas::vector<double>& sample)
{
    boost::array<boost::unique_future<som::thread_result>, 2> futures;
    unsigned length = ceil(map_size_ / 4.0);
    for (unsigned i = 0; i != 4; ++i)
    {
        unsigned offset = length * i;
        std::cout << "offset: " << offset << std::endl;
        if (offset + length > map_size_)
        {
            length = map_size_ - offset;
            std::cout << "length: " << length << std::endl;
        }
        boost::packaged_task<som::thread_result> pt(boost::bind(&calculate_distance, boost::ref(*this), boost::ref(sample), offset, length));
        futures.at(i) = pt.get_future();
        boost::thread t(boost::move(pt));
    }
    double min_distance = 2.0; // the euclidean norm between v1=(1,1,1) and v2=(1,1,1) is ||v1,v2||=1.732 < 2.0
    som::position p;
    for (boost::array<boost::unique_future<som::thread_result>, 2>::size_type i = 0; i != futures.size(); ++i)
    {
        futures.at(i).wait();
        double d = futures.at(i).get().min_distance;
        if (min_distance > d)
        {
            min_distance = d;
            p = futures.at(i).get().position;
        }
    }
    return p;
}

void
calculate_distance1(som::map& map, ublas::vector<double>& sample, unsigned offset, unsigned length, som::thread_result& r)
{
    som::distance_function distance_function = map.get_distance();
    r.min_distance = 2.0;
    for (unsigned i = offset; i != offset+length; ++i)
    {
        for (unsigned j = 0; j < map.get_size(); ++j)
        {
            double distance = distance_function(map(i,j)->get_weights(), sample);
            if (r.min_distance > distance)
            {
                r.min_distance = distance;
                r.position.x = i;
                r.position.y = j;
            }
        }
    }
    std::cout << r.position << std::endl;
}

som::position
som::map::get_bmu3(ublas::vector<double>& sample)
{
    boost::thread_group tg;
    std::vector<som::thread_result> results(2);
    unsigned length = ceil(map_size_ / 2.0);
    for (unsigned i = 0; i != 2; ++i)
    {
        unsigned offset = length * i;
        if (offset + length > map_size_)
        {
            length = map_size_ - offset;
        }
        tg.create_thread(boost::bind(&calculate_distance1, boost::ref(*this), boost::ref(sample), offset, length, boost::ref(results[i])));
    }
    tg.join_all();
    double min_distance = 2.0; // the euclidean norm between v1=(1,1,1) and v2=(1,1,1) is ||v1,v2||=1.732
    som::position p;
    for (std::vector<som::thread_result>::size_type i = 0; i != results.size(); ++i)
    {
        double d = results[i].min_distance;
        if (min_distance > d)
        {
            min_distance = d;
            p = results[i].position;
        }
    }
    return p;
}



