#ifndef SOM_HPP
#define SOM_HPP

/* work around some libpng errors */
#define png_infopp_NULL (png_infopp)NULL
#define int_p_NULL (int*)NULL

/* standard stuff */
#include <vector>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cmath>

/* vector operations */
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

/* array and multi_array storage types */
#include <boost/array.hpp>
#include <boost/multi_array.hpp>

/* random numbers */
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>

/* misc: save images, show progress */
#include <boost/progress.hpp>
#include <boost/gil/extension/io/png_io.hpp>

using namespace boost::numeric;

/* Using a global generator (we don't want to create a new one at every call */
namespace { boost::mt19937 twister; } // anonymous namespace - no one shall get their filthy hands on the twister!

namespace som {
/** \brief Random numbers uniformly distributed on the given interval.

    Underlying generator: Mersenne Twister (\c boost::mt19937)
 */
namespace random {
/** \brief Returns a double value \f$v\in[begin,end]\f$ (note the \e closed interval). */
double
double_range(double begin, double end)
{
    boost::uniform_real<> dist(begin, end);
    return dist(twister);
}

/** \brief Returns an integer value \f$v\in[begin,end]\f$ (again, \e closed interval). */
int
int_range(int begin, int end)
{
    boost::uniform_int<> dist(begin, end);
    return dist(twister);
}
} // namespace random

/** \brief Helper functions (read from file, etc) */
namespace util { namespace {
/*! \brief Split a space-separated string into tokens.

    Streams in C++ have the special ability of reading until a whitespace. The stringstream
    will use the output operator (>>) and put a string into \variable buf everytime a
    whitespace is met; buf is then used to push_back() into the vector.
    The vector tokens will contain all the words in str.

    \param str The string to split
 */
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
    return samples;
}
} // namespace util

/** \brief 3d point data type */
template <class T>
class vector3
{
    friend std::ostream& operator<< (std::ostream& os, const som::vector3<T>& p)
    {
        os << "(" << p.x << "," << p.y << "," << p.z << ")";
        return os;
    }
public:
    vector3() : x(0), y(0), z(0) {}
    vector3(T x, T y, T z) : x(x), y(y), z(z) {}
    T x; //! X coordinate
    T y; //! Y coordinate
    T z; //! Z coordinate
};

/** \brief The self-organizing map

    Holds its internal state and can self-organize by scaling its units to become
    more similar to the input samples (also called activity vectors).

    At each step k of the learning algorithm, the map learns an activity vector \f${AV}_k\f$ as follows:
        -# The best matching unit (BMU) is found, so that \f$d(BMU,{AV}_k)\f$ is minimal (where \f$d\f$ is the distance used)
        -# The BMU and its neighbours are changed to be more similar to \f${AV}_k\f$.\n\n
           The learning rate and the neighbourhood radius are given by two monotonically-decreasing functions:\n\n
           \f[ \eqalignno
               { L(k) &= L_0 \cdot \exp \left( { {-k} \cdot \frac{\ln \frac{L_0}{L_M} }{M} } \right)\ -\ \text{
learning rate at step }\bf k & (1) \cr
                 R(k) &= \left \lfloor R_0 \cdot \exp \left( -k \cdot \frac{\ln \left(2 \cdot R_0 \right)}{\frac{g}{100} \cdot M} \right)  \right \rceil\ -\ \text{radius at step }\bf k\ \textrm{(neighbourhood size around the BMU)} & (2)
               }
           \f]
           Where:
                - \f$M\f$ is the total number of training steps
                - \f$L_0\f$ and \f$L_M\f$ are the initial and final learning rates
                - \f$R_0\f$ is the initial radius of the neighbourhood and \f$g\f$ is the percentage of \f$M\f$
                  after which \f$R\f$ becomes 0 and only the BMU is modified.
                - \f$\lfloor x \rceil\f$ denotes the nearest integer function of real number \f$x\f$.
                .
\n
           Two constant variables can be extracted from the above relations, to ease the computation (hint: the logarithm parts).
  */
template < int N, int S >
class map
{
    friend std::ostream& operator<< (std::ostream& os, const som::map<N,3>& map)
    {
        for (int i = 0; i != N; ++i)
        {
            for (int j = 0; j != N; ++j)
            {
                for (int k = 0; j != N; ++j)
                {
                    os << map(i,j,k);
                }
            }
            os << std::endl;
        }
        return os;
    }

    typedef typename boost::multi_array<ublas::c_vector<double,S>, 3>::index index;

private:
    boost::multi_array<ublas::c_vector<double,S>, 3> grid3_; //!< 3d grid of nodes
    std::vector<ublas::c_vector<double,S> > samples_; //!< The input samples

    /* internal state - the following variables control the learning process */
    int step_;
    int learn_radius_;
    double learn_rate_;

    const int total_steps_;
    const int initial_learn_radius_;
    const double radius_gone_factor_;
    const double initial_learn_rate_;
    const double final_learn_rate_;

    double radius_decay_factor_;
    double learn_rate_decay_factor_;

    bool is_initialized_;

public:
    /** \brief Map constructor

        For now the map uses the euclidean distance (boost::numeric::ublas::norm_2)

        \todo Enable usage of custom metrics

        \tparam N Size (number of elements) for one dimension
        \tparam S Size of the sample vector
     */
    map<N,S>(int total_steps, double init_radius, double radius_gone_factor,
             double init_learn_rate, double final_learn_rate)
        : step_(0),
          total_steps_(total_steps),
          initial_learn_radius_(init_radius),
          radius_gone_factor_(radius_gone_factor),
          initial_learn_rate_(init_learn_rate),
          final_learn_rate_(final_learn_rate),
          radius_decay_factor_(0),
          learn_rate_decay_factor_(0),
          is_initialized_(false)
    {
    }

    /** \brief Initialization function

        The reasons this function exists are explained here:
        http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml#Doing_Work_in_Constructors
      */
    void init()
    {
        radius_decay_factor_ = std::exp(-std::log(initial_learn_radius_ * 2)/(radius_gone_factor_ * total_steps_));
        learn_rate_decay_factor_ = std::exp(-std::log(initial_learn_rate_/final_learn_rate_)/total_steps_);

        grid3_.resize(boost::extents[N][N][N]);

        std::cout << "Map dimensions: " << N << "x" << N << "x" << N << " (" << N*N*N << " elements)" << std::endl;

        for (index x = 0; x != N; ++x)
            for (index y = 0; y != N; ++y)
                for (index z = 0; z != N; ++z)
                {
                    for (int i = 0; i != S; ++i)
                    {
                        grid3_[x][y][z](i) = som::random::double_range(0.f,1.f);
                    }
                }
        is_initialized_ = true;
    }

    /*! \brief Accessor for the map's elements

        \returns The element at position \f$(i,j,k)\f$
    */
    const ublas::c_vector<double,S>&
    element(int i, int j, int k) const
    {
        return grid3_[i][j][k];
    }

    /*! \brief Get best matching unit

        Returns the position (in 3d coordinates) of the closest vector (in terms of distance to the sample vector, depending on what metric the map uses).

        By default the euclidean distance is used:
        \f[ d(p,q)=\sqrt{(p_1-q_1)^2+(p_2-q_2)^2+...+(p_i-q_i)^2+...+(p_n-q_n)^2} \f]

        \param sample The sample vector
        \returns A som::point3 representing the coordinates of the closest node
     */
    som::vector3<int>
    best_matching_unit(const ublas::c_vector<double,S>& sample) const
    {
        som::vector3<int> p;
        double min_distance = 2.0; // trick: the euclidean norm for a unit vector will never be greater than sqrt(3)
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

    /** \brief Load input samples from a som::util::samples object to the map's internal (and slightly different) sample vector

        This function basically transforms each \c std::vector (corresponding to an input sample read from file) into an \c ublas::c_vector (vector dimensions must match). The values in the original vector are normalized to be in \f$[0,1]\f$ by dividing the vector (element-wise) with it's magnitude: \f[ \hat{v}=\frac{v}{||v||}\f] so that internally the map will always operate with unit vectors.

        \todo Work out a way to make sure the dimensions will match.

        \param samples The input sample vector
     */
    void load_samples(const som::util::samples& samples)
    {
        if (!is_initialized_)
        {
            std::cerr << "Error: map not initialized. Initializing..";
            init();
        }

        samples_.resize(samples.size());
        for (typename std::vector<ublas::c_vector<double,S> >::size_type i = 0; i != samples_.size(); ++i)
        {
            for (typename ublas::c_vector<double,S>::size_type j = 0; j != samples[i].size(); ++j)
            {
                samples_[i](j) = samples[i][j];
            }
            samples_[i] /= ublas::norm_2(samples_[i]); // normalization
            std::cout << "\t" << samples_[i] << std::endl;
        }
    }

    /** \brief Attempts to map 3d (rgb) or 4d (rgba) points (when the sample vectors are of length 3 or 4) to a 2d plane

        Uses boost::gil to create a png file, however, loosing one or 2 dimensions makes the image look weird. But it's an
        interesting pattern nevertheless and it can show that the map really learned something.

        \param filename Filename to be written to disk
      */
    void
    save_image(const std::string& filename) const
    {
        const int dim = ceil(sqrt(pow(N,3)));
        const int dim2 = dim * dim;
        boost::uint8_t r[dim2];
        boost::uint8_t g[dim2];
        boost::uint8_t b[dim2];
        boost::uint8_t a[dim2];
        boost::uint64_t x = 0;
        for (index i = 0; i != N; ++i)
            for (index j = 0; j != N; ++j)
                for (index k = 0; k != N; ++k)
                {
                    ublas::c_vector<double,S> v = grid3_[i][j][k] * 255;
                    r[x] = v[0];
                    g[x] = v[1];
                    b[x] = v[2];
                    a[x] = v[3];
                    ++x;
                }
//        boost::gil::rgba8c_planar_view_t view = boost::gil::planar_rgba_view(dim, dim, r, g, b, a, dim);
        boost::gil::rgb8c_planar_view_t view = boost::gil::planar_rgb_view(dim, dim, r, g, b, dim);
        boost::gil::png_write_view(filename, view);
    }

    void
    learn()
    {
        if (!is_initialized_)
        {
            std::cerr << "Error: map not initialized. Initializing..";
            init();
        }
        boost::progress_display pt(total_steps_);
        for (step_ = 0; step_ != total_steps_; ++step_)
        {
            learn_radius_ = round(initial_learn_radius_ * pow(radius_decay_factor_, step_));
            learn_rate_ = initial_learn_rate_ * pow(learn_rate_decay_factor_, step_);

            for (unsigned i = 0; i != samples_.size(); ++i)
            {
                som::vector3<int> bmu = best_matching_unit(samples_[i]);
                scale_neighbours(bmu, samples_[i]);
            }
            ++pt;
        }
    }

private:
    /*! \brief Modify the best matching unit and its neighbours according to the current sample

        \param bmu The best matching unit
      */
    inline void
    scale_neighbours(const som::vector3<int>& bmu, const ublas::c_vector<double,S>& sample)
    {
        int xmin = (bmu.x > learn_radius_) ? (bmu.x - learn_radius_) : 0;
        int ymin = (bmu.y > learn_radius_) ? (bmu.y - learn_radius_) : 0;
        int zmin = (bmu.z > learn_radius_) ? (bmu.z - learn_radius_) : 0;
        int xmax = ( (bmu.x + learn_radius_) < N) ? (bmu.x + learn_radius_) : N;
        int ymax = ( (bmu.y + learn_radius_) < N) ? (bmu.y + learn_radius_) : N;
        int zmax = ( (bmu.z + learn_radius_) < N) ? (bmu.z + learn_radius_) : N;

        double radius2 = learn_radius_ * learn_radius_;
        /* double sigma = learn_radius_ / 3.0 therefore sigma2 = radius2 / 9.0 */
        double sigma2 = radius2 / 9.0;

        for (int i = xmin; i != xmax; ++i)
            for (int j = ymin; j != ymax; ++j)
                for (int k = zmin; k != zmax; ++k)
                {
                    int x = (i-bmu.x), y = (j-bmu.y), z = (k-bmu.z);
                    double dist2 = x*x + y*y + z*z;

                    if (dist2 < radius2)
                    {
                        double radius_func = learn_rate_ * exp(-dist2 / (2*sigma2));
                        grid3_[i][j][k] += (sample-grid3_[i][j][k]) * radius_func;
                    }
                }
    }
}; // end of map class definition

} // namespace som

#endif // SOM_HPP
