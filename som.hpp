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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

/* array and multi_array storage types */
#include <boost/multi_array.hpp>

/*  threadpool */
#include <boost/threadpool.hpp>

/* random numbers */
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>

/* misc: save images, show progress */
#include <boost/progress.hpp>
#include <boost/gil/extension/io/png_io.hpp>
#include <boost/lexical_cast.hpp>

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
    return dist(::twister);
}

/** \brief Returns an integer value \f$v\in[begin,end]\f$ (again, \e closed interval). */
int
int_range(int begin, int end)
{
    boost::uniform_int<> dist(begin, end);
    return dist(::twister);
}
} // namespace random

/* forward declaration of map */
template <int N, int S>
class map;

/** \brief Helper functions (read from file, save, etc) */
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

/** \brief Reads sample vectors line by line from a given filename

    \todo Maybe replace helper function \c split_by_spaces with \c boost::tokenizer something
*/
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

/** \brief Attempts to map 3d points (rgb or rgba) to a 2d plane

    This function uses the <a href="http://en.wikipedia.org/wiki/Isometric_projection">isometric projection</a> to draw the 3d lattice onto the planar surface.

    Uses \c boost::gil to create a png file.

    \param map The self-organizing map
    \param filename Filename to be written to disk
    \param slices Number of slices (cross sections of a given width along the x-axis) to divide the cube into (for easier representation). This parameter is optional.
  */
template <int N, int S>
void
save_to_image(const som::map<N,S>& map, const std::string& filename, int slices=0)
{
    typedef typename boost::multi_array<ublas::c_vector<double,S>, 3>::index index;

    /* offsets/padding */
    double ox = N;
    double oy = 50;

    /* image dimensions */
    unsigned d = 0.75*N; // horizontal distance between slices
    unsigned width;
    if (slices > 0)
        width = ox + slices * d + N/sqrt(2);
    else
        width = 2 * N;
    unsigned height = 2 * (N + oy);

    /* color channels - dynamically allocated to be able to create large images */
    boost::shared_array<boost::uint8_t> r(new boost::uint8_t[width*height]);
    boost::shared_array<boost::uint8_t> g(new boost::uint8_t[width*height]);
    boost::shared_array<boost::uint8_t> b(new boost::uint8_t[width*height]);
    boost::shared_array<boost::uint8_t> a(new boost::uint8_t[width*height]);

    /* black background, alpha set to 78% */
    for (unsigned i = 0; i != width*height; ++i)
    {
        r[i] = g[i] = b[i] = 0;
        a[i] = 200;
    }

    /* rotation angles for the isometric projection */
    double alpha = 0.615472907; // 35.264° = arcsin(tan(30°))
    double beta = 0.785398163; // 45°

    /*rotation matrixes */
    double r1_[3][3] = { {1,0,0}, {0,cos(alpha),sin(alpha)}, {0,-sin(alpha),cos(alpha)} };
    double r2_[3][3] = { {cos(beta),0,-sin(beta)}, {0,1,0}, {sin(beta), 0, cos(beta)} };

    ublas::matrix<double> r1(3,3), r2(3,3);
    for (unsigned i = 0; i != r1.size1(); ++i)
        for (unsigned j = 0; j != r1.size2(); ++j)
        {
            r1(i,j) = r1_[i][j];
            r2(i,j) = r2_[i][j];
        }

    boost::gil::rgb8_planar_view_t view = boost::gil::planar_rgb_view(width, height, r.get(), g.get(), b.get(), width);
//    boost::gil::rgba8_planar_view_t view = boost::gil::planar_rgba_view(width, height, r, g, b, a, width);
//
    for (index i = 0; i != N; ++i)
        for (index j = 0; j != N; ++j)
            for (index k = 0; k != N; ++k)
            {
                if (slices > 0)
                {
                    if (j == N-1 && k == N-1 && i % (N/slices) == 0 && i > 0)
                        ox += d;
                }
                unsigned v_[] = {i,j,k};

                ublas::c_vector<unsigned,3> v;
                std::copy(v_, v_+3, v.begin());

                ublas::c_vector<int,3> c = ublas::prod(ublas::prod(r1,r2), v);

                ublas::c_vector<double,S> m_vector = map.element(i,j,k) * 255; // model vector

                c[0] += ox;
                c[1] += oy;

                /* set the projected pixel */
                *view.at(c[0],c[1]) = boost::gil::rgb8_pixel_t(m_vector[0], m_vector[1], m_vector[2]);
            }

    boost::gil::png_write_view(filename, view);
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

template <int S>
class runnable
{
private:
    ublas::c_vector<double,S> sample_;
    int learn_radius_;
    double learn_rate_;
    som::vector3<int> bmu_;

public:
    runnable(int learn_radius, double learn_rate, som::vector3<int>& bmu)
        : learn_radius_(learn_radius),
          learn_rate_(learn_rate),
          bmu_(bmu)
    {}

    void
    run(boost::multi_array_ref<ublas::c_vector<double,S>, 3> grid3, ublas::c_vector<double,S>& sample,
        int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
    {
        double radius2 = learn_radius_ * learn_radius_;
        /* double sigma = learn_radius_ / 3.0 therefore sigma2 = radius2 / 9.0 */
        double sigma2 = radius2 / 9.0;

        for (int i = xmin; i != xmax; ++i)
            for (int j = ymin; j != ymax; ++j)
                for (int k = zmin; k != zmax; ++k)
                {
                    int x = (i-bmu_.x), y = (j-bmu_.y), z = (k-bmu_.z);
                    double dist2 = x*x + y*y + z*z;

                    if (dist2 < radius2)
                    {
                        double radius_func = learn_rate_ * exp(-dist2 / (2*sigma2));
                        grid3[i][j][k] += (sample - grid3[i][j][k]) * radius_func;
                    }
                }
    }
};

/** \brief The self-organizing map

    It learns (self-organizes) by scaling its units to become more similar to the input samples (also called activity vectors).

    At each step \f$k\f$ of the learning algorithm, the map learns an activity vector \f${AV}_k\f$ as follows:
        -# The best matching unit (\f$BMU\f$) is found, so that \f$d(BMU,{AV}_k)\f$ is minimal (where \f$d\f$ is the distance used)
        -# The \f$BMU\f$ and its neighbours are changed to be more similar to \f${AV}_k\f$.\n\n
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
           Two constants can be extracted from the above relations, to ease the computation (hint: the logarithm parts).\n\n
           The learning rate for model vectors within the \f$BMU\f$ neighbourhood (given by (2)) is scaled by a factor corresponding to a 3d Gaussian envelope with a standard deviation of \f$R(k)/3\f$:
           \f[ \eqalignno
                { {MV}_k[x,y,z] &= {MV}_{k-1}[z,y,z]+({AV}_k-{MV}_{k-1}[z,y,z]) \cdot L(k) \cdot
                  e^{ - \frac{1}{2} \cdot \left( \frac{\| BMU-{MV}_k \|}{R(k)/3} \right)^2 } & (3) }
           \f]
           Where \f$ \|BMU-{MV}_k\| = \sqrt {(x_{BMU}-x)^2+(y_{BMU}-y)^2+(z_{BMU}-z)^2} \f$ is the euclidean distance between the best matching unit and the current model vector within the neighbourhood.
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
    typedef typename boost::multi_array<ublas::c_vector<double,S>, 3>::iterator iterator3;
    typedef typename boost::multi_array<ublas::c_vector<double,S>, 3>::template subarray<2>::type::iterator iterator2;
    typedef typename boost::multi_array<ublas::c_vector<double,S>, 3>::template subarray<1>::type::iterator iterator1;

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

    /* thread pool */
    boost::threadpool::pool thread_pool_;

public:
    /** \brief Map constructor

        For now the map uses the euclidean distance (\c boost::numeric::ublas::norm_2)

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
        thread_pool_.size_controller().resize(boost::thread::hardware_concurrency());

        radius_decay_factor_ = std::exp(-std::log(initial_learn_radius_ * 2)/(radius_gone_factor_ * total_steps_));
        learn_rate_decay_factor_ = std::exp(-std::log(initial_learn_rate_/final_learn_rate_)/total_steps_);

        grid3_.resize(boost::extents[N][N][N]);

        std::cout << "Map dimensions: " << N << "x" << N << "x" << N << " (" << powl(N,3) << " elements)" << std::endl;

        for ( iterator3 it3 = grid3_.begin(); it3 != grid3_.end(); ++it3 )
            for ( iterator2 it2 = (*it3).begin(); it2 != (*it3).end(); ++it2 )
                for (iterator1 it1 = (*it2).begin(); it1 != (*it2).end(); ++it1 )
                    for (typename ublas::c_vector<double,S>::iterator it = (*it1).begin(); it != (*it1).end(); ++it)
                        *it = som::random::double_range(0.f,1.f);

        is_initialized_ = true;
    }

    /*! \brief Accessor for the map's elements

        \returns The model vector at position \f$(i,j,k)\f$ within the 3d lattice
    */
    const ublas::c_vector<double,S>&
    element(int i, int j, int k) const
    {
        return grid3_[i][j][k];
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
        for (typename std::vector<ublas::c_vector<double,S> >::size_type it1 = 0; it1 != samples_.size(); ++it1)
        {
            for (typename ublas::c_vector<double,S>::size_type it2 = 0; it2 != samples[it1].size(); ++it2)
            {
                samples_[it1](it2) = samples[it1][it2];
            }
            samples_[it1] /= ublas::norm_2(samples_[it1]); // normalization
        }
    }

    /** \brief Enter the learning loop

        Runs for a total number of epochs
        - at each epoch (discrete time step) the map scales the weight vectors of the \f$BMU\f$ and it's neighbours

        \todo Find a way to share work among threads (although this is overcooking it for a simple SOM)

      */
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
            if (step_ % 10 == 0)
            {
                std::string str = boost::lexical_cast<std::string>(step_/10);
                som::util::save_to_image(*this, str + ".png", 10);
            }

            learn_radius_ = round(initial_learn_radius_ * pow(radius_decay_factor_, step_));
            learn_rate_ = initial_learn_rate_ * pow(learn_rate_decay_factor_, step_);

            for (unsigned i = 0; i != samples_.size(); ++i)
            {
                som::vector3<int> bmu = best_matching_unit(samples_[i]);

                    int xmin = (bmu.x > learn_radius_) ? (bmu.x - learn_radius_) : 0;
                    int xmax = ( (bmu.x + learn_radius_) < N) ? (bmu.x + learn_radius_) : N;

                    int ymin = (bmu.y > learn_radius_) ? (bmu.y - learn_radius_) : 0;
                    int ymax = ( (bmu.y + learn_radius_) < N) ? (bmu.y + learn_radius_) : N;

                    int zmin = (bmu.z > learn_radius_) ? (bmu.z - learn_radius_) : 0;
                    int zmax = ( (bmu.z + learn_radius_) < N) ? (bmu.z + learn_radius_) : N;

                    som::runnable<S> run_me_over(learn_radius_, learn_rate_, bmu);

                    boost::threadpool::schedule(thread_pool_,
                                                boost::bind(&som::runnable<S>::run, &run_me_over,
                                                            boost::ref(grid3_), boost::ref(samples_[i]),
                                                            xmin, xmax, ymin, ymax, zmin, zmax));
                    thread_pool_.wait();
            }
            ++pt;
        }
    }

private:
    /*! \brief Modify the best matching unit and its neighbours according to the current sample

      \deprecated Will be completely replaced by instances of som::runnable

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

    /*! \brief Get best matching unit

        Returns the position (in 3d coordinates) of the closest vector (in terms of distance to the sample vector, depending on what metric the map uses).

        By default the euclidean distance is used:
        \f[ d(p,q)=\sqrt{(p_1-q_1)^2+(p_2-q_2)^2+...+(p_i-q_i)^2+...+(p_n-q_n)^2} \f]

        \param sample The sample vector
        \returns A som::point3 representing the coordinates of the closest node
     */
    inline som::vector3<int>
    best_matching_unit(const ublas::c_vector<double,S>& sample)
    {
        som::vector3<int> p;
        double min_distance = 2.0; // trick: the euclidean norm for a unit vector will never be greater than sqrt(3)
        for (index i = 0; i != N; ++i)
            for (index j = 0; j != N; ++j)
                for (index k = 0; k != N; ++k)
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
        return p;
    }
}; // end of map class definition

} // namespace som

#endif // SOM_HPP
