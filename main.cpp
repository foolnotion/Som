#include <SFML/System.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>

#include "som.hpp"

class som_window
{
public:
    som_window(unsigned w, unsigned h, const std::string& name)
        : width_(w), height_(h)
    {
        window_ptr = boost::make_shared<sf::RenderWindow>(sf::VideoMode(width_, height_), name);
    }

    void main_loop(som::map& map)
    {
        window_ptr->SetFramerateLimit(30);
        while (window_ptr->IsOpened())
        {
            sf::Event Event;
            while (window_ptr->GetEvent(Event))
            {
                // Close window : exit
                if (Event.Type == sf::Event::Closed)
                    window_ptr->Close();

                // Escape key : exit
                if ((Event.Type == sf::Event::KeyPressed) && (Event.Key.Code == sf::Key::Escape))
                    window_ptr->Close();
            }
            display_map(map);
            window_ptr->Display();
        }
    }

    void display_map(som::map& map)
    {
        unsigned map_size = map.get_size();
        unsigned dx = width_ / map_size;
        unsigned dy = height_ / map_size;

        for (unsigned i = 0; i < map_size; ++i)
        {
            for (unsigned j = 0; j < map_size; ++j)
            {
                boost::numeric::ublas::vector<double> weights = map(i,j)->get_weights();
                sf::Shape r = sf::Shape::Rectangle(i*dx, j*dy, i*dx+dx, j*dy+dy, sf::Color(weights[0]*255, weights[1]*255, weights[2]*255));
                window_ptr->Draw(r);
            }
        }
    }

private:
    unsigned width_, height_;
    boost::shared_ptr<sf::RenderWindow> window_ptr;
};

int main (void)
{
    unsigned map_size = 4;
    unsigned sample_size = 3; // rgb colors (3 components)

    som::map map(map_size, sample_size, boost::make_shared<som::euclidean_distance>());
    std::cout << map;

    std::cout << "Best matching unit for " << map(2,2)->get_weights() << std::endl;
    som::position p = map.get_bmu(map(2,2)->get_weights());
    std::cout << map(p.x,p.y)->get_weights() << std::endl;

//    som_window window(600, 600, "som");
//    window.main_loop(map);

    return EXIT_SUCCESS;
}
