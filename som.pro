HEADERS += \
    som.hpp

SOURCES += \
    main.cpp


LIBS += -lboost_thread -lsfml-graphics -lsfml-window -lsfml-system -lpng

QMAKE_CXXFLAGS = "-O3 -march=native -DNDEBUG -DBOOST_DISABLE_ASSERTS"
