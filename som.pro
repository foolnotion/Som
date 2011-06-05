HEADERS += \
    som.hpp

SOURCES += \
    main.cpp

LIBS += -lpng -lboost_thread

QMAKE_CXXFLAGS = "-O3 -march=native -ffast-math -DNDEBUG -DBOOST_DISABLE_ASSERTS"
#QMAKE_CXXFLAGS = "-g -march=native"
