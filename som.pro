HEADERS += \
    som.hpp

SOURCES += \
    main.cpp

LIBS += -lpng

#QMAKE_CXXFLAGS = "-O3 -march=native -DNDEBUG -DBOOST_DISABLE_ASSERTS"
QMAKE_CXXFLAGS = "-g -march=native"
