HEADERS += \
    som.hpp

SOURCES += \
    main.cpp


LIBS += -lboost_thread

QMAKE_CXXFLAGS = "-O3 -march=native"
