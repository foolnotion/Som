HEADERS += \
    som.hpp

SOURCES += \
    main.cpp \
    som.cpp

LIBS += -lboost_thread -lsfml-graphics -lsfml-window -lsfml-system

QMAKE_CXXFLAGS = "-O3 -march=native"
