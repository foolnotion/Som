HEADERS += \
    som.hpp

SOURCES += \
    main.cpp \
    som.cpp

#LIBS += -lboost_thread -lsfml-graphics -lsfml-window -lsfml-system
LIBS += -lboost_thread

QMAKE_CXXFLAGS = "-O3 -march=native"
