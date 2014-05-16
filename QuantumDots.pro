TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    adammath.cpp

HEADERS += \
    adammath.h
LIBS+=-larmadillo
