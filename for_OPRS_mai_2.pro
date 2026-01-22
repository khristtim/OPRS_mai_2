TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += D:/d/repositories/OPRS_mai

SOURCES += \
        main.cpp \
        custom.cpp \
        integrator.cpp \
        model.cpp
HEADERS += \
        custom.h \
        integrator.h \
    linalg.h \
        model.h \
