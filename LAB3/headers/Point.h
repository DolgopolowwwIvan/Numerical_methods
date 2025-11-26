#pragma once
#ifndef POINT_H
#define POINT_H

namespace Com_Methods {
    class Point {
    private:
        double X, Y, Z;
    public:
        Point(double x = 0, double y = 0, double z = 0);
        double x() const;
        double y() const;
        double z() const;
    };
}

#endif