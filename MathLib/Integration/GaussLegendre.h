
#pragma once

#include <vector>

namespace MathLib
{

class GaussLegendre
{
public:
    static double getPoint(long grd, long pkt) {
        switch (grd)
        {
        case 1:
            return 0.0;
        case 2:
            switch (pkt)
            {
            case 0:
                return 0.577350269189626;
            case 1:
                return -0.577350269189626;
            }
            break;
        case 3:
            switch (pkt)
            {
            case 0:
                return 0.774596669241483;
            case 1:
                return 0.0;
            case 2:
                return -0.774596669241483;
            }
            break;
        case 4:
            switch (pkt)
            {
            case 0:
                return 0.861136311594053;
            case 1:
                return 0.339981043584856;
            case 2:
                return -0.339981043584856;
            case 3:
                return -0.861136311594053;
            }
            break;
        }
        return 0.0;
    }

    static double getWeight(long grd, long pkt)
    {
        switch (grd)
        {
        case 1:
            return 2.0;
        case 2:
            switch (pkt)
            {
            case 0:
                return 1.0;
            case 1:
                return 1.0;
            }
            break;
        case 3:
            switch (pkt)
            {
            case 0:
                return 0.555555555555556;
            case 1:
                return 0.888888888888889;
            case 2:
                return 0.555555555555556;
            }
            break;
        case 4:
            switch (pkt)
            {
            case 0:
                return 0.347854845137454;
            case 1:
                return 0.652145154862546;
            case 2:
                return 0.652145154862546;
            case 3:
                return 0.347854845137454;
            }
            break;
        }
        return 0.0;
    }

    double integrate(double (*fun)(double), int sampling_size) {
        this->setPointAndWeight(sampling_size);
        double val = .0;
        for (int i=0; i<sampling_size; i++) {
            val += weights[i]*fun(points[i]);
        }
        return val;
    }

private:
    std::vector<double> points;
    std::vector<double> weights;

    void setPointAndWeight(int nQ) {
        points.resize(nQ);
        weights.resize(nQ);

        for (int i=0; i<nQ; i++) {
            points[i] = getPoint(nQ, i);
            weights[i] = getWeight(nQ, i);
        }
    }
};


}
