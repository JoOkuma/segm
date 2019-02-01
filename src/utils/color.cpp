
#include "color.h"
#include <cmath>
#include <stdexcept>

#include <iostream>

inline double gammaCorrection(double u) {
    return ((u <= 0.0031308) ? 12.92 * u : 1.055 * pow(u, 1/2.4) - 0.055);
}


inline double invGammaCorrection(double u) {
    return ((u <= 0.04045) ? u / 12.92 : pow((u + 0.055)/1.055, 2.4));
}


inline double between0and1(double u) {
    return ((u > 1.0) ? 1.0 : ((u < 0.0) ? 0.0 : u));
}


inline double infSlopeCorrection(double u) {
    double delta_3 = 0.00885645;
    double delta_2_times_3 = 0.12841855;
    return ((u > delta_3) ? pow(u, 1.0/3.0) : u/delta_2_times_3 + 0.1379310345);
}


inline double invInfSlopeCorrection(double u) {
    double delta = 0.2068965517;
    double delta_2_times_3 = 0.12841855;
    return ((u > delta) ? pow(u, 3.0) : delta_2_times_3 * (u - 0.1379310345));
}


inline int argmax(const double *color) {
    int idx = ((color[0] >= color[1]) ? 0 : 1);
    return ((color[idx] >= color[2]) ? idx : 2);
}


inline int argmin(const double *color) {
    int idx = ((color[0] <= color[1]) ? 0 : 1);
    return ((color[idx] <= color[2]) ? idx : 2);
}


void rgb2xyz(const double *in_color, double *out_color)
{
    double r_lin = invGammaCorrection(in_color[0]);
    double g_lin = invGammaCorrection(in_color[1]);
    double b_lin = invGammaCorrection(in_color[2]);

    out_color[0] = 0.4123950 * r_lin + 0.3575957 * g_lin + 0.1804978 * b_lin;
    out_color[1] = 0.2125974 * r_lin + 0.7151977 * g_lin + 0.0721989 * b_lin;
    out_color[2] = 0.0192998 * r_lin + 0.1191998 * g_lin + 0.9504999 * b_lin;
}


void xyz2rgb(const double *in_color, double *out_color)
{
    double r_lin =  3.2406255 * in_color[0] - 1.537208  * in_color[1] - 0.4986286 * in_color[2];
    double g_lin = -0.9689307 * in_color[0] + 1.8757561 * in_color[1] + 0.0415175 * in_color[2];
    double b_lin =  0.0557101 * in_color[0] - 0.2040211 * in_color[1] + 1.0569959 * in_color[2];

    out_color[0] = gammaCorrection(r_lin);
    out_color[1] = gammaCorrection(g_lin);
    out_color[2] = gammaCorrection(b_lin);

    out_color[0] = between0and1(out_color[0]);
    out_color[1] = between0and1(out_color[1]);
    out_color[2] = between0and1(out_color[2]);
}


void xyz2lab(const double *in_color, double *out_color)
{
    const double xn = 0.95047;
    const double yn = 1.0;
    const double zn = 1.08883;

    double fx = infSlopeCorrection(in_color[0] / xn);
    double fy = infSlopeCorrection(in_color[1] / yn);
    double fz = infSlopeCorrection(in_color[2] / zn);

    out_color[0] = 116.0 * fy - 16.0;
    out_color[1] = 500.0 * (fx - fy);
    out_color[2] = 200.0 * (fy - fz);
}


void lab2xyz(const double *in_color, double *out_color)
{
    const double xn = 0.95047;
    const double yn = 1.0;
    const double zn = 1.08883;

    double fL = (in_color[0] + 16.0) / 116.0;
    out_color[0] = xn * invInfSlopeCorrection(fL + (in_color[1] / 500.0));
    out_color[1] = yn * invInfSlopeCorrection(fL);
    out_color[2] = zn * invInfSlopeCorrection(fL - (in_color[2] / 200.0));
}


void rgb2lab(const double *in_color, double *out_color)
{
    double xyz[3];
    rgb2xyz(in_color, xyz);
    xyz2lab(xyz, out_color);
}


void lab2rgb(const double *in_color, double *out_color)
{
    double xyz[3];
    lab2xyz(in_color, xyz);
    xyz2rgb(xyz, out_color);
}


void rgb2ypbpr(const double *in_color, double *out_color)
{
    const double Kr = 0.2627;
    const double Kg = 0.6780;
    const double Kb = 0.0593;

    // should the rgb values be gamma corrected?
    out_color[0] = Kr * in_color[0] + Kg * in_color[1] + Kb * in_color[2];
    out_color[1] = 0.5 * (in_color[2] - out_color[0]) / (1.0 - Kb);
    out_color[2] = 0.5 * (in_color[0] - out_color[0]) / (1.0 - Kr);
}


void ypbpr2rgb(const double *in_color, double *out_color)
{
    const double Kr = 0.2627;
    const double Kg = 0.6780;
    const double Kb = 0.0593;

    // should the rgb values be gamma corrected?
    out_color[0] = 2.0 * in_color[2] * (1.0 - Kr) + in_color[0];
    out_color[2] = 2.0 * in_color[1] * (1.0 - Kb) + in_color[0];
    out_color[1] = (in_color[0] - Kr * out_color[0] - Kb * out_color[2]) / Kg;
}


void rgb2hsv(const double *in_color, double *out_color)
{
    const double epsilon = 1e-7;

    int max = argmax(in_color);
    int min = argmin(in_color);
    double delta = in_color[max] - in_color[min];
    if (max != min)
    {
        if (max == 0) { // max = R
            out_color[0] = (in_color[1] - in_color[2]) / delta; // TODO check if mod6 is necessary
        } else if (max == 1) { // max = G
            out_color[0] = 2 + (in_color[2] - in_color[0]) / delta;
        } else { // max = B
            out_color[0] = 4 + (in_color[0] - in_color[1]) / delta;
        }
        out_color[0] *= 60;
        if (out_color[0] < 0.0)
            out_color[0] += 360.0;
    } else {
        out_color[0] = 0.0;
    }
    out_color[1] = (in_color[max] < epsilon) ? 0.0 : delta / in_color[max];
    out_color[2] = in_color[max];
}


void hsv2rgb(const double *in_color, double *out_color)
{
    const double epsilon = 1e-7;

    if (in_color[1] < epsilon) {    // s
        out_color[0] = in_color[2]; // v
        out_color[1] = in_color[2]; // v
        out_color[2] = in_color[2]; // v
    } else {
        double h = ((in_color[0] >= 360.0) ? 0.0 : in_color[0] / 60.0); // h
        long integer = (long) h;
        double diff = h - integer;
        double p = in_color[2] * (1.0 - in_color[1]); // v, s
        double q = in_color[2] * (1.0 - (in_color[1] * diff)); // v, s
        double t = in_color[2] * (1.0 - (in_color[1] * (1.0 - diff))); // v, s

        switch (integer)
        {
            case 0:
                out_color[0] = in_color[2]; // v
                out_color[1] = t;
                out_color[2] = p;
                break;
            case 1:
                out_color[0] = q;
                out_color[1] = in_color[2]; // v
                out_color[2] = p;
                break;
            case 2:
                out_color[0] = p;
                out_color[1] = in_color[2]; // v
                out_color[2] = t;
                break;
            case 3:
                out_color[0] = p;
                out_color[1] = q;
                out_color[2] = in_color[2]; // v
                break;
            case 4:
                out_color[0] = t;
                out_color[1] = p;
                out_color[2] = in_color[2]; // v
                break;
            case 5:
                out_color[0] = in_color[2]; // v
                out_color[1] = p;
                out_color[2] = q;
                break;
            default:
                throw std::runtime_error("Something is wrong in hsv2rgb function");
        }
    }
}


void gray2rgb(const double *in_color, double *out_color)
{
    out_color[0] = in_color[0];
    out_color[1] = in_color[0];
    out_color[2] = in_color[0];
}


void rgb2gray(const double *in_color, double *out_color)
{
    double lab[3];
    rgb2lab(in_color, lab);
    out_color[0] = lab[0] / 100.0;
}


void rgb2rgchroma(const double *in_color, double *out_color)
{
    double sum = in_color[0] + in_color[1] + in_color[2];
    out_color[0] = in_color[0] / sum;
    out_color[1] = in_color[1] / sum;
    out_color[2] = in_color[2] / sum;
}