
#ifndef SEGM_COLOR_H
#define SEGM_COLOR_H

// TODO
//  - check when it is really sRGB versus RGB linear

/**
 * @param in_color      sRGB color components between (0.0, 1.0)
 * @param out_color     CIE XYZ color components  (Illuminant D65)
 * @details             https://en.wikipedia.org/wiki/SRGB#The_sRGB_transfer_function_.28.22gamma.22.29
 */

void rgb2xyz(const double *in_color, double *out_color);

/**
 * @param in_color      CIE XYZ color components  (Illuminant D65)
 * @param out_color     sRGB color components between (0.0, 1.0)
 * @details             https://en.wikipedia.org/wiki/SRGB#The_sRGB_transfer_function_.28.22gamma.22.29
 */

void xyz2rgb(const double *in_color, double *out_color);

/**
 * @param in_color      CIE XYZ color components
 * @param out_color     CIE LAB color components (Illuminant D65)
 * @details             https://en.wikipedia.org/wiki/CIELAB_color_space
 */

void xyz2lab(const double *in_color, double *out_color);

/**
 * @param in_color      CIE LAB color components (Illuminant D65)
 * @param out_color     CIE XYZ color components
 * @details             https://en.wikipedia.org/wiki/CIELAB_color_space
 */

void lab2xyz(const double *in_color, double *out_color);

/**
 * @param in_color      RGB color components
 * @param out_color     CIE LAB color components
 * @details             transformation computed with rgb2xyz and xyz2lab
 */

void rgb2lab(const double *in_color, double *out_color);

/**
 * @param in_color      CIE LAB color components
 * @param out_color     RGB color components
 * @details             transformation computed with rgb2xyz and xyz2lab
 */

void lab2rgb(const double *in_color, double *out_color);

/**
 * @param in_color      RGB color components
 * @param out_color     YPbPr color components, Y between (0, 1.0), Pb and Pr between (-0.5, 0.5)
 * @details             https://en.wikipedia.org/wiki/Rec._2020
 */

void rgb2ypbpr(const double *in_color, double *out_color);

/**
 * @param in_color      YPbPr color components, Y between (0, 1.0), Pb and Pr between (-0.5, 0.5)
 * @param out_color     RGB color components
 * @details             https://en.wikipedia.org/wiki/Rec._2020
 */

void ypbpr2rgb(const double *in_color, double *out_color);

/**
 * @param in_color      RGB color components
 * @param out_color     HSV (HSB) color components
 * @details             https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
 */

void rgb2hsv(const double *in_color, double *out_color);

/**
 * @param in_color      HSV (HSB) color components
 * @param out_color     RGB color components
 * @details             https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
 */
void hsv2rgb(const double *in_color, double *out_color);

/**
 * @attention in_color is a pointer just to keep the standard, it must 1 dimensional
 * @param in_color      Scalar between (0, 1)
 * @param out_color     RGB color components
 * @details
 */
void gray2rgb(const double *in_color, double *out_color);

/**
 * @attention out_color is a pointer just to keep the standard, it is 1 dimensional
 * @param in_color      RGB color components
 * @param out_color     Scalar, the Luminance from L*a*b color space
 * @details
 */
void rgb2gray(const double *in_color, double *out_color);

/**
 * @param in_color      RGB color components
 * @param out_color     RG Chromaticity color components
 * @details             https://en.wikipedia.org/wiki/Rg_chromaticity
 */
void rgb2rgchroma(const double *in_color, double *out_color);

#endif //SEGM_COLOR_H
