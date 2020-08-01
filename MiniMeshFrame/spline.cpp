#include "spline.h"

#if SPLINE_NORMALIZE
#else
const float Spline::sumOfWeights[]={1,4,6,2,6};
#endif

const void (*Spline::D0WeightFunctions[])(const float,float&,float&,float&,float&)={LinearD0Weights,QuadraticD0Weights,CubicD0Weights,CatmullRomD0Weights,CubicBSplineD0Weights};
const void (*Spline::D1WeightFunctions[])(const float,float&,float&,float&,float&)={LinearD1Weights,QuadraticD1Weights,CubicD1Weights,CatmullRomD1Weights,CubicBSplineD1Weights};
const void (*Spline::D2WeightFunctions[])(const float,float&,float&,float&,float&)={LinearD2Weights,QuadraticD2Weights,CubicD2Weights,CatmullRomD2Weights,CubicBSplineD2Weights};
const void (*Spline::D3WeightFunctions[])(const float,float&,float&,float&,float&)={LinearD3Weights,QuadraticD3Weights,CubicD3Weights,CatmullRomD3Weights,CubicBSplineD3Weights};
