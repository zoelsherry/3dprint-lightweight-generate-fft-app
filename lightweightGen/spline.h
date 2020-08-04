#ifndef SPLINE_INCLUDED
#define SPLINE_INCLUDED

#define SPLINE_NORMALIZE 0
class Spline{
public:
	const enum{
		LINEAR,
		QUADRATIC,
		CUBIC,
		CATMULL_ROM,
		CUBIC_B_SPLINE,
		SPLINE_TYPES
	};
#if SPLINE_NORMALIZE
#else
	const static float sumOfWeights[SPLINE_TYPES];
#endif
	const static void (*D0WeightFunctions[SPLINE_TYPES])(const float,float&,float&,float&,float&);
	const static void (*D1WeightFunctions[SPLINE_TYPES])(const float,float&,float&,float&,float&);
	const static void (*D2WeightFunctions[SPLINE_TYPES])(const float,float&,float&,float&,float&);
	const static void (*D3WeightFunctions[SPLINE_TYPES])(const float,float&,float&,float&,float&);

	static const void LinearD0Weights		(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void QuadraticD0Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CubicD0Weights		(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CatmullRomD0Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CubicBSplineD0Weights	(const float t,float &w1,float& w2,float &w3,float &w4);

	static const void LinearD1Weights		(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void QuadraticD1Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CubicD1Weights		(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CatmullRomD1Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CubicBSplineD1Weights	(const float t,float &w1,float& w2,float &w3,float &w4);

	static const void LinearD2Weights		(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void QuadraticD2Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CubicD2Weights		(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CatmullRomD2Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CubicBSplineD2Weights	(const float t,float &w1,float& w2,float &w3,float &w4);

	static const void LinearD3Weights		(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void QuadraticD3Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CubicD3Weights		(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CatmullRomD3Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
	static const void CubicBSplineD3Weights	(const float t,float &w1,float& w2,float &w3,float &w4);
};


////////////////////////////////////
// Inline Spline Weight Functions //
////////////////////////////////////
const inline void Spline::LinearD0Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	w2=1.f-t;
	w3=t;
	w1=w4=0;
}
const inline void Spline::QuadraticD0Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	float t2;
	t2=t*t;
#if SPLINE_NORMALIZE
	w1=( t2-t    )/4;
	w2=(-t2-3*t+4)/4;
	w3=(-t2+5*t  )/4;
	w4=( t2-  t  )/4;
#else
	w1=( t2-t    );
	w2=(-t2-3*t+4);
	w3=(-t2+5*t  );
	w4=( t2-  t  );
#endif
}
const inline void Spline::CubicD0Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	float t2,t3;
	t2=t*t;
	t3=t2*t;
#if SPLINE_NORMALIZE
	w1=(-  t3+3*t2-2*t  )/6;
	w2=( 3*t3-6*t2-3*t+6)/6;
	w3=(-3*t3+3*t2+6*t  )/6;
	w4=(   t3-       t  )/6;
#else
	w1=(-  t3+3*t2-2*t  );
	w2=( 3*t3-6*t2-3*t+6);
	w3=(-3*t3+3*t2+6*t  );
	w4=(   t3-       t  );
#endif
}
const inline void Spline::CatmullRomD0Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	float t2,t3;
	t2=t*t;
	t3=t2*t;
#if SPLINE_NORMALIZE
	w1=(-  t3+2*t2-  t  )/2;
	w2=( 3*t3-5*t2    +2)/2;
	w3=(-3*t3+4*t2+  t  )/2;
	w4=(   t3-  t2      )/2;
#else
	w1=(-  t3+2*t2-  t  );
	w2=( 3*t3-5*t2    +2);
	w3=(-3*t3+4*t2+  t  );
	w4=(   t3-  t2      );
#endif
}
const inline void Spline::CubicBSplineD0Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	float t2,t3;
	t2=t*t;
	t3=t2*t;
#if SPLINE_NORMALIZE
	w1=(-  t3+3*t2-3*t+1)/6;
	w2=( 3*t3-6*t2    +4)/6;
	w3=(-3*t3+3*t2+3*t+1)/6;
	w4=(   t3           )/6;
#else
	w1=(-  t3+3*t2-3*t+1);
	w2=( 3*t3-6*t2    +4);
	w3=(-3*t3+3*t2+3*t+1);
	w4=(   t3           );
#endif
}
const inline void Spline::LinearD1Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	w2=-1;
	w3= 1;
	w1=w4=0;
}
const inline void Spline::QuadraticD1Weights(const float t,float& w1,float& w2,float& w3,float& w4){
#if SPLINE_NORMALIZE
	w1=( 2*t-1)/4;
	w2=(-2*t-3)/4;
	w3=(-2*t+5)/4;
	w4=( 2*t-1)/4;
#else
	w1=( 2*t-1);
	w2=(-2*t-3);
	w3=(-2*t+5);
	w4=( 2*t-1);
#endif
}
const inline void Spline::CubicD1Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	float t2=t*t;
#if SPLINE_NORMALIZE
	w1=(-3*t2+ 6*t-2)/6;
	w2=( 9*t2-12*t-3)/6;
	w3=(-9*t2+ 6*t+6)/6;
	w4=( 3*t2-     1)/6;
#else
	w1=(-3*t2+ 6*t-2);
	w2=( 9*t2-12*t-3);
	w3=(-9*t2+ 6*t+6);
	w4=( 3*t2-     1);
#endif
}
const inline void Spline::CatmullRomD1Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	float t2=t*t;
#if SPLINE_NORMALIZE
	w1=(-3*t2+ 4*t-1)/2;
	w2=( 9*t2-10*t  )/2;
	w3=(-9*t2+ 8*t+1)/2;
	w4=( 3*t2- 2*t  )/2;
#else
	w1=(-3*t2+ 4*t-1);
	w2=( 9*t2-10*t  );
	w3=(-9*t2+ 8*t+1);
	w4=( 3*t2- 2*t  );
#endif
}
const inline void Spline::CubicBSplineD1Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	float t2=t*t;
#if SPLINE_NORMALIZE
	w1=(-3*t2+ 6*t-3)/6;
	w2=( 9*t2-12*t  )/6;
	w3=(-9*t2+ 6*t+3)/6;
	w4=( 3*t2       )/6;
#else
	w1=(-3*t2+ 6*t-3);
	w2=( 9*t2-12*t  );
	w3=(-9*t2+ 6*t+3);
	w4=( 3*t2       );
#endif
}
const inline void Spline::LinearD2Weights(const float t,float& w1,float& w2,float& w3,float& w4){
	w1=w2=w3=w4=0;
}
const inline void Spline::QuadraticD2Weights(const float t,float& w1,float& w2,float& w3,float& w4){
#if SPLINE_NORMALIZE
	w1=w4= 0.5f;
	w2=w3=-0.5f;
#else
	w1=w4= 2;
	w2=w3=-2;
#endif
}
const inline void Spline::CubicD2Weights(const float t,float& w1,float& w2,float& w3,float& w4){
#if SPLINE_NORMALIZE
	w1=(-  t+1);
	w2=( 3*t-2);
	w3=(-3*t+1);
	w4=(   t  );
#else
	w1=(- 6*t+ 6);
	w2=( 18*t-12);
	w3=(-18*t+ 6);
	w4=(  6*t   );
#endif
}
const inline void Spline::CatmullRomD2Weights(const float t,float& w1,float& w2,float& w3,float& w4){
#if SPLINE_NORMALIZE
	w1=(-3*t+2);
	w2=( 9*t-5);
	w3=(-9*t+4);
	w4=( 3*t-1);
#else
	w1=(- 6*t+ 4);
	w2=( 18*t-10);
	w3=(-18*t+ 8);
	w4=(  6*t- 2);
#endif
}
const inline void Spline::CubicBSplineD2Weights(const float t,float& w1,float& w2,float& w3,float& w4){
#if SPLINE_NORMALIZE
	w1=(-  t+1);
	w2=( 3*t-2);
	w3=(-3*t+1);
	w4=(   t  );
#else
	w1=(- 6*t+ 6);
	w2=( 18*t-12);
	w3=(-18*t+ 6);
	w4=(  6*t   );
#endif
}
const inline void Spline::LinearD3Weights(const float t,float& w1,float& w2,float& w3,float& w4){w1=w2=w3=w4=0;}
const inline void Spline::QuadraticD3Weights(const float t,float& w1,float& w2,float& w3,float& w4){w1=w2=w3=w4=0;}

const inline void Spline::CubicD3Weights(const float t,float& w1,float& w2,float& w3,float& w4){
#if SPLINE_NORMALIZE
	w1=-1;
	w2= 3;
	w3=-3;
	w4= 1;
#else
	w1=-6;
	w2= 18;
	w3=-18;
	w4= 6;
#endif
}
const inline void Spline::CatmullRomD3Weights(const float t,float& w1,float& w2,float& w3,float& w4){
#if SPLINE_NORMALIZE
	w1=-3;
	w2= 9;
	w3=-9;
	w4= 3;
#else
	w1=-6;
	w2= 18;
	w3=-18;
	w4= 6;
#endif
}
const inline void Spline::CubicBSplineD3Weights(const float t,float& w1,float& w2,float& w3,float& w4){
#if SPLINE_NORMALIZE
	w1=-1;
	w2= 3;
	w3=-3;
	w4= 1;
#else
	w1=-6;
	w2= 18;
	w3=-18;
	w4= 6;
#endif
}

#endif // SPLINE_INCLUDED

