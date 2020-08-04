/*
Copyright (c) 2005, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <fftw3.h>
#include "spline.h"
#include "voxelData.h"
#include "splinterp.h"

#define PI     3.1415926535897932384
#define SQRT_2 1.4142135623709504880
#define SQRT_3 1.7320508075688772935

#define INLINE_MARCHING_CUBES 1

class IntIndex{
public:
	int x;
	int value;
	static inline int Compare(const void* v1,const void* v2){return ((IntIndex*)v1)->value-((IntIndex*)v2)->value;}
};

////////////////////////
// Factorization Code //
////////////////////////
inline int Factor(double a1,double a0,double roots[1][2]){
	if(a1==0){return 0;}
	roots[0][1]=-a0/a1;
	roots[0][0]=0;
	return 1;
}
inline int Factor(double a2,double a1,double a0,double roots[2][2]){
	double d;
	if(a2==0){return Factor(a1,a0,roots);}

	a1/=a2;
	a0/=a2;
	d=a1*a1-4*a0;
	a1/=2;
	if(d<0){
		d=sqrt(-d)/2;
		roots[0][0]=roots[1][0]=-a1;
		roots[0][1]=-d;
		roots[1][1]= d;
	}
	else{
		d=sqrt(d)/2;
		roots[0][1]=roots[1][1]=0;
		roots[0][0]=-a1-d;
		roots[1][0]=-a1+d;
	}
	return 2;
}
/////////////////////////////////////////////////////////////////////////
// Solution taken from: http://mathworld.wolfram.com/CubicFormula.html //
/////////////////////////////////////////////////////////////////////////
inline int Factor(double a3,double a2,double a1,double a0,double roots[3][2]){
	double Q,R,R2,Q3,D,S,SqrQ,SqrD;

	if(a3==0){return Factor(a2,a1,a0,roots);}
	a2/=a3;
	a1/=a3;
	a0/=a3;

	Q = (3*a1-a2*a2)/9;
	R = (9*a2*a1-27*a0-2*a2*a2*a2)/54;
	R2= R*R;
	Q3= Q*Q*Q;
	D = Q3+R2;

	a2/=3;

	if(D<=0){
		double theta,cTheta,sTheta;
		SqrQ  = sqrt(-Q);
		theta = acos(R/(SqrQ*Q));
		cTheta= cos(theta/3)*SqrQ;
		sTheta= sin(theta/3)*SqrQ*SQRT_3;

		roots[0][1]=roots[1][1]=roots[2][1]=0;
		roots[0][0]=-2*cTheta-a2;
		roots[1][0]=(cTheta+sTheta)-a2;
		roots[2][0]=(cTheta-sTheta)-a2;
	}
	else{
		double temp;
		SqrD=sqrt(D);
		temp=R+SqrD;
		if(temp<0){S=-pow(-temp,1.0/3);}
		else{S=pow(temp,1.0/3);}
		temp=R-SqrD;
		if(temp<0){S-=pow(-temp,1.0/3);}
		else{S+=pow(temp,1.0/3);}

		roots[0][1]=0;
		roots[0][0]=S-a2;
		S/=2;
		roots[1][0]= roots[2][0]=-S-a2;
		roots[1][1]= SQRT_3*S;
		roots[2][1]=-roots[1][1];
	}
	return 3;
}

////////////////////
// AdjacencyTable //
////////////////////

AdjacencyTable::AdjacencyTable(void){xTable=yTable=zTable=NULL;}
AdjacencyTable::~AdjacencyTable(void){
	if(xTable){delete[] xTable;}
	if(yTable){delete[] yTable;}
	if(zTable){delete[] zTable;}
	xTable=yTable=zTable=NULL;
}
int AdjacencyTable::set(const int res,const int padded,const int toroidal){
	int r,r2;

	if(padded){r=res+2;}
	else{r=res;}
	r2=r*res;

	if(xTable){delete[] xTable;}
	if(yTable){delete[] yTable;}
	if(zTable){delete[] zTable;}
	xTable=new int[res+4];
	yTable=new int[res+4];
	zTable=new int[res+4];
	if(!xTable || !yTable || !zTable){
		if(xTable){delete[] xTable;}
		if(yTable){delete[] yTable;}
		if(zTable){delete[] zTable;}
		xTable=yTable=zTable=NULL;
		return 0;
	}
	if(toroidal){
		xTable[0]=(res-1)*r2;
		yTable[0]=(res-1)*r;
		zTable[0]=(res-1);
	}
	else{xTable[0]=yTable[0]=zTable[0]=0;}
	for(int i=1;i<res+4;i++){
		zTable[i]=i-1;
		if(zTable[i]>res-1){
			if(toroidal){zTable[i]%=res;}
			else{zTable[i]=res-1;}
		}
		yTable[i]=zTable[i]*r;
		xTable[i]=zTable[i]*r2;
	}

	return 1;
}


///////////////////////
// VoxelEdgeFunction //
///////////////////////
const int VoxelEdgeFunction::getIndex(const int idx,const int e){
	switch(e){
	case 0:	return idx*3;
	case 1:	return (idx+res*res)*3+1;
	case 2:	return (idx+res)*3;
	case 3: return idx*3+1;

	case 4: return (idx+1)*3;
	case 5: return (idx+res*res+1)*3+1;
	case 6: return (idx+res+1)*3;
	case 7: return (idx+1)*3+1;

	case 8: return idx*3+2;
	case 9: return (idx+res*res)*3+2;
	case 10:return (idx+res*res+res)*3+2;
	case 11:return (idx+res)*3+2;
	}
	return -1;
}

///////////////
// VoxelData //
///////////////
VoxelData::VoxelData(void){
	array=NULL;
	resolution=0;
	padded=0;
}
VoxelData::~VoxelData(void){
	if(array){delete[] array;}
	array=NULL;
	resolution=0;
	padded=0;
}
void VoxelData::clear(void){
	int r=resolution*resolution*resolution;
	if(padded){r+=2*resolution*resolution;}
	memset(array,0,sizeof(float)*r);
}
int VoxelData::set(const int res,const int pad){
	if(array){delete[] array;}
	array=NULL;
	resolution=0;
	padded=0;
	if(res<=0){return 1;}
	int r=res*res*res;
	if(pad){r+=2*res*res;}
	array=new float[r];
	if(!array){return 0;}
	padded=pad;
	resolution=res;
	clear();
	return 1;
}
const int VoxelData::write(const char* fileName){
	FILE* fp=fopen(fileName,"wb");
	int res2,res=resolution;
	if(padded){res2=res+2;}
	else{res2=res;}
	if(!fp){return 0;}
	fwrite(&resolution,sizeof(int),1,fp);
	for(int i=0;i<res;i++){
		for(int j=0;j<res;j++){
			fwrite(&array[i*res*res2+j*res2],sizeof(float),res,fp);
		}
	}
	fclose(fp);
	return 1;
}
const float VoxelData::VertexValue(const float p[3],const AdjacencyTable& t,const int intType){
	float v[4][4],vv[4];
	int idx0,idx1,idx2,idx3;
	int i,j;
	float w1,w2,w3,w4;
	int x=int(p[0]);
	int y=int(p[1]);
	int z=int(p[2]);
	float dx=float(p[0]-x);
	float dy=float(p[1]-y);
	float dz=float(p[2]-z);

	// Intepolate across the x-axis
	Spline::D0WeightFunctions[intType](dx,w1,w2,w3,w4);
	for(i=0;i<4;i++){
		idx0=t.xTable[x  ]+t.yTable[i+y];
		idx1=t.xTable[x+1]+t.yTable[i+y];
		idx2=t.xTable[x+2]+t.yTable[i+y];
		idx3=t.xTable[x+3]+t.yTable[i+y];
		for(j=0;j<4;j++){
			v[i][j]=
				array[idx0+t.zTable[j+z]]*w1+
				array[idx1+t.zTable[j+z]]*w2+
				array[idx2+t.zTable[j+z]]*w3+
				array[idx3+t.zTable[j+z]]*w4;
		}
	}
	
	Spline::D0WeightFunctions[intType](dy,w1,w2,w3,w4);
	// Intepolate across the y-axis
	for(i=0;i<4;i++){
		vv[i]=
			v[0][i]*w1+
			v[1][i]*w2+
			v[2][i]*w3+
			v[3][i]*w4;
	}
	// Intepolate across the z-axis
	Spline::D0WeightFunctions[intType](dz,w1,w2,w3,w4);
#if SPLINE_NORMALIZE
	return (vv[0]*w1+vv[1]*w2+vv[2]*w3+vv[3]*w4);
#else
	return (vv[0]*w1+vv[1]*w2+vv[2]*w3+vv[3]*w4)/(Spline::sumOfWeights[intType]*Spline::sumOfWeights[intType]*Spline::sumOfWeights[intType]);
#endif
}

const void VoxelData::VertexGradient(const float p[3],float grad[3],const AdjacencyTable& t,const int intType){
	float v[4][4],vv[4];
	int idx0,idx1,idx2,idx3;
	int i,j;
	float w1,w2,w3,w4;
	int x=int(p[0]);
	int y=int(p[1]);
	int z=int(p[2]);
	float dx=float(p[0]-x);
	float dy=float(p[1]-y);
	float dz=float(p[2]-z);

	for(int d=0;d<3;d++){
		// Intepolate across the x-axis
		if(d==0){Spline::D1WeightFunctions[intType](dx,w1,w2,w3,w4);}
		else if(d==1){Spline::D0WeightFunctions[intType](dx,w1,w2,w3,w4);}
		if(d!=2){
			for(i=0;i<4;i++){
				idx0=t.xTable[x  ]+t.yTable[i+y];
				idx1=t.xTable[x+1]+t.yTable[i+y];
				idx2=t.xTable[x+2]+t.yTable[i+y];
				idx3=t.xTable[x+3]+t.yTable[i+y];
				for(j=0;j<4;j++){
					v[i][j]=
						array[idx0+t.zTable[j+z]]*w1+
						array[idx1+t.zTable[j+z]]*w2+
						array[idx2+t.zTable[j+z]]*w3+
						array[idx3+t.zTable[j+z]]*w4;
				}
			}
		}
		
		if(d==1){Spline::D1WeightFunctions[intType](dy,w1,w2,w3,w4);}
		else{Spline::D0WeightFunctions[intType](dy,w1,w2,w3,w4);}
		// Intepolate across the y-axis
		for(i=0;i<4;i++){
			vv[i]=
				v[0][i]*w1+
				v[1][i]*w2+
				v[2][i]*w3+
				v[3][i]*w4;
		}
		// Intepolate across the z-axis
		if(d==2){Spline::D1WeightFunctions[intType](dz,w1,w2,w3,w4);}
		else{Spline::D0WeightFunctions[intType](dz,w1,w2,w3,w4);}
#if SPLINE_NORMALIZE
		grad[d]=(vv[0]*w1+vv[1]*w2+vv[2]*w3+vv[3]*w4);
#else
		grad[d]=(vv[0]*w1+vv[1]*w2+vv[2]*w3+vv[3]*w4)/(Spline::sumOfWeights[intType]*Spline::sumOfWeights[intType]*Spline::sumOfWeights[intType]);
#endif
	}
}

const void VoxelData::VertexInterp(const float iso,const int edge,float p[3],const AdjacencyTable& t,const int intType){
	int cc,cnt,i,j,k,e=edge%3;
	double a,b,c,d;
	double v0,v1,v2,v3,mu;
	double roots[3][2];
	int e2=edge/3;

	k=e2%resolution;
	j=(e2/resolution)%resolution;
	i=(e2/(resolution*resolution))%resolution;
	p[0]=float(i);
	p[1]=float(j);
	p[2]=float(k);

	if(intType==Spline::LINEAR){
		switch(e){
		case 0:
			v1=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+1]];
			v2=array[t.xTable[i+2]+t.yTable[j+1]+t.zTable[k+1]];
			break;
		case 1:
			v1=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+1]];
			v2=array[t.xTable[i+1]+t.yTable[j+2]+t.zTable[k+1]];
			break;
		case 2:
			v1=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+1]];
			v2=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+2]];
			break;
		}
		mu = (iso-v1) / (v2 - v1);
		p[e]+=float(mu);
		return;
	}
	if(intType!=Spline::CUBIC_B_SPLINE){
		switch(e){
		case 0:
			v0=array[t.xTable[i  ]+t.yTable[j+1]+t.zTable[k+1]];
			v1=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+1]];
			v2=array[t.xTable[i+2]+t.yTable[j+1]+t.zTable[k+1]];
			v3=array[t.xTable[i+3]+t.yTable[j+1]+t.zTable[k+1]];
			break;
		case 1:
			v0=array[t.xTable[i+1]+t.yTable[j  ]+t.zTable[k+1]];
			v1=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+1]];
			v2=array[t.xTable[i+1]+t.yTable[j+2]+t.zTable[k+1]];
			v3=array[t.xTable[i+1]+t.yTable[j+3]+t.zTable[k+1]];
			break;
		case 2:
			v0=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k  ]];
			v1=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+1]];
			v2=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+2]];
			v3=array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+3]];
			break;
		}
	}

	if(intType==Spline::QUADRATIC){
		a=( v0-  v1-  v2+v3)/4;
		b=(-v0-3*v1+5*v2-v3)/4;
		c=v1-iso;
		cc=Factor(a,b,c,roots);
		if(cc==1){mu=roots[0][0];}
		else if(cc==2){
			cnt=0;
			mu=0;
			for(i=0;i<cc;i++){
				if(fabs(roots[i][1])<10e-4 && roots[i][0]>=0 && roots[i][0]<=1){
					cnt++;
					mu+=roots[i][0];
				}
			}
			if(!cnt){mu=0;}
			else{mu/=cnt;}
		}
	}
	else{
		if(intType==Spline::CUBIC){
			a =  -v0/6+v1/2-v2/2+v3/6;
			b = (-v1+(v0+v2)/2);
			c = (-v0/3-v1/2+v2-v3/6);
			d = ( v1-iso);
		}
		else if(intType==Spline::CATMULL_ROM){
			a=-0.5*v0+1.5*v1-1.5*v2+0.5*v3;
			b=     v0-2.5*v1+  2*v2-0.5*v3;
			c=-0.5*v0+       0.5*v2;
			d=            v1-              iso;
		}
		else if(intType==Spline::CUBIC_B_SPLINE){
			double t1,t2,t3;
			double vv[4];
			int ii;
			switch(e){
			case 0:
				for(ii=0;ii<4;ii++){
					t1=	array[t.xTable[i+ii]+t.yTable[j  ]+t.zTable[k  ]]+
						array[t.xTable[i+ii]+t.yTable[j+1]+t.zTable[k  ]]*4+
						array[t.xTable[i+ii]+t.yTable[j+2]+t.zTable[k  ]];
					t2=	array[t.xTable[i+ii]+t.yTable[j  ]+t.zTable[k+1]]+
		 				array[t.xTable[i+ii]+t.yTable[j+1]+t.zTable[k+1]]*4+
						array[t.xTable[i+ii]+t.yTable[j+2]+t.zTable[k+1]];
					t3=	array[t.xTable[i+ii]+t.yTable[j  ]+t.zTable[k+2]]+
						array[t.xTable[i+ii]+t.yTable[j+1]+t.zTable[k+2]]*4+
						array[t.xTable[i+ii]+t.yTable[j+2]+t.zTable[k+2]];
					vv[ii]=(t1+t2*4+t3);
				}
				break;
			case 1:
				for(ii=0;ii<4;ii++){
					t1=	array[t.xTable[i  ]+t.yTable[j+ii]+t.zTable[k  ]]+
						array[t.xTable[i  ]+t.yTable[j+ii]+t.zTable[k+1]]*4+
						array[t.xTable[i  ]+t.yTable[j+ii]+t.zTable[k+2]];
					t2=	array[t.xTable[i+1]+t.yTable[j+ii]+t.zTable[k  ]]+
						array[t.xTable[i+1]+t.yTable[j+ii]+t.zTable[k+1]]*4+
						array[t.xTable[i+1]+t.yTable[j+ii]+t.zTable[k+2]];
					t3=	array[t.xTable[i+2]+t.yTable[j+ii]+t.zTable[k  ]]+
						array[t.xTable[i+2]+t.yTable[j+ii]+t.zTable[k+1]]*4+
						array[t.xTable[i+2]+t.yTable[j+ii]+t.zTable[k+2]];
					vv[ii]=(t1+t2*4+t3);
				}
				break;
			case 2:
				for(ii=0;ii<4;ii++){
					t1=	array[t.xTable[i  ]+t.yTable[j  ]+t.zTable[k+ii]]+
						array[t.xTable[i  ]+t.yTable[j+1]+t.zTable[k+ii]]*4+
						array[t.xTable[i  ]+t.yTable[j+2]+t.zTable[k+ii]];
					t2=	array[t.xTable[i+1]+t.yTable[j  ]+t.zTable[k+ii]]+
						array[t.xTable[i+1]+t.yTable[j+1]+t.zTable[k+ii]]*4+
						array[t.xTable[i+1]+t.yTable[j+2]+t.zTable[k+ii]];
					t3=	array[t.xTable[i+2]+t.yTable[j  ]+t.zTable[k+ii]]+
						array[t.xTable[i+2]+t.yTable[j+1]+t.zTable[k+ii]]*4+
						array[t.xTable[i+2]+t.yTable[j+2]+t.zTable[k+ii]];
					vv[ii]=(t1+t2*4+t3);
				}
				break;
			}

			a=(-vv[0]+3*vv[1]-3*vv[2]+vv[3]);
			b=( vv[0]-2*vv[1]+  vv[2]      )*3;
			c=(-vv[0]+          vv[2]      )*3;
			d=( vv[0]+4*vv[1]+  vv[2]      )-iso*216;
		}
		
		cc=Factor(a,b,c,d,roots);
		mu=0;
		cnt=0;
		for(i=0;i<cc;i++){
			if(fabs(roots[i][1])<10e-4 && roots[i][0]>=0 && roots[i][0]<=1){
				cnt++;
				mu+=roots[i][0];
			}
		}
		if(!cnt){mu=0;}
		else{mu/=cnt;}
	}
	p[e]+=float(mu);
}
const int VoxelData::setCubicBSplineValues(VoxelData& v){
	int i,j,k,r,r2;
	int idx;
	float temp;

	if(padded){r=resolution+2;}
	else{r=resolution;}
	r2=r*resolution;

	
	if(!v.set(resolution,padded)){return 0;}
	for(k=0;k<resolution;k++){
		for(j=0;j<resolution;j++){
			idx=j*r+k;
			for(i=idx+r2;i<idx+(resolution-1)*r2;i+=r2){
				v.array[i]=array[i-r2]+array[i]*4+array[i+r2];
			}
		}
	}
	
	for(i=0;i<resolution;i++){
		for(k=0;k<resolution;k++){
			idx=i*r2+k;
			temp=v.array[idx];
			for(j=idx+r;j<idx+(resolution-1)*r;j+=r){
				float tt=v.array[j];
				v.array[j]=temp+tt*4+v.array[j+r];
				temp=tt;
			}
		}
	}
	for(i=0;i<resolution;i++){
		for(j=0;j<resolution;j++){
			idx=i*r2+j*r;
			temp=v.array[idx];
			for(k=idx+1;k<idx+r-1;k++){
				float tt=v.array[k];
				v.array[k]=(temp+tt*4+v.array[k+1])/216;
				temp=tt;
			}
		}
	}
	return 1;
}
const int VoxelData::setCubicBSplineValues(VoxelData& v,int fftwType){
	int i,j,k,r,r2,bw=resolution/2;
	int idx1,idx2;
	float c,c1,c2;
	float *fIn;
	fftwf_complex *out,*fOut;
	fftwf_plan plan,iplan,fplan;

	if(padded){r=resolution+2;}
	else{r=resolution;}
	r2=r*resolution;

	if(!v.set(resolution,1)){return 0;}
	out = (fftwf_complex*)v.array;
	fOut = (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*(bw+1));
	fIn  = (float*)fftw_malloc(sizeof(float)*resolution);
	if(!fIn || !fOut){
		if(fIn){fftwf_free(fIn);}
		if(fOut){fftwf_free(fOut);}
		return 0;
	}
	
	for(i=0;i<resolution;i++){fIn[i]=0;}
	fIn[0]=4.f/6;
	fIn[1]=fIn[resolution-1]=1.f/6;
	c=0;
	fplan = fftwf_plan_dft_r2c_1d(resolution,fIn,fOut,fftwType);
	fftwf_execute(fplan);
	fftw_free(fIn);
	fftwf_destroy_plan(fplan);
	
	for(i=0;i<resolution;i++){
		for(j=0;j<resolution;j++){
			idx1=(i*resolution+j)*(resolution+2);
			idx2=i*r2+j*r;
			memcpy(&v.array[idx1],&array[idx2],resolution*sizeof(float));
		}
	}
	
	plan  = fftwf_plan_dft_r2c_3d(resolution,resolution,resolution,v.array,out,fftwType);
	fftwf_execute(plan);
	idx1=0;
	for(i=0;i<resolution;i++){
		if(i>bw){c1=fOut[resolution-i][0];}
		else{c1=fOut[i][0];}
		for(j=0;j<resolution;j++){
			if(j>bw){c2=c1*fOut[resolution-j][0];}
			else{c2=c1*fOut[j][0];}
			for(k=0;k<=bw;k++){
				c=c2*fOut[k][0];
				out[idx1][0]*=c;
				out[idx1][1]*=c;
				idx1++;
			}
		}
	}
	iplan = fftwf_plan_dft_c2r_3d(resolution,resolution,resolution,out,v.array,fftwType);
	fftwf_execute(iplan);
	fftwf_free(fOut);
	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(iplan);

	return 1;
}

const void VoxelData::isoSurface(vector<TriangleIndex>& triangles,const float isoValue,const int start,const int stop){
	int i,j,k,r,r2;
#if INLINE_MARCHING_CUBES
#else
	double z[2][2][2];
	double (*zz1)[2],(*zz2)[2],(*tempZ)[2];
#endif
	int idx;
	VoxelEdgeFunction f;

	if(padded){r=resolution+2;}
	else{r=resolution;}
	r2=r*resolution;

	f.res=resolution;

	triangles.clear();
	for(i=start;i<stop;i++){
		for(j=start;j<stop;j++){
			int idx11=(i  )*r2+(j  )*r;
			int idx12=(i  )*r2+(j+1)*r;
			int idx21=(i+1)*r2+(j  )*r;
			int idx22=(i+1)*r2+(j+1)*r;
#if INLINE_MARCHING_CUBES
			idx=0;
			if (array[idx11+start] < isoValue) idx |=  16;  // {0,0,0}
			if (array[idx21+start] < isoValue) idx |=  32;  // {1,0,0}
			if (array[idx22+start] < isoValue) idx |=  64;  // {1,1,0}
			if (array[idx12+start] < isoValue) idx |= 128;  // {0,1,0}
#else
			z[0][0][0]=array[idx11+start];
			z[0][1][0]=array[idx21+start];
			z[0][0][1]=array[idx12+start];
			z[0][1][1]=array[idx22+start];
			zz1=z[0];
			zz2=z[1];
#endif
			for(k=start+1;k<=stop;k++){
#if INLINE_MARCHING_CUBES
#else
				zz2[0][0]=array[idx11+k];
				zz2[1][0]=array[idx21+k];
				zz2[0][1]=array[idx12+k];
				zz2[1][1]=array[idx22+k];
#endif

#if INLINE_MARCHING_CUBES
				// Marching cubes code pulled to here
				idx>>=4;
				if (array[idx11+k] < isoValue) idx |=  16;  // {0,0,1}
				if (array[idx21+k] < isoValue) idx |=  32;  // {1,0,1}
				if (array[idx22+k] < isoValue) idx |=  64;  // {1,1,1}
				if (array[idx12+k] < isoValue) idx |= 128;  // {0,1,1}
				
				// Cube is entirely in/out of the surface
				if (!MarchingCubes::edges[idx]){continue;}
				
				// Find the vertices where the surface intersects the cube
				int c,ii=1;
				int vIndex=(i*resolution+j)*resolution+k-1;
				for(c=0;c<12;c++){
					if(MarchingCubes::edges[idx] & ii){MarchingCubes::vertexIndexList[c] = f.getIndex(vIndex,c);}
					ii<<=1;
				}
				// Create the triangle
				c=0;
				while(MarchingCubes::triangles[idx][c]!=-1){
					TriangleIndex tri;
					tri.idx[0]=MarchingCubes::vertexIndexList[MarchingCubes::triangles[idx][c++]];
					tri.idx[1]=MarchingCubes::vertexIndexList[MarchingCubes::triangles[idx][c++]];
					tri.idx[2]=MarchingCubes::vertexIndexList[MarchingCubes::triangles[idx][c++]];
					triangles.push_back(tri);
				}
#else
				MarchingCubes::AddTriangles(zz1,zz2,isoValue,&f,(i*resolution+j)*resolution+k-1,triangles,idx,k-start-1);
				tempZ=zz1;
				zz1=zz2;
				zz2=tempZ;
#endif
			}
		}
	}
}
const void VoxelData::isoSurface(vector<TriangleIndex>& triangles,vector<Vertex>& vertices,const float isoValue,const AdjacencyTable& t,const int intType,const int fftwType){
	int i,j,count,tNum;
	int* table;
	IntIndex* indices;

	if(intType==Spline::CUBIC_B_SPLINE){
		VoxelData temp;
		setCubicBSplineValues(temp,fftwType);
		temp.isoSurface(triangles,isoValue,1,resolution-2);
	}
	else if(intType==Spline::LINEAR){isoSurface(triangles,isoValue,0,resolution-1);}
	else{isoSurface(triangles,isoValue,1,resolution-2);}

	tNum=triangles.size();
	indices=new IntIndex[tNum*3];
	table=new int[tNum*3];
	for(i=0;i<tNum;i++){
		for(j=0;j<3;j++){
			indices[3*i+j].x=3*i+j;
			indices[3*i+j].value=triangles[i].idx[j];
		}
	}
	qsort(indices,tNum*3,sizeof(IntIndex),IntIndex::Compare);
	count=0;
	for(i=0;i<tNum*3;i++){if(!i || indices[i].value!=indices[i-1].value){count++;}}
	vertices.clear();
	vertices.resize(count);
	count=0;
	for(i=0;i<tNum*3;i++){
		if(!i || indices[i].value!=indices[i-1].value){
			VertexInterp(isoValue,indices[i].value,vertices[count].v,t,intType);
			table[indices[i].x]=count;
			count++;
		}
		else{table[indices[i].x]=count-1;}
	}
	for(i=0;i<tNum;i++){for(j=0;j<3;j++){triangles[i].idx[j]=table[3*i+j];}}
	delete[] table;
	delete[] indices;
}
const void VoxelData::isoSurfaceCubicBSpline(VoxelData& cubicBSplineValues,vector<TriangleIndex>& triangles,vector<Vertex>& vertices,const float isoValue,const AdjacencyTable& t,const int fftwType){
	int i,j,count,tNum;
	int* table;
	IntIndex* indices;

	cubicBSplineValues.isoSurface(triangles,isoValue,1,resolution-2);
	tNum=triangles.size();
	indices=new IntIndex[tNum*3];
	table=new int[tNum*3];
	for(i=0;i<tNum;i++){
		for(j=0;j<3;j++){
			indices[3*i+j].x=3*i+j;
			indices[3*i+j].value=triangles[i].idx[j];
		}
	}
	qsort(indices,tNum*3,sizeof(IntIndex),IntIndex::Compare);
	count=0;
	for(i=0;i<tNum*3;i++){if(!i || indices[i].value!=indices[i-1].value){count++;}}
	vertices.clear();
	vertices.resize(count);
	count=0;
	for(i=0;i<tNum*3;i++){
		if(!i || indices[i].value!=indices[i-1].value){
			VertexInterp(isoValue,indices[i].value,vertices[count].v,t,Spline::CUBIC_B_SPLINE);
			table[indices[i].x]=count;
			count++;
		}
		else{table[indices[i].x]=count-1;}
	}
	for(i=0;i<tNum;i++){for(j=0;j<3;j++){triangles[i].idx[j]=table[3*i+j];}}
	delete[] table;
	delete[] indices;
}



void VoxelData::interpCell(vector<float>& oldPhi, vector<float>& newPhi, int interpNum)
{
	double ms_ori = exp(log(oldPhi.size()) / 3); // can not use int because lose accuracy
	//
	//cout << ceil(ms_ori) << endl;
	//
	ms_ori = ceil(ms_ori); // transform the accurate int

	int len = ms_ori * ms_ori * ms_ori;
	float *data = new float[len];
	for (int i = 0; i < oldPhi.size(); i++)
	{
		data[i] = oldPhi[i];
	}

	int len_tar = interpNum * interpNum * interpNum;
	float* linex = new float[interpNum];
	float* liney = new float[interpNum];
	float* linez = new float[interpNum];
	lineSpace<float>(linex, 0, ms_ori - 1, interpNum);
	lineSpace<float>(liney, 0, ms_ori - 1, interpNum);
	lineSpace<float>(linez, 0, ms_ori - 1, interpNum);
	float* xx = new float[len_tar];
	float* yy = new float[len_tar];
	float* zz = new float[len_tar];
	meshGrid<float>(xx, yy, zz, linex, liney, linez, interpNum, interpNum, interpNum);

	float *temp = new float[len_tar];
	splinterp::parallel_interp3(splinterp::interp3_F<float>, data, ms_ori, ms_ori, ms_ori, xx, yy, zz, len_tar, temp, 0);

	for (int i = 0; i < len_tar; i++)
	{
		newPhi.push_back(temp[i]);
	}

	delete[] data;
	delete[] linex;
	delete[] liney;
	delete[] linez;

	delete[] xx;
	delete[] yy;
	delete[] zz;
	delete[] temp;
}

void VoxelData::fillUniformStructure_interior(const char* fileName, const float isoValue)
{
	vector<float> cell;

	ifstream fp;
	fp.exceptions(ifstream::badbit); // delete the `failbit` otherwise error
	try
	{
		fp.open(fileName, ifstream::in);
		double a;
		while (fp >> a)
			cell.push_back(a);

	}
	catch (ifstream::failure e)
	{
		std::cerr << "ERROR::CELL::Caught an exception: " << e.what() << endl;
	}
	
	for (auto &p : cell)
	{
		p = p * 2 - 1;
	}

	vector<float> newCell;
	int maxBox = 10;
	interpCell(cell, newCell, maxBox);

	int r, r2;

	if (padded) { r = resolution + 2; }
	else { r = resolution; }
	r2 = r * resolution;

	int start = 1;
	//int stop = resolution - 2;
	int stop = resolution - 1;


	//int np = floor(stop / (maxBox - 1));

	for (int i = start; i < stop; i++)
	{
		for (int j = start; j < stop; j++)
		{
			// the size of cell reduce x0y x0z y0z
			for (int k = start; k < stop; k++)
			{
				array[(i)*r2 + (j)*r + k] = std::min(array[(i)*r2 + (j)*r + k],
					newCell[((i - 1) % (maxBox - 1))*maxBox*maxBox + ((j - 1) % (maxBox - 1))*maxBox + (k - 1) % (maxBox - 1)]);

			}
		/*	array[idx0 + start] = std::min(array[idx0 + start], 
				newCell[((i - 1) % (maxBox - 1))*maxBox*maxBox + ((j - 1) % (maxBox - 1))*maxBox]);

			for (int k = start + 1; k <= stop; k++)
			{
				array[idx0 + k] = std::min(array[idx0 + k],
					newCell[((i - 1) % (maxBox - 1))*maxBox*maxBox + ((j - 1) % (maxBox - 1))*maxBox + (k - 1) % (maxBox - 1)]);

			}*/
		}
	}

#if (0)
	cout << "output the debug data..." << endl;
	std::ofstream fout;
	fout.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		fout.open("./outdata/array_l.txt");
		for (int i = start; i < stop; i++)
		{
			for (int j = start; j < stop; j++)
			{
				int idx0 = (i)*r2 + (j)*r;

				for (int k = start + 1; k <= stop; k++)
				{
					fout << array[idx0 + k] << endl;

				}


			}
		}

		fout.close();
	}
	catch (std::ofstream::failure e)
	{
		std::cout << "ERROR::SIGNEDDIST::FILE FAIL TO WRITE!" << std::endl;
	}
#endif

}

void VoxelData::fillUniformStructure(const char * fileName, const float isoValue)
{
	vector<float> cell;

	ifstream fp;
	fp.exceptions(ifstream::badbit); // delete the `failbit` otherwise error
	try
	{
		fp.open(fileName, ifstream::in);
		double a;
		while (fp >> a)
			cell.push_back(a);

	}
	catch (ifstream::failure e)
	{
		std::cerr << "ERROR::CELL::Caught an exception: " << e.what() << endl;
	}

	for (auto &p : cell)
	{
		p = p * 2 - 1;
	}

	vector<float> newCell;
	int maxBox = 10;
	interpCell(cell, newCell, maxBox);

	int r, r2;

	if (padded) { r = resolution + 2; }
	else { r = resolution; }
	r2 = r * resolution;

	int start = 1;
	//int stop = resolution - 2;
	int stop = resolution - 1;


	//int np = floor(stop / (maxBox - 1));

	for (int i = start; i < stop; i++)
	{
		for (int j = start; j < stop; j++)
		{
			// the size of cell reduce x0y x0z y0z
			for (int k = start; k < stop; k++)
			{
				array[(i)*r2 + (j)*r + k] = std::min(array[(i)*r2 + (j)*r + k],
					newCell[((i - 1) % (maxBox - 1))*maxBox*maxBox + ((j - 1) % (maxBox - 1))*maxBox + (k - 1) % (maxBox - 1)]);

			}

		}
	}
}

void VoxelData::fillAdaptiveStructure_interior(const char * fileName, const float isoValue, const int depth)
{
	vector<float> cell;

	ifstream fp;
	fp.exceptions(ifstream::badbit); // delete the `failbit` otherwise error
	try
	{
		fp.open(fileName, ifstream::in);
		double a;
		while (fp >> a)
			cell.push_back(a);

	}
	catch (ifstream::failure e)
	{
		std::cerr << "ERROR::CELL::Caught an exception: " << e.what() << endl;
	}

	for (auto &p : cell)
	{
		p = p * 2 - 1;
	}

	vector<vector<float>> newCell;
	int maxBox = pow(2, 6);
	int tempBox = maxBox;
	for (int i = 0; i < depth; i++)
	{
		vector<float> tempCell;
		interpCell(cell, tempCell, tempBox);
		newCell.push_back(tempCell);
		tempBox = maxBox / 2;
	}

	int r, r2;

	if (padded) { r = resolution + 2; }
	else { r = resolution; }
	r2 = r * resolution;

	int start = 1;
	//int stop = resolution - 2;
	int stop = resolution - 1;


	queue<float*> q;
	q.push(array);

	while (!q.empty())
	{
		float* tempArray = q.front();
		
	}
	

}

void VoxelData::fillAdaptiveStructure_recursive_interior(const char * fileName, const float isoValue, const int depth)
{
	vector<float> cell;

	ifstream fp;
	fp.exceptions(ifstream::badbit); // delete the `failbit` otherwise error
	try
	{
		fp.open(fileName, ifstream::in);
		double a;
		while (fp >> a)
			cell.push_back(a);

	}
	catch (ifstream::failure e)
	{
		std::cerr << "ERROR::CELL::Caught an exception: " << e.what() << endl;
	}

	for (auto &p : cell)
	{
		p = p * 2 - 1;
	}

	vector<vector<float>> newCell;
	int maxBox = pow(2, 6);
	int tempBox = maxBox;
	for (int i = 0; i < depth; i++)
	{
		vector<float> tempCell;
		interpCell(cell, tempCell, tempBox);
		newCell.push_back(tempCell);
		tempBox = maxBox / 2;
	}

	int r, r2;

	if (padded) { r = resolution + 2; }
	else { r = resolution; }
	r2 = r * resolution;


	int cellLength = maxBox;

	int start = 1;
	//int stop = resolution - 2;
	int stop = resolution - 1;


	int np = floor(stop / (maxBox - 1));
	int cellNum = pow(cellLength, 3);
	for (int ix = 0; ix < np; ix++)
	{
		for (int jy = 0; jy < np; jy++)
		{
			// the size of cell reduce x0y x0z y0z
			for (int kz = 0; kz < np; kz++)
			{

				int idx = (ix)*np*np + (jy)*np + kz;
				int count = 0;
				for (int i = 0; i < cellLength; i++)
				{
					for (int j = 0; j < cellLength; j++)
					{
						for (int k = 0; k < cellLength; k++)
						{
							if (array[idx + i*cellLength*cellLength + j * cellLength + k] > 0)
								count++;
						}
					}
				}

				if (count == cellNum)
				{
					for (int i = 0; i < cellLength; i++)
					{
						for (int j = 0; j < cellLength; j++)
						{
							for (int k = 0; k < cellLength; k++)
							{
								array[idx + i * cellLength*cellLength + j * cellLength + k] 
									= newCell[0][i*cellLength*cellLength + j * cellLength + k];
							}
						}
					}
				}
				else
				{

				}

			}

		}
	}
}
