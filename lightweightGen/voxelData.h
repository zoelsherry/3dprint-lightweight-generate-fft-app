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
#ifndef VOXEL_DATA_INCLUDED
#define VOXEL_DATA_INCLUDED

#define USE_ADJACENT_TABLE 1
#include <vector>
#include <math.h>
#include "marchingCubes.h"

using std::vector;

int Factor(                    double a1,double a0,double roots[1][2]);
int Factor(          double a2,double a1,double a0,double roots[2][2]);
int Factor(double a3,double a2,double a1,double a0,double roots[3][2]);

class Vertex{
public:
	float v[3];
};

class AdjacencyTable{
public:
	int *xTable;
	int *yTable;
	int *zTable;
	AdjacencyTable(void);
	~AdjacencyTable(void);
	int set(const int res,const int padded,const int toroidal=0);
};
class VoxelEdgeFunction : public EdgeFunction{
public:
	int res;
	const int getIndex(const int idx,const int e);
};
class VoxelData{
public:
	int padded;
	int resolution;
	float* array;

	VoxelData(void);
	~VoxelData(void);
	int set(const int res,const int pad=0);
	void clear(void);
	const int write(const char* fileName);

	const float index(const int* idx,const float* dx);
	const static void Position(const int res,const int padded,const float x,const float y,const float z,int* idx,float* dx);

	const void  VertexInterp	(const float iso,const int edge,float v[3],const AdjacencyTable& table,const int intType);
	const float VertexValue		(const float v[3],const AdjacencyTable& table,const int intType);
	const void  VertexGradient	(const float v[3],float grad[3],const AdjacencyTable& table,const int intType);

	const void isoSurface(vector<TriangleIndex>& triangles,const float isoValue,const int start,const int stop);
	const void isoSurface(vector<TriangleIndex>& triangles,vector<Vertex>& vertices,const float isoValue,const AdjacencyTable& table,const int intType,const int fftwType);
	const void isoSurfaceCubicBSpline(VoxelData& cubicBSplineValues,vector<TriangleIndex>& triangles,vector<Vertex>& vertices,const float isoValue,const AdjacencyTable& table,const int fftwType);
	const int setCubicBSplineValues(VoxelData& v);
	const int setCubicBSplineValues(VoxelData& v,const int fftwType);

	void interpCell(vector<float>& oldPhi, vector<float>& newPhi, int interpNum);
	void fillUniformStructure_interior(const char* fileName, const float isoValue);
	void fillUniformStructure(const char* fileName, const float isoValue); // TODO: this should modify latterly
	void fillAdaptiveStructure_interior(const char* fileName, const float isoValue, const int depth); // TODO: realize latterly
	void fillAdaptiveStructure_recursive_interior(const char* fileName, const float isoValue, const int depth);


};

////////////////////////////
// Inline VoxelData Stuff //
////////////////////////////
const inline float VoxelData::index(const int* idx,const float* dx){
	float temp;
	temp =array[idx[0]]*dx[0];
	temp+=array[idx[1]]*dx[1];
	temp+=array[idx[2]]*dx[2];
	temp+=array[idx[3]]*dx[3];
	temp+=array[idx[4]]*dx[4];
	temp+=array[idx[5]]*dx[5];
	temp+=array[idx[6]]*dx[6];
	temp+=array[idx[7]]*dx[7];
	return temp;
}
const inline void VoxelData::Position(const int res,const int padded,const float x,const float y,const float z,int* idx,float* dx){
	int xx=int(x);
	int yy=int(y);
	int zz=int(z);
	float ex=float(x-xx);
	float ey=float(y-yy);
	float ez=float(z-zz);
	int res2;
	if(padded){res2=res+2;}
	else{res2=res;}

	if(xx  >=0 && xx  <res && yy  >=0 && yy  <res && zz  >=0 && zz  <res){idx[0]=(xx  )*res*res2+(yy  )*res2+(zz  );dx[0]=(float)((1-ex)*(1-ey)*(1-ez));}
	else{idx[0]=0;dx[0]=0;}
	if(xx  >=0 && xx  <res && yy  >=0 && yy  <res && zz+1>=0 && zz+1<res){idx[1]=(xx  )*res*res2+(yy  )*res2+(zz+1);dx[1]=(float)((1-ex)*(1-ey)*(  ez));}
	else{idx[1]=0;dx[1]=0;}
	if(xx  >=0 && xx  <res && yy+1>=0 && yy+1<res && zz+1>=0 && zz+1<res){idx[2]=(xx  )*res*res2+(yy+1)*res2+(zz+1);dx[2]=(float)((1-ex)*(  ey)*(  ez));}
	else{idx[2]=0;dx[2]=0;}
	if(xx  >=0 && xx  <res && yy+1>=0 && yy+1<res && zz  >=0 && zz  <res){idx[3]=(xx  )*res*res2+(yy+1)*res2+(zz  );dx[3]=(float)((1-ex)*(  ey)*(1-ez));}
	else{idx[3]=0;dx[3]=0;}
	if(xx+1>=0 && xx+1<res && yy  >=0 && yy  <res && zz  >=0 && zz  <res){idx[4]=(xx+1)*res*res2+(yy  )*res2+(zz  );dx[4]=(float)((  ex)*(1-ey)*(1-ez));}
	else{idx[4]=0;dx[4]=0;}
	if(xx+1>=0 && xx+1<res && yy  >=0 && yy  <res && zz+1>=0 && zz+1<res){idx[5]=(xx+1)*res*res2+(yy  )*res2+(zz+1);dx[5]=(float)((  ex)*(1-ey)*(  ez));}
	else{idx[5]=0;dx[5]=0;}
	if(xx+1>=0 && xx+1<res && yy+1>=0 && yy+1<res && zz+1>=0 && zz+1<res){idx[6]=(xx+1)*res*res2+(yy+1)*res2+(zz+1);dx[6]=(float)((  ex)*(  ey)*(  ez));}
	else{idx[6]=0;dx[6]=0;}
	if(xx+1>=0 && xx+1<res && yy+1>=0 && yy+1<res && zz  >=0 && zz  <res){idx[7]=(xx+1)*res*res2+(yy+1)*res2+(zz  );dx[7]=(float)((  ex)*(  ey)*(1-ez));}
	else{idx[7]=0;dx[7]=0;}
}
#endif
