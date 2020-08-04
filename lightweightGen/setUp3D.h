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
#ifndef SET_UP_INCLUDED
#define SET_UP_INCLUDED
#include <vector>
#include <fftw3.h>
#include "voxelData.h"

using namespace std;


class OrientedVertex{
public:
	float v[3];
	float n[3];
};

double SquareLength(const float v[3]);
double Length(const float v[3]);

int ReadWisdom(char* fileName,int bw,std::vector<int>& supportedBW);
int WriteWisdom(char* fileName,int bw,std::vector<int>& supportedBW);
int InitWisdom(int bw,int mult,int fftwType);

int SetPositions(int res,int padded,vector<OrientedVertex>& points,int* &idx,float* &dx);

void SetUpPoints(vector<OrientedVertex>& points,int bw,float inScale,float translate[3],float& outScale,float& volume);

int SetUpWeightedPoints					(vector<OrientedVertex>& points,int* idx,float* dx,int bw,float weight,int fftwType);
int SetUpInvertedWeightedPoints			(vector<OrientedVertex>& points,int* idx,float* dx,int bw,float weight,int fftwType);
int SetUpNormalWeightedPoints			(vector<OrientedVertex>& points,int* idx,float* dx,int bw,float weight,int fftwType);
int SetUpNormalInvertedWeightedPoints	(vector<OrientedVertex>& points,int* idx,float* dx,int bw,float weight,int fftwType);

void SplatPositions(vector<OrientedVertex>& points,int* idx,float* dx,VoxelData &g);
void SplatNormalsX (vector<OrientedVertex>& points,int* idx,float* dx,VoxelData &g);
void SplatNormalsY (vector<OrientedVertex>& points,int* idx,float* dx,VoxelData &g);
void SplatNormalsZ (vector<OrientedVertex>& points,int* idx,float* dx,VoxelData &g);
void SplatValues   (vector<OrientedVertex>& points,int* idx,float* dx,VoxelData &g,vector<float>& values);


int GetInterior(vector<OrientedVertex>& points,int* idx,float* dx,int bw,VoxelData& interior,int fftwType,float smooth=0);
int GetInterior(vector<OrientedVertex>& points,int* idx,float* dx,int bw,VoxelData& interior,VoxelData& bSplineValues,int fftwType,float smooth=0);

double AverageValue(vector<OrientedVertex>& points,VoxelData& v,AdjacencyTable& table,int intType);

int GetWeightedAverage(vector<OrientedVertex>& points,vector<float>& sampleValues,int* idx,float* dx,int bw,VoxelData& interior,int fftwType,float smooth);

///////////////////
// Inline Length //
///////////////////
inline double SquareLength(const float v[3]){return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];}
inline double Length(const float v[3]){return sqrt(SquareLength(v));}


#endif
