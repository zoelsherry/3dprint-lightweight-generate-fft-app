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
#ifndef MARCHING_CUBES_INCLUDED
#define MARCHING_CUBES_INCLUDED
#include <vector>
using namespace std;

class EdgeFunction{
public:
	virtual const int getIndex(const int idx,const int e)=0;
};

class TriangleIndex{
public:
	int idx[3];
};
class MarchingCubes{
public:
	static const int edges[256];
	static const int triangles[256][16];
	static int vertexIndexList[12];

	static int AddTriangles(
		double (*z1)[2],
		double (*z2)[2],
		double isoValue,
		EdgeFunction* f,
		int idx,
		vector<TriangleIndex>& triangles,
		int& zIdx,
		int useZFront
		);
};
class MarchingCubesWriter{
public:
	static int vertices1[8][3];
	static int vertices2[2][2][2];
	static int edges1[12][2];
	static int edges2[8][8];
	static void writeEdges(FILE* fp);
	static void writeTriangles(FILE* fp);
	
	static int Flip(int in);

	static void RotateX(int in[3],int out[3]);
	static void RotateY(int in[3],int out[3]);
	static void RotateZ(int in[3],int out[3]);
	static int  RotateX(int in);
	static int  RotateY(int in);
	static int  RotateZ(int in);
	static int  Antipode(int in);

	static int  Transform (int idx,int x,int y,int z,int a);
	static int  ITransform(int idx,int x,int y,int z,int a);
	static int  ITransformEdge(int e,int x,int y,int z,int a);
	static void ITransformTriangle(int* tIdx,int x,int y,int z,int a);
	static void ITransformTriangle(int* tIdx,int x,int y,int z,int a,int f);

	static int MinIndex(int in,int& x,int& y,int& z,int& a,int &f);
	static int MinIndex(int in,int& x,int& y,int& z,int &f);
	static int GetTriangles(int in,int* triangleIndices);

};
#endif //MARCHING_CUBES_INCLUDED