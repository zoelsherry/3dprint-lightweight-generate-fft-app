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
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <sys/timeb.h>
#include <fstream>
#include <iostream>
#ifndef WIN32
#include <sys/time.h>
#endif
#include "setUp3D.h"
#include "cmdLineParser.h"
#include "spline.h"
#include "ply.h"

using std::cout;
using std::endl;

double GetTime(void){
#ifdef WIN32
	struct _timeb t;
	_ftime(&t);
	return t.time+t.millitm/1000.0;
#else
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+(double)t.tv_usec/1000000;
#endif
}
void ShowUsage(char* ex){
	printf("Usage %s:\n",ex);
	printf("\t--in  <input points>\n");
	printf("\t--out <ouput triangle mesh>\n");
	printf("\t--res <reconstruction resolution>: should be a power of 2)\n");
	printf("\t[--smooth <width>]: the width, in voxels, of the smoothing Gaussian\n");
	printf("\t[--weight <width>]: the width, in voxels, of the re-weighting Gaussian\n"); 
	printf("\t[--noNormalize]: specifies that the sizes of the input normals\n");
	printf("\t\tcorrespond to the sampling density\n");
	printf("\t[--scale <scale factor>]: specifies the factor of the bounding cube\n");
	printf("\t\tthat the input samples should fit into\n");
	printf("\t[--order <interpolation type>]: specifies the type of function that\n");
	printf("\t\tshould be used to interpolate/approximate the voxel data:\n");
	printf("\t\t   [1] -- Linear Fitting\n");
	printf("\t\t   [2] -- Quadratic Fitting\n");
	printf("\t\t   [3] -- Cubic  Fitting\n");
	printf("\t\t   [4] -- Catmull-Rom Fitting\n");
	printf("\t\t   [5] -- Uniform Cubic B-Spline Fitting\n");
	printf("\t[--normal]: specifies that surface normals should be written out with\n");
	printf("\t\twith the triangle mesh\n");
	printf("\t[--fftwWisdom <wisdom file>]: the name of the file into which FFTW\n");
	printf("\t\tshould write in the wisdom it accumulates as it sets up the\n");
	printf("\t\tplans. (The first time this argument is used, the program will\n");
	printf("\t\ttake a long time to run. Also, you can run the program without\n");
	printf("\t\tspecifying input and output files, as long as a resolution is\n");
	printf("\t\tspecifed.)\n");
}

int WriteTriangles(char* fileName,vector<TriangleIndex>& triangles,vector<Vertex>& vertices,vector<Vertex>& normals){
	int i,j,ret;
	FILE* fp;
	char* extension;
	
	extension=GetFileExtension(fileName);
	if(!strcasecmp(extension, "ply")){
		ret=PlyWriteTriangles(fileName,triangles,vertices,normals,PlyDefaultFileType());
	}
	else{		
		fp=fopen(fileName,"w");
		if(!fp){ret=0;}
		else{
			for(i=0;i<triangles.size();i++){
				// Write the triangle vertices
				for(j=0;j<3;j++){
					fprintf(fp,"%f %f %f",vertices[triangles[i].idx[j]].v[0],vertices[triangles[i].idx[j]].v[1],vertices[triangles[i].idx[j]].v[2]);
					if(j<2){fprintf(fp,"\t");}
					else{fprintf(fp,"\n");}
				}
				for(j=0;j<3;j++){
					fprintf(fp,"%f %f %f",normals[triangles[i].idx[j]].v[0],normals[triangles[i].idx[j]].v[1],normals[triangles[i].idx[j]].v[2]);
					if(j<2){fprintf(fp,"\t");}
					else{fprintf(fp,"\n");}
				}
			}
			fclose(fp);
		}
	}
	delete[] extension;
	return ret;
}
int WriteTriangles(char* fileName,vector<TriangleIndex>& triangles,vector<Vertex>& vertices){
	int ret;
	FILE* fp;
	char* extension;
	
	extension=GetFileExtension(fileName);
	if(!strcasecmp(extension,"ply")){
		ret=PlyWriteTriangles(fileName,triangles,vertices,PlyDefaultFileType());
	}
	else{		
		fp=fopen(fileName,"w");
		if(!fp){ret=0;}
		else{
			for(int i=0;i<triangles.size();i++){
				// Write the triangle vertices
				for(int j=0;j<3;j++){
					fprintf(fp,"%f %f %f",vertices[triangles[i].idx[j]].v[0],vertices[triangles[i].idx[j]].v[1],vertices[triangles[i].idx[j]].v[2]);
					if(j<2){fprintf(fp,"\t");}
					else{fprintf(fp,"\n");}
				}
			}
			fclose(fp);
		}
	}
	delete[] extension;
	return ret;
}
int ReadOrientedPoints(char* fileName,vector<OrientedVertex>& points){
	char* extension;
	FILE* fp;
	OrientedVertex v;
	int ret;

	points.clear();
	extension=GetFileExtension(fileName);
	if(!strcasecmp(extension,"ply")){
		int ft;
		ret=PlyReadOrientedPoints(fileName,points,ft);
	}
	else{
		fp=fopen(fileName,"r");
		if(!fp){ret=0;}
		else{
			while(fscanf(fp," %f %f %f %f %f %f",&v.v[0],&v.v[1],&v.v[2],&v.n[0],&v.n[1],&v.n[2])==6){points.push_back(v);}
			fclose(fp);
		}
	}
	delete[] extension;
	return ret;
}

int GetCount(int dim,int order){
	int cnt=0;
	if(!order || dim==1){return 1;}
	if(order==1){return dim;}
	for(int i=0;i<dim;i++){cnt+=GetCount(dim-i,order-1);}
	return cnt;
}

int main(int argc,char* argv[]){
	std::vector<int> supportedBW;
	float translate[3],outScale;
	VoxelData v,cubicBSplineValues;
	double t,tt,avg;
	int i,bw,fftwType,order;
	AdjacencyTable table;
	vector<OrientedVertex> pts;
	vector<TriangleIndex> triangles;
	vector<Vertex> vertices;
	vector<Vertex> normals;
	float inScale,smooth=-1;
	float volume;
	int* idx;
	float* dx;
	char* paramNames[]={"res","in","iso","out","scale","weight","smooth","noNormalize","order","noInvert","volume","normal","fftwWisdom","cell"};
	cmdLineReadable* params[14];
	cmdLineInt Res,Order;
	cmdLineString In,Iso,Out,FFTWWisdom,Cell;
	cmdLineFloat Scale,Weight,Smooth,Volume;
	cmdLineReadable NoNormalize,NoInvert,Normal;

	params[0] = &Res;
	params[1] = &In;
	params[2] = &Iso;
	params[3] = &Out;
	params[4] = &Scale;
	params[5] = &Weight;
	params[6] = &Smooth;
	params[7] = &NoNormalize;
	params[8] = &Order;
	params[9] = &NoInvert;
	params[10] = &Volume;
	params[11] = &Normal;
	params[12] = &FFTWWisdom;
	params[13] = &Cell;
	cmdLineParse(argc-1,&argv[1],paramNames,14,params,0);

	if (Cell.set)
	{
		printf("Enable lightweight structrue production!\n");
	}


	if(!Res.set){
		printf("No reconstruction resolution specified\n");
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}
	bw=Res.value/2;

	fftwType=FFTW_ESTIMATE;
	if(FFTWWisdom.set){
		fftwType=FFTW_PATIENT;
		tt=GetTime();
		if(!ReadWisdom(FFTWWisdom.value,bw,supportedBW)){
			InitWisdom(bw,2,fftwType);
			if(!WriteWisdom(FFTWWisdom.value,bw,supportedBW)){printf("Could not export FFTW wisdom to: %s\n",FFTWWisdom.value);}
		}
		printf("FFTW Wisdom acquisition time: %f\n",GetTime()-tt);
	}

	if(!In.set){
		printf("No input file specified\n");
		ShowUsage(argv[0]);
		return EXIT_FAILURE;
	}

	if(!ReadOrientedPoints(In.value,pts)){
		printf("Failed to read: %s\n",In.value);
		return EXIT_FAILURE;
	}
	if(Order.set){
		switch(Order.value){
		case 1:
			order=Spline::LINEAR;
			break;
		case 2:
			order=Spline::QUADRATIC;
			break;
		case 3:
			order=Spline::CUBIC;
			break;
		case 4:
			order=Spline::CATMULL_ROM;
			break;
		case 5:
			order=Spline::CUBIC_B_SPLINE;
			break;
		default:
			printf("Only linear, quadratic, and cubic interpolation are supported\n");
			ShowUsage(argv[0]);
			return EXIT_FAILURE;
		}
	}
	else{order=Spline::QUADRATIC;}

	if(Smooth.set){smooth=Smooth.value;}
	else{smooth=0;}
	if(Scale.set){inScale=Scale.value;}
	else{inScale=1;}
	if(Volume.set){volume=Volume.value;}
	else{volume=1;}
	if(!NoNormalize.set && !Volume.set){
		int cnt=0;
		for(i=0;i<pts.size();i++){
			double l=(float)Length(pts[i].n);
			if(l<10e-8 || l>10e8){
				pts[i]=pts.back();
				pts.pop_back();
			}
			else{
				pts[i].n[0]/=l;
				pts[i].n[1]/=l;
				pts[i].n[2]/=l;
				cnt++;
			}
		}
		for(i=0;i<pts.size();i++){
			pts[i].n[0]/=cnt;
			pts[i].n[1]/=cnt;
			pts[i].n[2]/=cnt;
		}
	}
	t=GetTime();
	SetUpPoints(pts,bw,inScale,translate,outScale,volume);
	if(!SetPositions(2*bw,1,pts,idx,dx)){
		printf("Failed to set positions\n");
		return EXIT_FAILURE;
	}
	if(Weight.set){
		if(Weight.value>0){
			if(NoInvert.set){SetUpWeightedPoints(pts,idx,dx,bw,Weight.value,fftwType);}
			else{SetUpInvertedWeightedPoints(pts,idx,dx,bw,Weight.value,fftwType);}
		}
		else{
			if(NoInvert.set){SetUpNormalWeightedPoints(pts,idx,dx,bw,-Weight.value,fftwType);}
			else{SetUpNormalInvertedWeightedPoints(pts,idx,dx,bw,-Weight.value,fftwType);}
		}
	}

	tt=GetTime()-t;
	printf("Set up time: %f\n",GetTime()-t);

	t=GetTime();

	if(order==Spline::CUBIC_B_SPLINE){
		if(!GetInterior(pts,idx,dx,bw,v,cubicBSplineValues,fftwType,smooth)){
			printf("Failed to compute interior\n");
			return EXIT_FAILURE;
		}
	}
	else{
		if(!GetInterior(pts,idx,dx,bw,v,fftwType,smooth)){
			printf("Failed to compute interior\n");
			return EXIT_FAILURE;
		}
	}
	if(v.padded){
		if(!table.set(Res.value,1,0)){
			printf("Could not allocated memory for adjacency table\n");
			return EXIT_FAILURE;
		}
	}
	else{
		if(!table.set(Res.value,0,0)){
			printf("Could not allocated memory for adjacency table\n");
			return EXIT_FAILURE;
		}
	}
	if(Volume.set){avg=-volume/(8*bw*bw*bw);}
	else{avg=AverageValue(pts,v,table,order);}
	delete[] idx;
	delete[] dx;


	int r;
	r=8*bw*bw*bw;
	if(v.padded){r+=8*bw*bw;}
	for(i=0;i<r;i++){v.array[i]-=(float)avg;}
	if(order==Spline::CUBIC_B_SPLINE){
		r=8*bw*bw*bw;
		if(cubicBSplineValues.padded){r+=8*bw*bw;}
		for(i=0;i<r;i++){cubicBSplineValues.array[i]-=(float)avg;}
	}

	tt+=GetTime()-t;
	printf("Iso-Function Time: %f\n",GetTime()-t);
#if (0)
	cout << "output the debug data..." << endl;
	std::ofstream fout;
	fout.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	try
	{
		fout.open("./outdata/array.txt");

		for (int i = 0; i < r; i++)
			fout << v.array[i] << endl;

		fout.close();
	}
	catch (std::ofstream::failure e)
	{
		std::cout << "ERROR::SIGNEDDIST::FILE FAIL TO WRITE!" << std::endl;
	}
#endif
	//if(Iso.set){v.write(Iso.value);}
	//if(Out.set){
	//	t=GetTime();
	//	if(order==Spline::CUBIC_B_SPLINE){v.isoSurfaceCubicBSpline(cubicBSplineValues,triangles,vertices,0,table,fftwType);}
	//	else{v.isoSurface(triangles,vertices,0,table,order,fftwType);} // extract isosurface
	//	
	//	if(Normal.set){
	//		normals.resize(vertices.size());
	//		for(i=0;i<vertices.size();i++){v.VertexGradient(vertices[i].v,normals[i].v,table,order);}
	//	}
	//	tt+=GetTime()-t;
	//	printf("Iso-Extraction Time: %f\n",GetTime()-t);
	//	
	//	if(Normal.set){
	//		for(i=0;i<vertices.size();i++){
	//			normals[i].v[0]=-normals[i].v[0]/outScale;
	//			normals[i].v[1]=-normals[i].v[1]/outScale;
	//			normals[i].v[2]=-normals[i].v[2]/outScale;
	//			vertices[i].v[0]=(vertices[i].v[0]-translate[0])/outScale;
	//			vertices[i].v[1]=(vertices[i].v[1]-translate[1])/outScale;
	//			vertices[i].v[2]=(vertices[i].v[2]-translate[2])/outScale;
	//		}
	//		for(i=0;i<triangles.size();i++){ // TODO: reserve the normal?
	//			int temp=triangles[i].idx[0];
	//			triangles[i].idx[0]=triangles[i].idx[2];
	//			triangles[i].idx[2]=temp;
	//		}
	//		printf("Triangle Count: %d\n",triangles.size());
	//		WriteTriangles(Out.value,triangles,vertices,normals); // write .ply file
	//	}
	//	else{
	//		//for(i=0;i<vertices.size();i++){
	//		//	vertices[i].v[0]=(vertices[i].v[0]-translate[0])/outScale;
	//		//	vertices[i].v[1]=(vertices[i].v[1]-translate[1])/outScale;
	//		//	vertices[i].v[2]=(vertices[i].v[2]-translate[2])/outScale;
	//		//}
	//		//for(i=0;i<triangles.size();i++){ // TODO: reserve the normal?
	//		//	int temp=triangles[i].idx[0];
	//		//	triangles[i].idx[0]=triangles[i].idx[2];
	//		//	triangles[i].idx[2]=temp;
	//		//}
	//		//printf("Triangle Count: %d\n",triangles.size());
	//		//WriteTriangles(Out.value,triangles,vertices); // write .ply file
	//	}
	//}
	////printf("Total Time: %f\n",tt);

	if(FFTWWisdom.set){if(!WriteWisdom(FFTWWisdom.value,bw,supportedBW)){printf("Could not export FFTW wisdom to: %s\n",FFTWWisdom.value);}}

	
	
	////////// fill lightweight structure ///////////
	// TODO: read the voxel data from v.arrary

	t = GetTime();
	v.fillUniformStructure_interior(Cell.value, 0);
	tt += GetTime() - t;
	printf("Fill Structure Time: %f\n", GetTime() - t);

	t = GetTime();
	v.isoSurface(triangles, vertices, 0, table, order, fftwType);
	tt += GetTime() - t;
	printf("Iso-Extraction Time: %f\n", GetTime() - t);


	for (i = 0; i < vertices.size(); i++) {
		vertices[i].v[0] = (vertices[i].v[0] - translate[0]) / outScale;
		vertices[i].v[1] = (vertices[i].v[1] - translate[1]) / outScale;
		vertices[i].v[2] = (vertices[i].v[2] - translate[2]) / outScale;
	}
	for (i = 0; i < triangles.size(); i++) { // TODO: reserve the normal?
		int temp = triangles[i].idx[0];
		triangles[i].idx[0] = triangles[i].idx[2];
		triangles[i].idx[2] = temp;
	}
	printf("Triangle Count: %d\n", triangles.size());
	WriteTriangles(Out.value, triangles, vertices); // write .ply file

	printf("Total Time: %f\n",tt);


	return EXIT_SUCCESS;
}
