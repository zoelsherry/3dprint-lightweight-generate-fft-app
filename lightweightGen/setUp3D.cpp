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
#include <string.h>
#include <math.h>
#include "setUp3D.h"

#define PI 3.1415926535897932384

////////////////////
// Splatting Code //
////////////////////
void SplatValues(vector<OrientedVertex>& points,int* idx,float* dx,VoxelData &g,vector<float>& values){
	g.clear();
	int i=0,c=points.size()<<3;
	while(i<c){
		float v=Length(points[i>>3].n)*values[i>>3];
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
	}
}
void SplatPositions(vector<OrientedVertex>& points,int* idx,float* dx,VoxelData& g){
	g.clear();
	int i=0,c=points.size()<<3;
	while(i<c){
		float v=Length(points[i>>3].n);
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
		g.array[idx[i]]+=dx[i++]*v;
	}
}
void SplatNormalsX(vector<OrientedVertex>& points,int* idx,float* dx,VoxelData& X){
	X.clear();
	int i=0,c=points.size()<<3;
	while(i<c){
		float v=points[i>>3].n[0];
		X.array[idx[i]]+=dx[i++]*v;
		X.array[idx[i]]+=dx[i++]*v;
		X.array[idx[i]]+=dx[i++]*v;
		X.array[idx[i]]+=dx[i++]*v;
		X.array[idx[i]]+=dx[i++]*v;
		X.array[idx[i]]+=dx[i++]*v;
		X.array[idx[i]]+=dx[i++]*v;
		X.array[idx[i]]+=dx[i++]*v;
	}
}
void SplatNormalsY(vector<OrientedVertex>& points,int* idx,float* dx,VoxelData& Y){
	Y.clear();
	int i=0,c=points.size()<<3;
	while(i<c){
		float v=points[i>>3].n[1];
		Y.array[idx[i]]+=dx[i++]*v;
		Y.array[idx[i]]+=dx[i++]*v;
		Y.array[idx[i]]+=dx[i++]*v;
		Y.array[idx[i]]+=dx[i++]*v;
		Y.array[idx[i]]+=dx[i++]*v;
		Y.array[idx[i]]+=dx[i++]*v;
		Y.array[idx[i]]+=dx[i++]*v;
		Y.array[idx[i]]+=dx[i++]*v;
	}

}
void SplatNormalsZ(vector<OrientedVertex>& points,int* idx,float* dx,VoxelData& Z){
	Z.clear();
	int i=0,c=points.size()<<3;
	while(i<c){
		float v=points[i>>3].n[2];
		Z.array[idx[i]]+=dx[i++]*v;
		Z.array[idx[i]]+=dx[i++]*v;
		Z.array[idx[i]]+=dx[i++]*v;
		Z.array[idx[i]]+=dx[i++]*v;
		Z.array[idx[i]]+=dx[i++]*v;
		Z.array[idx[i]]+=dx[i++]*v;
		Z.array[idx[i]]+=dx[i++]*v;
		Z.array[idx[i]]+=dx[i++]*v;
	}
}

/////////////////
// Wisdom Code //
/////////////////
int ReadWisdom(char* fileName,int bw,std::vector<int>& supportedBW){
	FILE* fp=fopen(fileName,"r");
	if(!fp){return 0;}
	int i,cnt;
	fscanf(fp," %d ",&cnt);
	if(!cnt){
		fclose(fp);
		return 0;
	}
	supportedBW.resize(cnt);
	for(i=0;i<cnt;i++){fscanf(fp," %d ",&supportedBW[i]);}
	if(!fftwf_import_wisdom_from_file(fp)){
		fclose(fp);
		return 0;
	}
	fclose(fp);
	for(i=0;i<supportedBW.size();i++){if(supportedBW[i]==bw){return 1;}}
	return 0;
}
int WriteWisdom(char* fileName,int bw,std::vector<int>& supportedBW){
	int i;
	FILE* fp=fopen(fileName,"w");
	if(!fp){return 0;}
	for(i=0;i<supportedBW.size();i++){if(supportedBW[i]==bw){break;}}
	if(i==supportedBW.size()){supportedBW.push_back(bw);}
	fprintf(fp,"%d\n",supportedBW.size());
	for(i=0;i<supportedBW.size();i++){fprintf(fp," %d",supportedBW[i]);}
	fprintf(fp,"\n");
	fftwf_export_wisdom_to_file(fp);
	fclose(fp);
	return 1;
}

int InitWisdom(int bw,int mult,int fftwType){
	VoxelData v,filter;
	int r=2*bw;
	float *fIn1,*fIn2;
	fftwf_complex *out1,*out2,*fOut1,*fOut2;
	fftwf_plan plan,iplan1,iplan2,fplan1,fplan2;

	if(!v.set(r,1) || !filter.set(bw+1)){return 0;}
	out1 = (fftwf_complex*)v.array;
	out2 = (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*r*r*(bw+1));
	fOut1 = (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*(bw+1));
	fOut2= (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*(mult*bw+1));
	fIn1  = (float*)fftw_malloc(sizeof(float)*r);
	fIn2 = (float*)fftw_malloc(sizeof(float)*mult*r);
	if(!fIn1 || !fOut1 || !out2 || !fIn2 || !fOut2){
		if(fIn1){fftwf_free(fIn1);}
		if(fIn2){fftwf_free(fIn2);}
		if(fOut1){fftwf_free(fOut1);}
		if(fOut2){fftwf_free(fOut2);}
		if(out2){fftwf_free(out2);}
		return 0;
	}

	plan   = fftwf_plan_dft_r2c_3d(r,r,r,v.array,out1,fftwType);
	fplan1 = fftwf_plan_dft_r2c_1d(r,fIn1,fOut1,fftwType);
	fplan2 = fftwf_plan_dft_r2c_1d(mult*r,fIn2,fOut2,fftwType);
	iplan1 = fftwf_plan_dft_c2r_3d(r,r,r,out1,v.array,fftwType);
	iplan2 = fftwf_plan_dft_c2r_3d(r,r,r,out2,v.array,fftwType);

	fftwf_free(out2);
	fftwf_free(fIn1);
	fftwf_free(fIn2);
	fftwf_free(fOut1);
	fftwf_free(fOut2);

	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(fplan1);
	fftwf_destroy_plan(fplan2);
	fftwf_destroy_plan(iplan1);
	fftwf_destroy_plan(iplan2);
	return 1;
}

/////////////////
// Set Up Code //
/////////////////
int SetPositions(int res,int padded,vector<OrientedVertex>& points,int* &idx,float* &dx){
	int pNum=points.size();
	idx=new int[pNum*8];
	dx=new float[pNum*8];

	if(!idx || !dx){
		if(idx){delete[] idx;}
		if(dx){delete[] dx;}
		return 0;
	}
	int res2;
	if(padded){res2=res+2;}
	else{res2=res;}
	for(int i=0;i<pNum;i++){VoxelData::Position(res,padded,points[i].v[0],points[i].v[1],points[i].v[2],&idx[i*8],&dx[i*8]);}

	return 1;
}
void SetUpPoints(vector<OrientedVertex>& points,int bw,float inScale,float translate[3],float& outScale,float& volume){
	int i;
	float bRadius;
	Vertex center;
	int pNum=points.size();


	// Compute the center
	center.v[0]=center.v[1]=center.v[2]=0;
	for(i=0;i<pNum;i++){
		center.v[0]+=points[i].v[0];
		center.v[1]+=points[i].v[1];
		center.v[2]+=points[i].v[2];
	}
	center.v[0]/=pNum;
	center.v[1]/=pNum;
	center.v[2]/=pNum;

	for(i=0;i<pNum;i++){
		Vertex temp1;
		double temp2;
		temp1.v[0]=center.v[0]-points[i].v[0];
		temp1.v[1]=center.v[1]-points[i].v[1];
		temp1.v[2]=center.v[2]-points[i].v[2];
		temp2=SquareLength(temp1.v);
		if(!i || temp2>bRadius){bRadius=float(temp2);}
	}
	bRadius=float(sqrt(bRadius));


	outScale=inScale/bRadius*bw;
	translate[0]=-center.v[0]*outScale+bw;
	translate[1]=-center.v[1]*outScale+bw;
	translate[2]=-center.v[2]*outScale+bw;
	for(i=pNum-1;i>=0;i--){
		points[i].v[0]=points[i].v[0]*outScale+translate[0];
		points[i].v[1]=points[i].v[1]*outScale+translate[1];
		points[i].v[2]=points[i].v[2]*outScale+translate[2];
		if(points[i].v[0]<0 || points[i].v[0]>=2*bw || points[i].v[1]<0 || points[i].v[1]>2*bw || points[i].v[2]<0 || points[i].v[2]>2*bw){
			points[i]=points[points.size()-1];
			points.pop_back();
			continue;
		}
		points[i].n[0]*=outScale*outScale;
		points[i].n[1]*=outScale*outScale;
		points[i].n[2]*=outScale*outScale;
	}
	volume*=outScale*outScale*outScale;
}

int SetUpWeightedPoints(vector<OrientedVertex>& points,int* id,float* dx,int bw,float weight,int invertWeight,int fftwType){
	VoxelData v;
	float c,c1,c2;
	int idx,i,j,k,r=2*bw;
	float *fIn;
	fftwf_complex *out,*fOut;
	fftwf_plan plan,iplan,fplan;

	if(!v.set(r,1)){return 0;}
	out = (fftwf_complex*)v.array;
	fOut = (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*(bw+1));
	fIn  = (float*)fftw_malloc(sizeof(float)*r);
	if(!fIn || !fOut){
		if(fIn){fftwf_free(fIn);}
		if(fOut){fftwf_free(fOut);}
		return 0;
	}

	for(i=0;i<=bw;i++){
		fIn[i]=(float)exp(-i*i/((double)(2*weight*weight)));
		if(i){fIn[r-i]=fIn[i];}
	}
	c=0;
	for(i=0;i<r;i++){c+=fIn[i];}
	for(i=0;i<r;i++){fIn[i]/=c*r;}

	fplan = fftwf_plan_dft_r2c_1d(r,fIn,fOut,fftwType);
	fftwf_execute(fplan);
	fftw_free(fIn);
	fftwf_destroy_plan(fplan);

	SplatPositions(points,id,dx,v);

	plan  = fftwf_plan_dft_r2c_3d(r,r,r,v.array,out,fftwType);
	fftwf_execute(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){c1=fOut[r-i][0];}
		else{c1=fOut[i][0];}
		for(j=0;j<r;j++){
			if(j>bw){c2=c1*fOut[r-j][0];}
			else{c2=c1*fOut[j][0];}
			for(k=0;k<=bw;k++){
				c=c2*fOut[k][0];
				out[idx][0]*=c;
				out[idx][1]*=c;
				idx++;
			}
		}
	}
	iplan = fftwf_plan_dft_c2r_3d(r,r,r,out,v.array,fftwType);
	fftwf_execute(iplan);
	fftwf_free(fOut);
	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(iplan);

	c1=c2=0;
	for(i=0;i<points.size();i++){
		double temp=0;
		if(points[i].v[0]<0 || points[i].v[0]>=2*bw || points[i].v[1]<0 || points[i].v[1]>2*bw || points[i].v[2]<0 || points[i].v[2]>2*bw){;}
		else{
			c1+=(float)Length(points[i].n);
			temp=v.index(&id[i*8],&dx[i*8]);
			if(invertWeight){temp=1.0/temp;}
		}
		points[i].n[0]*=(float)temp;
		points[i].n[1]*=(float)temp;
		points[i].n[2]*=(float)temp;
		c2+=(float)Length(points[i].n);
	}
	c1/=c2;
	for(i=0;i<points.size();i++){
		points[i].n[0]*=c1;
		points[i].n[1]*=c1;
		points[i].n[2]*=c1;
	}
	return 1;
}


int SetUpWeightedPoints(vector<OrientedVertex>& points,int* idx,float* dx,int bw,float weight,int fftwType){return SetUpWeightedPoints(points,idx,dx,bw,weight,0,fftwType);}
int SetUpInvertedWeightedPoints(vector<OrientedVertex>& points,int* idx,float* dx,int bw,float weight,int fftwType){return SetUpWeightedPoints(points,idx,dx,bw,weight,1,fftwType);}

int SetUpNormalWeightedPoints(vector<OrientedVertex>& points,int* id,float* dx,int bw,float weight,int invertWeight,int fftwType){
	VoxelData v,temp,filter;
	float c,c1,c2;
	int idx,i,j,k,x,y,z,r=2*bw;
	float *fIn;
	fftwf_complex *out,*fOut;
	fftwf_plan plan,iplan,fplan;

	if(!v.set(r,1) || !temp.set(r,1) || !filter.set(bw+1)){return 0;}
	out    = (fftwf_complex*)temp.array;
	fIn    = (float*)fftwf_malloc(sizeof(float)*r);
	fOut   = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*(bw+1));
	if(!fIn || !fOut){
		if(fIn){fftwf_free(fIn);}
		if(fOut){fftwf_free(fOut);}
		return 0;
	}
	for(i=0;i<=bw;i++){
		fIn[i]=(float)exp(-i*i/((double)(2*weight*weight)));
		if(i){fIn[r-i]=fIn[i];}
	}
	c=0;
	for(i=0;i<r;i++){c+=fIn[i];}
	for(i=0;i<r;i++){fIn[i]/=c*r;}
	fplan = fftwf_plan_dft_r2c_1d(r,fIn,fOut,fftwType);
	fftwf_execute(fplan);
	fftw_free(fIn);
	fftwf_destroy_plan(fplan);
	for(i=0;i<=bw;i++){
		x=i*(bw+1)*(bw+1);
		for(j=0;j<=bw;j++){
			y=j*(bw+1);
			for(k=0;k<=bw;k++){
				z=k;
				filter.array[x+y+z]=fOut[i][0]*fOut[j][0]*fOut[k][0];
			}
		}
	}
	fftwf_free(fOut);

	plan  = fftwf_plan_dft_r2c_3d(r,r,r,temp.array,out,fftwType);
	iplan = fftwf_plan_dft_c2r_3d(r,r,r,out,temp.array,fftwType);

	// Next update with the value of the normal convolution
	// First do the X component
	SplatNormalsX(points,id,dx,temp);
	fftwf_execute(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){x=r-i;}
		else{x=i;}
		x*=(bw+1)*(bw+1);
		for(j=0;j<r;j++){
			if(j>bw){y=r-j;}
			else{y=j;}
			y*=(bw+1);
			for(k=0;k<=bw;k++){
				z=k;
				out[idx][0]*=filter.array[x+y+z];
				out[idx][1]*=filter.array[x+y+z];
				idx++;
			}
		}
	}
	fftwf_execute(iplan);
	for(i=0;i<(r+2)*r*r;i++){v.array[i]=temp.array[i]*temp.array[i];}

	// Next the Y component
	SplatNormalsY(points,id,dx,temp);
	fftwf_execute(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){x=r-i;}
		else{x=i;}
		x*=(bw+1)*(bw+1);
		for(j=0;j<r;j++){
			if(j>bw){y=r-j;}
			else{y=j;}
			y*=(bw+1);
			for(k=0;k<=bw;k++){
				z=k;
				out[idx][0]*=filter.array[x+y+z];
				out[idx][1]*=filter.array[x+y+z];
				idx++;
			}
		}
	}
	fftwf_execute(iplan);
	for(i=0;i<(r+2)*r*r;i++){v.array[i]+=temp.array[i]*temp.array[i];}

	// And finally the X component
	SplatNormalsZ(points,id,dx,temp);
	fftwf_execute(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){x=r-i;}
		else{x=i;}
		x*=(bw+1)*(bw+1);
		for(j=0;j<r;j++){
			if(j>bw){y=r-j;}
			else{y=j;}
			y*=(bw+1);
			for(k=0;k<=bw;k++){
				z=k;
				out[idx][0]*=filter.array[x+y+z];
				out[idx][1]*=filter.array[x+y+z];
				idx++;
			}
		}
	}
	fftwf_execute(iplan);
	for(i=0;i<(r+2)*r*r;i++){v.array[i]=(float)sqrt(temp.array[i]*temp.array[i]+v.array[i]);}

	fftwf_destroy_plan(plan);
	fftwf_destroy_plan(iplan);

	c1=c2=0;
	for(i=0;i<points.size();i++){
		double tmp=0;
		tmp=v.index(&id[i*8],&dx[i*8]);
		if(invertWeight){tmp=1.0/tmp;}
		c1+=(float)Length(points[i].n);
		points[i].n[0]*=float(tmp);
		points[i].n[1]*=float(tmp);
		points[i].n[2]*=float(tmp);
		c2+=(float)Length(points[i].n);
	}
	c1/=c2;
	for(i=0;i<points.size();i++){
		points[i].n[0]*=c1;
		points[i].n[1]*=c1;
		points[i].n[2]*=c1;
	}
	return 1;
}
int SetUpNormalWeightedPoints(vector<OrientedVertex>& points,int* idx,float* dx,int bw,float weight,int fftwType){return SetUpNormalWeightedPoints(points,idx,dx,bw,weight,0,fftwType);}
int SetUpNormalInvertedWeightedPoints(vector<OrientedVertex>& points,int* idx,float* dx,int bw,float weight,int fftwType){return SetUpNormalWeightedPoints(points,idx,dx,bw,weight,1,fftwType);}


//////////////////////////////////////////
// Compute Characteristic Function Code //
//////////////////////////////////////////
int MyGetInterior(vector<OrientedVertex>& points,int* id,float* dx,int bw,VoxelData& interior,int fftwType,float smooth,VoxelData* fCoefficients=NULL){
	VoxelData filter;
	int radius,norm,idx,i,j,k,x,xx,x1,y,yy,y1,z,zz,z1,r=2*bw;
	fftwf_complex *temp,*out,*fOut,*fOut1;
	fftwf_plan plan,iplan,fplan;
	float *fIn,*fIn1;
	float c;

	if(!interior.set(r,1) || !filter.set(bw+1)){return 0;}
	if(fCoefficients){
		if(!fCoefficients->set(r,1)){return 0;}
		out=(fftwf_complex*)fCoefficients->array;
	}
	else{out = (fftwf_complex*)fftw_malloc(sizeof(fftw_complex)*r*r*(bw+1));}

	temp = (fftwf_complex*)interior.array;
	fIn  = (float*)fftw_malloc(sizeof(float)*2*r);
	fOut = (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*(2*bw+1));
	fIn1 = (float*)fftw_malloc(sizeof(float)*r);
	fOut1= (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*(bw+1));

	if(!out || !fIn || !fOut || !fIn1 || !fOut1){
		if(out && !fCoefficients){fftw_free(out);}
		if(fIn){fftw_free(fIn);}
		if(fOut){fftw_free(fOut);}
		if(fIn1){fftw_free(fIn1);}
		if(fOut1){fftw_free(fOut1);}
		return 0;
	}
	if(smooth!=0){
		for(i=0;i<=bw;i++){
			fIn1[i]=(float)exp(-(float)(i*i)/(2.0*smooth*smooth));
			if(i){fIn1[r-i]=fIn1[i];}
		}
		// Make the filter have constant norm 1
		c=0;
		for(i=0;i<r;i++){c+=fIn1[i];}
		for(i=0;i<r;i++){fIn1[i]/=c;}
	}
	else{
		for(i=0;i<r;i++){fIn1[i]=0;}
		fIn1[0]=1;
	}
	fplan = fftwf_plan_dft_r2c_1d(r,fIn1,fOut1,fftwType);
	fftwf_execute(fplan);
	if(smooth<0){for(i=0;i<=bw;i++){fOut1[i][0]=1.f/fOut1[i][0];}}
	fftw_free(fIn1);
	fftwf_destroy_plan(fplan);

#if 1
	for(i=0;i<=bw*2;i++){
		fIn[i]=(float)exp(-i*i/((double)(4)));
		if(i){fIn[2*r-i]=fIn[i];}
	}
#else 
	for(i=0;i<r*2;i++){fIn[i]=0;}
	for(i=0;i<2;i++){
		fIn[i]=(float)(2-i);
		if(i){fIn[r*2-i]=fIn[i];}
	}
#endif
	// Make the filter have constant norm 1
	c=0;
	for(i=0;i<2*r;i++){c+=fIn[i];}
	for(i=0;i<2*r;i++){fIn[i]/=c;}
	fplan = fftwf_plan_dft_r2c_1d(2*r,fIn,fOut,fftwType);
	fftwf_execute(fplan);
	fftw_free(fIn);
	fftwf_destroy_plan(fplan);

	radius=4*bw*bw;
	float scale=1.0/(2.0*PI/r)/(r*r*r);
	for(x=0;x<=bw;x++){
		x1=x-r;
		for(y=0;y<=bw;y++){
			y1=y-r;
			float f1=fOut1[x][0]*fOut1[y][0];
			int norm11=x*x+y*y;
			int norm12=x*x+y1*y1;
			int norm21=x1*x1+y*y;
			int norm22=x1*x1+y1*y1;
			float f11=fOut[ x ][0]*fOut[ y ][0];
			float f12=fOut[ x ][0]*fOut[-y1][0];
			float f21=fOut[-x1][0]*fOut[ y ][0];
			float f22=fOut[-x1][0]*fOut[-y1][0];
			int idx1=x*(bw+1)*(bw+1)+y*(bw+1);
			int idx2=x*(bw+1)*(bw+1)+y;
			for(z=0;z<=y;z++){
				z1=z-r;
				if(x || y || z){
					double t1=0,t2=0;
					
					if(norm11<=radius){
						norm=norm11+z*z;
						if(norm<=radius){t1+=f11*fOut[ z ][0]/norm;}
						norm=norm11+z1*z1;
						if(norm<=radius){t1+=f11*fOut[-z1][0]/norm;}
					}
					if(norm12<=radius){
						norm=norm12+z1*z1;
						if(norm<=radius){t1+=f12*fOut[-z1][0]/norm;}
						norm=norm12+z*z;
						if(norm<=radius){t1+=f12*fOut[ z ][0]/norm;}
					}
					
					if(norm21<=radius){
						norm=norm21+z*z;
						if(norm<=radius){t2+=f21*fOut[ z ][0]/norm;}
						norm=norm21+z1*z1;
						if(norm<=radius){t2+=f21*fOut[-z1][0]/norm;}
					}
					if(norm22<=radius){
						norm=norm22+z1*z1;
						if(norm<=radius){t2+=f22*fOut[-z1][0]/norm;}
						norm=norm22+z*z;
						if(norm<=radius){t2+=f22*fOut[ z ][0]/norm;}
					}
					t1*=x;
					t2*=x1;
					
					filter.array[idx1+z]=filter.array[idx2+z*(bw+1)]=(t1+t2)*scale*f1*fOut1[z][0];
				}
			}
		}
	}
	fftwf_free(fOut);
	fftwf_free(fOut1);

	plan  = fftwf_plan_dft_r2c_3d(r,r,r,interior.array,temp,fftwType);
	SplatNormalsX(points,id,dx,interior);
	fftwf_execute(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){xx=r-i;x=-r+i;}
		else{xx=x=i;}
		xx*=(bw+1)*(bw+1);
		for(j=0;j<r;j++){
			if(j>bw){yy=r-j;}
			else{yy=j;}
			yy*=bw+1;
			for(k=0;k<=bw;k++){
				zz=k;
				float im;
				if(x<0){im=-filter.array[xx+yy+zz];}
				else{im=filter.array[xx+yy+zz];}
				out[idx][0]=-temp[idx][1]*im;
				out[idx][1]= temp[idx][0]*im;
				idx++;
			}
		}
	}
	SplatNormalsY(points,id,dx,interior);
	fftwf_execute(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){xx=r-i;}
		else{xx=i;}
		xx*=bw+1;
		for(j=0;j<r;j++){
			if(j>bw){yy=r-j;y=-r+j;}
			else{yy=y=j;}
			yy*=(bw+1)*(bw+1);
			for(k=0;k<=bw;k++){
				zz=k;
				float im;
				if(y<0){im=-filter.array[yy+xx+zz];}
				else{im=filter.array[yy+xx+zz];}
				out[idx][0]+=-temp[idx][1]*im;
				out[idx][1]+= temp[idx][0]*im;
				idx++;
			}
		}
	}
	SplatNormalsZ(points,id,dx,interior);
	fftwf_execute(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){xx=r-i;}
		else{xx=i;}
		xx*=bw+1;
		for(j=0;j<r;j++){
			if(j>bw){yy=r-j;}
			else{yy=j;}
			for(k=0;k<=bw;k++){
				zz=k;
				zz*=(bw+1)*(bw+1);
				float im;
				im=filter.array[zz+xx+yy];
				out[idx][0]+=-temp[idx][1]*im;
				out[idx][1]+= temp[idx][0]*im;
				idx++;
			}
		}
	}
	fftwf_destroy_plan(plan);
	if(fCoefficients){
		memcpy(interior.array,out,sizeof(fftwf_complex)*r*r*(bw+1));
		iplan = fftwf_plan_dft_c2r_3d(r,r,r,(fftwf_complex*)interior.array,interior.array,fftwType);
		fftwf_execute(iplan);
	}
	else{
		iplan = fftwf_plan_dft_c2r_3d(r,r,r,out,interior.array,fftwType);
		fftwf_execute(iplan);
		fftw_free(out);
		interior.padded=0;
	}
	fftwf_destroy_plan(iplan);
	return 1;
}
int GetInterior(vector<OrientedVertex>& points,int* id,float* dx,int bw,VoxelData& interior,int fftwType,float smooth){
	return MyGetInterior(points,id,dx,bw,interior,fftwType,smooth);
}
int GetInterior(vector<OrientedVertex>& points,int* id,float* dx,int bw,VoxelData& interior,VoxelData& bSplineValues,int fftwType,float smooth){
	int idx,i,j,k;
	float c,c1,c2;
	float *fIn;
	fftwf_complex *out,*fOut;
	fftwf_plan iplan,fplan;

	if(!MyGetInterior(points,id,dx,bw,interior,fftwType,smooth,&bSplineValues)){return 0;}
	out  = (fftwf_complex*)bSplineValues.array;

	fIn  = (float*)fftw_malloc(sizeof(float)*2*bw);
	fOut = (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*(bw+1));
	if(!fIn || !fOut){
		if(fIn){fftwf_free(fIn);}
		if(fOut){fftwf_free(fOut);}
		return 0;
	}
	
	for(i=0;i<2*bw;i++){fIn[i]=0;}
	fIn[0]=4.f/6;
	fIn[1]=fIn[2*bw-1]=1.f/6;
	fplan = fftwf_plan_dft_r2c_1d(2*bw,fIn,fOut,fftwType);
	fftwf_execute(fplan);
	fftw_free(fIn);
	fftwf_destroy_plan(fplan);
	idx=0;
	for(i=0;i<2*bw;i++){
		if(i>bw){c1=fOut[2*bw-i][0];}
		else{c1=fOut[i][0];}
		for(j=0;j<2*bw;j++){
			if(j>bw){c2=c1*fOut[2*bw-j][0];}
			else{c2=c1*fOut[j][0];}
			for(k=0;k<=bw;k++){
				c=c2*fOut[k][0];
				out[idx][0]*=c;
				out[idx][1]*=c;
				idx++;
			}
		}
	}
	fftwf_free(fOut);

	iplan = fftwf_plan_dft_c2r_3d(2*bw,2*bw,2*bw,out,bSplineValues.array,fftwType);
	fftwf_execute(iplan);
	fftwf_destroy_plan(iplan);
	return 1;

}
int GetWeightedAverage(vector<OrientedVertex>& points,vector<float>& sampleValues,int* id,float* dx,int bw,VoxelData& interior,int fftwType,float smooth){
	VoxelData v;
	float c,c1,c2;
	int idx,i,j,k,r=2*bw;
	float *fIn;
	fftwf_complex *outI,*outV,*fOut;
	fftwf_plan plan,iplan,fplan;

	if(!v.set(r,1) || !interior.set(r,1)){return 0;}
	outI = (fftwf_complex*)interior.array;
	outV = (fftwf_complex*)v.array;
	fOut = (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex)*(bw+1));
	fIn  = (float*)fftw_malloc(sizeof(float)*r);
	if(!fIn || !fOut){
		if(fIn){fftwf_free(fIn);}
		if(fOut){fftwf_free(fOut);}
		return 0;
	}

	for(i=0;i<=bw;i++){
		fIn[i]=(float)exp(-i*i/((double)(2*smooth*smooth)))+smooth/r;
		if(i){fIn[r-i]=fIn[i];}
	}
	fplan = fftwf_plan_dft_r2c_1d(r,fIn,fOut,fftwType);
	fftwf_execute(fplan);
	fftwf_destroy_plan(fplan);
	fftw_free(fIn);

	// Get the weighted sum of the sample values
	SplatValues(points,id,dx,interior,sampleValues);
	plan  = fftwf_plan_dft_r2c_3d(r,r,r,interior.array,outI,fftwType);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){c1=fOut[r-i][0];}
		else{c1=fOut[i][0];}
		for(j=0;j<r;j++){
			if(j>bw){c2=c1*fOut[r-j][0];}
			else{c2=c1*fOut[j][0];}
			for(k=0;k<=bw;k++){
				c=c2*fOut[k][0];
				outI[idx][0]*=c;
				outI[idx][1]*=c;
				idx++;
			}
		}
	}
	iplan = fftwf_plan_dft_c2r_3d(r,r,r,outI,interior.array,fftwType);
	fftwf_execute(iplan);
	fftwf_destroy_plan(iplan);

	// Get the sum of the weights used
	SplatPositions(points,id,dx,v);
	plan  = fftwf_plan_dft_r2c_3d(r,r,r,v.array,outV,fftwType);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	idx=0;
	for(i=0;i<r;i++){
		if(i>bw){c1=fOut[r-i][0];}
		else{c1=fOut[i][0];}
		for(j=0;j<r;j++){
			if(j>bw){c2=c1*fOut[r-j][0];}
			else{c2=c1*fOut[j][0];}
			for(k=0;k<=bw;k++){
				c=c2*fOut[k][0];
				outV[idx][0]*=c;
				outV[idx][1]*=c;
				idx++;
			}
		}
	}
	iplan = fftwf_plan_dft_c2r_3d(r,r,r,outV,v.array,fftwType);
	fftwf_execute(iplan);
	fftwf_destroy_plan(iplan);

	fftwf_free(fOut);
	for(i=0;i<r;i++){for(j=0;j<r;j++){for(k=0;k<r;k++){
		if(v.array[(i*r+j)*(r+2)+k]<=0){
			printf("Argggg1: %d %d %d\n",i,j,k);
		}
		if(interior.array[(i*r+j)*(r+2)+k]<=0){
			printf("Argggg2: %d %d %d\n",i,j,k);
		}
		interior.array[(i*r+j)*(r+2)+k]/=v.array[(i*r+j)*(r+2)+k];
		if(interior.array[(i*r+j)*(r+2)+k]<=0){
			printf("Argggg3: %d %d %d\n",i,j,k);
		}
	}}}
	return 1;
}

double AverageValue(vector<OrientedVertex>& points,VoxelData& v,AdjacencyTable& table,int intType){
	double ll,l=0,a=0;
	for(int i=0;i<points.size();i++){
		ll=Length(points[i].n);
		a+=v.VertexValue(points[i].v,table,intType)*ll;
		l+=ll;
	}
	return a/l;
}
