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

#include "ply.h"

//
// PLY data structures
//
char *elem_names[] = { "vertex", "face" };

typedef struct PlyVertex {
	float x, y, z;
} PlyVertex;

typedef struct PlyOrientedVertex {
	float x, y, z , nx, ny, nz;
} PlyOrientedVertex;

typedef struct PlyFace {
	unsigned char nr_vertices;
	int *vertices;
	int segment;
} PlyFace;

static PlyProperty vert_props[] = {
	{"x", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,x), 0, 0, 0, 0},
	{"y", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,y), 0, 0, 0, 0},
	{"z", PLY_FLOAT, PLY_FLOAT, offsetof(PlyVertex,z), 0, 0, 0, 0}
};
static PlyProperty oriented_vert_props[] = {
	{"x",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,x ), 0, 0, 0, 0},
	{"y",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,y ), 0, 0, 0, 0},
	{"z",  PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,z ), 0, 0, 0, 0},
	{"nx", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,nx), 0, 0, 0, 0},
	{"ny", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,ny), 0, 0, 0, 0},
	{"nz", PLY_FLOAT, PLY_FLOAT, offsetof(PlyOrientedVertex,nz), 0, 0, 0, 0}
};

// List of property information for a vertex
static PlyProperty face_props[] = {
	{"vertex_indices", PLY_INT, PLY_INT, offsetof(PlyFace,vertices),
		1, PLY_UCHAR, PLY_UCHAR, offsetof(PlyFace,nr_vertices)},
};

int PlyDefaultFileType(void){return PLY_ASCII;}

int PlyReadOrientedPoints(char* fileName,vector<OrientedVertex>& points,int& file_type){
	int nr_elems;
	char **elist;
	float version;
	int i,j;
	PlyFile* ply;
	char* elem_name;
	int num_elems;
	int nr_props;
	PlyProperty** plist;
	PlyOrientedVertex ply_vertex;

	ply = ply_open_for_reading(fileName, &nr_elems, &elist, &file_type, &version);
	if(!ply){return 0;}
	
	for (i=0; i < nr_elems; i++) {
		elem_name = elist[i];
		plist = ply_get_element_description(ply, elem_name, &num_elems, &nr_props);
		if(!plist){
			for(i=0;i<nr_elems;i++){
				free(ply->elems[i]->name);
				free(ply->elems[i]->store_prop);
				for(j=0;j<ply->elems[i]->nprops;j++){
					free(ply->elems[i]->props[j]->name);
					free(ply->elems[i]->props[j]);
				}
				free(ply->elems[i]->props);
			}
			for(i=0;i<nr_elems;i++){free(ply->elems[i]);}
			free(ply->elems);
			for(i=0;i<ply->num_comments;i++){free(ply->comments[i]);}
			free(ply->comments);
			for(i=0;i<ply->num_obj_info;i++){free(ply->obj_info[i]);}
			free(ply->obj_info);
			ply_free_other_elements (ply->other_elems);
			
			
			for(i=0;i<nr_elems;i++){free(elist[i]);}
			free(elist);
			ply_close(ply);
			return 0;
		}
		
		if (equal_strings("vertex", elem_name)) {
			ply_get_property (ply, elem_name, &oriented_vert_props[0]);
			ply_get_property (ply, elem_name, &oriented_vert_props[1]);
			ply_get_property (ply, elem_name, &oriented_vert_props[2]);
			ply_get_property (ply, elem_name, &oriented_vert_props[3]);
			ply_get_property (ply, elem_name, &oriented_vert_props[4]);
			ply_get_property (ply, elem_name, &oriented_vert_props[5]);
			for (j=0; j < num_elems; j++) {
				ply_get_element (ply, (void *) &ply_vertex);
				OrientedVertex v;
				v.v[0]=ply_vertex.x;
				v.v[1]=ply_vertex.y;
				v.v[2]=ply_vertex.z;
				v.n[0]=ply_vertex.nx;
				v.n[1]=ply_vertex.ny;
				v.n[2]=ply_vertex.nz;
				points.push_back(v);
			}  // for, read vertices
		}  // if vertex
		else{ply_get_other_element (ply, elem_name, num_elems);}

		for(j=0;j<nr_props;j++){
			free(plist[j]->name);
			free(plist[j]);
		}
		free(plist);
	}  // for each type of element
	
	for(i=0;i<nr_elems;i++){
		free(ply->elems[i]->name);
		free(ply->elems[i]->store_prop);
		for(j=0;j<ply->elems[i]->nprops;j++){
			free(ply->elems[i]->props[j]->name);
			free(ply->elems[i]->props[j]);
		}
		if(ply->elems[i]->props && ply->elems[i]->nprops){free(ply->elems[i]->props);}
	}
	for(i=0;i<nr_elems;i++){free(ply->elems[i]);}
	free(ply->elems);
	for(i=0;i<ply->num_comments;i++){free(ply->comments[i]);}
	free(ply->comments);
	for(i=0;i<ply->num_obj_info;i++){free(ply->obj_info[i]);}
	free(ply->obj_info);
	ply_free_other_elements (ply->other_elems);
	
	
	for(i=0;i<nr_elems;i++){free(elist[i]);}
	free(elist);
	ply_close(ply);
	return 1;
}

int PlyWriteTriangles(char* fileName,std::vector<TriangleIndex>& triangles,std::vector<Vertex>& vertices,std::vector<Vertex>& normals,int file_type){
	int i;
	int nr_vertices = vertices.size();
	int nr_faces = triangles.size();
	float version;
	PlyFile *ply = ply_open_for_writing(fileName, 2, elem_names, file_type, &version);
	if (!ply){return 0;}
	
	//
	// describe vertex and face properties
	//
	ply_element_count(ply, "vertex", nr_vertices);
	ply_describe_property(ply, "vertex", &oriented_vert_props[0]);
	ply_describe_property(ply, "vertex", &oriented_vert_props[1]);
	ply_describe_property(ply, "vertex", &oriented_vert_props[2]);
	ply_describe_property(ply, "vertex", &oriented_vert_props[3]);
	ply_describe_property(ply, "vertex", &oriented_vert_props[4]);
	ply_describe_property(ply, "vertex", &oriented_vert_props[5]);
	
	ply_element_count(ply, "face", nr_faces);
	ply_describe_property(ply, "face", &face_props[0]);
	
	ply_header_complete(ply);
	
	// write vertices
	ply_put_element_setup(ply, "vertex");
	for (i=0; i < nr_vertices; i++) {
		PlyOrientedVertex ply_vertex;
		ply_vertex.x  = vertices[i].v[0];
		ply_vertex.y  = vertices[i].v[1];
		ply_vertex.z  = vertices[i].v[2];
		ply_vertex.nx = normals[i].v[0];
		ply_vertex.ny = normals[i].v[1];
		ply_vertex.nz = normals[i].v[2];
		ply_put_element(ply, (void *) &ply_vertex);
		
	}  // for, write vertices
	
	// write faces
	ply_put_element_setup(ply, "face");
	for (i=0; i < nr_faces; i++) {
		//
		// create and fill a struct that the ply code can handle
		//
		PlyFace ply_face;
		ply_face.nr_vertices = 3;
		ply_face.vertices = new int[3];
		for(int j=0; j < 3; j++){ply_face.vertices[j] = triangles[i].idx[j];}
		ply_put_element(ply, (void *) &ply_face);
		delete[] ply_face.vertices;
		
	}  // for, write faces
	
	ply_close(ply);
	return 1;
}

int PlyWriteTriangles(char* fileName,vector<TriangleIndex>& triangles,vector<Vertex>& vertices,int file_type){
	int i;
	int nr_vertices = vertices.size();
	int nr_faces = triangles.size();
	
	float version;
	PlyFile *ply = ply_open_for_writing(fileName, 2, elem_names, file_type, &version);
	if (!ply){return 0;}
	
	//
	// describe vertex and face properties
	//
	ply_element_count(ply, "vertex", nr_vertices);
	ply_describe_property(ply, "vertex", &vert_props[0]);
	ply_describe_property(ply, "vertex", &vert_props[1]);
	ply_describe_property(ply, "vertex", &vert_props[2]);
	
	ply_element_count(ply, "face", nr_faces);
	ply_describe_property(ply, "face", &face_props[0]);
	
	ply_header_complete(ply);
	
	// write vertices
	ply_put_element_setup(ply, "vertex");
	for (i=0; i < nr_vertices; i++) {
		PlyVertex ply_vertex;
		ply_vertex.x = vertices[i].v[0];
		ply_vertex.y = vertices[i].v[1];
		ply_vertex.z = vertices[i].v[2];
		ply_put_element(ply, (void *) &ply_vertex);
		
	}  // for, write vertices
	
	// write faces
	ply_put_element_setup(ply, "face");
	for (i=0; i < nr_faces; i++) {
		//
		// create and fill a struct that the ply code can handle
		//
		PlyFace ply_face;
		ply_face.nr_vertices = 3;
		ply_face.vertices = new int[3];
		for(int j=0; j < 3; j++){ply_face.vertices[j] = triangles[i].idx[j];}
		ply_put_element(ply, (void *) &ply_face);
		delete[] ply_face.vertices;
		
	}  // for, write faces
	
	ply_close(ply);
	return 1;
}
