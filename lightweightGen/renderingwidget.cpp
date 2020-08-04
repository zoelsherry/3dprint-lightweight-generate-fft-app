#include "renderingwidget.h"
#include <QKeyEvent>
#include <QColorDialog>
#include <QFileDialog>
#include <iostream>
#include <QtWidgets/QMenu>
#include <QtWidgets/QAction>
#include <QTextCodec>
#include <qelapsedtimer.h>
#include <qdebug.h>
#include <gl/GLU.h>
#pragma comment(lib, "GlU32.Lib")       /* link with SGI OpenGL Utility lib */
//#include <GL/glut.h> // Modified

#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <sys/timeb.h>
#include <fstream>

#include "gl_helper.h"
#include <algorithm>
#include "mainwindow.h"
#include "ArcBall.h"
#include "globalFunctions.h"
#include "HE_mesh/Mesh3D.h"
#include "setUp3D.h"
#include "spline.h"


RenderingWidget::RenderingWidget(QWidget *parent, MainWindow* mainwindow)
: QGLWidget(parent), ptr_mainwindow_(mainwindow), eye_distance_(5.0),
has_lighting_(false), is_draw_point_(true), is_draw_edge_(false), is_draw_face_(false), is_draw_texture_(false)
{
	ptr_arcball_ = new CArcBall(width(), height());
	ptr_mesh_ = new Mesh3D();

	is_load_texture_ = false;
	is_draw_axes_ = false;

	eye_goal_[0] = eye_goal_[1] = eye_goal_[2] = 0.0;
	eye_direction_[0] = eye_direction_[1] = 0.0;
	eye_direction_[2] = 1.0;
}

RenderingWidget::~RenderingWidget()
{
	SafeDelete(ptr_arcball_);
	SafeDelete(ptr_mesh_);
}

void RenderingWidget::initializeGL()
{
	glClearColor(0.3, 0.3, 0.3, 0.0);
	glShadeModel(GL_SMOOTH);

	glEnable(GL_DOUBLEBUFFER);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1);

	SetLight();

}

void RenderingWidget::resizeGL(int w, int h)
{
	h = (h == 0) ? 1 : h;

	ptr_arcball_->reSetBound(w, h);

	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(45.0, GLdouble(w) / GLdouble(h), 0.001, 1000);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void RenderingWidget::paintGL()
{
	glShadeModel(GL_SMOOTH);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (has_lighting_)
	{
		SetLight();
	}
	else
	{
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
	}

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	register vec eyepos = eye_distance_*eye_direction_;
	gluLookAt(eyepos[0], eyepos[1], eyepos[2],
		eye_goal_[0], eye_goal_[1], eye_goal_[2],
		0.0, 1.0, 0.0);
	glPushMatrix();

	glMultMatrixf(ptr_arcball_->GetBallMatrix());

	Render();
	glPopMatrix();
}

void RenderingWidget::timerEvent(QTimerEvent * e)
{
	updateGL();
}

void RenderingWidget::mousePressEvent(QMouseEvent *e)
{
	switch (e->button())
	{
	case Qt::LeftButton:
		ptr_arcball_->MouseDown(e->pos());
		break;
	case Qt::MidButton:
		current_position_ = e->pos();
		break;
	default:
		break;
	}

	updateGL();
}
void RenderingWidget::mouseMoveEvent(QMouseEvent *e)
{
	switch (e->buttons())
	{
		setCursor(Qt::ClosedHandCursor);
	case Qt::LeftButton:
		ptr_arcball_->MouseMove(e->pos());
		break;
	case Qt::MidButton:
		eye_goal_[0] -= 4.0*GLfloat(e->x() - current_position_.x()) / GLfloat(width());
		eye_goal_[1] += 4.0*GLfloat(e->y() - current_position_.y()) / GLfloat(height());
		current_position_ = e->pos();
		break;
	default:
		break;
	}

	updateGL();
}
void RenderingWidget::mouseDoubleClickEvent(QMouseEvent *e)
{
	switch (e->button())
	{
	case Qt::LeftButton:
		break;
	default:
		break;
	}
	updateGL();
}
void RenderingWidget::mouseReleaseEvent(QMouseEvent *e)
{
	switch (e->button())
	{
	case Qt::LeftButton:
		ptr_arcball_->MouseUp(e->pos());
		setCursor(Qt::ArrowCursor);
		break;

	case Qt::RightButton:
		break;
	default:
		break;
	}
}

void RenderingWidget::wheelEvent(QWheelEvent *e)
{
	eye_distance_ += e->delta()*0.001;
	eye_distance_ = eye_distance_ < 0 ? 0 : eye_distance_;

	updateGL();
}

void RenderingWidget::keyPressEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_A:
		break;
	default:
		break;
	}
}

void RenderingWidget::keyReleaseEvent(QKeyEvent *e)
{
	switch (e->key())
	{
	case Qt::Key_A:
		break;
	default:
		break;
	}
}

void RenderingWidget::Render()
{
	DrawAxes(is_draw_axes_);

	DrawPoints(is_draw_point_);
	DrawEdge(is_draw_edge_);
	DrawFace(is_draw_face_);
	DrawTexture(is_draw_texture_);
}

void RenderingWidget::SetLight()
{
	static GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	static GLfloat mat_shininess[] = { 50.0 };
	static GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
	static GLfloat white_light[] = { 0.8, 0.8, 0.8, 1.0 };
	static GLfloat lmodel_ambient[] = { 0.3, 0.3, 0.3, 1.0 };

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
	glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
}

void RenderingWidget::SetBackground()
{
	QColor color = QColorDialog::getColor(Qt::white, this, tr("background color"));
	GLfloat r = (color.red()) / 255.0f;
	GLfloat g = (color.green()) / 255.0f;
	GLfloat b = (color.blue()) / 255.0f;
	GLfloat alpha = color.alpha() / 255.0f;
	glClearColor(r, g, b, alpha);
	updateGL();
}


void RenderingWidget::New()
{
	SafeDelete(ptr_mesh_);
	ptr_mesh_ = new Mesh3D();

	is_load_texture_ = false;
	is_draw_axes_ = false;

	eye_goal_[0] = eye_goal_[1] = eye_goal_[2] = 0.0;
	eye_direction_[0] = eye_direction_[1] = 0.0;
	eye_direction_[2] = 1.0;

	updateGL();
}

void RenderingWidget::ReadMesh()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Read Mesh"), "..", tr("Meshes (*.obj)"));

	if (filename.isEmpty())
	{
		emit(operatorInfo(QString("Read Mesh Failed!")));
		return;
	}
	//support chinese path
	QTextCodec *code = QTextCodec::codecForName("gd18030");
	QTextCodec::setCodecForLocale(code);

	QByteArray byfilename = filename.toLocal8Bit();
	//ptr_mesh_->LoadFromOBJFile(byfilename.data()); 
	ptr_mesh_->LoadFromOBJFile(byfilename.data());    // TODO: some matter


	//	m_pMesh->LoadFromOBJFile(filename.toLatin1().data());
	emit(operatorInfo(QString("Read Mesh from") + filename + QString(" Done")));
	emit(meshInfo(ptr_mesh_->num_of_vertex_list(), ptr_mesh_->num_of_edge_list(), ptr_mesh_->num_of_face_list()));
	updateGL();
}

void RenderingWidget::WriteMesh()
{
	if (ptr_mesh_->num_of_vertex_list() == 0)
	{
		emit(QString("The Mesh is Empty !"));
		return;
	}
	QString filename = QFileDialog::
		getSaveFileName(this, tr("Write Mesh"),
		"..", tr("Meshes (*.obj)"));

	if (filename.isEmpty())
		return;

	ptr_mesh_->WriteToOBJFile(filename.toLatin1().data());

	emit(operatorInfo(QString("Write Mesh to ") + filename + QString(" Done")));
}

void RenderingWidget::LoadTexture()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Load Texture"),
		"..", tr("Images(*.bmp *.jpg *.png *.jpeg)"));
	if (filename.isEmpty())
	{
		emit(operatorInfo(QString("Load Texture Failed!")));
		return;
	}


	glGenTextures(1, &texture_[0]);
	QImage tex1, buf;
	if (!buf.load(filename))
	{
		//        QMessageBox::warning(this, tr("Load Fialed!"), tr("Cannot Load Image %1").arg(filenames.at(0)));
		emit(operatorInfo(QString("Load Texture Failed!")));
		return;
		/*
		QImage dummy(128, 128, QImage::Format_ARGB32);
		dummy.fill(Qt::green);
		buf = dummy;
		*/
	}

	tex1 = QGLWidget::convertToGLFormat(buf);
	glBindTexture(GL_TEXTURE_2D, texture_[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_NEAREST);
	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, tex1.width(), tex1.height(),
		GL_RGBA, GL_UNSIGNED_BYTE, tex1.bits());

	is_load_texture_ = true;
	emit(operatorInfo(QString("Load Texture from ") + filename + QString(" Done")));
}

void RenderingWidget::FillStructure()
{
	std::vector<int> supportedBW;
	float translate[3], outScale;
	VoxelData v, cubicBSplineValues;
	double tt, avg;
	int i, bw, fftwType, order;
	AdjacencyTable table;
	vector<OrientedVertex> pts;
	vector<TriangleIndex> triangles;
	vector<Vertex> vertices;
	vector<Vertex> normals;
	float inScale, smooth = -1;
	float volume;
	int* idx;
	float* dx;

	int Res = 256;

	QElapsedTimer timer;
	qDebug("lightweight structure fabrication start...");
	bw = Res / 2;

	fftwType = FFTW_ESTIMATE;

	timer.start();
	const std::vector<HE_vert*>& verts = *(ptr_mesh_->get_vertex_list());
	
	OrientedVertex ver;
	for (size_t i = 0; i != ptr_mesh_->num_of_vertex_list(); ++i)
	{
		ver.v[0] = verts[i]->position().data()[0];
		ver.v[1] = verts[i]->position().data()[1];
		ver.v[2] = verts[i]->position().data()[2];

		ver.n[0] = verts[i]->normal().data()[0];
		ver.n[1] = verts[i]->normal().data()[1];
		ver.n[2] = verts[i]->normal().data()[2];

		pts.push_back(ver);
	}
	qDebug() << "Convert to half-edge Time: " << timer.elapsed() / 1000.0;
	

	order = Spline::QUADRATIC;

	smooth = 0;
	inScale = 1;
	volume = 1; 
	
	int cnt = 0;
	for (i = 0; i < pts.size(); i++) {
		double l = (float)Length(pts[i].n);
		if (l<10e-8 || l>10e8) {
			pts[i] = pts.back();
			pts.pop_back();
		}
		else {
			pts[i].n[0] /= l;
			pts[i].n[1] /= l;
			pts[i].n[2] /= l;
			cnt++;
		}
	}
	for (i = 0; i < pts.size(); i++) {
		pts[i].n[0] /= cnt;
		pts[i].n[1] /= cnt;
		pts[i].n[2] /= cnt;
	}
	
	timer.start();
	SetUpPoints(pts, bw, inScale, translate, outScale, volume);
	if (!SetPositions(2 * bw, 1, pts, idx, dx)) {
		printf("Failed to set positions\n");
		return;
	}

	tt = timer.elapsed() / 1000.0;
	qDebug() << "Set up time: " << timer.elapsed()/1000.0;

	timer.start();

	if (!GetInterior(pts, idx, dx, bw, v, fftwType, smooth)) {
		printf("Failed to compute interior\n");
		return;
	}
	

	if (!table.set(Res, 0, 0)) {
		printf("Could not allocated memory for adjacency table\n");
		return;
	}
	

	avg = AverageValue(pts, v, table, order);
	delete[] idx;
	delete[] dx;


	int r;
	r = 8 * bw*bw*bw;
	for (i = 0; i < r; i++) { v.array[i] -= (float)avg; }

	tt += timer.elapsed() / 1000.0;
	qDebug() << "Iso-Function Time: " << timer.elapsed() / 1000.0;

	////////// fill lightweight structure ///////////
	// TODO: read the voxel data from v.arrary

	timer.start();
	char* cell = "various_surface/P50.m";
	v.fillUniformStructure_interior(cell, 0);
	tt += timer.elapsed() / 1000.0;
	qDebug() << "Fill Structure Time: " << timer.elapsed() / 1000.0;

	timer.start();
	v.isoSurface(triangles, vertices, 0, table, order, fftwType);
	tt += timer.elapsed() / 1000.0;
	qDebug() << "Iso-Extraction Time: " << timer.elapsed() / 1000.0;

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

	qDebug() << "Total Time: " << tt;


	//////////////////////////convert to half-edge
	timer.start();
	SafeDelete(ptr_mesh_);
	ptr_mesh_ = new Mesh3D();

	for (int i = 0; i < vertices.size(); i++)
	{
		Vec3f nvv;
		nvv[0] = vertices[i].v[0];
		nvv[1] = vertices[i].v[1];
		nvv[2] = vertices[i].v[2];
		ptr_mesh_->InsertVertex(nvv);
	}
	
	for (i = 0; i < triangles.size(); i++) 
	{
		std::vector<HE_vert* > s_faceid;
		for (int j = 0; j < 3; j++)
		{
			HE_vert* hv = ptr_mesh_->get_vertex(triangles[i].idx[j]);
			if (hv != NULL)
			{
				s_faceid.push_back(hv);
			}
		}
		ptr_mesh_->InsertFace(s_faceid);
	}
	qDebug() << "Convert to half-edge Time: " << timer.elapsed() / 1000.0;

	ptr_mesh_->UpdateMesh();

	emit(operatorInfo(QString("Fill structure... Done")));
	emit(meshInfo(ptr_mesh_->num_of_vertex_list(), ptr_mesh_->num_of_edge_list(), ptr_mesh_->num_of_face_list()));
	updateGL();
	//Unify(2.f);

}

void RenderingWidget::CheckDrawPoint(bool bv)
{
	is_draw_point_ = bv;
	updateGL();
}
void RenderingWidget::CheckDrawEdge(bool bv)
{
	is_draw_edge_ = bv;
	updateGL();
}
void RenderingWidget::CheckDrawFace(bool bv)
{
	is_draw_face_ = bv;
	updateGL();
}
void RenderingWidget::CheckLight(bool bv)
{
	has_lighting_ = bv;
	updateGL();
}
void RenderingWidget::CheckDrawTexture(bool bv)
{
	is_draw_texture_ = bv;
	if (is_draw_texture_)
		glEnable(GL_TEXTURE_2D);
	else
		glDisable(GL_TEXTURE_2D);

	updateGL();
}
void RenderingWidget::CheckDrawAxes(bool bV)
{
	is_draw_axes_ = bV;
	updateGL();
}

void RenderingWidget::DrawAxes(bool bV)
{
	if (!bV)
		return;
	//x axis
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0.7, 0.0, 0.0);
	glEnd();
	glPushMatrix();
	glTranslatef(0.7, 0, 0);
	glRotatef(90, 0.0, 1.0, 0.0);
	glutSolidCone(0.02, 0.06, 20, 10);
	glPopMatrix();

	//y axis
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0.0, 0.7, 0.0);
	glEnd();

	glPushMatrix();
	glTranslatef(0.0, 0.7, 0);
	glRotatef(90, -1.0, 0.0, 0.0);
	glutSolidCone(0.02, 0.06, 20, 10);
	glPopMatrix();

	//z axis
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0.0, 0.0, 0.7);
	glEnd();
	glPushMatrix();
	glTranslatef(0.0, 0, 0.7);
	glutSolidCone(0.02, 0.06, 20, 10);
	glPopMatrix();

	glColor3f(1.0, 1.0, 1.0);
}

void RenderingWidget::DrawPoints(bool bv)
{
	if (!bv || ptr_mesh_ == NULL)
		return;
	if (ptr_mesh_->num_of_vertex_list() == 0)
	{
		return;
	}

	const std::vector<HE_vert*>& verts = *(ptr_mesh_->get_vertex_list());
	glBegin(GL_POINTS);
	for (size_t i = 0; i != ptr_mesh_->num_of_vertex_list(); ++i)
	{
		glNormal3fv(verts[i]->normal().data());
		glVertex3fv(verts[i]->position().data());
	}
	glEnd();
}

void RenderingWidget::DrawEdge(bool bv)
{
	if (!bv || ptr_mesh_ == NULL)
		return;

	if (ptr_mesh_->num_of_face_list() == 0)
	{
		return;
	}

	const std::vector<HE_face *>& faces = *(ptr_mesh_->get_faces_list());
	for (size_t i = 0; i != faces.size(); ++i)
	{
		glBegin(GL_LINE_LOOP);
		HE_edge *pedge(faces.at(i)->pedge_);
		do
		{
			glNormal3fv(pedge->pvert_->normal().data());
			glVertex3fv(pedge->pvert_->position().data());

			pedge = pedge->pnext_;

		} while (pedge != faces.at(i)->pedge_);
		glEnd();
	}
}

void RenderingWidget::DrawFace(bool bv)
{
	if (!bv || ptr_mesh_ == NULL)
		return;

	if (ptr_mesh_->num_of_face_list() == 0)
	{
		return;
	}

	const std::vector<HE_face *>& faces = *(ptr_mesh_->get_faces_list());

	glBegin(GL_TRIANGLES);

	for (size_t i = 0; i != faces.size(); ++i)
	{
		HE_edge *pedge(faces.at(i)->pedge_);
		do
		{
			glNormal3fv(pedge->pvert_->normal().data());
			glVertex3fv(pedge->pvert_->position().data());

			pedge = pedge->pnext_;

		} while (pedge != faces.at(i)->pedge_);
	}

	glEnd();
}

void RenderingWidget::DrawTexture(bool bv)
{
	if (!bv)
		return;
	if (ptr_mesh_->num_of_face_list() == 0 || !is_load_texture_)
		return;

	//默认使用球面纹理映射，效果不好
	ptr_mesh_->SphereTex();

	const std::vector<HE_face *>& faces = *(ptr_mesh_->get_faces_list());

	glBindTexture(GL_TEXTURE_2D, texture_[0]);
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i != faces.size(); ++i)
	{
		HE_edge *pedge(faces.at(i)->pedge_);
		do
		{
			/*请在此处绘制纹理，添加纹理坐标即可*/
			glTexCoord2fv(pedge->pvert_->texCoord_.data());
			glNormal3fv(pedge->pvert_->normal().data());
			glVertex3fv(pedge->pvert_->position().data());

			pedge = pedge->pnext_;

		} while (pedge != faces.at(i)->pedge_);
	}

	glEnd();
}
