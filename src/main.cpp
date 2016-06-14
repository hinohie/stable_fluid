#include<stdlib.h>
#include<GL/glew.h>
#include<glut.h>

#include<stdio.h>
#include<algorithm>
#include<vector>

#include<time.h>
#include<Windows.h>
/*
#include<FreeImage.h>
*/
#define FRAME_TIME 32
using namespace std;

#include "stablefluid.h"

static int width;
static int height;

double wide = 50;
vector3f eye,look;
StableFluid *f;
static int resolution = 50;
static long t1,t2;

double mZoom;
double mRotateX, mRotateY, mRotateZ, mAngle;
double mTransX, mTransY, mTransZ, mTransSize;
double mScaleX, mScaleY, mScaleZ;
void reshape(int w,int h)
{
        glViewport(0, 0, w, h);
		width = w;
		height = h;
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
//		glOrtho(-wide,wide,-wide,wide, 0.5,1000);
		gluPerspective(60, 1.0, 0.1, resolution * 4);
}
/*
static int frame=0;
void capture()
{
	if(frame > 500) exit(0);
	// Make the BYTE array, factor of 3 because it's RBG.
	BYTE* pixels = new BYTE[ 3 * width * height];

	glReadPixels(0, 0, width, height, GL_BGR, GL_UNSIGNED_BYTE, pixels);

	// Convert to FreeImage format & save to file
	FIBITMAP* image = FreeImage_ConvertFromRawBits(pixels, width, height, 3 * width, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);
	char filename[200];
	sprintf(filename,"F:/Result_Seen/Stable_fluid/Image/%05d.png",frame);
	frame++;
	FreeImage_Save(FIF_PNG, image, filename, 0);

	// Free resources
	FreeImage_Unload(image);
	delete [] pixels;
}
*/
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(eye.x, eye.y, eye.z,
		look.x, look.y, look.z,
		0, 1, 0);
	
	glTranslatef(look.x, look.y, look.z);
	glTranslatef(0, 0, -mZoom);
	glTranslatef(mTransX, mTransY, mTransZ);
	glRotatef(mRotateX, 1.0, 0.0, 0.0);
	glRotatef(mRotateY, 0.0, 1.0, 0.0);
	glRotatef(mRotateZ, 0.0, 0.0, 1.0);	
	glTranslatef(-look.x, -look.y, -look.z);


	int i, j;

	f->pipeline();
	f->draw();
//	capture();

	glutSwapBuffers();

	t2 = clock();
	if(FRAME_TIME > (t2 - t1))
	{
		Sleep(FRAME_TIME - (t2 - t1));
	}
	printf("%ld : ms/frame\n",t2-t1);
	t1 = clock();
	glutPostRedisplay();
}
void camera_update(vector3f e,vector3f l)
{
//	printf("%lf %lf %lf  %lf %lf %lf\n",e.x,e.y,e.z, l.x,l.y,l.z);
	eye = e;
	look = l;
	glutPostRedisplay();
}
void keyboard(unsigned char key, int x,int y)
{
	
	switch(key)
	{
	case 27:
	case 'q':
	case 'Q':
		exit(0);
		break;
	case 'r':
	case 'R':
		f->reset();
		break;
	case 'v':
	case 'V':
		f->show_velo ^= true;
		break;
		/*
	case 'w':
	case 'W':
		camera_update(vector3f(eye.x, eye.y+10 , eye.z), vector3f(look.x, look.y+10, look.z));
		break;
	case 's':
	case 'S':
		camera_update(vector3f(eye.x, eye.y-10 , eye.z), vector3f(look.x, look.y-10, look.z));
		break;
	case 'a':
	case 'A':
		camera_update(vector3f(eye.x+10, eye.y , eye.z), vector3f(look.x+10, look.y, look.z));
		break;
	case 'd':
	case 'D':
		camera_update(vector3f(eye.x-10, eye.y , eye.z), vector3f(look.x-10, look.y, look.z));
		break;
	case '+':
		camera_update(vector3f(eye.x, eye.y , eye.z+0.5), vector3f(look.x, look.y, look.z+0.5));
		break;
	case '-':
		camera_update(vector3f(eye.x, eye.y , eye.z-0.5), vector3f(look.x, look.y, look.z-0.5));
		break;
		*/
	}
}

static bool is_click;
static bool is_crick;
static bool is_cwick;
vector2f cur_mouse_pos;
void Cmouse(int button, int state, int x,int y)
{
	double cx, cy;
	int i, j;

	cx = x / (double)width;
	cy = y / (double)height;
	if(state == GLUT_DOWN)
	{
		if(button == GLUT_LEFT_BUTTON)
		{
			is_click = true;
			cur_mouse_pos = vector2f(cx,cy);
		}
		if(button == GLUT_RIGHT_BUTTON)
		{
			is_crick = true;
			cur_mouse_pos = vector2f(cx,cy);
		}
		if(button == GLUT_MIDDLE_BUTTON)
		{
			is_cwick = true;
			cur_mouse_pos = vector2f(cx,cy);
		}
	}
	else
	{
		if(button == GLUT_LEFT_BUTTON)
		{
			is_click = false;
		}
		if(button == GLUT_RIGHT_BUTTON)
		{
			is_crick = false;
		}
		if(button == GLUT_MIDDLE_BUTTON)
		{
			is_cwick = false;
		}
	}
}
void rotate_eye(vector3f axis, double f)
{
	int i, j, k;
	vector3f diff = eye - look;
	double hei = diff * axis;
	vector3f dy = hei * axis;
	vector3f dx = diff - dy;

	vector3f dz = vector3f(axis.y*dx.z - axis.z*dx.y, axis.z*dx.x - axis.x*dx.z, axis.x*dx.y - axis.y*dx.x);

	vector3f res = cos(f) * dx + sin(f) * dz;
	res = res + dy;
	camera_update(look + res,look);
}
void Mmouse(int x,int y)
{
	double cx = x / (double)width;
	double cy = y / (double)height;

	vector2f new_mouse_pos = vector2f(cx,cy);
	vector2f diff = new_mouse_pos - cur_mouse_pos;
	vector3f eye_diff = eye - look;
	if(is_cwick)
	{
		mZoom -= 0.5 * resolution * diff.y;
	}
	if(is_click)
	{
		mRotateX += 1.0 * resolution * diff.y;
		mRotateY += 1.0 * resolution * diff.x;	
	}
	if(is_crick)
	{
		camera_update(eye - vector3f(diff.x,0,diff.y)*0.5*resolution, look - vector3f(diff.x,0,diff.y)*0.5*resolution);
	}
	cur_mouse_pos = new_mouse_pos;
}


void data_init()
{
	resolution = 25;
	mZoom = 5.0f;
//	mRotateX = 0, mRotateY = 0, mRotateZ = 0, mAngle = 45.0;
	mRotateX = 14.20833, mRotateY = -18.8750, mRotateZ = 0, mAngle = 45.0;
	mTransX = 0, mTransY = 0, mTransZ = 0, mTransSize = 0.5;
	mScaleX = 1.0, mScaleY = 1.0, mScaleZ = 1.0;

	eye = vector3f(resolution/2.0, resolution/2.0, resolution + resolution );
	look = vector3f(resolution/2.0, resolution/2.0, resolution/2.0);
	
//	f = new StableFluid(resolution, resolution);				// 2차원  15ms/frame
	f = new StableFluid(resolution, resolution, resolution);	// 3차원 300ms/frame
//	f->mian();
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(600, 600);
	glutCreateWindow("Stable");
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glutDisplayFunc( display );
	glutReshapeFunc( reshape );
	glutKeyboardFunc( keyboard );
	glutMouseFunc( Cmouse );
	glutMotionFunc( Mmouse );
	data_init();

	glutMainLoop();
	return 0;
}