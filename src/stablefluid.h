#ifndef __STABLE_FLUID_H__
#define __STABLE_FLUID_H__

#include "vector3f.h"
class StableFluid
{
public:
	double ***a;
	double ***u;
	double ***v;
	double ***z;
	double ***tu;
	double ***tv;
	double ***tz;
	double ***ta;

	vector3f ***f;
	int width, height, depth;
	double dx, dy, dz;
	double visc;
	double dt;
	double diff;

	bool show_velo;

	StableFluid(int hh, int ww,int dd);
	StableFluid(int hh, int ww);
	~StableFluid();
	void mian();
	void reset();
	void pipeline();
	void draw();
	void draw_box(double sx, double sy, double sz, double ex,double ey, double ez,double alpha);
	void draw_velo(double sx, double sy, double sz, double vx,double vy, double vz);
	void draw_fog(int h, int w, int d, double ***a,double dx,double dy,double dz);
};

#endif