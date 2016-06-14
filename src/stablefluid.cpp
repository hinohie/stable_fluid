#include"stablefluid.h"
#include"vector3f.h"
#include<glut.h>
#include<omp.h>
bool open_dookung = false;
void set_bnd(int h,int w,int d,int b, double ***v)
{
	int i, j, k, l;
	
	#pragma omp parallel for
	for(l=0; l<h*d;l++)
	{
		int i, k;
		i = l/d+1;
		k = l%d+1;
		v[i][0][k] = b==2 ? v[i][1][k] * (-1) : v[i][1][k];
		v[i][w+1][k] = b==2 ? v[i][w][k] * (-1) : v[i][w][k];
	}
	#pragma omp parallel for
	for(l=0; l<w*d;l++)
	{
		int i=l/d+1;
		int k=l%d+1;
		v[0][i][k] = b==1 ? v[1][i][k] * (-1) : v[1][i][k];
		v[h+1][i][k] = b==1 ? v[h][i][k] * (-1) : v[h][i][k];
		if(open_dookung) v[h+1][i][k] = v[h][i][k] = 0;
	}
	#pragma omp parallel for
	for(l=0; l<w*h;l++)
	{
		int i, j;
		i = l/w+1;
		j = l%w+1;
		v[i][j][0] = b==3? v[i][j][1] * (-1) : v[i][j][1];
		v[i][j][d+1] = b==3 ? v[i][j][d] * (-1) : v[i][j][d];
	}
	
	#pragma omp parallel for
	for(i=1;i<=h;i++)
	{
		v[i][0][0] = (v[i][0][1] + v[i][1][0])/2;
		v[i][w+1][0] = (v[i][w+1][1] + v[i][w][0])/2;
		v[i][w+1][d+1] = (v[i][w+1][d] + v[i][w][d+1])/2;
		v[i][0][d+1] = (v[i][0][d] + v[i][1][d+1])/2;
	}
	#pragma omp parallel for
	for(i=1;i<=w;i++)
	{
		v[0][i][0] = (v[0][i][1] + v[1][i][0])/2;
		v[h+1][i][0] = (v[w+1][i][1] + v[w][i][0])/2;
		v[h+1][i][d+1] = (v[w+1][i][d] + v[w][i][d+1])/2;
		v[0][i][d+1] = (v[0][i][d] + v[1][i][d+1])/2;
	}
	#pragma omp parallel for
	for(i=1;i<=d;i++)
	{
		v[0][0][i] = (v[0][1][i] + v[1][0][i])/2;
		v[h+1][0][i] = (v[h+1][1][i] + v[h][0][i])/2;
		v[h+1][w+1][i] = (v[h+1][w][i] + v[h][w+1][i])/2;
		v[0][w+1][i] = (v[0][w][i] + v[1][w+1][i])/2;
	}
	v[0][0][0] = (v[0][0][1] + v[0][1][0] + v[1][0][0])/3;
	v[0][0][d+1] = (v[0][0][d] + v[0][1][d+1] + v[1][0][d+1])/3;
	v[0][w+1][0] = (v[0][w+1][1] + v[0][w][0] + v[1][w+1][0])/3;
	v[0][w+1][d+1] = (v[0][w+1][d] + v[0][w][d+1] + v[1][w+1][d+1])/3;
	v[h+1][0][0] = (v[h+1][0][1] + v[h+1][1][0] + v[h][0][0])/3;
	v[h+1][0][d+1] = (v[h+1][0][d] + v[h+1][1][d+1] + v[h][0][d+1])/3;
	v[h+1][w+1][0] = (v[h+1][w+1][1] + v[h+1][w][0] + v[h][w+1][0])/3;
	v[h+1][w+1][d+1] = (v[h+1][w+1][d] + v[h+1][w][d+1] + v[h][w+1][d+1])/3;
	if(open_dookung) v[h+1][0][0] = 0;
	if(open_dookung) v[h+1][0][d+1] = 0;
	if(open_dookung) v[h+1][w+1][0] = 0;
	if(open_dookung) v[h+1][w+1][d+1] = 0;
}
#include<iostream>
void diffuse(int h,int w,int d,int b, double ***tv, double ***v, double diff, double dt)
{
	int i, j, k, iter;
	double da = dt * diff * h*w;
	
	/*
	int aaa = (h+2)*(w+2)*(d+2);
	SparseMatrix<double> A(aaa,aaa);
	BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double> > solver;
	VectorXd bb(aaa),x(aaa);
	vector<Triplet<double> > Tr;
	bb.setZero();
	x.setZero();
	for(i=0;i<=h+1;i++)
		for(j=0;j<=w+1;j++)
			for(k=0;k<=d+1;k++)
			{
				int pid,qid;
				int spid,si,sj,sk;
				pid = (i + 0)*(w+2)*(d+2) + (j + 0)*(d+2) + (k + 0);
				si=i;sj=j;sk=k;
				if(si<=0)si=1;if(si>h)si=h;
				if(sj<=0)sj=1;if(sj>w)sj=w;
				if(sk<=0)sk=1;if(sk>d)sk=d;
				spid = (si + 0)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk + 0);
				{
					Tr.push_back(Triplet<double>(pid,spid,1+6*da));
					qid = (si + 1)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk + 0);
					Tr.push_back(Triplet<double>(pid,qid,-da));
					qid = (si - 1)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk + 0);
					Tr.push_back(Triplet<double>(pid,qid,-da));
					qid = (si - 0)*(w+2)*(d+2) + (sj + 1)*(d+2) + (sk + 0);
					Tr.push_back(Triplet<double>(pid,qid,-da));
					qid = (si - 0)*(w+2)*(d+2) + (sj - 1)*(d+2) + (sk + 0);
					Tr.push_back(Triplet<double>(pid,qid,-da));
					qid = (si - 0)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk + 1);
					Tr.push_back(Triplet<double>(pid,qid,-da));
					qid = (si - 0)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk - 1);
					Tr.push_back(Triplet<double>(pid,qid,-da));
				}
				bb[pid] = v[si][sj][sk];
				x[pid] = v[si][sj][sk];
			}
	A.setFromTriplets(Tr.begin(),Tr.end());
	solver.setTolerance(0.00000001);
	solver.compute(A);
	int cnt=0;
	do{
		x = solver.solveWithGuess(bb,x);
	for(i=1;i<=h;i++)
		for(j=1;j<=w;j++)
			for(k=1;k<=d;k++)
			{
				int pid;
				pid = (i + 0)*(w+2)*(d+2) + (j + 0)*(d+2) + (k + 0);
				tv[i][j][k] = x[pid];
			}
	set_bnd(h, w, d, b, tv);
	cnt++;
	}while(cnt<2);
//	printf("wow %d %d %.10lf\n",solver.iterations(),solver.info(), solver.error());
	*/

	//gause-saidel
	for(iter=0;iter<20;iter++)
	{
		int l;
	#pragma omp parallel for
	for(l=0;l<h*w*d;l++)
		{
			int i,j,k;
			i = l/w/d + 1;
			j = (l/d)%w+1;
			k = l%d+1;
				tv[i][j][k] = (v[i][j][k] + (tv[i-1][j][k] + tv[i][j-1][k] + tv[i+1][j][k] + tv[i][j+1][k] + tv[i][j][k-1] + tv[i][j][k+1])*da)/(1+6*da);
			} 
		set_bnd(h, w, d, b, tv);
	}
}
void advect(int h, int w, int d, int b, double ***ta, double ***a, double ***u, double ***v, double ***z,double dt)
{
	int l;
	double dt0 = dt * 10;
	#pragma omp parallel for
	for(l=0;l<h*w*d;l++)
		{
			int i,j,k;
			i = l/w/d + 1;
			j = (l/d)%w+1;
			k = l%d+1;
	int pi, pj, pk, qi, qj, qk;
	double ti,tj,tk, s0, s1, t0, t1, f0, f1;
			ti = i - dt0 * u[i][j][k];
			tj = j - dt0 * v[i][j][k];
			tk = k - dt0 * z[i][j][k];
			if(ti<0.5)ti = 0.5; if(tj<0.5) tj = 0.5; if(tk<0.5) tk = 0.5;
			if(ti>h+0.5) ti = h+0.5; if(tj>w+0.5) tj=w+0.5; if(tk > d+0.5) tk = d + 0.5;
			pi = (int)ti;
			pj = (int)tj;
			pk = (int)tk;
			qi = pi + 1;
			qj = pj + 1;
			qk = pk + 1;
			s1 = ti - pi; s0 = 1 - s1;
			t1 = tj - pj; t0 = 1 - t1;
			f1 = tk - pk; f0 = 1 - f1;
			ta[i][j][k] = ((a[pi][pj][pk] * t0 + a[pi][qj][pk] * t1) * s0 + (a[qi][pj][pk] * t0 + a[qi][qj][pk] * t1 ) * s1) * f0 +
					((a[pi][pj][qk] * t0 + a[pi][qj][qk] * t1) * s0 + (a[qi][pj][qk] * t0 + a[qi][qj][qk] * t1 ) * s1) * f1;
		}
	set_bnd(h, w, d, b, ta);
}
void project(int h,int w,int d,double ***tu, double ***tv, double ***tz, double ***div, double ***p)
{
	int i, j, k, l;
	int iter;
	double da = 1.0/max(h,max(w,d));
	
	#pragma omp parallel for
	for(l=0;l<h*w*d;l++)
		{
			int i,j,k;
			i = l/w/d + 1;
			j = (l/d)%w+1;
			k = l%d+1;
			div[i][j][k] = -0.5 * (tu[i+1][j][k] - tu[i-1][j][k] + tv[i][j+1][k] - tv[i][j-1][k] + tz[i][j][k+1] - tz[i][j][k-1])*da;
			p[i][j][k] = 0;
		}
	set_bnd(h, w, d, 0, div);
	set_bnd(h, w, d, 0, p);
	/*
	int aaa = (h+2)*(w+2)*(d+2);
	SparseMatrix<double> A(aaa,aaa);
	BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double> > solver;
	VectorXd bb(aaa),x(aaa);
	vector<Triplet<double> > Tr;
	bb.setZero();
	x.setZero();
	for(i=0;i<=h+1;i++)
		for(j=0;j<=w+1;j++)
			for(k=0;k<=d+1;k++)
			{
				int pid,qid;
				int spid,si,sj,sk;
				pid = (i + 0)*(w+2)*(d+2) + (j + 0)*(d+2) + (k + 0);
				si=i;sj=j;sk=k;
				if(si<=0)si=1;if(si>h)si=h;
				if(sj<=0)sj=1;if(sj>w)sj=w;
				if(sk<=0)sk=1;if(sk>d)sk=d;
				spid = (si + 0)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk + 0);
				{
					Tr.push_back(Triplet<double>(pid,spid,6));
					qid = (si + 1)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk + 0);
					Tr.push_back(Triplet<double>(pid,qid,-1));
					qid = (si - 1)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk + 0);
					Tr.push_back(Triplet<double>(pid,qid,-1));
					qid = (si - 0)*(w+2)*(d+2) + (sj + 1)*(d+2) + (sk + 0);
					Tr.push_back(Triplet<double>(pid,qid,-1));
					qid = (si - 0)*(w+2)*(d+2) + (sj - 1)*(d+2) + (sk + 0);
					Tr.push_back(Triplet<double>(pid,qid,-1));
					qid = (si - 0)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk + 1);
					Tr.push_back(Triplet<double>(pid,qid,-1));
					qid = (si - 0)*(w+2)*(d+2) + (sj + 0)*(d+2) + (sk - 1);
					Tr.push_back(Triplet<double>(pid,qid,-1));
				}
				bb[pid] = div[si][sj][sk];
				x[pid] = div[si][sj][sk];
			}
	A.setFromTriplets(Tr.begin(),Tr.end());
	solver.setTolerance(0.00000001);
	solver.compute(A);
	int cnt=0;
	do{
		x = solver.solveWithGuess(bb,x);
	for(i=1;i<=h;i++)
		for(j=1;j<=w;j++)
			for(k=1;k<=d;k++)
			{
				int pid;
				pid = (i + 0)*(w+2)*(d+2) + (j + 0)*(d+2) + (k + 0);
				p[i][j][k] = x[pid];
			}
	set_bnd(h, w, d, 0, p);
	cnt++;
	}while(cnt<2);
//	printf("pro %d %d %.10lf\n",solver.iterations(),solver.info(), solver.error());
	*/
	
	// gause-saidel
	for(iter=0;iter<20;iter++)
	{
	#pragma omp parallel for
	for(l=0;l<h*w*d;l++)
		{
			int i,j,k;
			i = l/w/d + 1;
			j = (l/d)%w+1;
			k = l%d+1;
				p[i][j][k] = (div[i][j][k] + p[i+1][j][k]  + p[i-1][j][k] + p[i][j+1][k] + p[i][j-1][k] + p[i][j][k+1] + p[i][j][k-1])/6;
			}
		set_bnd(h, w, d, 0, p);
	}

	#pragma omp parallel for
	for(l=0;l<h*w*d;l++)
		{
			int i,j,k;
			i = l/w/d + 1;
			j = (l/d)%w+1;
			k = l%d+1;
			tu[i][j][k] -= 0.5 * (p[i+1][j][k] - p[i-1][j][k]) / da;
			tv[i][j][k] -= 0.5 * (p[i][j+1][k] - p[i][j-1][k]) / da;
			tz[i][j][k] -= 0.5 * (p[i][j][k+1] - p[i][j][k-1]) / da;
		}
	set_bnd(h, w, d, 1, tu);
	set_bnd(h, w, d, 2, tv);
	set_bnd(h, w, d, 3, tz);
}
/////////////////////////////////////////////////////////////////////////////
static int h;
static int w;
static int d;

StableFluid::StableFluid(int hh, int ww)
{
	int i, j;
	h = hh;
	w = ww;
	d = 2;
	
	width = ww;
	height = hh;
	depth = 2;
	
	a = new double**[height+2];
	u = new double**[height+2];
	v = new double**[height+2];
	z = new double**[height+2];
	f = new vector3f**[height+2];
	tu = new double**[h+2];
	tv = new double**[h+2];
	ta = new double**[h+2];
	tz = new double**[h+2];
	for(int i=0;i<height+2;i++)
	{
		a[i] = new double*[width+2];
		u[i] = new double*[width+2];
		v[i] = new double*[width+2];
		z[i] = new double*[width+2];
		f[i] = new vector3f*[width+2];
		tu[i] = new double*[w+2];
		tv[i] = new double*[w+2];
		ta[i] = new double*[w+2];
		tz[i] = new double*[w+2];
		for(j=0;j<width+2;j++)
		{
			a[i][j] = new double[depth+2];
			u[i][j] = new double[depth+2];
			v[i][j] = new double[depth+2];
			z[i][j] = new double[depth+2];
			f[i][j] = new vector3f[depth+2];
			tu[i][j] = new double[d+2];
			tv[i][j] = new double[d+2];
			tz[i][j] = new double[d+2];
			ta[i][j] = new double[d+2];
		}
	}

	reset();
}
StableFluid::StableFluid(int hh, int ww,int dd)
{
	omp_set_num_threads(4);
	int i, j;
	h = hh;
	w = ww;
	d = dd;
	
	width = ww;
	height = hh;
	depth = dd;
	
	a = new double**[height+2];
	u = new double**[height+2];
	v = new double**[height+2];
	z = new double**[height+2];
	f = new vector3f**[height+2];
	tu = new double**[h+2];
	tv = new double**[h+2];
	ta = new double**[h+2];
	tz = new double**[h+2];
	for(int i=0;i<height+2;i++)
	{
		a[i] = new double*[width+2];
		u[i] = new double*[width+2];
		v[i] = new double*[width+2];
		z[i] = new double*[width+2];
		f[i] = new vector3f*[width+2];
		tu[i] = new double*[w+2];
		tv[i] = new double*[w+2];
		ta[i] = new double*[w+2];
		tz[i] = new double*[w+2];
		for(j=0;j<width+2;j++)
		{
			a[i][j] = new double[depth+2];
			u[i][j] = new double[depth+2];
			v[i][j] = new double[depth+2];
			z[i][j] = new double[depth+2];
			f[i][j] = new vector3f[depth+2];
			tu[i][j] = new double[d+2];
			tv[i][j] = new double[d+2];
			tz[i][j] = new double[d+2];
			ta[i][j] = new double[d+2];
		}
	}

	show_velo = false;
	reset();

}
StableFluid::~StableFluid()
{
	int i, j;
	for(i=0;i<height+2;i++)
	{
		for(j=0;j<width+2;j++)
		{
			delete a[i][j];
			delete v[i][j];  
			delete u[i][j];
			delete z[i][j];
			delete f[i][j];
			delete ta[i][j];
			delete tv[i][j];
			delete tz[i][j];
			delete tu[i][j];
		}
		delete a[i];
		delete v[i];
		delete f[i];
		delete u[i];
		delete z[i];
		delete ta[i];
		delete tv[i];
		delete tu[i];
		delete tz[i];
	}
	delete a;
	delete v;
	delete f;
	delete u;
	delete z;
	delete ta;
	delete tv;
	delete tu;
	delete tz;
}
static int timer = 0;
void StableFluid::reset()
{
	int i, j, k;
	
	for(i=0;i<height+2;i++)
		for(j=0;j<width+2;j++)
			for(k=0;k<depth+2;k++)
		{
			a[i][j][k] = 0; 
			v[i][j][k] = 0;
			u[i][j][k] = 0;
			z[i][j][k] = 0;
			f[i][j][k] = vector3f();
			tu[i][j][k] = 0;
			tv[i][j][k] = 0;
			ta[i][j][k] = 0;
			tz[i][j][k] = 0;
		}
	dx = 1.0;
	dy = 1.0;
	dz = 1.0;
	visc = 0;
	diff = 0.0001;
	dt = 0.1;
	timer = 0;

}
vector3f wahaha(int h,int w, int d, int si, int sj,int sk)
{
	/*
	double px = rand()%101 - 50;
	double py = rand()%101 - 50;
	px /= 400;
	py /= 400;
	return vector2f(px,py);

	vector3f res = vector3f();
	if(timer%40 < 20)
	res += vector3f(0,(d/2 - sk)/100.0,(sj - w/2)/100.0);
	else
	res += vector3f(0,(d/2 - sk)/100.0,(sj - w/2)/100.0);
//	if(abs(atan2((sj-w/2),(si-h/2)) - timer/10.0) < 0.4) res = vector3f(cos(timer/10.0)*(sk*2.0)/d/5,-sin(timer/10.0)*(sk*2.0)/d/5,0);
	res.normalize();
	res /= 4.0;
	return res;
	*/
	vector3f res = vector3f();
	
	if( ((sk>(float)d/4.0f) && (sk<(float)d/4.0f*3.0f)) && 
		((sj>(float)w/4.0f) && (sj<(float)w/4.0f*3.0f)) && 
		(si<(float)h/7.0f))
		res = vector3f(2, 0, 0);


	return res;

}
void StableFluid::pipeline()
{
	int i, j, k;
	timer++;
	for(i=1;i<=h;i++)
		for(j=1;j<=w;j++)
			for(k=1;k<=d;k++)
		{
		if(timer % 20 < 10)
				f[i][j][k] += wahaha(h,w,d,i,j,k);
			u[i][j][k] += f[i][j][k].x;
			v[i][j][k] += f[i][j][k].y;
			z[i][j][k] += f[i][j][k].z;
			if(0.0 != f[i][j][k].dist())
				a[i][j][k] += 0.1f;
			f[i][j][k] = vector3f();
		}
	// set velocity
	diffuse(h, w, d, 1, tu, u, visc, dt);
	diffuse(h, w, d, 2, tv, v, visc, dt);
	diffuse(h, w, d, 3, tz, z, visc, dt);
	project(h, w, d, tu, tv, tz, u, v);
	advect(h, w, d, 1, u, tu, tu, tv, tz, dt);
	advect(h, w, d, 2, v, tv, tu, tv, tz, dt);
	advect(h, w, d, 3, z, tz, tu, tv, tz, dt);
	project(h, w, d, u, v, z, tu, tv);

	// set ipja
	diffuse(h, w, d, 0, ta, a, diff, dt);
	advect(h, w, d, 0, a, ta, u, v, z, dt);
	
	/*
	for(i=-h/25;i<h/25; i++)
		for(j=-w/25;j<w/25;j++)
		{
			a[h/2+i][w/2+j] = 1.0;
		}
	{
		if(timer % 20 < 10)
		{
			for(i=1;i<=h;i++)
				a[i][w/2][d/2] = 1;
			for(i=1;i<=w;i++)
				a[h/2][i][d/2] = 1;
			for(i=1;i<=d;i++)
				a[h/2][w/2][i] = 1;
		}
	}
		*/

}

void draw_boundarybox(int h, int w,int d, double dx, double dy, double dz)
{
	
	glColor3d(1, 1, 1);
	glLineWidth(2);
	double hei = (h+2) * dy;
	double wid = (w+2) * dx;
	double dep = (d+2) * dz;
	glBegin(GL_LINES);
	// draw boundary box
	glColor3d(1, 1, 1);
	glVertex3d(0.0, 0.0, 0.0);
	glColor3d(0, 1, 1);
	glVertex3d(0.0, 0.0, dep);
	glColor3d(1, 1, 1);
	glVertex3d(0.0, 0.0, 0.0);
	glColor3d(0, 1, 1);
	glVertex3d(0.0, hei, 0.0);
	glColor3d(1, 1, 1);
	glVertex3d(0.0, 0.0, 0.0);
	glColor3d(0, 1, 1);
	glVertex3d(wid, 0.0, 0.0);
	glColor3d(1, 1, 1);
	glVertex3d(wid, 0.0, 0.0);
	glColor3d(1, 0, 1);
	glVertex3d(wid, hei, 0.0);
	glColor3d(1, 1, 1);
	glVertex3d(wid, 0.0, 0.0);
	glVertex3d(wid, 0.0, dep);
	glVertex3d(0.0, hei, 0.0);
	glVertex3d(wid, hei, 0.0);
	glVertex3d(0.0, hei, 0.0);
	glVertex3d(0.0, hei, dep);
	glVertex3d(0.0, 0.0, dep);
	glVertex3d(wid, 0.0, dep);
	glVertex3d(0.0, 0.0, dep);
	glVertex3d(0.0, hei, dep);
	glColor3d(1, 1, 1);
	glVertex3d(wid, hei, 0.0);
	glColor3d(1, 1, 0);
	glVertex3d(wid, hei, dep);
	glColor3d(1, 1, 1);
	glVertex3d(wid, 0.0, dep);
	glColor3d(1, 1, 0);
	glVertex3d(wid, hei, dep);
	glColor3d(1, 1, 1);
	glVertex3d(0.0, hei, dep);
	glColor3d(1, 1, 0);
	glVertex3d(wid, hei, dep);
	glEnd();
	glColor3d(1, 1, 1);
}
void StableFluid::draw_box(double sx, double sy, double sz, double ex,double ey, double ez,double alpha)
{
	if(alpha > 1.0) alpha = 1;
	if(alpha < 0) alpha = 0;
	glPushMatrix();
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glColor4d(alpha, alpha*alpha, alpha*alpha*alpha*alpha, alpha);
	glBegin(GL_TRIANGLES);
	glVertex3d(sx,sy,sz);	glVertex3d(sx,sy,ez);	glVertex3d(sx,ey,sz);
	glVertex3d(sx,ey,sz);	glVertex3d(sx,sy,ez);	glVertex3d(sx,ey,ez);
	glVertex3d(sx,sy,ez);	glVertex3d(ex,sy,ez);	glVertex3d(sx,ey,ez);
	glVertex3d(sx,ey,ez);	glVertex3d(ex,sy,ez);	glVertex3d(ex,ey,ez);
	glVertex3d(ex,sy,ez);	glVertex3d(ex,sy,sz);	glVertex3d(ex,ey,ez);
	glVertex3d(ex,ey,ez);	glVertex3d(ex,sy,sz);	glVertex3d(ex,ey,sz);
	glVertex3d(ex,sy,sz);	glVertex3d(sx,sy,sz);	glVertex3d(ex,ey,sz);
	glVertex3d(ex,ey,sz);	glVertex3d(sx,sy,sz);	glVertex3d(sx,ey,sz);

	
	glVertex3d(sx,ey,sz);	glVertex3d(sx,ey,ez);	glVertex3d(ex,ey,sz);
	glVertex3d(ex,ey,sz);	glVertex3d(sx,ey,ez);	glVertex3d(ex,ey,ez);
	glVertex3d(sx,sy,sz);	glVertex3d(ex,sy,sz);	glVertex3d(sx,sy,ez);
	glVertex3d(sx,sy,ez);	glVertex3d(ex,sy,sz);	glVertex3d(ex,sy,ez);
	glEnd();
	glColor4d(1.0, 1.0, 1.0, 1.0);
	glPopMatrix();
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
}
void StableFluid::draw_velo(double sx, double sy, double sz, double vx,double vy, double vz)
{
	if(vy > 0)
	glColor4d(0, 1, 1, 1);
	else
	glColor4d(1, 1, 0, 1);
	glBegin(GL_LINES);
	glVertex3d(sx, sy, sz);
	glVertex3d(sx, sy+vy, sz);
	glEnd();
	glColor4d(1.0, 1.0, 1.0, 1.0);
}

int mmx, mmy, mmz;
double *vv;
void StableFluid::mian()
{
	FILE *fp;

	int mx, my, mz;
	int nx, ny, nz;
	int nw = max(nx, max(ny, nz));
	int nn;
	mx = h;
	my = w;
	mz = d;
	mmx = mx;
	mmy = my;
	mmz = mz;
	vv = new double[mx*my*mz];

	fp = fopen("test.sdf", "rt");
	fscanf(fp, "%d %d %d", &nx, &ny, &nz);
	double offset_x, offset_y, offset_z;
	fscanf(fp, "%lf %lf %lf", &offset_x, &offset_y, &offset_z);
	double dx;
	fscanf(fp, "%lf", &dx);
	nn = nx*ny*nz;
	double *y_sdf = new double[nn];
	auto yyy = [nx, ny, nz](int i, int j, int k)->int{
		int id = i*ny*nz + j*nz + k;
		return id;
	};
	int i, j, k;

	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			for (k = 0; k < nz; k++)
			{
				int id = yyy(i, j, k);
				double v;
				fscanf(fp, "%lf", &v);
				y_sdf[id] = v;
			}

	fclose(fp);

	auto xxx = [&yyy, &y_sdf, nx, ny, nz, mx, my, mz](double i, double j, double k) ->double {

		float xx = ((float)i / mx) * nx;
		float yy = ((float)j / my) * ny;
		float zz = ((float)k / mz) * nz;

		int sx, sy, sz;
		int tx, ty, tz;

		sx = (int)xx;
		tx = sx + 1;
		sy = (int)yy;
		ty = sy + 1;
		sz = (int)zz;
		tz = sz + 1;

		if (sx < 0) sx = 0; if (sx >= nx) sx = nx - 1;
		if (sy < 0) sy = 0; if (sy >= ny) sy = ny - 1;
		if (sz < 0) sz = 0; if (sz >= nz) sz = nz - 1;
		if (tx < 0) tx = 0; if (tx >= nx) tx = nx - 1;
		if (ty < 0) ty = 0; if (ty >= ny) ty = ny - 1;
		if (tz < 0) tz = 0; if (tz >= nz) tz = nz - 1;

		float pp, pq, qp, qq;
		float rp, rq;
		pp = xx - sx; pq = 1 - pp;
		qp = yy - sy; qq = 1 - qp;
		rp = zz - sz; rq = 1 - rp;

		double res = ((y_sdf[yyy(sx, sy, sz)] * pq + y_sdf[yyy(tx, sy, sz)] * pp) * qq +
			(y_sdf[yyy(sx, ty, sz)] * pq + y_sdf[yyy(tx, ty, sz)] * pp) * qp) * rq +
			((y_sdf[yyy(sx, sy, tz)] * pq + y_sdf[yyy(tx, sy, tz)] * pp) * qq +
			(y_sdf[yyy(sx, ty, tz)] * pq + y_sdf[yyy(tx, ty, tz)] * pp) * qp) * rp;
//			printf("xxx - %lf %lf %lf  %d %d %d  %lf %lf %lf   %lf %lf\n", xx, yy, zz, sx, sy, sz, pp, qp, rp, y_sdf[yyy(sx,sy,sz)], res);
		return res;

	};
	for (k = 0; k < mz; k++)
	{
		for (i = 0; i < mx; i++)
		{
			for (j = 0; j < my; j++)
			{
				int id = i * my * mz + j * mz + k;

				double v = xxx(i, j, k);
				vv[id] = -v * 10;
				//(v + 1.0) / 2.0;
			}
		}
	}
	delete[] y_sdf;
}
void StableFluid::draw_fog(int h, int w, int d, double ***a,double dx,double dy,double dz)
{
	int i, j, k;
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	//NOTE : i is y axis, j is x axis
	if(show_velo)
	{
		for(i=1;i<=h;i++)
		{
			for(j=1;j<=w;j++)
				for(k=1;k<=d;k++)
			{
				draw_velo(j*dx, i*dy, k*dz, v[i][j][k], u[i][j][k], z[i][j][k]);
			}
		}
	}
	else
	{
		for(i=1;i<=h;i++)
		{
//			printf("%d\n",i);
			for(j=1;j<=w;j++)
				for(k=1;k<=d;k++)
			{
				draw_box(j*dx, i*dy, k*dz, j*dx + dx, i*dy + dy, k*dz + dz, a[i][j][k]);
//				draw_box(j*dx, i*dy, k*dz, j*dx + dx, i*dy + dy, k*dz + dz, vv[(i-1)*w*d + (j-1)*d + (k-1)]);
			}
		}
	}

}
void StableFluid::draw()
{
	draw_boundarybox(h, w, d, dx, dy,dz);
	draw_fog(h,w,d, a,dx,dy,dz);
}