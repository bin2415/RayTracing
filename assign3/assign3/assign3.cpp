/*
CSCI 480
Assignment 3 Raytracer

Name: <庞成宾>
*/
#define _CRT_SECURE_NO_WARNINGS
#define GLUT_DISABLE_ATEXIT_HACK
#include <pic.h>

#include <windows.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <vector>
#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10


char *filename=0;

#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

#define WIDTH 640
#define HEIGHT 480

#define fov 60.0
#define PI 3.1415926


int MaxStep = 10;
int Steps = 0;
//求出投射的屏幕的x和y的最大范围的坐标
double yMax = tan((double)PI*fov / (2 * 180));
double xMax = yMax*((double)WIDTH) / ((double)HEIGHT);

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];    //漫射
  double color_specular[3];	  //反射
  double normal[3];
  double shininess;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

//表示点的结构
struct point
{
	double x;
	double y;
	double z;
};
//点在里面还是在外面
struct isIn
{
	bool in;		//当在三角形外面时为0，当在里面时为1
	double bary[3];
};

//交点的结构
struct intexPoint
{
	point p;
	double t;
	int tID;	//物体在数组中的标号
	int tObj;  //如果是1表示点在圆上，如果为2表示在三角形上
	isIn iO;	
};
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

struct point cam;

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

point reflect(intexPoint p, point dir);
point Render(point p, point dir);
//两个点进行相减
point minusPoint(point A, point B)
{
	point C;
	C.x = A.x - B.x;
	C.y = A.y - B.y;
	C.z = A.z - B.z;

	return C;
}

//向量除以一个数
point DivConst(point A, double a)
{
	point B;

	B.x = 0.0;
	B.y = 0.0;
	B.z = 0.0;

	if (abs(a) > 1e-10)
	{
		B.x = A.x / a;
		B.y = A.y / a;
		B.z = A.z / a;
	}

	return B;
}
//计算向量的大小
double caculateSize(point A)
{
	double size;
	size = sqrt(pow(A.x, 2) + pow(A.y, 2) + pow(A.z, 2));
	return size;
}

//单位化向量
point unitize(point A)
{
	point uni;
	double size;
	size = caculateSize(A);
	uni = DivConst(A, size);

	return uni;
}
//点乘
double dot(point A, point B)
{
	double C;

	C = (A.x*B.x + A.y*B.y + A.z*B.z);
	return C;
}
//叉乘
point cross(point A, point B)
{
	point C;
	C.x = (A.y*B.z - B.y*A.z);
	C.y = (B.x*A.z - A.x*B.z);
	C.z = (A.x*B.y - A.y*B.x);

	return C;

}

//求出距离src，方向为dir，长度为t的点
point caluPoint(point src, point dir, double t)
{
	point p;

	p.x = src.x + t*(dir.x);
	p.y = src.y + t*(dir.y);
	p.z = src.z + t*(dir.z);

	return p;
}

//与球的交点
double intersectSphere(Sphere sphere, point src, point dir)
{
	double b, c, t, t1, t2;
	t1 = 0;
	t2 = 0;

	c = pow((src.x - sphere.position[0]), 2) + pow((src.y - sphere.position[1]), 2)
		+ pow((src.z - sphere.position[2]), 2) - pow(sphere.radius, 2);
	b = 2 * (dir.x*(src.x - sphere.position[0])
		+ dir.y*(src.y - sphere.position[1])
		+ dir.z*(src.z - sphere.position[2])
		);

	//检查判别式是否大于0
	if ((pow(b, 2) - 4 * c) > 0)
	{
		t1 = (((-1)*b) + sqrt(pow(b, 2) - 4 * c)) / 2;
		t2 = (((-1)*b) - sqrt(pow(b, 2) - 4 * c)) / 2;

		if (t1 <= t2)
			t = t1;
		else
			t = t2;
		if (t < 0)
			t = -1;
		else if (t<1e-10)
		{
			if (t1<1e-15 && t2>1e-15)
				t = t2;
			else if (t2<1e-15 && t1>1e-15)
				t = t1;
			else t = -1; 
		}
	}
	else t = -1; 

	return t;
}

//获得两个顶点之间的边
point getSide(Vertex v1, Vertex v2)
{
	point c;
	c.x = v1.position[0] - v2.position[0];
	c.y = v1.position[1] - v2.position[1];
	c.z = v1.position[2] - v2.position[2];

	return c;
}
//判断两个点是否相等
bool checkEqual(point A, point B)
{
	bool equ;
	if ((abs(A.x - B.x)<1e-10) && (abs(A.y - B.y)<1e-10) && (abs(A.z - B.z)<1e-10))
		equ = 1;
	else
		equ = 0;
	return equ;
}

double intersectTriangle(Triangle triangle, point src, point dir,isIn* iO)
{
	point AB, AC, DirxAC;
	float u, v, t;
	//isIn iO;
	AB = getSide(triangle.v[1], triangle.v[0]);
	AC = getSide(triangle.v[2], triangle.v[0]);
	DirxAC = cross(dir, AC);
	float det = dot(AB, DirxAC);
	point T;
	point p1;
	p1.x = triangle.v[0].position[0];
	p1.y = triangle.v[0].position[1];
	p1.z = triangle.v[0].position[2];
	if (det >0)
   {
	   T = minusPoint(src, p1);
}
	else
	{
		T = minusPoint(p1, src);
	    det = -det;
 }
	if (det < 1e-10)
	{
		iO->in = -1;
		return -1;
	}

	u = dot(T, DirxAC);
	if (u < 0.0f || u > det)
	{
		iO->in = -1;
		return -1;

	}
	point Q = cross(T, AB);
	v = dot(dir, Q);
	if (v < 0.0f || u + v > det)
	{
		iO->in = -1;
		return -1;
	}

	t = dot(AC, Q);
	//t = -t;
	float fInvDet = 1.0f / det;
	t *= fInvDet;
	u *= fInvDet;
	v *= fInvDet;

	iO->in = 1;
	iO->bary[0] = u;
	iO->bary[1] = v;
	iO->bary[2] = 1 - u - v;
	return t;
}

//找到p点在圆上的的法向量，tID为spheres数组中的标号
point findSphereNormal(point p, int tID)
{
	point n;

	// based on the equation
	n.x = (p.x - spheres[tID].position[0]) / spheres[tID].radius;
	n.y = (p.y - spheres[tID].position[1]) / spheres[tID].radius;
	n.z = (p.z - spheres[tID].position[2]) / spheres[tID].radius;
	return n;
}

//三角形的线性插值,如果ID为0则表示法线插值,如果ID为1则表示漫反射，如果为2则表示镜面反射
point chazhi(Triangle triangle, isIn iO, int ID)
{
	point P;
	if (ID == 0)
	{
		P.x = iO.bary[0] * triangle.v[0].normal[0]
			+ iO.bary[1] * triangle.v[1].normal[1]
			+ iO.bary[2] * triangle.v[2].normal[2];

		P.y = iO.bary[0] * triangle.v[0].normal[1]
			+ iO.bary[1] * triangle.v[1].normal[1]
			+ iO.bary[2] * triangle.v[2].normal[1];

		P.z = iO.bary[0] * triangle.v[0].normal[2]
			+ iO.bary[1] * triangle.v[1].normal[2]
			+ iO.bary[2] * triangle.v[2].normal[2];
	}
	else if (ID == 1)
	{
		P.x = iO.bary[0] * triangle.v[0].color_diffuse[0]
			+ iO.bary[1] * triangle.v[1].color_diffuse[0]
			+ iO.bary[2] * triangle.v[2].color_diffuse[0];

		P.y = iO.bary[0] * triangle.v[0].color_diffuse[1]
			+ iO.bary[1] * triangle.v[1].color_diffuse[1]
			+ iO.bary[2] * triangle.v[2].color_diffuse[1];

		P.z = iO.bary[0] * triangle.v[0].color_diffuse[2]
			+ iO.bary[1] * triangle.v[1].color_diffuse[2]
			+ iO.bary[2] * triangle.v[2].color_diffuse[2];
	}

	else if (ID == 2)
	{
		P.x = iO.bary[0] * triangle.v[0].color_specular[0]
			+ iO.bary[1] * triangle.v[1].color_specular[0]
			+ iO.bary[2] * triangle.v[2].color_specular[0];

		P.y = iO.bary[0] * triangle.v[0].color_specular[1]
			+ iO.bary[1] * triangle.v[1].color_specular[1]
			+ iO.bary[2] * triangle.v[2].color_specular[1];

		P.z = iO.bary[0] * triangle.v[0].color_specular[2]
			+ iO.bary[1] * triangle.v[1].color_specular[2]
			+ iO.bary[2] * triangle.v[2].color_specular[2];
	}

	return P;
}

//公式模型I_spec = k_z * I_l(V • ((2N • L)N - L))^n_s
//其中R = (2N • L)N - L
point phong(point p, int id, int Obj, isIn iO, Light light, point camera)
{
	point n, l, v, r, kd, ks;
	point po;
	double lDotN, rDotV, n_s;

	//用point表示light的位置
	l.x = light.position[0];
	l.y = light.position[1];
	l.z = light.position[2];

	//l = l - p;
	//入射光线
	l = unitize(minusPoint(l, p));

	//v = camera - p;
	//观察者到p点的射线
	v = unitize(minusPoint(camera, p));

	if (Obj == 1)
	{
		n = findSphereNormal(p, id);

		kd.x = spheres[id].color_diffuse[0];
		kd.y = spheres[id].color_diffuse[1];
		kd.z = spheres[id].color_diffuse[2];

		ks.x = spheres[id].color_specular[0];
		ks.y = spheres[id].color_specular[1];
		ks.z = spheres[id].color_specular[2];

		n_s = spheres[id].shininess;
	}

	else
		if (Obj == 2)
		{
			n = unitize(chazhi(triangles[id], iO, 0));
			kd = chazhi(triangles[id], iO, 1);
			ks = chazhi(triangles[id], iO, 2);

			n_s = iO.bary[0] * triangles[id].v[0].shininess
				+ iO.bary[1] * triangles[id].v[1].shininess
				+ iO.bary[2] * triangles[id].v[2].shininess;
		}


	lDotN = dot(l, n);
	if (lDotN<0)
		lDotN = 0;
	else if (lDotN>1.f)
		lDotN = 1.f;

	//R = (2N • L)N - L
	r.x = 2 * lDotN*n.x - l.x;
	r.y = 2 * lDotN*n.y - l.y;
	r.z = 2 * lDotN*n.z - l.z;

	
	rDotV = dot(r, v);
	if (rDotV<0)
		rDotV = 0;
	else if (rDotV>1.f)
		rDotV = 1.f;

	//计算该点的颜色r,g,b
	po.x = light.color[0] * ((kd.x)*lDotN + ((ks.x)*pow((rDotV), (n_s)))); // r
	po.y = light.color[1] * ((kd.y)*lDotN + ((ks.y)*pow((rDotV), (n_s)))); // g 
	po.z = light.color[2] * ((kd.z)*lDotN + ((ks.z)*pow((rDotV), (n_s)))); // b
	return po;
}

/*
与物体相交
*/
intexPoint intersectObjects(point p1, point p2, point dir, int flag)
{
	point p, q,raySrc, pixPoint;
	intexPoint intxObj;
	isIn iO;
	iO.in = -1;
	double t, t1, t2, tS, tT;
	int id, Obj;

	q = p1;
	t1 = 0;
	id = -1;
	Obj = -1;

	
	raySrc = p1;
	

	//找到最近的交点
	for (int i = 0; i < num_spheres; i++)
	{
		tS = intersectSphere(spheres[i], raySrc, dir);
		if (t1 == 0 && tS > 1e-10)
		{
			t1 = tS;
			id = i;
			Obj = 1;
		}
		else if (tS <= t1 && tS > 1e-10)
		{
			t1 = tS;
			id = i;
			Obj = 1;
		}
	}

	for (int i = 0; i < num_triangles; i++)
	{
		tT = intersectTriangle(triangles[i], raySrc, dir, &iO);
		p = caluPoint(raySrc, dir, tT); //找到交点
									   //iO = isInTest(triangles[i], p); //判断p点是否在三角形内

		if (iO.in == 1)
		{
			if (t1 == 0 && tT > 1e-5)
			{
				t1 = tT;
				id = i;
				Obj = 2;
				if (flag == 0)
					q = p;
				intxObj.iO.bary[0] = iO.bary[0];
				intxObj.iO.bary[1] = iO.bary[1];
				intxObj.iO.bary[2] = iO.bary[2];

			}
			else if (tT<t1 && tT>1e-5)
			{
				t1 = tT;
				id = i;
				Obj = 2;
				if (flag == 0)
					q = p;
				intxObj.iO.bary[0] = iO.bary[0];
				intxObj.iO.bary[1] = iO.bary[1];
				intxObj.iO.bary[2] = iO.bary[2];
			}
		}
	}

	if (flag == 1)
	{
		//如果为1则计算从light到p点的距离t2
		if (dir.x != 0)
		{
			t2 = (p2.x - raySrc.x) / dir.x;
		}
		else if (dir.y != 0)
		{
			t2 = (p2.y - raySrc.y) / dir.y;
		}
		else if (dir.z != 0)
		{
			t2 = (p2.z - raySrc.z) / dir.z;
		}
		else t2 = 0;

		//t2和t1进行比较，如果t1小于t2则说明没有物体挡住光线	
		if (t1 >= t2)
		{
			Obj = -1;
			id = -1;
		}
	}
	else if ((t1 >= 0) && (Obj == 1))
		q = caluPoint(raySrc, dir, t1);

	intxObj.p = q;
	intxObj.t = t1;
	intxObj.tID = id;
	intxObj.tObj = Obj;
	return intxObj;
}

point findColor(int x, int y)
{
	//point p, q, dir, light, lightS;
	//point black, pixColor, temp, tempN;
	//intexPoint intxObj, intxFlag;

	point black, pixColor, p;
	black.x = 0.0;
	black.y = 0.0;
	black.z = 0.0;
	pixColor = black;

	//将像素点转换为世界坐标
	p.x = (((double)x / (double)WIDTH) * 2 * xMax) - xMax;
	p.y = (((double)y / (double)HEIGHT) * 2 * yMax) - yMax;
	p.z = -1;

	
	point dir1, p1 = p;
	dir1 = minusPoint(p, cam);
	dir1 = unitize(dir1);
	pixColor = Render(p, dir1);
	if (pixColor.x > 1) pixColor.x = 1.f;
	if (pixColor.y > 1) pixColor.y = 1.f;
	if (pixColor.z > 1) pixColor.z = 1.f;
	return pixColor;
}

//迭代渲染
point Render(point p, point dir)
{
	Steps++;
	point blackColor, pixColor, q, light, dir1, temp, temp1;
	point reflect_ray;
	intexPoint intxObj, intxFlag;
	blackColor.x = 0.0;
	blackColor.y = 0.0;
	blackColor.z = 0.0;
	pixColor = blackColor;
	if (Steps > MaxStep)
	{
		Steps = 0;
		return blackColor;
	}
	point p1 = p;
	intxObj = intersectObjects(p, p1, dir, 0);
	//如果和一个物体有交点
		if (intxObj.tID != -1)
		{
			q = intxObj.p;

			reflect_ray = reflect(intxObj, dir);
			pixColor.x += ambient_light[0];
			pixColor.y += ambient_light[1];
			pixColor.z += ambient_light[2];

			for (int h = 0; h < num_lights; h++)
			{
				light.x = lights[h].position[0];
				light.y = lights[h].position[1];
				light.z = lights[h].position[2];

				dir1 = minusPoint(light, q);
				dir1 = unitize(dir1);
				intxFlag = intersectObjects(q, light, dir1, 1);

				//如果没有物体遮挡
				if (intxFlag.tID == -1)
				{
					//phong模型求出颜色
					temp = phong(q, intxObj.tID, intxObj.tObj, intxObj.iO, lights[h], cam);
					pixColor.x += temp.x;
					pixColor.y += temp.y;
					pixColor.z += temp.z;
				}
			}
			temp1 = Render(q, reflect_ray);
			int id = intxObj.tID;
			if (intxObj.tObj == 1)
			{
				pixColor.x += temp1.x * spheres[id].color_specular[0];
				pixColor.y += temp1.y * spheres[id].color_specular[1];
				pixColor.z += temp1.z * spheres[id].color_specular[2];
			}
			else
			{
				point tp = chazhi(triangles[id], intxObj.iO, 2);
				pixColor.x += temp1.x * tp.x;
				pixColor.y += temp1.y * tp.y;
				pixColor.z += temp1.z * tp.z;
			}
			
		}
		else
		{
			pixColor = blackColor;
			Steps = 0;
		}
		return pixColor;
}
//反射的光线的方向
point reflect(intexPoint intx, point dir)
{
	point result;
	point n;
	if (intx.tObj == 1)
	{
		n = findSphereNormal(intx.p, intx.tID);
	}
	else
	{
		n = unitize(chazhi(triangles[intx.tID], intx.iO, 0));
	}
	
	dir.x = -dir.x;
	dir.y = -dir.y;
	dir.z = -dir.z;
	double r1 = dot(n, dir);
	point n2;
	n2.x = n.x * 2 * r1;
	n2.y = n.y * 2 * r1;
	n2.z = n.z * 2 * r1;
	
	result = n2;
	return result;
}
void drawColor()
{
	unsigned int x, y;
	point pixColor;

	for (x = 0; x < WIDTH; x++)
	{
		for (y = 0; y < HEIGHT; y++)
		{
			pixColor = findColor(x, y);

			plot_pixel_jpeg(x, y, abs(pixColor.x) * 255, abs(pixColor.y) * 255, abs(pixColor.z) * 255);
		}
	}
}
//MODIFY THIS FUNCTION
void draw_scene()
{
  unsigned int x,y;

  //glPointSize(2.0);
  //glBegin(GL_POINTS);
  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
     // plot_pixel(x,y,x%256,y%256,(x+y)%256);
		plot_pixel_display(x, y, buffer[HEIGHT - y - 1][x][0], buffer[HEIGHT - y - 1][x][1], buffer[HEIGHT - y - 1][x][2]);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}



/*void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}*/

void parse_check(char *expected,char *found)
{
  if(stricmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(stricmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(stricmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(stricmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{
	cam.x = 0.0;
	cam.y = 0.0;
	cam.z = 0.0;
	glLoadIdentity();

	drawColor();
	draw_scene();
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if (!once)
  {
	  draw_scene();
//	  if (mode == MODE_JPEG)
//		  save_jpg();
  }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
