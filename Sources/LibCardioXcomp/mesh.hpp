#ifndef __CardioXcomp__mesh__
#define __CardioXcomp__mesh__

#include <iostream>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include "cassert"
#include <fstream>

#include "platform.h"

extern "C" {
#include <libmesh5.h>
#undef call
  // call conflicts with boost/regex.hpp and will be remove from libmesh5.h
  // in future versions.

  extern int GmfOpenMesh(char *, int, ...);
  extern int GmfCloseMesh(int);
  extern int GmfStatKwd(int, int, ...);
  extern int GmfGotoKwd(int, int);
  extern int GmfSetKwd(int, int, ...);
  extern void GmfGetLin(int, int, ...);
  extern void GmfSetLin(int, int, ...);

}


namespace cardioxcomp {
  
  template<class T> inline T Max (const T &a,const T & b,const T & c){return std::max(std::max(a,b),c);}
  template<class T> inline T Min (const T &a,const T & b,const T & c){return std::max(std::max(a,b),c);}
  
  
  //=============================================================================
  //
  // Class Point
  //
  //=============================================================================

  class Point{
    friend std::ostream& operator <<(std::ostream& f, const Point & P )
    { f << '(' << P.x << ',' << P.y << ')'  ; return f; }
    friend std::istream& operator >>(std::istream& f,  Point & P)
    { f >>  P.x >>  P.y  ; return f; }
    
  public:
    
    double x,y;
    Point ():x(0),y(0){}
    Point (double a,double b):x(a),y(b){}
    Point (const Point& A):x(A.x),y(A.y){}
    Point (Point A,Point B):x(B.x-A.x),y(B.y-A.y) {}
    Point   operator+(Point P)const   {return Point(x+P.x,y+P.y);}
    Point   operator+=(Point P)  {x += P.x;y += P.y;return *this;}
    Point   operator-(Point P)const   {return Point(x-P.x,y-P.y);}
    Point   operator-=(Point P) {x -= P.x;y -= P.y;return *this;}
    Point   operator-()const  {return Point(-x,-y);}
    Point   operator+()const  {return *this;}
    double   operator,(Point P)const  {return  x*P.x+y*P.y;}
    double   operator^(Point P)const {return  x*P.y-y*P.x;}
    Point   operator*(double c)const {return Point(x*c,y*c);}
    Point   operator/(double c)const {return Point(x/c,y/c);}
    double  &  operator[](int i){ return (&x)[i];}
    Point   perp() {return Point(-y,x);}
    friend Point operator*(double c,Point P) {return P*c;}
    friend Point operator/(double c,Point P) {return P/c;}
  };
  
  
  inline double area2(const Point A,const Point B,const Point C){return (B-A)^(C-A);}
  inline double norm2_2(const Point & A){ return (A,A);}
  inline double norm2(const Point & A){ return sqrt((A,A));}
  inline double norm_infty(const Point & A){return std::max(fabs(A.x),fabs(A.y));}
  inline double theta(Point P){ return atan2(P.y,P.x);}
  
  const short VerticesOfTriangularEdge[3][2] = {{1,2},{2,0},{0,1}};
  const short EdgesVertexTriangle[3][2] = {{1,2},{2,0},{0,1}};
  const short OppositeVertex[3] = {0,1,2};
  const short OppositeEdge[3] =  {0,1,2};
  const short NextEdge[3] = {1,2,0};
  const short PreviousEdge[3] = {2,0,1};
  const short NextVertex[3] = {1,2,0};
  const short PreviousVertex[3] = {2,0,1};
  const int onWhatIsEdge[3][7] = { {0,1,3, 2,0,0, 0}, // edge 0
    {3,0,1, 0,2,0, 0},
    {1,3,0, 0,0,2, 0}};
  
  const Point TriangleHat[3]= { Point(0,0),Point(1,0),Point(0,1) } ;
  
  //=============================================================================
  //
  // Class Label
  //
  //=============================================================================

  class Label {  // reference number for the physics
  public:
    int ref;
    Label(int r=0):ref(r){}
    Label(const Label& lab):ref(lab.ref){}
    int operator!() const{return !ref;}
  };
  
  inline std::ostream& operator <<(std::ostream& f,const Label & r  )
  { f <<  r.ref ; return f; }
  inline std::istream& operator >>(std::istream& f, Label & r  )
  { f >>  r.ref ; return f; }
  
  //=============================================================================
  //
  // Class Vertex
  //
  //=============================================================================

  class Vertex :
  public Point,
  public Label {
  public:
    Vertex() : Point(),Label(){}
    Vertex(Point P,int r=0): Point(P),Label(r){}
    Vertex(const Vertex& v): Point( (Point) v), Label( (Label) v){}
  };
  
  inline std::ostream& operator <<(std::ostream& f, const Vertex & v )
  { f << (Point) v << ' ' << (Label) v   ; return f; }
  
  inline std::istream& operator >> (std::istream& f,  Vertex & v )
  { f >> (Point&) v >> (Label&) v ; return f; }
  
  //=============================================================================
  //
  // Class Triangle
  //
  //=============================================================================

  
  class Triangle:
  public Label
  {
  public:
    friend std::ostream& operator <<(std::ostream& f, const Triangle & t );
    double area;
    Vertex *vertices[3]; // an array of 3 pointer to vertex
    std::vector<int> numNeu;// global number of vertex on Triangle

    Triangle(){}              // constructor empty for array
    Vertex & operator[](int i) const// to see triangle as a array of vertex
    {return *vertices[i];}
    
    Triangle(Vertex * v0,int i0,int i1,int i2,int r,double a=0.0):
    Label(r)
    {
      Point A=*(vertices[0]=v0+i0);
      Point B=*(vertices[1]=v0+i1);
      Point C=*(vertices[2]=v0+i2);
      area = a ==0 ? (( B-A)^(C-A))*0.5 : a;
      assert(area>=0);

      numNeu.push_back(i0);
      numNeu.push_back(i1);
      numNeu.push_back(i2);

    }
    
    Point Edge(int i) const // opposite edge vertex i
    {return (Point) *vertices[(i+2)%3]-(Point) *vertices[(i+1)%3];}
    
    Vertex & Edge(int j,int i) const // Vertex j of edge i
    {assert(j==0 || j==1 );return  *vertices[(i+j+1)%3];}
    
    Point n(int i) const //  unit exterior normal
    {Point E=Edge(i);return Point(E.y,-E.x)/norm2(E);}
    
    Point H(int i)  const  // heigth ($\nabla \lambda_i$)
    {Point E=Edge(i);return Point(-E.y,E.x)/(2*area);}
    
    double lenEdge(int i) const {Point E=Edge(i);return sqrt((E,E));}
    double h() const { return Max(lenEdge(0),lenEdge(1),lenEdge(2));}
    
    void Renum(Vertex   *v0, long * r)  {
      for (int i=0;i<3;i++)
        vertices[i]=v0+r[vertices[i]-v0];}
    
    Vertex & VerticeOfEdge(int i,int j) const  // vertex j of edge i
    {return  *vertices[(i+1+j)%3];}  // vertex j of edge i
    
    void SetVertex(int j,Vertex *v){vertices[j]=v;}
    
    void computeArea(){
      Point A=*(vertices[0]);
      Point B=*(vertices[1]);
      Point C=*(vertices[2]);
      area = (( B-A)^(C-A))*0.5;
      assert(area>=0);
    }
    Point operator()  (const Point & P) const{ // local to Global in triangle
      return    (const Point &) *vertices[0] * (1-P.x-P.y)
      +  (const Point &) *vertices[1] * (P.x)
      +  (const Point &) *vertices[2] * (P.y) ;}
    
  };
  
  inline std::ostream& operator <<(std::ostream& f, const Triangle & t )
  { f << (Point) *(t.vertices[0]) << ' ' <<(Point) *(t.vertices[1]) << ' ' <<(Point) *(t.vertices[2]) << ' ' << (Label) t; return f; }

  //=============================================================================
  //
  // Class BoundaryEdge
  //
  //=============================================================================

  class BoundaryEdge: public Label {
  public:
    Vertex *vertices[2];
    BoundaryEdge(Vertex * v0,int i0,int i1,int r): Label(r)
    { vertices[0]=v0+i0; vertices[1]=v0+i1; }
    bool in(const Vertex * pv) const {return pv == vertices[0] || pv == vertices[1];}
    BoundaryEdge(){} // constructor empty for array
    Vertex & operator[](int i) const {return *vertices[i];}
    void Renum(Vertex   *v0, long * r) {
      for (int i=0;i<2;i++)
        vertices[i]=v0+r[vertices[i]-v0];}
  };
  
  //=============================================================================
  //
  // Class Mesh
  //
  //=============================================================================

  
  
  class Mesh
  {
  public:
    int nt,nv,neb;
    double area;
    Vertex *vertices;
    Triangle *triangles;
    BoundaryEdge  *bedges;

    Triangle & operator[](int i) const {assert(i>=0 && i<nt);return triangles[i];}
    Vertex & operator()(int i) const {assert(i>=0 && i<nv);return vertices[i];}
    Mesh(const char * filename,int rankProc) {read(filename,rankProc);} // read on a file
    Mesh(const std::string s,int rankProc) {read(s.c_str(),rankProc);}
    Mesh(const Mesh& global_mesh,const std::vector<int>& eltPartition,int rank,std::vector<PetscInt>& loc2globElem);
    ~Mesh();
    int number(const Triangle & t) const {return &t - triangles;}
    int number(const Triangle * t) const {return t  - triangles;}
    int number(const Vertex & v)   const {return &v - vertices;}
    int number(const Vertex * v)   const {return v  - vertices;}
    int operator()(const Triangle & t) const {return &t - triangles;}
    int operator()(const Triangle * t) const {return t  - triangles;}
    int operator()(const Vertex & v)   const {return &v - vertices;}
    int operator()(const Vertex * v)   const {return v  - vertices;}
    int operator()(int it,int j) const {return number(triangles[it][j]);}
    void BoundingBox(Point & Pmin,Point &Pmax) const;
    int renum();
    int gibbsv (long* ptvoi,long* vois,long* lvois,long* w,long* v);
    int TriangleAdj(int it,int &j) const
    {int i=TheAdjacencesLink[3*it+j];j=i%3;return i/3;}
    void VerticesNumberOfEdge(Triangle & T,int j,int & j0,int & j1) const
    {j0 =  number(T[(j+1)%3]),j1=  number(T[(j+ 2)%3]);}
    int BoundaryTriangle(int be,int & edgeInT) const {
      int i= BoundaryEdgeHeadLink[be];
      edgeInT = i%3;
      return i/3;
    }
    Triangle * Find(const Point & P) const ;
    
    BoundaryEdge * TheBoundaryEdge(int i,int j)  const
    {
      int ii=std::min(i,j);
      for (int p=BoundaryAdjacencesHead[ii];p>=0;p=BoundaryAdjacencesLink[p])
        if ( bedges->in(vertices+ii) )   return bedges+p;
      return 0;
    }
    friend std::ostream& operator <<(std::ostream& f, const Mesh & m );
    void print(int rankProc) const;
    void initVerticesPerLabel(int rankProc);
    std::set<int> labelOfVertex;
    std::map< int, std::vector<int> > vertexPerLabel;
  private:
    void read(const char * filename,int rankProc);
    // to construct the adj triangle
    int *TheAdjacencesLink;
    int *BoundaryEdgeHeadLink;
    void ConsAdjacence();
    int *BoundaryAdjacencesHead;
    int *BoundaryAdjacencesLink;

  };
  
}

#endif
