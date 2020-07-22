#include "stdafx.h"

#include <iostream>
#include <fstream>
#include <petscsys.h>

#include "mesh.hpp"

#undef VERBOSE_DEBUG

using namespace std;
namespace cardioxcomp
{
  
  
  /*!
   \brief Create the local mesh from the global mesh
   \param[in] meshGlobal The global mesh
   \param[in] eltPartition The partitioning of element.
   \param[in] rankProc The rank of the process.
   \param[out] loc2globElem Local to global mapping of the elements.
   */
  Mesh::Mesh(const Mesh& global_mesh,const std::vector<int>& eltPartition,int rankProc,std::vector<PetscInt>& loc2globElem)
  {
    nv = global_mesh.nv;
    vertices = new Vertex[nv];
    for (int i=0; i<nv; i++) {
      vertices[i] = global_mesh.vertices[i];
    }
    nt = 0;
    for (int it=0; it<global_mesh.nt; it++) {
      if (eltPartition[it]==rankProc) {
        nt++;
      }
    }
    loc2globElem.resize(nt);
    triangles = new Triangle[nt];
    int it_local=0;
    for (int it=0; it<global_mesh.nt; it++) {
      if (eltPartition[it]==rankProc) {
        Triangle& locTria = triangles[it_local];
        Triangle& globTria = global_mesh.triangles[it];
        loc2globElem[it_local] = it;
        //        std::cout << "Rank #"<< rankProc << " loc2globElem[" << it_local << "] = " << it << std::endl;
        for(int i=0;i<3;i++ )locTria.SetVertex(i, &vertices[global_mesh.number(globTria[i])]);
        locTria.computeArea();
        locTria.ref = globTria.ref;
        it_local++;
      }
    }
    neb = 0;
    // to do: boundary edges
    //
    initVerticesPerLabel(rankProc);
  }

  void Mesh::initVerticesPerLabel(int rankProc){
    int label;
    for (int i=0; i<nv; i++) {
      label = vertices[i].ref;
      labelOfVertex.insert(label);
      vertexPerLabel[label].push_back(i);

      //For the MEA's corners
      if(vertices[i].x == 1 && vertices[i].y == 1){
        vertexPerLabel[2].push_back(i);
      }

      if(vertices[i].x == 0 && vertices[i].y == 0){
        vertexPerLabel[1].push_back(i);
      }

      if(vertices[i].x == 1 && vertices[i].y == 0){
        vertexPerLabel[1].push_back(i);
      }

      if(vertices[i].x == 0 && vertices[i].y == 1){
        vertexPerLabel[3].push_back(i);
      }

    }
#ifdef VERBOSE_DEBUG
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n========  List of vertices per label on Processor #%d  ========\n",rankProc);
 //   PetscSynchronizedPrintf(PETSC_COMM_WORLD,"numDofProc = %d\n",numDofProc);
   /* for (int i=0; i<nv; i++) {
      label = vertices[i].ref;
      verticesPerLabel[label].push_back(i);
    }*/
    for(set<int>::const_iterator ilabel=labelOfVertex.begin(); ilabel !=labelOfVertex.end(); ilabel++){
      PetscSynchronizedPrintf(PETSC_COMM_WORLD," Label #%d :", *ilabel);
      for (int i=0; i<vertexPerLabel[*ilabel].size(); i++) {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD," %d", vertexPerLabel[*ilabel][i]);
      }
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"........ End of vertices per label on Processor #%d  ........\n",rankProc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif
  }
  
  Mesh::~Mesh()
  {
    std::cout << "Delete mesh " << this << std::endl;
    delete [] triangles;
    delete [] vertices;
    delete [] bedges;
    delete [] TheAdjacencesLink;
    delete [] BoundaryEdgeHeadLink;
    delete [] BoundaryAdjacencesHead;
    delete [] BoundaryAdjacencesLink;
  }
  
  void  Mesh::BoundingBox(Point &Pmin,Point &Pmax) const
  {
    assert(nv);
    Pmin=Pmax=vertices[0];
    for (int i=0;i<nv;i++)
    {
      const Point & P=vertices[i];
      Pmin.x = std::min(Pmin.x,P.x);
      Pmin.y = std::min(Pmin.y,P.y);
      Pmax.x = std::max(Pmax.x,P.x);
      Pmax.y = std::max(Pmax.y,P.y);
    }
    std::cout << " Bounding Box = " << Pmin <<" " <<  Pmax << std::endl;

  }
  
  
  void Mesh::read(const char * filename,int rankProc){ // read the mesh
    
    //std::cout << "Reading from file: "<< filename << std::endl;
      
    //------------------------
    // Sizes:
    //------------------------
    
    int InpMsh = -1;
    int meshVersion;
    int dimension;


    BoundaryAdjacencesHead=0;
    BoundaryAdjacencesLink=0;
    BoundaryEdgeHeadLink=0;
    TheAdjacencesLink =0;
    area = 0.;
    
    string name = filename;  

    char* file;
    file = (char*)name.c_str();
    
    InpMsh = GmfOpenMesh(file, GmfRead, &meshVersion, &dimension);
    nv = GmfStatKwd(InpMsh, GmfVertices);//Number of vertices
    nt = GmfStatKwd(InpMsh, GmfTriangles);//Number of triangles
    neb = GmfStatKwd(InpMsh, GmfEdges);//Number of edges

    vertices = new Vertex[nv];     
    triangles = new Triangle [nt];      
    bedges = new BoundaryEdge[neb];      
    
    assert(triangles && vertices && bedges);
    
    //------------------------
    //------------------------
    
    //------------------------
    // Vertices:
    //------------------------
    
    GmfGotoKwd(InpMsh, GmfVertices);
    
    if(meshVersion==1){
      float xcoor=0., ycoor=0., zcoor=0.;
      int readlabel=0;
      for(int iv=0 ; iv<nv ; iv++) {
        if(dimension == 2) GmfGetLin(InpMsh, GmfVertices, &xcoor, &ycoor, &readlabel);
        if(dimension == 3) GmfGetLin(InpMsh, GmfVertices, &xcoor, &ycoor, &zcoor, &readlabel);
        Point P(xcoor,ycoor);
        vertices[iv] = Vertex(P,readlabel);
  
      }
      
    }
    else{ // mesh version = 2 or 3 : double!
      double xcoor=0., ycoor=0., zcoor=0.;
      int readlabel=0;
      for(int iv=0; iv< nv ; iv++) {
        if(dimension == 2) GmfGetLin(InpMsh, GmfVertices, &xcoor, &ycoor, &readlabel);
        if(dimension == 3) GmfGetLin(InpMsh, GmfVertices, &xcoor, &ycoor, &zcoor, &readlabel);
        Point P(xcoor,ycoor);
        vertices[iv] = Vertex(P,readlabel);                    
        }
      
    }
    
    //------------------------
    //------------------------      
    
    //------------------------
    // Triangles:
    //------------------------
    
    GmfGotoKwd(InpMsh, GmfTriangles);
    
    for(int it=0 ; it<nt ; it++){
      
      int v1,v2,v3,readlabel;
      GmfGetLin(InpMsh, GmfTriangles, &v1, &v2, &v3, &readlabel);        

      triangles[it] = Triangle(vertices,v1-1,v2-1,v3-1,readlabel);
      if(vertices[v1-1].ref == vertices[v2-1].ref && vertices[v3-1].ref == vertices[v2-1].ref && vertices[v1-1].ref >= 10){
        triangles[it].ref = vertices[v1-1].ref;
      }


      //TEST NEW MEA
      /*
      if(vertices[v1-1].ref == vertices[v2-1].ref && vertices[v3-1].ref == vertices[v2-1].ref && vertices[v1-1].ref == 5){
        triangles[it].ref = vertices[v1-1].ref;
      }
      */
      area += triangles[it].area;//Compute domain area

    }

    //------------------------
    //------------------------
    
    //------------------------
    // Edges
    //------------------------
    
    GmfGotoKwd(InpMsh, GmfEdges);
    
    for(int ieb=0 ; ieb<neb ; ieb++){
      
      int i1,i2,readlabel;
      GmfGetLin(InpMsh, GmfEdges, &i1, &i2, &readlabel);      
      bedges[ieb] = BoundaryEdge(vertices,i1-1,i2-1,readlabel);
      
    }
    
    //------------------------
    //------------------------
    
    
    //------------------------
    // Close the mesh
    //------------------------
    
    GmfCloseMesh(InpMsh);
    
    //------------------------
    //------------------------


    initVerticesPerLabel(rankProc);
    ConsAdjacence();
    //PetscPrintf(PETSC_COMM_WORLD,"    Mesh area: %f\n", area);
    
  }
  
  
  void Mesh::ConsAdjacence()
  {
    //  -----------
    int NbCollision=0,NbOfEdges=0,NbOfBEdges=0,NbOfMEdges=0;
    const char MaskEdge[]={1,2,4};

    if(TheAdjacencesLink) return; //
    TheAdjacencesLink = new int[3*nt];
    const int NbCode = 2*nv;
    char * TonBoundary = new char [nt]; // the edge is 1 2 4   AddMortar = 8
    
    {
      int * Head = new int[NbCode];
      //  make the list
      int i,j,k,n,j0,j1;
      for ( i=0; i<NbCode; i++ ) Head[i]=-1; // empty list
      n=0; // make all the link
      for (i=0;i<nt;i++)
      {
        Triangle & T=triangles[i];
        for( j=0; j<3; j++,n++ )
        {
          VerticesNumberOfEdge(T,j,j0,j1);
          k = j0+j1;
          TheAdjacencesLink[n]=Head[k];
          Head[k]=n; //
        }
      }
      //
      for (int k=0;k<nt;k++) TonBoundary[k]=0;
      
      BoundaryEdgeHeadLink = new int[neb];
      for (i=0;i<neb;i++)
      {
        BoundaryEdge & be(bedges[i]);
        int n;
        int i0=number(be.vertices[0]);
        int i1=number(be.vertices[1]);
        assert(i0 >=0 && i0 < nv);
        assert(i1 >=0 && i1 < nv);
        int im=std::min(i0,i1);
        BoundaryEdgeHeadLink[i]=-1;
        for ( n=Head[i0+i1]; n>=0; n=TheAdjacencesLink[n])
        {
          int jj=n%3,ii=n/3, jj0,jj1;
          VerticesNumberOfEdge(triangles[ii],jj,jj0,jj1);
          if(im==std::min(jj0,jj1)) // same edge
          {
            TonBoundary[n/3] += MaskEdge[n%3];
            BoundaryEdgeHeadLink[i]=n;
            break;
          }
        }
        if ( BoundaryEdgeHeadLink[i] <0)
          std::cout << "WARNING: boundary edge " << i
          << " is not in the mesh " <<i0 << " " << i1 <<  std::endl;
      }
      
      for (i=0;i<nt;i++)
      {
        Triangle & T=triangles[i];
        for( j=0; j<3; j++,n++ )
        {
          VerticesNumberOfEdge(T,j,j0,j1);
          k = j0+j1; // code of current edge
          int jm = std::min(j0,j1), NbAdj=0, He=-1;
          int *pm=Head+k;
          while (*pm>=0) // be carefull
          {
            int m=*pm,jj=m%3,ii=m/3, jj0,jj1;
            VerticesNumberOfEdge(triangles[ii],jj,jj0,jj1);
            if(jm==std::min(jj0,jj1)) // same edge
            {
              NbAdj++;
              // remove from  the liste
              *pm=TheAdjacencesLink[m];
              TheAdjacencesLink[m]=He;  // link to He
              He = m;
            }
            else
            {
              NbCollision++;
              pm=TheAdjacencesLink+*pm; // next
            }
          }
          //  make a circular link
          if (NbAdj>0)
          {
            int m=He;
            while(TheAdjacencesLink[m]>=0)
              m=TheAdjacencesLink[m]; // find end of list
            // close the List of common edges
            TheAdjacencesLink[m] = He;
          }
          if (NbAdj >2)
          {
            int m=He;
            do {
              m=TheAdjacencesLink[m];
            } while(TheAdjacencesLink[m]!=He);
          }
          
          if (NbAdj) NbOfEdges++;
          if(NbAdj==1){
            if (! (TonBoundary[i]& MaskEdge[j]))
            {
              NbOfMEdges++;
              //TonBoundary[i]+= AddMortar[j];
            }
            else {
              NbOfBEdges++;
            }
          }
        }
      }
     // std::cout << " Nb of edges on Boundary = " << NbOfBEdges << ", neb = " << neb <<  std::endl;
      delete [] Head; // cleanning memory
    }
    
    delete [] TonBoundary;
    /*
    std::cout << " Number of Edges                 = " << NbOfEdges << std::endl;
    std::cout << " Number of MEdges        = " << NbOfMEdges << std::endl;
    std::cout << " Number of Boundary Edges        = " << NbOfBEdges << std::endl;
    std::cout << " Euler Number nt- NbOfEdges + nv = "
    << nt - NbOfEdges + nv << " = Nb of Connected Componant - Nb Of Hole "
    <<std::endl;
    */
  }
  
  std::ostream& operator <<(std::ostream& f, const Mesh & m )
  {
    f << "\n========= Mesh ===========" << std::endl;
    f << "Vertices:" << std::endl;
    for (int i=0;i<m.nv;i++)
      f << m.vertices[i]  << std::endl;
    f << "Triangles:" << std::endl;
    for (int i=0;i<m.nt;i++){
      Triangle& t = m.triangles[i];
      f << " Vertices:" << m(t[0])  << ' ' << m(t[1]) << ' ' << m(t[2]) << " Label:" << (Label) t <<std::endl;
    }
    f << "Boundary edges:" << std::endl;
    for (int i=0;i<m.neb; i++){
      BoundaryEdge& be = m.bedges[i];
      f << " Vertices:" << m(be[0]) << ' ' << m(be[1]) <<" Label:" << (Label) be <<std::endl;
    }
    f << "....... End of mesh ......." << std::endl;
    return f;
    
  }

  void Mesh::print(int rankProc) const{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n======== Print mesh on Processor #%d  ========\n",rankProc);
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Vertices:\n");
    //   for (int i=0;i<m.nv;i++) f << m.vertices[i]  << std::endl;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d triangles:\n",nt);
    for (int i=0;i<nt;i++){
      Triangle& t = triangles[i];
      PetscSynchronizedPrintf(PETSC_COMM_WORLD," Vertices: %d %d %d, Label: %d \n",operator()(t[0]),operator()(t[1]),operator()(t[2]),t.ref);
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"........ End of mesh on Processor #%d  .........\n",rankProc);

#ifdef PETSC_3_6
    PetscSynchronizedFlush(PETSC_COMM_WORLD,NULL);
#else
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif

   /*
    f << "Boundary edges:" << std::endl;
    for (int i=0;i<m.neb; i++){
      BoundaryEdge& be = m.bedges[i];
      f << " Vertices:" << m(be[0]) << ' ' << m(be[1]) <<" Label:" << (Label) be <<std::endl;
    }
    f << "======= End of mesh =======" << std::endl;
    return f;
    */
  }

}
