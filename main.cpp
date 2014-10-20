#include <iostream>
#include <cassert>
#include <vector>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include "data.h"
#include "matrix.h"
#include "readpng.h"
#include "trimesh.h"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

#ifdef MACOSX
#include "OpenGL/gl.h"
#include "OpenGl/glu.h"
#include "GLUT/glut.h"
#else
#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glut.h"
#endif

#define RAD_TO_DEG 180/M_PI

using namespace std;

static int lastX, lastY;
static int downX, downY;

// to remove: vec.erase(vec.begin() + index);
int xRes, yRes;

int userTransformType;

int drawingType;

Eigen::SparseMatrix<double> eigenD0;
Eigen::SparseMatrix<double> eigenD1;
Eigen::SparseMatrix<double> eigenH0;
Eigen::SparseMatrix<double> eigenH1;
Eigen::SparseMatrix<double> eigenH2;
Eigen::SparseMatrix<double> Delta;
Eigen::SparseMatrix<double> DeltaInterior;
vector <Vertex> PrimalVertices;
vector <Edge> PrimalEdges;
vector <Triangle> PrimalTriangles;

vector <Vertex> DualVertices;
vector <DualEdge> DualEdges;
vector <vector <Vertex> > DualFaces;
vector <vector <int> > DualFacesIndices;

vector <Vertex> Particles;

int currentSelection = -1;
int currentDualFace = -1;
int gDebug = 0;
int gFrame = 0;
int gRun = 0;
int gDrawVelocity = 0;
int gDrawMesh = 0;
int gDrawDualMesh = 0;
int gDrawBacktracked = 0;
int gMeshWidth = 0;
int gMeshHeight = 0;
GLuint gTexNames[4];
Vertex gLoc1(0.15, 0);
Vertex gLoc2(-0.15, 0);
double gSumVorticity = 0.0;

/** PROTOTYPES **/
void redraw();
void idle();
void initGL();
void resize(GLint w, GLint h);
void keyfunc(GLubyte key, GLint x, GLint y);
void updateEverything();
void setTexture(const std::string &fileName, int textureNumber);

// given a value of the pixel that we clicked on, return the value in
// object space (which is on a scale from -1 to 1). 
float getNDCCoord(int pixel, int min, int max, int resolution) 
{

    float smallScale = max - min;
    float result = smallScale * (pixel - resolution/smallScale) / resolution;

    //std::cout << "getNDCCoord " << result << std::endl;
    return result;
}

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

bool isLeft(Vec3 a, Vec3 b, Vec3 c){
    return (((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)) >= 0);
}

int findClosestDual(double x, double y) {

    size_t closest;
    double closestDistance = 1000;
    for (size_t i = 0; i < PrimalVertices.size(); i++) {
        Vertex &vertex = PrimalVertices[i];
        double x_squared = (x - vertex.x) * (x - vertex.x);
        double y_sqared = (y - vertex.y) * (y - vertex.y);
        double distance = sqrt(x_squared + y_sqared);
        if (distance < closestDistance) {
            closestDistance = distance;
            closest = i;
        }
    }
    return (int)closest;
}

int findTriangle(double x, double y) 
{
    Vertex point(x,y);
    for (size_t i = 0; i < PrimalTriangles.size(); i++) {
        const Triangle &triangle= PrimalTriangles[i];
        if (triangle.isPointInside(point)) {
            return (int)i;
        }
    }
    std::cout << "COULD NOT FIND" << std::endl;
    return -1;
}

int findClosestTriangle(double x, double y) {

    double closest = 1000;
    double closestTriangle;
    for (size_t i = 0; i < PrimalTriangles.size(); i++) {
        double xSq = (PrimalTriangles[i].x - x) * (PrimalTriangles[i].x - x);
        double ySq = (PrimalTriangles[i].y - y) * (PrimalTriangles[i].y - y);
        double d = sqrt(xSq + ySq);
        if (d < closest) {
            closest = d;
            closestTriangle = i;
        }
    }
    return int(closestTriangle);
}

void drawDualMesh() {
    glColor3d(1,0,1);
    for (size_t i = 0; i < DualFaces.size(); i++) {
        if (i != 1399) {
            //continue;
        }
        vector <Vertex> &faceVertices = DualFaces[i];
        //glColor3d(0,1,0);
        glBegin(GL_LINE_LOOP);
        for (size_t j = 0; j < faceVertices.size(); j++) {
            Vertex &vertex = faceVertices[j];
            glVertex2d(vertex.x, vertex.y);
        }
        glEnd();
    }
}

void drawBacktracked() {
    glColor3d(1,1,0);
    for (size_t i = 0; i < PrimalVertices.size(); i++) {
        vector <Vertex> &vertices = DualFaces[i];
        vector <int> &triangles = DualFacesIndices[i];
        if (gDebug) std::cout << "vorticity here: " << PrimalVertices[i].vorticity << std::endl;
        if (gDebug) std::cout << "Vertices: " << std::endl;
            
        if (PrimalVertices[i].boundary == 1) {
            continue;
        }
            
        glBegin(GL_LINE_LOOP);
            
        // circ += (1/2)(vi + vj)(ci + cj)
        vector <Vec3> backtrackedPoints;
        vector <Vec3> backtrackedVelocities;
        //std::cout << "dual size: " << vertices.size() << std::endl;
        for (size_t j = 0; j < vertices.size(); j++) {
             //std::cout << vertices[j].x << " " << vertices[j].y << std::endl;
             int tri = vertices[j].triIndex;
             Triangle &t = PrimalTriangles[tri];
             //std::cout << "triangle: " << tri << std::endl;
            
             Vec3 backPoint;
             Vec3 backVelocity;
             double backTrackScale = 100.0;
             backPoint.x = vertices[j].x - t.vx/backTrackScale;
             backPoint.y = vertices[j].y - t.vy/backTrackScale;
             glVertex2d(vertices[j].x - t.vx/backTrackScale, vertices[j].y - t.vy/backTrackScale);

         }
         glEnd();

    }
}

double signedArea(Vertex &v1, Vertex &v2, Vertex &v3) {

    double term1 = (v1.x * v2.y) - (v1.x * v3.y); // x1*y2 - x1*y3
    double term2 = (v3.x * v1.y) - (v2.x * v1.y); // x3*y1 - x2*y1
    double term3 = (v2.x * v3.y) - (v3.x * v2.y); // x2*y3 - x3*y2

    double area = 0.5 * (term1 + term2 + term3);
    return area;
}

double distance(double x1, double y1, double x2, double y2) {
    double xSq = (x2 - x1) * (x2 - x1);
    double ySq = (y2 - y1) * (y2 - y1);
    double d = sqrt(xSq + ySq);
    return d;
}

#if 0
Vec3 interpolateVelocity(double x, double y) {
    //std::cout << "int " << x << " " << y << std::endl;
    int currentDualFace = findClosestDual(x, y);
    //std::cout << "Interpolate " << currentDualFace << std::endl;
    vector <int> &triangles = DualFacesIndices[currentDualFace];
    double totalDistance = 0;
    vector <double> inverseDistances;
    vector <Vec3> velocities;
    for (size_t i = 0; i < triangles.size(); i++) {
        //std::cout << triangles[i] << std::endl;
        Triangle &t = PrimalTriangles[triangles[i]];
        //std::cout << t.x << " " << t.y << std::endl;
        double d = distance(x, y, t.x, t.y); // get distance to each dual vertex
        //std::cout << "distance: " << d << std::endl;
        totalDistance += d;
        inverseDistances.push_back(d);
        Vec3 v(t.vx, t.vy, 0);
        velocities.push_back(v);
    }
    //std::cout << "TOTAL: " << totalDistance << std::endl;
    double total_x = 0;
    double total_y = 0;

    double totalInverse = 0;
    for (size_t i = 0; i < inverseDistances.size(); i++) {
        double d = inverseDistances[i];
        if (d == 0) {
            continue;
        }
        //std::cout << "D: " << d << std::endl;
        double scale = (totalDistance / d);
        totalInverse += scale;
    }
    for (size_t i = 0; i < inverseDistances.size(); i++) {
        double d = inverseDistances[i];
        if (d == 0) {
            continue;
        }
        Vec3 &velocity = velocities[i];
        double scale = (totalDistance / d) / totalInverse;
        //std::cout << "sc: " << scale << std::endl;
        total_x += scale * velocity.x;
        total_y += scale * velocity.y;
    } 

    Vec3 a(total_x, total_y, 0);
    return a;
}

#else
Vec3 interpolateVelocity(double x, double y) {
    //std::cout << "int " << x << " " << y << std::endl;
    int currentDualFace = findClosestDual(x, y);
    //std::cout << "Interpolate " << currentDualFace << std::endl;
    vector <int> &triangles = DualFacesIndices[currentDualFace];
    double totalDistance = 0;
    vector <double> inverseDistances;
    vector <Vec3> velocities;
    double sumVx = 0.0;
    double sumVy = 0.0;
    for (size_t i = 0; i < triangles.size(); i++) {
        //std::cout << triangles[i] << std::endl;
        Triangle &t = PrimalTriangles[triangles[i]];

        sumVx += t.vx;
        sumVy += t.vy;
    }

    sumVx /= triangles.size();
    sumVy /= triangles.size();

    Vec3 a(sumVx, sumVy, 0);
    return a;
}
#endif

void backtrackDualVertices(int n) {

    for (size_t j = 0; j < PrimalTriangles.size(); j++) {
        Triangle &t = PrimalTriangles[j];
        PrimalTriangles[j].backtrackX = t.x; // reset first before calculating
        PrimalTriangles[j].backtrackY = t.y;
        }

    for (int i = 0; i < n; i++) {
        for (size_t j = 0; j < PrimalTriangles.size(); j++) {
            Triangle &t = PrimalTriangles[j];
            PrimalTriangles[j].backtrackX -= t.vx/10.0;
            PrimalTriangles[j].backtrackY -= t.vy/10.0;
            
            /*
            if (j == 0) {
                interpolateVelocity(PrimalTriangles[j].backtrackX, PrimalTriangles[j].backtrackY);
                std::cout << PrimalTriangles[j].backtrackX << " " << PrimalTriangles[j].backtrackY << std::endl;
            }
            */
            
        }
    }
}

void drawBacktrackedVertices() {
    glColor3d(0.5,0,0.55);
    glPointSize(3);
    glBegin(GL_POINTS);
    for (size_t j = 0; j < PrimalTriangles.size(); j++) {
        //Triangle &t = PrimalTriangles[j];
        //PrimalTriangles[j].backtrackX -= t.vx;
        //PrimalTriangles[j].backtrackY -= t.vy;
        glVertex2d(PrimalTriangles[j].backtrackX, PrimalTriangles[j].backtrackY);
    }
    glEnd();
}

void calculateVorticity() {
    if (gDebug) std::cout << "calculateVorticity start" << std::endl;
    for (size_t i = 0; i < PrimalVertices.size(); i++) {
        vector <Vertex> &vertices = DualFaces[i];
        vector <int> &triangles = DualFacesIndices[i];
        if (gDebug) std::cout << "vorticity here: " << PrimalVertices[i].vorticity << std::endl;
        if (gDebug) std::cout << "Vertices: " << std::endl;
            
        if (PrimalVertices[i].boundary == 1) {
            continue;
        }
            
        //glColor3d(0,1,0);
        //glBegin(GL_LINE_LOOP);
            
        // circ += (1/2)(vi + vj)(ci + cj)
        vector <Vec3> backtrackedPoints;
        vector <Vec3> backtrackedVelocities;
        //std::cout << "dual size: " << vertices.size() << std::endl;
        for (size_t j = 0; j < vertices.size(); j++) {
             //std::cout << vertices[j].x << " " << vertices[j].y << std::endl;
             int tri = vertices[j].triIndex;
             Triangle &t = PrimalTriangles[tri];
             //std::cout << "triangle: " << tri << std::endl;
            
             Vec3 backPoint;
             Vec3 backVelocity;
             double backTrackScale = 100.0;
             backPoint.x = vertices[j].x - t.vx/backTrackScale;
             backPoint.y = vertices[j].y - t.vy/backTrackScale;
             //glVertex2d(vertices[j].x - t.vx/backTrackScale, vertices[j].y - t.vy/backTrackScale);
             //std::cout << "interpolating..." << t.vx << " " << t.vy << " " << backPoint.x << " " << backPoint.y << std::endl;
             Vec3 velocity = interpolateVelocity(backPoint.x, backPoint.y);
             //std::cout << velocity.x << " " << velocity.y << std::endl;
             backVelocity.x = velocity.x;//t.vx;
             backVelocity.y = velocity.y;//t.vy;
             //backVelocity.x = t.vx;
             //backVelocity.y = t.vy;
             backtrackedPoints.push_back(backPoint);
             backtrackedVelocities.push_back(backVelocity);
         }
         //glEnd();

         // add back the first vertex to the end of the vector
         Vec3 firstPoint = backtrackedPoints[0];
         Vec3 firstVelocity = backtrackedVelocities[0];
         backtrackedPoints.push_back(firstPoint);
         backtrackedVelocities.push_back(firstVelocity);

         double circulation = 0;
         for (size_t k = 1; k < backtrackedPoints.size(); k++) {
             //std::cout << k << " " << circulation << std::endl;
             Vec3 c0 = backtrackedPoints[k-1];
             Vec3 c1 = backtrackedPoints[k];
             Vec3 v0 = backtrackedVelocities[k-1];
             Vec3 v1 = backtrackedVelocities[k];

             double avg_vx = 0.5 * (v0.x + v1.x);
             double avg_vy = 0.5 * (v0.y + v1.y);
             //std::cout << avg_vx << " " << avg_vy << std::endl;
             double diff_cx = c0.x - c1.x;
             double diff_cy = c0.y - c1.y;
             double dot = (avg_vx * diff_cx) + (avg_vy * diff_cy);
             circulation += dot;
         }
         if (gDebug) std::cout << "Final circulation: " << -circulation << " for vertex " << i << std::endl;
         PrimalVertices[i].vorticity = -circulation;
         /*
         for (size_t j = 0; j < triangles.size(); j++) {
             std::cout << triangles[j] << std::endl;
             Triangle t = PrimalTriangles[triangles[j]];
             glVertex2d(t.backtrackX, t.backtrackY);
         }
         */
         //glEnd();
    }

    double totalW = 0.0;
    for (size_t i = 0; i < PrimalVertices.size(); i++) {
        if (PrimalVertices[i].boundary != 1) {
            totalW += PrimalVertices[i].vorticity;
        }
    }

    double divisor = DeltaInterior.rows();
    for (size_t i = 0; i < PrimalVertices.size(); i++) {
        if (PrimalVertices[i].boundary == 1) {
            continue;
        }
        PrimalVertices[i].vorticity -= (totalW - gSumVorticity)/divisor;
    }

    double fullW = 0;
    totalW = 0; // do another totalW summation
    for (int i = 0; i < PrimalVertices.size(); i++) {
        fullW += PrimalVertices[i].vorticity;
        if (PrimalVertices[i].boundary != 1) {
            totalW += PrimalVertices[i].vorticity;
        }
    }
    if (gDebug) std::cout << "Total w: " << totalW << ", All w: " << fullW << std::endl;
} 

// updates the particles for 1 time step
void moveParticles() {

    double scale = 50.0;
    for (size_t i = 0; i < Particles.size(); i++) {
        Vertex &vertex = Particles[i];
        Triangle &t = PrimalTriangles[vertex.triIndex];

        //std::cout << "t: " << t.index << " " << t.x << " " << t.y << std::endl;
        Vec3 v = interpolateVelocity(vertex.x, vertex.y);
        //std::cout << "Vel: " << v.x << " " << v.y << std::endl;
        
        Vertex newVertex(vertex.x + v.x/scale, vertex.y + v.y/scale);
        
        int tri = findTriangle(newVertex.x, newVertex.y);
        if (tri == -1) {
            /* could potential just restart particles at origin.. */
            newVertex.x = 0;
            newVertex.y = 0;
            tri = findTriangle(newVertex.x, newVertex.y);
        }
        vertex.x = newVertex.x;
        vertex.y = newVertex.y;
        vertex.triIndex = tri;
    }
}

void drawParticles() {
    glEnable(GL_BLEND);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_POINT_SPRITE_ARB);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_CONSTANT_COLOR);
    glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glDepthMask(GL_FALSE);

    int halfParticles = Particles.size() / 2;
    glColor4d(0,0,0.75, 0.5);
    glPointSize(20);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < Particles.size(); i++) {
        if (i == halfParticles) {
            glColor4d(0.75, 0,0,0.5);
        }
        glBindTexture(GL_TEXTURE_2D, gTexNames[0]);
        glVertex2d(Particles[i].x, Particles[i].y);
    }
    glEnd();
    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_POINT_SPRITE_ARB);
    glDepthMask(GL_TRUE);
}

void calculateVelocity() {

    for (size_t i = 0; i < PrimalTriangles.size(); i++) {
        Triangle &t = PrimalTriangles[i];
        Edge &edge1 = PrimalEdges[t.e1-1];
        Edge &edge2 = PrimalEdges[t.e2-1];
        Edge &edge3 = PrimalEdges[t.e3-1];
        double total = 0;
        if (edge1.start == t.v1) {
            total += edge1.flux;
        } else {
            total -= edge1.flux;
        }
        if (edge2.start == t.v2) {
            total += edge2.flux;
        } else {
            total -= edge2.flux;
        }
        if (edge3.start == t.v3) {
            total += edge3.flux;
        } else {
            total -= edge3.flux;
        }

        //std::cout << edge1.start << " " << edge1.end << std::endl;
        //std::cout << "Tot: " << total << "  indices: " << edge1.index << " " << edge2.index << " " << edge3.index << std::endl;

        Vertex &vert1a = PrimalVertices[edge1.start-1];
        Vertex &vert1b = PrimalVertices[edge1.end-1];
        double e1x = vert1b.x - vert1a.x;
        double e1y = vert1b.y - vert1a.y;

        //std::cout << "e1: " << e1x << " " << e1y << std::endl;
        Vertex &vert2a = PrimalVertices[edge2.start-1];
        Vertex &vert2b = PrimalVertices[edge2.end-1];
        double e2x = vert2b.x - vert2a.x;
        double e2y = vert2b.y - vert2a.y;

        //std::cout << "e2: " << e2x << " " << e2y << std::endl;

        Vertex &vert3a = PrimalVertices[edge3.start-1];
        Vertex &vert3b = PrimalVertices[edge3.end-1];
        double e3x = vert3b.x - vert3a.x;
        double e3y = vert3b.y - vert3a.y;

        //std::cout << "e3: " << e3x << " " << e3y << std::endl;

        double f1 = edge1.flux;
        double f2 = edge2.flux;
        double f3 = edge3.flux;

        //f1 *= -1;

        double vy = (f2 - ((e2y/e1y) * f1)) / (((e1x*e2y) / e1y) - e2x);
        double vx = (f1 + vy*e1x) / e1y;

        if (e1y == 0) { // avoid dividing by 0, use different edge
            vy = (f3 - ((e3y/e2y) * f2)) / (((e2x*e3y) / e2y) - e3x);
            vx = (f2 + vy*e2x) / e2y;
        }

        PrimalTriangles[i].vx = vx; // update the triangles velocity components for backtracking
        PrimalTriangles[i].vy = vy;

        //std::cout << i << " " <<
        //std::cout << "vx: " << vx << "  vy: " << vy << std::endl;

        // flux = vx*ey - vy*ex
        //double flux1 = (vx * e1y) - (vy * e1x);
        //double flux2 = (vx * e2y) - (vy * e2x);
        //double flux3 = (vx * e3y) - (vy * e3x);
        //std::cout << "flux1: " << flux1  << "   " << edge1.flux << std::endl;
        //std::cout << "flux2: " << flux2 << "   " << edge2.flux << std::endl;
        //std::cout << "flux3: " << flux3 << "   " << edge3.flux << std::endl;
    }
}

void drawFluxes() {

    glColor3d(0,1.0,0);
    for (size_t i = 0; i < PrimalTriangles.size(); i++) {
        Triangle &t = PrimalTriangles[i];
        double vx = PrimalTriangles[i].vx;
        double vy = PrimalTriangles[i].vy;
        glBegin(GL_LINES);  // drawing triangle meshes
        glVertex2d(t.x, t.y);
        glVertex2d(t.x + vx/10.0, t.y + vy/10.0);
        glEnd();
    }
}

/* Draws the current object, based on gMaterial propereries and the current
 * gSeperator block.
*/
void drawObject() 
{

    glPushMatrix();  // save the current matrix
    
    glDisable(GL_LIGHTING);

    if (gDrawMesh) {
        glPointSize(2);
        glColor3d (1,0,0);
        glBegin(GL_POINTS);  // drawing triangle meshes
        for (size_t i = 0; i < PrimalVertices.size(); i++) {
            const Vertex &vert = PrimalVertices[i];
            glVertex2d(vert.x, vert.y);
        }

        glColor3d(0,0,1);
        for (size_t i = 0; i < PrimalTriangles.size(); i++) {
            const Triangle &triangle = PrimalTriangles[i];
            glVertex2d(triangle.x, triangle.y);
        }
        glEnd();

        glColor3d (1,1,1); // make lines white
        for (size_t i = 0; i < PrimalEdges.size(); i++) {
            const Edge &edge = PrimalEdges[i];
            const Vertex &start = PrimalVertices[edge.start-1];
            const Vertex &end = PrimalVertices[edge.end-1];
            glBegin(GL_LINES);  // drawing triangle meshes
            glVertex2d(start.x, start.y);
            glVertex2d(end.x, end.y);
            //glVertex2d(vertC.x, vertC.y);
            glEnd();
        }
    }

    if (gDrawDualMesh) {
        drawDualMesh();
    }

    if (gDrawVelocity) {
        drawFluxes();
    }

    if (gDrawBacktracked) {
        drawBacktracked();
    }

    drawParticles();
   
    if (currentSelection != -1) { // draw selected triangle (if any)
        Triangle &triangle = PrimalTriangles[currentSelection];

        const Vertex &vertA = triangle.mVert1;
        const Vertex &vertB = triangle.mVert2;
        const Vertex &vertC = triangle.mVert3;
        glColor3d(1,1,0); 
        glBegin(GL_TRIANGLES);  // drawing triangle meshes
        glVertex2d(vertA.x, vertA.y);
        glVertex2d(vertB.x, vertB.y);
        glVertex2d(vertC.x, vertC.y);
        glEnd();

        if (currentDualFace != -1) {
            const vector <Vertex> &faceIndices = DualFaces[currentDualFace];
            glColor3d(0,1,0);
            glBegin(GL_POLYGON);
            for (size_t j = 0; j < faceIndices.size(); j++) {
                const Vertex &vertex = faceIndices[j];
                glVertex2d(vertex.x, vertex.y);
            }
            glEnd();
        }        
    }
    
    glEnable(GL_LIGHTING);
    
    glPopMatrix();  // restore old matrix.

}

struct keeper {
    keeper() {};
    // return true to keep entry, false to remove entry
    bool operator() (const int& row, 
                     const int& col, const double& value) const 
    { 
      int index = row * gMeshWidth + col;
      Vertex &v1 = PrimalVertices[row];
      Vertex &v2 = PrimalVertices[col];
      bool isV1Boundary = (v1.boundary == 1);
      bool isV2Boundary = (v2.boundary == 1);
      bool keep = !(isV1Boundary || isV2Boundary);
      if (gDebug) {
        std::cout << "keep: " << keep << " " << row << " " << col << std::endl;
      }

      return keep;
    }
};

void setupDeltaInterior() {
    if (gDebug) std::cout << "setupDeltaInterior start" << std::endl;
    if (gDebug) {
        std::cout << "Delta: " << Delta.rows() << std::endl;
        std::cout << Delta << std::endl;
    }

    Delta.prune(keeper());

    int interiorSize = (gMeshWidth - 2) * (gMeshHeight - 2);
    Eigen::SparseMatrix<double> subMat(interiorSize, interiorSize);

    for (int k=0; k< Delta.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(Delta,k); it; ++it)
    {
        int newRow = it.row()-(gMeshWidth - 1) - 2 * (int)(it.row()/(gMeshWidth));
        int newCol = it.col()-(gMeshWidth - 1) - 2 * (int)(it.col()/(gMeshHeight));

        if (newRow < 0 || newCol < 0) {
            continue;
        }

        subMat.insert(newRow, newCol) = it.value();
    }

    if (gDebug) {
        std::cout << "subMat size " << subMat.rows() << " " << subMat.cols() << std::endl;
        std::cout << subMat << std::endl;
    }

    DeltaInterior = subMat;
    if (gDebug) std::cout << "setupDeltaInterior end " << std::endl;
}

/* Updates flux values on each edge, given vorticities */
void calculateFlux() {

    if (gDebug) std::cout << "calculateFlux start" << std::endl;
    int size = PrimalVertices.size();
    int counter = DeltaInterior.rows();

    //int n = PrimalVertices.size();
    Eigen::VectorXd X(counter), B(counter);
    int counterFill = 0;
    for (int i = 0; i < size; i++) {
        if (PrimalVertices[i].boundary == 0) {
            B[counterFill] = PrimalVertices[i].vorticity;
            if (gDebug) std::cout << i << " " << PrimalVertices[i].vorticity << std::endl;
            counterFill++;
        }
    }
    if (gDebug) {
        std::cout << "B" << std::endl;
        std::cout << B << std::endl;
    }

    
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
    cg.compute(DeltaInterior);
    X = cg.solve(B);
    if (gDebug) {
        std::cout << "X: " << std::endl;
        std::cout << X << std::endl;
    }

    Eigen::SparseMatrix<double> Theta(PrimalVertices.size(),1);
    int counterInsert = 0;
    for (int i = 0; i < PrimalVertices.size(); i++) {
        //std::cout << i << std::endl;
        if (PrimalVertices[i].boundary == 1) {
            Theta.insert(i,0) = 0;
        } else {
            Theta.insert(i,0) = X[counterInsert];
            counterInsert++;
        }
    }

    if (gDebug) {
        std::cout << "Theta: " << std::endl;
        std::cout << Theta << std::endl;
    }

    Eigen::SparseMatrix<double> Fluxes = eigenD0.transpose() * Theta;
    /*
    std::cout << "Fluxes: " << std::endl;
    std::cout << Fluxes << std::endl;
    */
    
    //zero all out first
    
    for (int i = 0; i < PrimalEdges.size(); i++) {
        PrimalEdges[i].flux = 0;
    }

    double sumFluxes = 0.0;
    int divisor = 0;
    for (int k=0; k<Fluxes.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Fluxes,k); it; ++it){
            sumFluxes += it.value();
            divisor++;
        }
    }

    for (int k=0; k<Fluxes.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Fluxes,k); it; ++it){
            PrimalEdges[it.index()].flux = it.value() - sumFluxes/divisor;
        }
    }
    
    /*
    for (int k=0; k<Fluxes.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Fluxes,k); it; ++it){
            PrimalEdges[it.index()].flux = it.value();
        }
    }
    */

    // another sum after adjustment
    sumFluxes = 0.0;
    for (int k=0; k<Fluxes.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Fluxes,k); it; ++it){
            sumFluxes += PrimalEdges[it.index()].flux;
        }
    }
    if (gDebug) std::cout << "sumFluxes after " << sumFluxes << std::endl;
}


/** GLUT callback functions **/

/*
 * This function gets called every time the window needs to be updated
 * i.e. after being hidden by another window and brought back into view,
 * or when the window's been resized.
 * You should never call this directly, but use glutPostRedisply() to tell
 * GLUT to do it instead.
 */
void redraw()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    drawObject();
    glutSwapBuffers();
}

void idle()
{
    if (gRun) {
        updateEverything();
        glutPostRedisplay();
    }
}

/**
 * GLUT calls this function when the window is resized.
 * All we do here is change the OpenGL viewport so it will always draw in the
 * largest square that can fit in our window..
 */
void resize(GLint w, GLint h)
{
    if (h == 0)
        h = 1;

    // ensure that we are always square (even if whole window not used)
    if (w > h)
        w = h;
    else
        h = w;

    // Reset the current viewport and perspective transformation
    glViewport(0, 0, w, h);
    xRes = w;
    yRes = h;

    // Tell GLUT to call redraw()
    glutPostRedisplay();
}

/*
 * GLUT calls this function when any key is pressed while our window has
 * focus.  Here, we just quit if any appropriate key is pressed.  You can
 * do a lot more cool stuff with this here.
 */
void keyfunc(GLubyte key, GLint x, GLint y)
{
    if (key == 'f') {
        glShadeModel(GL_FLAT);  // set to Flat shading
        drawingType = 1;  // set up to draw actual triangles
        glEnable(GL_LIGHTING);
        glutPostRedisplay();
    } else if (key == 'g') {
        glShadeModel(GL_SMOOTH);  // set to Gouraud shading
        drawingType = 1; // set up to draw actual triangles
        glEnable(GL_LIGHTING);
        glutPostRedisplay();
    } else if (key == 'w') {
        drawingType = 0;  // set up for wireframe
        glDisable(GL_LIGHTING); // turn off lighting for wireframe
        glutPostRedisplay();

    } else if (key == 27 || key == 'q' || key =='Q') {
        exit(0); // escape or q or Q to exit the program
    } else if (key == 'm') {
        updateEverything();
        glutPostRedisplay();
    } else if (key == 't') {
        if (gDrawMesh == 0) {
            gDrawMesh = 1;
        } else {
            gDrawMesh = 0;
        }
        glutPostRedisplay();
    } else if (key == 'd') {
        if (gDrawDualMesh == 0) {
            gDrawDualMesh = 1;
        } else {
            gDrawDualMesh = 0;
        }
        glutPostRedisplay();
    } else if (key == 'b') {
        if (gDrawBacktracked == 0) {
            gDrawBacktracked = 1;
        } else {
            gDrawBacktracked = 0;
        }
        glutPostRedisplay();
    } else if (key == 'v') {
        if (gDrawVelocity == 0) {
            gDrawVelocity = 1;
        } else {
            gDrawVelocity = 0;
        }
        glutPostRedisplay();
    }
}

void updateEverything() {

    std::cout << "Frame: " << ++gFrame << std::endl;
    calculateVorticity();
    calculateFlux();
    calculateVelocity();
    moveParticles();
    /*
    for (size_t i = 0; i < PrimalTriangles.size(); i++) {
        Triangle t = PrimalTriangles[i];
        std::cout << "Velocity: " << i << "   " << t.vx << " " << t.vy << std::endl;
    }*/
    glutPostRedisplay();
    if (gDebug) std::cout << "\n" << std::endl;
}

// call this to see if there's any GL error, if there is,
// print out the error
void checkGLError()
{
    GLenum error = glGetError();
    if (error != GL_NO_ERROR) {
        std::cerr << "glError " << gluErrorString(error) << std::endl;
    }
}

/**
 * Set up OpenGL state.  This does everything so when we draw we only need to
 * actually draw the sphere, and OpenGL remembers all of our other settings.
 */
void initGL()
{
    // Tell openGL to use gouraud shading:
    glShadeModel(GL_SMOOTH);
    drawingType = 1;  // start at Gouraud setting
    
    // Enable back-face culling:
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Enable depth-buffer test.
    glEnable(GL_DEPTH_TEST);
    
    // Set up projection and modelview matrices ("camera" settings) 
    // Look up these functions to see what they're doing.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    checkGLError();

    // this is the only operation we are doing on GL_PROJECTION matrix

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glGenTextures(4, gTexNames);
    setTexture("particle.png", 0);
}

static void mouseButton(int button, int state, int x, int y) 
{
    double xNDC = getNDCCoord(x, -1, 1, xRes);
    double yNDC = getNDCCoord(y, -1, 1, yRes);
    std::cout << "Clicked: " << xNDC << " " << -yNDC << std::endl;
    if (state == GLUT_DOWN) {
        currentSelection = findTriangle(xNDC, -yNDC);
        currentDualFace = findClosestDual(xNDC, -yNDC);
        std::cout << "currentSelection: " << currentSelection << std::endl;
        std::cout << "dualFace: " << currentDualFace << std::endl;
        std::cout << "vorticity = " << PrimalVertices[currentDualFace].vorticity << std::endl;
        glutPostRedisplay();
    } else {
        currentSelection = -1;
        currentDualFace = -1;
        glutPostRedisplay();
    }
    switch(button) {
        case GLUT_LEFT_BUTTON:
            if (state == GLUT_DOWN) {
                lastX = x;  // note where the drag began
                lastY = y;
                if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
                    userTransformType = 1;  // "1" indicates ZOOM! // TODO: MOVE TO MIDDLE MOUSE         
                } else {
                    userTransformType = 0; // "0" indicates TRANSLATION! // TODO: MOVE TO MIDDLE MOUSE
                }
            } else if (state == GLUT_UP) {
                userTransformType = -1;  // reset
            }
            break;
        case GLUT_MIDDLE_BUTTON:
            if (state == GLUT_DOWN) {
                lastX = x;
                lastY = y;
            } else if (state == GLUT_UP) {
                userTransformType = -1; // reset
            }
            break;
        case GLUT_RIGHT_BUTTON:
            if(state == GLUT_DOWN) {
                lastX = x;
                lastY = y; 
                userTransformType = 2;
            } else if (state == GLUT_UP) {
                userTransformType = -1; // reset
            }
            break;
    }
    downX = x;
    downY = y;
}

void sortVertices2(vector <Vertex> &dualFace) { 
    double min_y = 1.0;
    int minVertex;
    for (size_t i = 0; i < dualFace.size(); i++) { // find lowest y coordinate
        if (dualFace[i].y < min_y) {
            min_y = dualFace[i].y;
            minVertex = i;
        }
    }

    Vertex minimumVertex = dualFace[minVertex];

    for (size_t i = 0; i < dualFace.size(); i++) {
        Vertex &v = dualFace[i];
        double x = v.x - minimumVertex.x;
        double y = v.y - minimumVertex.y;
        double result = atan2 (y,x);
        v.z = result;
    }

    while (true) {
        bool done = false;
        if (dualFace.size() == 0) {
            break;
        }
        for (size_t j = 0; j < (dualFace.size() - 1); j++) {
            if (dualFace[j].z > dualFace[j+1].z) {
                Vertex temp = dualFace[j+1];
                dualFace[j+1] = dualFace[j];
                dualFace[j] = temp;
                done = true;
            }
        }
        if (done == false) {
            break;
        }
    }

    // correction for boundary dual faces, where arctan might be the same for 2 
    // dual vertices.
    for (size_t j = 0; j < (dualFace.size() - 2); j++) {
        double x = dualFace[j].x;
        double y = dualFace[j].y;
        if (dualFace[j+1].x == x && dualFace[j+2].x == x && x > 0) {
            if (dualFace[j+1].y > dualFace[j+2].y) {
                Vertex temp = dualFace[j+2];
                dualFace[j+2] = dualFace[j+1];
                dualFace[j+1] = temp;
            }
        }
        if (dualFace[j+1].y == y && dualFace[j+2].y == y && y < 0) {
            if (dualFace[j+1].x > dualFace[j+2].x) {
                Vertex temp = dualFace[j+2];
                dualFace[j+2] = dualFace[j+1];
                dualFace[j+1] = temp;
            }    
        }
    }

}      

double polygonArea(vector <double> &X, vector <double> &Y, int numPoints) {
    double area = 0;         // Accumulates area in the loop
    double j = numPoints-1;  // The last vertex is the 'previous' one to the first

    for (size_t i=0; i<numPoints; i++) {
        area = area + (X[j]+X[i]) * (Y[j]-Y[i]); 
        j = i;
    }
    return area/2;
}


void setupEigenD0() {

    int numVertices = PrimalVertices.size();
    int numEdges = PrimalEdges.size();
    eigenD0.resize(numVertices,numEdges);
    for (size_t i = 0; i < PrimalEdges.size(); i++) {
        Edge edge = PrimalEdges[i];
        //std::cout << "Edge: " << edge.x << " " << edge.y << std::endl;
        eigenD0.insert(edge.start-1,i) = -1;
        eigenD0.insert(edge.end-1,i) = 1;
    }

    if (gDebug) {
        std::cout << "EIGEN d0: " << eigenD0.rows() << " " << eigenD0.cols() << std::endl;
        std::cout << eigenD0 << std::endl;
    }
}


bool edgeOrientation(Triangle triangle, Edge edge) {
    if (triangle.v1 == edge.start && triangle.v2 == edge.end) {
        return true;
    }
    if (triangle.v2 == edge.start && triangle.v3 == edge.end) {
        return true;
    }
    if (triangle.v3 == edge.start && triangle.v1 == edge.end) {
        return true;
    }
    return false;
}

void setupEigenD1() {
    int numEdges = PrimalEdges.size();
    int numTriangles = PrimalTriangles.size();
    eigenD1.resize(numEdges, numTriangles);
    for (size_t i = 0; i < PrimalTriangles.size(); i++) {
        Triangle &triangle = PrimalTriangles[i];

        int e1 = triangle.e1;
        int e2 = triangle.e2;
        int e3 = triangle.e3;

        Edge edge1 = PrimalEdges[e1-1];
        Edge edge2 = PrimalEdges[e2-1];
        Edge edge3 = PrimalEdges[e3-1];

        if (edgeOrientation(triangle, edge1)) {
            eigenD1.insert(e1-1, i) = 1;
        } else {
            eigenD1.insert(e1-1, i) = -1;
        }

        if (edgeOrientation(triangle, edge2)) {
            eigenD1.insert(e2-1, i) = 1;
        } else {
            eigenD1.insert(e2-1, i) = -1;
        }

        if (edgeOrientation(triangle, edge3)) {
            eigenD1.insert(e3-1, i) = 1;
        } else {
            eigenD1.insert(e3-1, i) = -1;
        }
    }

    if (gDebug) {
        std::cout << "EIGEN d1: " << std::endl;
        std::cout << eigenD1 << std::endl;
    }
}

void setupEigenH0() {
    int numVertices = PrimalVertices.size();
    eigenH0.resize(numVertices,numVertices);

    double total = 0;
    for (size_t i = 0; i < PrimalVertices.size(); i++) {
        if (gDebug) std::cout << "i" << std::endl;
        const vector <Vertex> &dualFace = DualFaces[i];
        if (dualFace.size() == 0) {
            eigenH0.insert(i,i) = 0;
        } else {
            vector <double> Xcoords;
            vector <double> Ycoords;
            int numPoints = dualFace.size();
            for (size_t j = 0; j < dualFace.size(); j++) {
                //int index = dualFace[j];
                //std::cout << "ID: " << index << std::endl;
                //std::cout << index << " ";
                Xcoords.push_back(dualFace[j].x);
                Ycoords.push_back(dualFace[j].y);
            }
            //std::cout << std::endl;
            double areaDual = polygonArea(Xcoords, Ycoords, numPoints);
            eigenH0.insert(i,i) = -areaDual; // For Hodge, dualFace/primalVertex = area/1 = area
            total += -areaDual;
        }
    }

    if (gDebug) {
        std::cout << "EIGEN h0: " << total << std::endl;
        std::cout << eigenH0 << std::endl;
    }
}

void setupEigenH1() {
    int numEdges = PrimalEdges.size();
    eigenH1.resize(numEdges, numEdges);

    for (size_t i = 0; i < PrimalEdges.size(); i++) {
        Edge &edge = PrimalEdges[i];

        Vertex &v1 = PrimalVertices[edge.start-1]; // account for 0-indexing here!
        Vertex &v2 = PrimalVertices[edge.end-1]; // account for 0-indexing here!

        double xSquared = (v2.x - v1.x) * (v2.x - v1.x);
        double ySquared = (v2.y - v1.y) * (v2.y - v1.y);
        double length = sqrt(xSquared + ySquared);
        //std::cout << "Length " << length << std::endl;
        if (edge.triangle1 >= 0 && edge.triangle2 >= 0) {
            Triangle triangle1 = PrimalTriangles[edge.triangle1-1];
            Triangle triangle2 = PrimalTriangles[edge.triangle2-1];
            double xSq = (triangle2.x-triangle1.x) * (triangle2.x-triangle1.x);
            double ySq = (triangle2.y-triangle1.y) * (triangle2.y-triangle1.y);
            double d = sqrt(xSq + ySq);
            eigenH1.insert(i,i) = d/length; // For Hodge, dualEdge/primalEdge = length/length
        } else {
            if (edge.triangle1 >= 0) {
                Triangle triangle1 = PrimalTriangles[edge.triangle1-1];
                double xSq = (edge.midX-triangle1.x) * (edge.midX-triangle1.x);
                double ySq = (edge.midY-triangle1.y) * (edge.midY-triangle1.y);
                double d = sqrt(xSq + ySq);
                eigenH1.insert(i,i) = d/length; // For Hodge, dualVertex/primalFace = 1/area
            } else {
                eigenH1.insert(i,i) = 0; // For Hodge, dualEdge/primalEdge = length/length
            }
            //double d = edge.dualEdge1 + edge.dualEdge2;
            //std::cout << "Half D: " << d  << " "<< edge.dualEdge1 << " " << edge.dualEdge2 << std::endl;
            //eigenH1.insert(i,i) = d/length; // For Hodge, dualVertex/primalFace = 1/area
        }
    }
    if (gDebug) {
        std::cout << "EIGEN h1: " << std::endl;
        std::cout << eigenH1 << std::endl;
    }
}

void setupEigenH2() {
    int numTriangles = PrimalTriangles.size();
    //eigenH0.resize(numVertices,numVertices);
    eigenH2.resize(numTriangles, numTriangles);

    double total = 0;
    for (size_t i = 0; i < PrimalTriangles.size(); i++) {
        //std::cout << "i" << std::endl;
        //vector <Vertex> dualFace = DualFaces[i];
        Triangle &triangle = PrimalTriangles[i];

        if (gDebug) std::cout << triangle.v1 << triangle.v2 << triangle.v3 << std::endl;
        Vertex &vertA = triangle.mVert1; 
        Vertex &vertB = triangle.mVert2;
        Vertex &vertC = triangle.mVert3;

        //std::cout << "xcords" << std::endl;

        vector <double> Xcoords;
        vector <double> Ycoords;
        Xcoords.push_back(vertA.x);
        Xcoords.push_back(vertB.x);
        Xcoords.push_back(vertC.x);
        Ycoords.push_back(vertA.y);
        Ycoords.push_back(vertB.y);
        Ycoords.push_back(vertC.y);

        //std::cout << "area" << std::endl;

        double areaFace = polygonArea(Xcoords, Ycoords, 3);
        //std::cout << areaFace << std::endl;
        total += -areaFace;
        eigenH2.insert(i,i) = -1/areaFace; // For Hodge, dualVertex/primalFace = 1/area
    }

    if (gDebug) {
        std::cout << "EIGEN h2: " << total << std::endl;
        std::cout << eigenH2 << std::endl;
    }
}

void setupInitialParticles() {
    double radius = 0.1;
    double numParticles = 500;

    /*
    double scale = 100.0;
    for (int i = 0; i < scale; i++) {

        Vec3 position;
        position.x = ((0.46 - 0.15) * i / scale) + 0.15;
        position.y = ((0.01 - (-0.335)) * i / scale) + (-0.335);
        int tri = findTriangle(position.x, position.y);
        position.z = tri;
        Particles.push_back(position);
    }
    */

    for (int l = 0; l < 2; l++) {
        Vertex loc = (l == 0) ? gLoc1 : gLoc2;
        for (int i = 0; i < numParticles; i++) {
            Vertex position;
            double rand1 = fRand(0, radius);
            double rand2 = M_PI * fRand(-1.0, 1.0);

            position.x = loc.x + rand1 * cos(rand2);
            position.y = loc.y + rand1 * sin(rand2);
            position.triIndex = findTriangle(position.x, position.y);
            Particles.push_back(position);
        }
    }
}

void setupInitialFlux() {

    gSumVorticity = 0.0;
    
    double vorticity = 0.0;

    for (int l = 0; l < 2; l++) {
        Vertex loc = (l == 0) ? gLoc1 : gLoc2;
        int row = 0.5 * (loc.y + 1) * (gMeshHeight - 1);
        int col = 0.5 * (loc.x + 1) * (gMeshWidth - 1);

        for(int i = -1; i < 2; i++) {  
            for(int k = -1; k < 2; k++) {
                int index = (row+i) * gMeshWidth + (col+k);
                if ((i == 0) && (k==0)) {
                    vorticity = 3.0/9.0;
                    //vorticity = 0.3;
                } else {
                    vorticity = 0.5/9.0;
                    //vorticity = 0.1;
                }
                PrimalVertices[index].vorticity = vorticity;
                gSumVorticity += vorticity;
            }
        }
    }

    if (gDebug) std::cout << "gSumVorticity " << gSumVorticity << std::endl;

    /* center point */
#if 0
    int row = 0.5 * (gMeshHeight - 1);
    int col = 0.5 * (gMeshWidth - 1);
    int index = (row) * gMeshWidth + (col);
    vorticity = 9.0/9.0;
    PrimalVertices[index].vorticity = vorticity;
    gSumVorticity += vorticity;
#endif

}

void solveLaplacian() {

    if (gDebug) std::cout << "solveLaplacian start" << std::endl;

    //Eigen::SparseMatrix<double> H1timesD0 = eigenH1 * eigenD0;

    Eigen::SparseMatrix<double> D0Transpose = eigenD0.transpose();
    if (gDebug) std::cout << D0Transpose << std::endl;

    Eigen::SparseMatrix<double> H1timesD0Trans = eigenH1 * D0Transpose;
    if (gDebug) std::cout << "*1 * d0T" << std::endl;
    if (gDebug) std::cout << H1timesD0Trans << std::endl;

    Eigen::SparseMatrix<double> D0timesH1timesD0Trans = eigenD0 * H1timesD0Trans;
    if (gDebug) std::cout << "d0 x *1 * d0T" << std::endl;
    if (gDebug) std::cout << D0timesH1timesD0Trans << std::endl;

    int numVertices = PrimalVertices.size();
    //eigenH0.resize(numVertices,numVertices);
    //Delta.resize(numVertices, numVertices);
    Delta = D0timesH1timesD0Trans; 

    if (gDebug) std::cout << "skipping initial solve" << std::endl;
#if 0
    int n = PrimalVertices.size();
    Eigen::VectorXd X(n), B(n);
    for (int i = 0; i < n; i++) {
        B[i] = PrimalVertices[i].vorticity;
    }
    if (gDebug) std::cout << "B" << std::endl;
    if (gDebug) std::cout << B << std::endl;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
    cg.compute(D0timesH1timesD0Trans);
    std::cout << "solve start" << std::endl;
    X = cg.solve(B);
    std::cout << "solve end" << std::endl;
    if (gDebug) std::cout << "X: " << std::endl;
    if (gDebug) std::cout << X << std::endl;

    Eigen::SparseMatrix<double> Theta(PrimalVertices.size(),1);
    for (int i = 0; i < PrimalVertices.size(); i++) {
        //std::cout << i << std::endl;
        Theta.insert(i,0) = X[i];
    }
    if (gDebug) std::cout << "X0: " << X[0] << std::endl;
    if (gDebug) std::cout << Theta << std::endl;

    Eigen::SparseMatrix<double> Fluxes = D0Transpose * Theta;
    if (gDebug) std::cout << "Fluxes: " << std::endl;
    if (gDebug) std::cout << Fluxes << std::endl;

    if (gDebug) std::cout << "iterator: " << Fluxes.outerSize()<<  std::endl;
    for (int k=0; k<Fluxes.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Fluxes,k); it; ++it){
            PrimalEdges[it.index()].flux = it.value(); // assign Flux!
        }
    }
#endif
    if (gDebug) std::cout << "solveLaplacian end" << std::endl;
}

// sets the texture stored in fileName
void setTexture(const std::string &fileName, int textureNumber) {

    int w, h;
    png_bytepp p = readpng(fileName.c_str(), &w, &h);

    unsigned char *data = new unsigned char [w * h * 3];

    for(int y=0; y<h; y++) {
        png_bytep r = p[h-1-y]; // get the row
        for (int x = 0; x < 3*w; x += 3) {
            int index = x + (3 * w * y);
            data[index] = r[x];
            data[index + 1] = r[x+1];
            data[index + 2] = r[x+2];
        }
    }

    glBindTexture(GL_TEXTURE_2D, gTexNames[textureNumber]);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // use texture colors as is
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // how pixels interpolated
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, data);

}

/**
 * Main entrance point, obviously.
 * Sets up some stuff then passes control to glutMainLoop() which never
 * returns.
 */
int main(int argc, char* argv[])
{
    if (argc > 4 || argc < 3) {
        std::cerr << "USAGE: oglRenderer [xRes] [yRes] [Run]" << std::endl;
        exit(1);
    }

    if (argc == 4) {
        gRun = atoi(argv[3]);
    }

    xRes = atoi(argv[1]);
    yRes = atoi(argv[2]);

    assert(xRes >= 0);
    assert(yRes >= 0);

    //std::string prefix = "Mesh2.1";
    //std::string prefix = "spiral.1";
    //std::string prefix = "simpleMesh.1";
    //std::string prefix = "square_circle_hole.1";
    std::string prefix = "gridMesh40.1";
    //std::string prefix = "gridMesh19.1";
    //std::string prefix = "gridMesh.1";

    gMeshWidth = 40;
    gMeshHeight = 40;

    std::string nodeFilename(prefix + ".node");
    std::string edgeFilename(prefix + ".edge");
    std::string eleFilename(prefix + ".ele");
    std::string vNodeFilename(prefix + ".v.node");
    std::string vEdgeFilename(prefix + ".v.edge");

    double maxScale = 1.0;

    std::ifstream nodeFile(nodeFilename.c_str());
    double a, b, c, d, e, f;
    nodeFile >> a >> b >> c >> d;
    int numVertices = a;
    int dimension = b;
    int numAttributes = c;
    int boundary = d;
    while (nodeFile >> a >> b >> c >> d)
    {
        //std::cout << a << " " << b << " " << c << " " << d << std::endl;

        Vertex vertex;
        vertex.index = a;
        vertex.x = b/maxScale;
        vertex.y = c/maxScale;
        vertex.boundary = d;
        /*
        if (b/maxScale < -0.93 || b/maxScale > 0.93) {
            vertex.boundary = 0;
        }
        */
        PrimalVertices.push_back(vertex);

    }

    std::ifstream edgeFile(edgeFilename.c_str());
    edgeFile >> a >> b;

    while (edgeFile >> a >> b >> c >> d)
    {
        //std::cout << a << " " << b << " " << c << " " << d  << " " << e << " " << f << std::endl;
        //Vec3 vert = Vec3(b,c,f); // 3rd coordinate will tell us if boundary or not...
        //std::cout << "Edge: " << a << "  vertex: " << b << " to vertex: " << c << "  boundary?" << d << std::endl;

        Edge edge(a, b, c, d);

        const Vertex &vertA = PrimalVertices[b-1]; // account for 0-indexing here!
        const Vertex &vertB = PrimalVertices[c-1];

        double avg_x = (vertA.x + vertB.x) / 2.0;
        double avg_y = (vertA.y + vertB.y) / 2.0;

        edge.midX = avg_x;
        edge.midY = avg_y;

        PrimalEdges.push_back(edge);
    }

    std::ifstream eleFile(eleFilename.c_str());
    eleFile >> a >> b >> c;
    //int numVertices = a;
    int nodesPerElem = b;
    //int numAttributes = c;
    while (eleFile >> a >> b >> c >> d)
    {
        //std::cout << "triangle: " << b << " " << c <<  " " << d << std::endl;

        Vertex &vertA = PrimalVertices[b-1]; // account for 0-indexing here!
        Vertex &vertB = PrimalVertices[c-1];
        Vertex &vertC = PrimalVertices[d-1];

        Triangle triangle(a, b, c, d, vertA, vertB, vertC);


        //std::cout << "A: " << vertA.x << " " << vertA.y << std::endl;
        //std::cout << "B: " << vertB.x << " " << vertB.y << std::endl;
        //std::cout << "C: " << vertC.x << " " << vertC.y << std::endl;
        double x1 = vertA.x;
        double x2 = vertB.x;
        double x3 = vertC.x;
        double y1 = vertA.y;
        double y2 = vertB.y;
        double y3 = vertC.y;

        double mid1X = (x1 + x2) / 2.0;
        double mid1Y = (y1 + y2) / 2.0;
        double mid2X = (x1 + x3) / 2.0;
        double mid2Y = (y1 + y3) / 2.0;

        double m1Inverse, m2Inverse;
        if ((x2 - x1) == 0) {
             m1Inverse = 9999.9;
        } else {
            double slope1 = (y2-y1) / (x2-x1);
             m1Inverse = -1 / slope1;
        }
        if ((x3-x1) == 0) {
             m2Inverse = 9999.9;
        } else {
            double slope2 = (y3-y1) / (x3-x1);
             m2Inverse = -1 / slope2;
        }
        //double m1 = (y2-y1) / (x2-x1);
        //double m2 = (y3-y1) / (x3-x1);
        double m1 = m1Inverse;
        double m2 = m2Inverse;

        double b1 = mid1Y - (m1 * mid1X);
        double b2 = mid2Y - (m2 * mid2X);

        double circumX = (b2-b1) / (m1-m2);
        double circumY = ((b1 * m2) - (b2 * m1))/ (m2 - m1);


        for (size_t i = 0; i < PrimalEdges.size(); i++) {
            const Edge &edge = PrimalEdges[i];
            int result = triangle.containEdge(edge);
            if (result == 0) {
                triangle.e1 = i+1;
            } else if (result == 1) {
                triangle.e2 = i+1;
            } else if (result == 2) {
                triangle.e3 = i+1;
            }

            if (result >= 0) {
                double xSq = (triangle.x - edge.midX) * (triangle.x - edge.midX);
                double ySq = (triangle.y - edge.midY) * (triangle.y - edge.midY);
                double distance = sqrt(xSq + ySq);
                if (edge.triangle1 < 0) {
                    PrimalEdges[i].triangle1 = a;
                    PrimalEdges[i].dualEdge1 = distance;
                } else {
                    PrimalEdges[i].triangle2 = a;
                    PrimalEdges[i].dualEdge2 = distance;
                }
            }
        }
        PrimalTriangles.push_back(triangle);
    }

    for (size_t i = 0; i < PrimalVertices.size(); i++) {
        //vector <int> sortedFace;
        Vertex vertex = PrimalVertices[i];
        vector <Vertex> dualFace;
        if (vertex.boundary == 1) { 
            //dualFaces.push_back(sortedFace); // find the boundary edges, add midpoints
            for (size_t k = 0; k < PrimalEdges.size(); k++) {
                Edge edge = PrimalEdges[k];
                if (edge.start == vertex.index || edge.end == vertex.index) {
                    if (edge.boundary == 0) {
                        continue;
                    }
                    Vertex mid(edge.midX, edge.midY);
                    mid.triIndex = findClosestTriangle(mid.x, mid.y);
                    dualFace.push_back(mid);
                }
            }
            Vertex actualVertex(vertex.x, vertex.y);
            actualVertex.triIndex = findClosestTriangle(actualVertex.x, 
                                                        actualVertex.y);
            dualFace.push_back(actualVertex);
            //continue;
        }
        int vertexNum = i + 1; // vertices start at 1
        for (size_t j = 0; j < PrimalTriangles.size(); j++) {
            Triangle &triangle = PrimalTriangles[j];
            if (triangle.v1 == vertexNum || triangle.v2 == vertexNum || triangle.v3 == vertexNum) {
                //std::cout << "MATCH: " << j << " " << vertexNum << std::endl;
                //std::cout << "Center: " << center.x << " " << center.y << std::endl;
                Vertex dualVertex(triangle.x, triangle.y);
                dualVertex.triIndex = j;
                dualFace.push_back(dualVertex);
            }
        }

        //std::cout << "Sorting " << dualFace.size() << std::endl;
        //std::cout << "Sorting..." << i << std::endl;
        sortVertices2(dualFace);
        //std::cout << "Sorted " << sortedDualFace.size() << std::endl;

        DualFaces.push_back(dualFace);
        vector <int> dualFaceByIndex;

        if (gDebug) std::cout << "NEW DUAL" << std::endl;
        for (size_t i = 0; i < dualFace.size(); i++) {
            Vertex v = dualFace[i];
            int t = v.triIndex;
            if (t == -1) {
                continue;
            }
            for (size_t j = 0; j < dualFaceByIndex.size(); j++) {
                if (dualFaceByIndex[j] == t) {
                    continue;
                }
            }
            if (gDebug) std::cout << t << std::endl;
            dualFaceByIndex.push_back(t);
        }
        DualFacesIndices.push_back(dualFaceByIndex);
    }

    setupEigenD0();
    setupEigenD1();

    if (gDebug) {
    Eigen::SparseMatrix<double> eigenDD = eigenD1.transpose() * eigenD0.transpose();//eigenD0 * eigenD1;
        std::cout << "eigenDD: (This should be 0!) " << std::endl;
        std::cout << eigenDD << std::endl;
    }
    setupEigenH0();
    setupEigenH1();
    setupEigenH2();

    setupInitialFlux();
    solveLaplacian();
    setupDeltaInterior();
    calculateFlux();

    setupInitialParticles();

    // this will set up initial velocity
    calculateVelocity();


    // OpenGL will take out any arguments intended for its use here.
    // Useful ones are -display and -gldebug.
    glutInit(&argc, argv);

    // Get a double-buffered, depth-buffer-enabled window, with an
    // alpha channel.
    // These options aren't really necessary but are here for examples.
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

    glutInitWindowSize(xRes, yRes);
    glutInitWindowPosition(300, 100);

    glutCreateWindow("CS176 Project");
    
    initGL();

    // set up GLUT callbacks.
    glutDisplayFunc(redraw);
    glutReshapeFunc(resize);
    glutKeyboardFunc(keyfunc);

    glutMouseFunc(mouseButton);
    if (gRun) {
        glutIdleFunc(idle);
    }


    // From here on, GLUT has control,
    glutMainLoop();

    // so we should never get to this point.
    return 1;
}

