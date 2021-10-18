#include "graph.h"
#include "program.h"
#include "mesh.h"
#include <cmath>
#include <chrono>

using namespace simit;

std::array<float, 3> hack(std::array<double, 3> k) {
  return {float(k[0]), float(k[1]), float(k[2])};
}

int main(int argc, char **argv)
{
  if (argc != 3) {
    std::cerr << "Usage: springs <path to simit code> <path to data>" << std::endl;
    return -1;
  }
  std::string codefile = argv[1];
  std::string datafile = argv[2];

  simit::init("gpu", sizeof(float));

  // Load mesh data using Simit's mesh loader.
  MeshVol mesh;
  mesh.loadTet(datafile+".node", datafile+".ele");
  mesh.loadTetEdge(datafile+".edge");
  mesh.makeTetSurf();

  // Move the data to a floor and normalize to 1
  std::array< float,3> mn = hack(mesh.v[0]), mx = hack(mesh.v[0]);
  for(size_t vi = 0; vi < mesh.v.size(); vi++) {
    for(int i = 0; i < 3; i++) {
      mn[i] = std::min(float(mesh.v[vi][i]), mn[i]);
      mx[i] = std::max(float(mesh.v[vi][i]), mx[i]);
    }
  }
  for(size_t vi = 0; vi<mesh.v.size(); vi++) {
    for(int i = 0; i<3; i++) {
      mesh.v[vi][i] = (mesh.v[vi][i] - mn[2]) / (mx[2] - mn[2]);
    }
  }

  // Create a graph and initialize it with mesh data
  Set points;
  Set springs(points, points);

  float stiffness = 1e3;
  float density   = 1e3;
  float radius    = 0.01;
  float pi        = 3.14159265358979;
  float zfloor    = 0.1;               // we fix everything below the floor

  // The fields of the points set 
  FieldRef<float,3> x     = points.addField<float,3>("x");
  FieldRef<float,3> v     = points.addField<float,3>("v");
  FieldRef<float>   m     = points.addField<float>("m");
  FieldRef<bool>     fixed = points.addField<bool>("fixed");

  // The fields of the springs set 
  FieldRef<float> k  = springs.addField<float>("k");
  FieldRef<float> l0 = springs.addField<float>("l0");

  std::vector<ElementRef> pointRefs;
  for(auto vertex : mesh.v) {
    ElementRef point = points.add();
    pointRefs.push_back(point);

    x.set(point, vertex);
    v.set(point, {0.0, 0.0, 0.0});
    fixed.set(point, vertex[2] < zfloor);
  }

  // Compute point masses
  std::vector<float> pointMasses(mesh.v.size(), 0.0);

  for(auto e : mesh.edges) {
    float dx[3];
    //float *x0 = &(mesh.v[e[0]][0]);
    //float *x1 = &(mesh.v[e[1]][0]);
    std::array<float, 3> x0 = hack(mesh.v[e[0]]), x1 = hack(mesh.v[e[1]]);
    dx[0] = x1[0] - x0[0];
    dx[1] = x1[1] - x0[1];
    dx[2] = x1[2] - x0[2];
    float l0_ = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    float vol = pi*radius*radius*l0_;
    float mass = vol*density;
    pointMasses[e[0]] += 0.5*mass;
    pointMasses[e[1]] += 0.5*mass;
    ElementRef spring = springs.add(pointRefs[e[0]], pointRefs[e[1]]);
    l0.set(spring, l0_);
    k.set(spring, stiffness);
  }
  for(size_t i = 0; i < mesh.v.size(); ++i) {
    ElementRef p = pointRefs[i];
    m.set(p, pointMasses[i]);
  }

  // Compile program and bind arguments
  Program program;
  program.loadFile(codefile);

  Function timestep = program.compile("timestep");
  timestep.bind("points",  &points);
  timestep.bind("springs", &springs);

  timestep.init();

  // Take 100 time steps
  for (int i = 1; i <= 1; ++i) {

    
    timestep.unmapArgs(); // Move data to compute memory space (e.g. GPU)
    auto start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < 1000; j++) {
    timestep.run();       // Run the timestep function
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration<double>(end - start).count() << "\n";
    timestep.mapArgs();   // Move data back to this memory space

    // Copy the x field to the mesh and save it to an obj file
    std::cout << "timestep " << i << std::endl;
    int vi = 0;
    for (auto &vert : points) {
      for(int ii = 0; ii < 3; ii++){
        mesh.v[vi][ii] = x.get(vert)(ii);
      }
      vi++;
    }
    mesh.updateSurfVert();
    mesh.saveTetObj(std::to_string(i)+".obj");
    float sum = 0;
    for(size_t vi = 0; vi<mesh.v.size(); vi++) {
      for(int i = 0; i<3; i++) {
        sum += float(mesh.v[vi][i]);
      }
    }
    printf("%f\n", sum / mesh.v.size() / 3);
  }
}
