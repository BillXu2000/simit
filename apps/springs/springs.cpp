#include "graph.h"
#include "program.h"
#include "mesh.h"
#include <cmath>
#include <chrono>
#include <fstream>

using namespace simit;

std::array<float, 3> hack(std::array<double, 3> k) {
  return {float(k[0]), float(k[1]), float(k[2])};
}

class float3 {
public:
  float x[3];

  float3() {}

  float3(float *k) {
    for (int i = 0; i < 3; i++) x[i] = k[i];
  }

  float3(std::array<double, 3> k) {
    for (int i = 0; i < 3; i++) x[i] = k[i];
  }

  float operator () (int i) {
    return x[i];
  }
};

float3 operator - (float3 a, float3 b) {
  float ans[] = {a(0) - b(0), a(1) - b(1), a(2) - b(2)};
  return float3(ans);
}

float determinant(float3 a, float3 b, float3 c) {
  float ans = 0;
  for (int i = 0; i < 3; i++) {
    ans += a(i) * b((i + 1) % 3) * c((i + 2) % 3);
    ans -= a(i) * b((i + 2) % 3) * c((i + 1) % 3);
  }
  return ans;
}

int main(int argc, char **argv)
{
  if (argc < 3) {
    std::cerr << "Usage: springs <path to simit code> <path to data>" << std::endl;
    return -1;
  }
  std::string codefile = argv[1];
  std::string datafile = argv[2];

  if (argc >= 4 && std::string(argv[3]) == "cpu") simit::init("cpu", sizeof(float));
  else simit::init("gpu", sizeof(float));

  // Load mesh data using Simit's mesh loader.
  MeshVol mesh;
  if (datafile.substr(datafile.size() - 3) == "obj") {
    std::fstream in(datafile, std::fstream::in);
    mesh.load_obj(in);
  }
  else {
    mesh.loadTet(datafile+".node", datafile+".ele");
    mesh.loadTetEdge(datafile+".edge");
  }
  //mesh.makeTetSurf();
  puts("hello");

  // Move the data to a floor and normalize to 1
  std::array< float,3> mn = hack(mesh.v[0]), mx = hack(mesh.v[0]);
  for(size_t vi = 0; vi < mesh.v.size(); vi++) {
    for(int i = 0; i < 3; i++) {
      mn[i] = std::min(float(mesh.v[vi][i]), mn[i]);
      mx[i] = std::max(float(mesh.v[vi][i]), mx[i]);
    }
  }
  int fix_axis = 1;
  for(size_t vi = 0; vi<mesh.v.size(); vi++) {
    for(int i = 0; i<3; i++) {
      mesh.v[vi][i] = (mesh.v[vi][i] - mn[i]) / (mx[fix_axis] - mn[fix_axis]);
    }
  }

  // Create a graph and initialize it with mesh data
  Set points;
  Set springs(points, points);

  float stiffness = 1e1;
  float density   = 1e3;
  float radius    = 0.01;
  float pi        = 3.14159265358979;
  float zfloor    = 0.1;               // we fix everything below the floor

  // The fields of the points set 
  FieldRef<float,3> x     = points.addField<float,3>("x");
  FieldRef<float,3> ox     = points.addField<float,3>("ox");
  FieldRef<float,3> v     = points.addField<float,3>("v");
  FieldRef<float>   m     = points.addField<float>("m");
  FieldRef<bool>     fixed = points.addField<bool>("fixed");

  // The fields of the springs set 
  FieldRef<float> k  = springs.addField<float>("k");
  //FieldRef<float> l0 = springs.addField<float>("l0");

  std::vector<ElementRef> pointRefs;
  for(auto vertex : mesh.v) {
    ElementRef point = points.add();
    pointRefs.push_back(point);

    x.set(point, vertex);
    ox.set(point, vertex);
    v.set(point, {0.0, 0.0, 0.0});
    fixed.set(point, vertex[fix_axis] < zfloor);
  }

  // Compute point masses
  std::vector<float> pointMasses(mesh.v.size(), 0.0);
  /*for (int i = 0; i < 10; i++) {
    printf("i = %d ", i);
    for (auto j : mesh.e[i]) {
      printf("%d ", j);
    }
    puts("");
  }*/
  /*for (auto e: mesh.e) {
    float3 k[3];
    for (int i = 0; i < 3; i++) {
      k[i] = float3(mesh.v[e[i]]) - float3(mesh.v[e[3]]);
    }
    float vol = abs(determinant(k[0], k[1], k[2]) / 6.0);
    for (int i = 0; i < 4; i++) {
      pointMasses[e[i]] += vol / 4.0;
    }
  }*/

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
    //pointMasses[e[0]] += 0.5*mass;
    //pointMasses[e[1]] += 0.5*mass;
    ElementRef spring = springs.add(pointRefs[e[0]], pointRefs[e[1]]);
    //l0.set(spring, l0_);
    //k.set(spring, stiffness * l0_);
    k.set(spring, stiffness);
  }
  for(size_t i = 0; i < mesh.v.size(); ++i) {
    ElementRef p = pointRefs[i];
    pointMasses[i] = 1.0 / mesh.v.size();
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
  timestep.unmapArgs(); // Move data to compute memory space (e.g. GPU)
  double timer = 0;
  int num_run = 2, num_step = 100;
  timestep.run();       // Run the timestep function
  puts("hello");
  for (int i = 1; i <= num_run; ++i) {

    
    auto start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < num_step; j++) {
    timestep.run();       // Run the timestep function
    }
    auto end = std::chrono::high_resolution_clock::now();
    timer += std::chrono::duration<double>(end - start).count();
    /*
    //std::cout << std::chrono::duration<double>(end - start).count() << "\n";
    timestep.mapArgs();   // Move data back to this memory space

    // Copy the x field to the mesh and save it to an obj file
    //std::cout << "timestep " << i << std::endl;
    int vi = 0;
    for (auto &vert : points) {
      for(int ii = 0; ii < 3; ii++){
        mesh.v[vi][ii] = x.get(vert)(ii);
      }
      vi++;
    }
    mesh.updateSurfVert();
    mesh.saveTetObj(std::to_string(i)+".obj");
    float sum[4];
    memset(sum, 0, sizeof(sum));
    for(size_t vi = 0; vi<mesh.v.size(); vi++) {
      for(int i = 0; i<3; i++) {
        float product = 1;
        for (int k = 0; k < 4; k++) {
          sum[k] += product;
          product *= mesh.v[vi][i];
        }
      }
    }
    printf("%f %f %f %f\n", sum[0], sum[1], sum[2], sum[3]);
    timestep.unmapArgs(); // Move data to compute memory space (e.g. GPU)*/
  }
  std::fstream out("timer.log", std::fstream::out);
  out << timer / (num_run * num_step) << "\n" << (num_run * num_step) << "\n";
  return 0;
}
