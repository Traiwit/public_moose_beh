#ifndef mesh_h_defined
#define mesh_h_defined

#include <fstream>
#include <string>


class MeshBaseType {
public:
  // pure-virtual function, i.e. subclass has to define it
  // and this class is an abstract class
  virtual unsigned int size() const = 0;
  virtual void readFromFile(std::ifstream& inFile) = 0;
};


class MeshStructuredPoints : public MeshBaseType {
public:
  unsigned int dims[3];  // resolution
  float org[3];    // origin
  float spc[3];    // spacing

  unsigned int size() const;
  void readFromFile(std::ifstream& inFile);

};

#endif
