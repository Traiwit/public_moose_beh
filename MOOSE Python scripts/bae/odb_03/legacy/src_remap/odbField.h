#ifndef odbField_h_defined
#define odbField_h_defined

#include <string>
#include <set>

#include <odb_API.h>

class OdbFieldReaderCommonData {
public:

  OdbFieldReaderCommonData(odb_Odb& odb);

  // instance
  odb_Instance* m_instance;

  // node and element sets required for the interpolation
  odb_Set * m_nset;
  odb_Set * m_elset;
  void initNset(const std::set<unsigned int>& nodes);
  void initElset(const std::set<unsigned int>& elems);

  // constant ids
  static const odb_String instanceName;
  static const odb_String nsetName;
  static const odb_String elsetName;

};


class OdbFieldDataMapBase
{
public:
  virtual void store(const odb_SequenceFieldValue& values) = 0;
};

class OdbFieldDataMapNodal : public OdbFieldDataMapBase,
  public std::map<unsigned int, const float*>
{
public:
  void store(const odb_SequenceFieldValue& values);
};

class OdbFieldDataMapIntegrationPt : public OdbFieldDataMapBase,
  public std::map<unsigned int, std::vector<const float*> >
{
public:
  void store(const odb_SequenceFieldValue& values);
};


class OdbFieldDataMapWholeElem : public OdbFieldDataMapBase,
  public std::map<unsigned int, const float*>
{
public:
  void store(const odb_SequenceFieldValue& values);
};


class OdbFieldReader {
public:
  // common data for all fields, pointer needs to be initialized externally
  static OdbFieldReaderCommonData* m_odbCommonData;

  std::string m_fieldName;
  odb_String * m_odbFieldName;
  std::string m_componentName;

  // nodal or integration pt, i.e.:
  // odb_Enum::NODAL, odb_Enum::INTEGRATION_POINT, odb_Enum::ELEMENT_NODAL,
  // odb_Enum::ELEMENT_FACE, odb_Enum::CENTROID, odb_Enum::WHOLE_ELEMENT
  odb_Enum::odb_ResultPositionEnum m_location;

  // scalar, vector, tensor, i.e.:
  // odb_Enum::SCALAR, odb_Enum::VECTOR or odb_Enum::TENSOR_3D_FULL
  odb_Enum::odb_DataTypeEnum m_dataType;

  // map to the actual data
  OdbFieldDataMapBase* m_map;

public:
  OdbFieldReader(const std::string& fieldName);
  virtual ~OdbFieldReader();

  // read field type information from the odb and fill the
  void getFieldFromFrame(odb_Frame& frame);

};

#endif
