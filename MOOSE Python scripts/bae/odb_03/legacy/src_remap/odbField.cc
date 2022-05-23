#include "odbField.h"

#include <sstream>

#include <stdlib.h>


/**************************************************************/
/**************************************************************/
/**************************************************************/

// constant ids
const odb_String OdbFieldReaderCommonData::instanceName("PART-1-1");
const odb_String OdbFieldReaderCommonData::nsetName("nodesForInterpolation");
const odb_String OdbFieldReaderCommonData::elsetName("elemsForInterpolation");

// constructor
OdbFieldReaderCommonData::OdbFieldReaderCommonData(odb_Odb& odb) {
  m_instance = &(odb.rootAssembly().instances()[instanceName]);
}

// create nset of nodes required for the interpolation
void OdbFieldReaderCommonData::initNset(const std::set<unsigned int>& nodes) {
    
  // initialize odb_SequenceInt object...
  odb_SequenceInt odbSeqNodes( nodes.size() );

  // copy to odb_SequenceInt object...
  for (std::set<unsigned int>::iterator it=nodes.begin();
       it != nodes.end();  ++it) {
    odbSeqNodes.append(*it);
  }
  
  // create odb_Set
  m_instance->NodeSet(nsetName, odbSeqNodes);
  m_nset = &(m_instance->nodeSets()[nsetName]);
}

// create nset of nodes required for the interpolation
void OdbFieldReaderCommonData::initElset(const std::set<unsigned int>& elems) {
    
  // initialize odb_SequenceInt object...
  odb_SequenceInt odbSeqElems( elems.size() );

  // copy to odb_SequenceInt object...
  for (std::set<unsigned int>::iterator it=elems.begin();
       it != elems.end();  ++it) {
    odbSeqElems.append(*it);
  }
  
  // create odb_Set
  m_instance->ElementSet(elsetName, odbSeqElems);
  m_elset = &(m_instance->elementSets()[elsetName]);
}


/**************************************************************/
/**************************************************************/
/**************************************************************/

void OdbFieldDataMapNodal::store(const odb_SequenceFieldValue& values)
{
  int nbComponents = 0; // will be supplied by the call to values[i].data(...)

  // Iterating over all values:
  for (int i=0; i<values.size(); i++) {
    const odb_FieldValue& val = values[i];  // just an abbreviation
    // Storing pointer to first component:
    (*this)[val.nodeLabel()] = val.data(nbComponents);
  }
}

void OdbFieldDataMapIntegrationPt::store(const odb_SequenceFieldValue& values)
{
  int nbComponents = 0; // will be supplied by the call to values[i].data(...)

  // Iterating over all values:
  for (int i=0; i<values.size(); i++) {
    const odb_FieldValue& val = values[i];  // just an abbreviation
    // Storing pointer to first component:
    (*this)[val.elementLabel()].push_back(val.data(nbComponents));
  }
}

void OdbFieldDataMapWholeElem::store(const odb_SequenceFieldValue& values)
{
  int nbComponents = 0; // will be supplied by the call to values[i].data(...)

  // Iterating over all values:
  for (int i=0; i<values.size(); i++) {
    const odb_FieldValue& val = values[i];  // just an abbreviation
    // Storing pointer to first component:
    (*this)[val.elementLabel()] = val.data(nbComponents);
  }
}

/**************************************************************/
/**************************************************************/
/**************************************************************/

// common data for all fields, pointer needs to be initialized externally
OdbFieldReaderCommonData* OdbFieldReader::m_odbCommonData = 0;

OdbFieldReader::OdbFieldReader(const std::string& fieldName)
  : m_map(0) {

  m_fieldName = fieldName;

  if (fieldName.substr(0, 3) == "SDV") {
    // SDVs are scalar and may contain a '_' in the odbFieldName
    m_odbFieldName = new odb_String(fieldName.c_str());
    m_componentName.clear();
  }

  else {

    std::istringstream sstr(fieldName);

    // extract odbFieldName from fieldName
    std::string x;
    std::getline(sstr, x, '_');
    m_odbFieldName = new odb_String(x.c_str());
  
    // store remainder in componentName (only if there is a '_')
    sstr.get();   // skip '_'
    if (sstr.good()) {
      sstr >> m_componentName;
    }
  }
}

OdbFieldReader::~OdbFieldReader() {
  if (m_odbFieldName) delete m_odbFieldName;
  if (m_map) delete m_map;
}

// read field type information from the odb and fill the
void OdbFieldReader::getFieldFromFrame(odb_Frame& frame) {

  if (!frame.fieldOutputs().isMember(*m_odbFieldName)) {
    cout << "Field " << m_odbFieldName->cStr() << " not found in frame "
	 << frame.incrementNumber() << "." << endl;
    exit(3);
  }
  odb_FieldOutput& odbField = frame.fieldOutputs()[*m_odbFieldName];

  // dataType: scalar, vector, tensor
  m_dataType = odbField.type();

  // location: nodal, integration point, ...
  if (odbField.locations().size() != 1) {
    // this field is not only defined on nodes or integration points but on
    // more than one locations at the same time...
    cout << "Odb field " << m_fieldName << " is defined on "
	 << odbField.locations().size() << " locations at the same time."
	 << " Can't handle that." << endl;
    exit(3);
  }
  m_location = odbField.locations(0).position();

  // depending on location (node, inegration point...)
  // ... get values (odb_SequenceFieldValue) of subset...
  odb_FieldOutput subField;

  if (m_location == odb_Enum::NODAL) {
    subField = odbField.getSubset(*(m_odbCommonData->m_nset));
    if (m_map) delete m_map;
    m_map = new OdbFieldDataMapNodal();
  }
  else if (m_location==odb_Enum::INTEGRATION_POINT) {
    subField = odbField.getSubset(*(m_odbCommonData->m_elset));
    if (m_map) delete m_map;
    m_map = new OdbFieldDataMapIntegrationPt();
  }
  else if (m_location==odb_Enum::WHOLE_ELEMENT) {
    subField = odbField.getSubset(*(m_odbCommonData->m_elset));
    if (m_map) delete m_map;
    m_map = new OdbFieldDataMapWholeElem();
  }
  else {
    cout << "Unrecognized field location " << m_location << " for field "
	 << m_fieldName << ". Stopping." << endl;
    exit(3);
  }

  // actually store the label -> value relation
  m_map->store(subField.values());
}
