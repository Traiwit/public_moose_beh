/* getMatNames
 *
 * This is to read the material assignment from the odb.
 *
 * If createFilterElset or setFilterElset has been called before only elements
 * from this set will be considered.
 *
 * input arguments: firstMatCode (int), nbMatNameFilters (int),
 *     matNameFilters (strings)
 *
 * firstMatCode determines the lowest (possible) material number in the result
 * corresponding to the first material name returned. It's typically =1.
 *
 * nbMatNameFilters is the number of items in the following matNameFilters
 * sequence of strings.
 * If you don't need this feature then pass zero for nbMatNameFilters and skip
 * the matNameFilters argument(s) alltogether.
 *
 * matNameFilters are strings. Only section / material names starting with one
 * of those strings will be considered in the output. Elements whose section /
 * material name does not begin with one of these strings will be ignored.
 *
 *
 * Output:
 *
 * The first output of getMatNames is an array of element labels.
 * The second output is an array of the material numbers.
 * The third output is a sequence of material names.
 *
 * Here is the data in detail that will be written to the pipe:
 * - nb of items =:nbElems
 * - integer array of element labels (nbElems ints)
 * - integer array material numbers (nbElems ints)
 * - nb of material names =:nbMatNames
 * - a sequence of nbMatNames strings giving the material (section) names from
 *   the odb.
 */

#include <algorithm>

int getMatNames() {
  // no matcodes allowed smaller or equal to this
  const int noMatCode(-9999999);
  
  /* get and check arguments */
  int firstMatCode = pipeReadInt();
  int nbMatNameFilters = pipeReadInt();
  std::vector<std::string> matNameFilters;
  for (int i=0; i<nbMatNameFilters; ++i) {
    matNameFilters.push_back(pipeReadString());
  }
  debug_print("This is getMatNames, firstMatCode %d, %d type filter strings.\n",
              firstMatCode, nbMatNameFilters);
  if (firstMatCode<=noMatCode) {
    error("getMatNames: firstMatCode must be greater than %d. Got %d.\n",
          noMatCode, firstMatCode);
  }
  
  /* initialize/get Abaqus instance */
  if (ptInstance == NULL) {
    setPtInstance();
  }

  /* determine applicable elements / apply filterset */
  std::set<int> filterset;
  if (filterElsetName && filterElsetName->length()>0) {
    odb_Set odbfilterset = ptInstance->elementSets()[*(filterElsetName)];
    const odb_SequenceElement& elements = odbfilterset.elements();
    for (int i=0; i<elements.size(); ++i) {
      filterset.insert(elements[i]);
    }
    debug_print("getMatNames: filter elset %s with %d elements.\n",
                odbfilterset.name().cStr(), filterset.size());
  }
  else {
    const odb_SequenceElement& elements = ptInstance.elements();
    for (int i=0; i<elements.size(); ++i) {
      filterset.insert(elements[i]);
    }
    debug_print("getMatNames: no filterset, considering all %d elements.\n",
                filterset.size());
  }

  /* examine section assignments */
  int nextNewMatCode = firstMatCode;
  int matCode;
  bool matNameIsInFilterTuple;
  
  // initialize result objects: arrays to store all data
  // because of it's (potential) size it's on the heap not the stack
  int nbElems = filterset.size();
  int* labels = new int[nbElems];
  int* matCodes = new int[nbElems];
  std::vector<std::string> matNamesList;
  const char* matName = NULL;
  
  odb_SequenceSectionAssignment sectionAssignmentSeq =
    ptInstance->sectionAssignments();
  for (int i = 0; i < sectionAssignmentSeq.size(); ++i) {
    odb_SectionAssignment sa = sectionAssignmentSeq[i];
    const odb_SequenceElement& elements = sa.region().elements();
    int nbElems = elements.size();
    const odb_Section& sect = sa.section();

    if (odb_isA(odb_CohesiveSection, sect)) {
      odb_CohesiveSection sect2 =
        odb_dynamicCast(odb_CohesiveSection, sect);
      matName = sect2.material().cStr();
      debug_print("getMatNames: found cohesive section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, nbElems);
    }
    else if (odb_isA(odb_HomogeneousSolidSection, sect)) {
      odb_HomogeneousSolidSection sect2 =
        odb_dynamicCast(odb_HomogeneousSolidSection, sect);
      matName = sect2.material().cStr();
      debug_print("getMatNames: found solid section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, nbElems);
    }
    else if (odb_isA(odb_BeamSection, sect)) {
      odb_BeamSection sect2 =
        odb_dynamicCast(odb_BeamSection, sect);
      matName = sect2.material().cStr();
      debug_print("getMatNames: found beam section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, nbElems);
    }
    else if (odb_isA(odb_HomogeneousShellSection, sect)) {
      odb_HomogeneousShellSection sect2 =
        odb_dynamicCast(odb_HomogeneousShellSection, sect);
      matName = sect2.material().cStr();
      debug_print("getMatNames: found shell section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, nbElems);
    }
    else if (odb_isA(odb_TrussSection, sect)) {
      odb_TrussSection sect2 =
        odb_dynamicCast(odb_TrussSection, sect);
      matName = sect2.material().cStr();
      debug_print("getMatNames: found truss section %s with material %s and"
                  " %d elements\n",
                  sa.sectionName().cStr(), matName, nbElems);
    }
    else {
      error("getMatNames: section type not implemented for section %s with"
            " %d elements.\n", sa.sectionName().cStr(), nbElems);
    }

    /* check that matName passes matNameFilters */
    matNameIsInFilterTuple = false;
    for (int j=0; j<matNameFilters.size(); ++j) {
      if (strncmp(matName, matNameFilters[j].c_str(),
                  matNameFilters[j].size()) == 0) {
        matNameIsInFilterTuple = true;
        break;
      }
    }

    if (matNameFilters.size() && !matNameIsInFilterTuple) {
      /* if filtertuple is specified and if current section type not in the
       * filter list then ignore this section assignment */
      debug_print("getMatNames: ignoring section assignment %s.\n",
                  sa.sectionName().cStr());
      // skip storing of new section assignment
      // and element label
      continue;
    }

    /* store matName in list */
    matNamesList.push_back(matName);

    /* get new matCode */
    matCode = nextNewMatCode;
    ++nextNewMatCode;

    /* loop over elements and update  */
    for (int i=0; i<nbElems; ++i) {
      elemNb = elements[i].label();
      if (std::find(filterset.begin(), filterset.end(), elemNb)
          != filterset.end()) {
        // store (only valid) matCode for element
        debug_print("getMatNames: storing element %d with matCode %d.\n",
                    elemNb, matCode);
        labels[i] = elemNb;
        matCodes[i] = matCode;
      }
    }
      }

  /* Send the result data through the pipe */
  debug_print("getMatNames: Will write %d elem numbers and matCodes.\n", nbElems);
  pipeWrite(nbElems);
  pipeWrite(labels, nbElems);
  pipeWrite(matCodes, nbElems);
  
  debug_print("getMatNames: Will write %d matNames.\n", sectionNamesList.size());
  pipeWrite((int) matNamesList.size());
  for (int i=0; i<matNamesList.size(); ++i) {
    pipeWrite(matNamesList[i]);
  }

  /* release data arrays */
  delete[] labels;
  delete[] matCodes;
  
  debug_print("getMatNames: Ending.\n");
  return 0;
}



/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
///////////////////////  old version     ////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


/* getMatNames
 *
 * This is to read the material assignment from the odb.
 *
 * If createFilterElset or setFilterElset has been called before only elements
 * from this set will be considered.
 *
 * input arguments: firstMatCode (int), nbSecTypeFilters (int),
 *     secTypeFilters (strings)
 *
 * firstMatCode determines the lowest (possible) material number in the result
 * corresponding to the first material name returned. It's typically =1.
 *
 * nbSecTypeFilters is the number of items in the following secTypeFilters
 * sequence of strings.
 * If you don't need this feature then pass zero for nbSecTypeFilters and skip
 * the secTypeFilters argument(s) alltogether.
 *
 * secTypeFilters are strings. Only section / material names starting with one
 * of those strings will be considered in the output. Elements whose section /
 * material name does not begin with one of these strings will be ignored.
 *
 *
 * Output:
 *
 * The first output of getMatNames is an array of element labels.
 * The second output is an array of the material numbers.
 * The third output is a sequence of material names.
 *
 * Here is the data in detail that will be written to the pipe:
 * - nb of items =:nbElems
 * - integer array of element labels (nbElems ints)
 * - integer array material numbers (nbElems ints)
 * - nb of material names =:nbMatNames
 * - a sequence of nbMatNames strings giving the material (section) names from
 *   the odb.
 *
 *
 * Suggestion to make this more economic:
 * Instead of iterating over elements and reading sectionCategory and then its
 * name: 
 */
int getMatNames() {
  // no matcodes allowed smaller or equal to this
  const int noMatCode(-9999999);
  
  /* get and check arguments */
  int firstMatCode = pipeReadInt();
  int nbSecTypeFilters = pipeReadInt();
  std::vector<std::string> secTypeFilters;
  for (int i=0; i<nbSecTypeFilters; ++i) {
    secTypeFilters.push_back(pipeReadString());
  }
  debug_print("This is getMatNames, firstMatCode %d, %d type filter strings.\n",
              firstMatCode, nbSecTypeFilters);
  if (firstMatCode<=noMatCode) {
    error("getMatNames: firstMatCode must be greater than %d. Got %d.\n",
          noMatCode, firstMatCode);
  }
  
  /* initialize/get Abaqus instance */
  if (ptInstance == NULL) {
    setPtInstance();
  }

  /* determine applicable elements / apply filterset */
  odb_Set* filterset = NULL;
  if (filterElsetName && filterElsetName->length()>0) {
    filterset =&(ptInstance->elementSets()[*(filterElsetName)]);
    debug_print("getMatNames: filter elset %s with %d elements.\n",
                filterset->name().cStr(), filterset->elements().size());
  }
  const odb_SequenceElement& elements =
    (filterset ? filterset->elements() : ptInstance->elements());

  /* initialize
   * - sectNameToMatCodeMap: map {sectName: matCode}
   * - sectionNamesList: */
  std::map<std::string, int> sectNameToMatCodeMap;
  std::map<std::string, int>::const_iterator sectNameToMatCodeSearch;
  int elemNb;
  int nextNewMatCode = firstMatCode;
  int matCode;
  bool secTypeIsInFilterTuple;
  
  // initialize result objects: arrays to store all data
  // because of it's (potential) size it's on the heap not the stack
  int nbElems = elements.size();
  int* labels = new int[nbElems];
  int* matCodes = new int[nbElems];
  std::vector<std::string> sectionNamesList;
  const char* sectName = NULL;
  const char* matName = NULL;
  
  /* loop over elements and update  */
  for (int i=0; i<nbElems; ++i) {
    elemNb = elements[i].label();
    sectName = elements[i].sectionCategory().name().cStr();

    sectNameToMatCodeSearch = sectNameToMatCodeMap.find(sectName);
    if (sectNameToMatCodeSearch == sectNameToMatCodeMap.end()) {
      // sectName not found in sectNameToMatCodeMap

      secTypeIsInFilterTuple = false;
      for (int j=0; j<secTypeFilters.size(); ++j) {
        if (strncmp(sectName, secTypeFilters[j].c_str(),
                    secTypeFilters[j].size()) == 0) {
          secTypeIsInFilterTuple = true;
          break;
        }
      }

      if (secTypeFilters.size() && !secTypeIsInFilterTuple) {
        /* if filtertuple is specified and if current section type not in the
         * filter list then ignore this section assignment */
        debug_print("getMatNames: ignoring section assignment %s.\n", sectName);
        sectNameToMatCodeMap[sectName] = noMatCode;
        // skip storing of new section assignment
        // and element label
        continue;
      }

      /* ... store new section name in sectionNamesList */
      debug_print("getMatNames: storing new section assignment %s.\n", sectName);
      // /////##### DEBUG

      // odb_SequenceSectionAssignment sectionAssignmentSeq = 
      //   ptInstance->sectionAssignments();  
      // int sects = sectionAssignmentSeq.size();
      // for (int s = 0; s < sects; ++s) {
      //   odb_SectionAssignment sa = sectionAssignmentSeq[s];
      //   odb_String sectionName = sa.sectionName();
      //   debug_print("getMatNames: sect %s \n", sectionName.cStr());
      //   odb_Set set = sa.region();
      //   const odb_SequenceElement& elements = set.elements();
      //   int size = elements.size();
      //   debug_print("getMatNames: elset name: %s \n", set.name().cStr());

      //   debug_print("getMatNames: section name: %s \n", sa.section().name().cStr());
      //   const odb_Section& sect = sa.section();
      //   if (odb_isA(odb_CohesiveSection, sect)) {
      //     odb_CohesiveSection cohesiveSection =
      //       odb_dynamicCast(odb_CohesiveSection, sect);
      //     debug_print("getMatNames: section material name: %s \n",
      //                 cohesiveSection.material().cStr());
      //   }
      //   debug_print("getMatNames: EOF section: %s \n", sa.section().name().cStr());
        
      // }
      
      // // if (strncmp(elements[i].sectionCategory().name().cStr(), "coh", 3)==0) {
      // //   debug_print("getMatNames: found cohesiove.\n");
      // //   sectionAssignments(i)
      // //   odb_section& .section

      // // }


      // /////##### DEBUG END
      sectionNamesList.push_back(sectName);

      /* get new matCode and store in sectNameToMatCodeMap */
      matCode = nextNewMatCode;
      ++nextNewMatCode;
      sectNameToMatCodeMap[sectName] = matCode;
    }
    else {
      // sectName found in sectNameToMatCodeMap
      matCode = sectNameToMatCodeSearch->second;
    }

    // store (only valid) matCode for element
    if (matCode != noMatCode) {
      debug_print("getMatNames: storing element %d with matCode %d.\n", elemNb, matCode);
      labels[i] = elemNb;
      matCodes[i] = matCode;
    }
  }

  /* Send the result data through the pipe */
  debug_print("getMatNames: Will write %d elem numbers and matCodes.\n", nbElems);
  pipeWrite(nbElems);
  pipeWrite(labels, nbElems);
  pipeWrite(matCodes, nbElems);
  
  debug_print("getMatNames: Will write %d matNames.\n", sectionNamesList.size());
  pipeWrite((int) sectionNamesList.size());
  for (int i=0; i<sectionNamesList.size(); ++i) {
    pipeWrite(sectionNamesList[i]);
  }

  /* release data arrays */
  delete[] labels;
  delete[] matCodes;
  
  debug_print("getMatNames: Ending.\n");
  return 0;
}
