/*
 * @Author: chaomy
 * @Date:   2017-12-31 16:04:15
 * @Last Modified by:  1mingfei
 * @Last Modified time: 2019-03-25
 */

#include "gbCnf.h"
#include "EMHome.h"

void EMHome::initParam() { readParam(); }

/**************************************************
 * read parameters
 **************************************************/
void EMHome::readParam() {
  ifstream fid(sparams["parfile"], std::ifstream::in);
  if (!fid.is_open()) cerr << " error opening " << sparams["parfile"] << endl;
  vector<string> segs;
  string buff;
  while (getline(fid, buff)) {
    segs.clear();
    split(buff, " ", segs);
    if (!segs[0].compare("Nconfigs"))
      iparams[segs[0]] = stoi(segs[1]);
    else if (!segs[0].compare("refFile"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("halfThick"))
      dparams[segs[0]] = stof(segs[1]);
    else if (!segs[0].compare("T"))
      dparams[segs[0]] = stof(segs[1]);
    else if (!segs[0].compare("N"))
      iparams[segs[0]] = stof(segs[1]);
    else if (!segs[0].compare("Rcut"))
      dparams[segs[0]] = stof(segs[1]);
    else if (!segs[0].compare("AlignFnx"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("elem0"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("elem1"))
      sparams[segs[0]] = segs[1];
  }
  fid.close();
}

/**************************************************
 * parse arguments
 **************************************************/
void EMHome::parseArgs(int argc, char* argv[]) {
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "--p") || !strcmp(argv[i], "-p"))
      sparams["parfile"] = string(argv[++i]);
    if (!strcmp(argv[i], "--i") || !strcmp(argv[i], "-i"))
      sparams["dmpfile"] = string(argv[++i]);
    if (!strcmp(argv[i], "--f") || !strcmp(argv[i], "-f"))
      sparams["potfile"] = string(argv[++i]);
  }
}
