/*
 * @Author: chaomy
 * @Date:   2018-07-07 16:58:27
 * @Last Modified by:  1mingfei 
 * @Last Modified time: 2019-4-17 05:32:23
 */

#include "gbCnf.h"
/**************************************************
 * write lammps data files
 **************************************************/
/*
  #lmp data config
  912 atoms
  3 atom types
  0.000000  31.772283 xlo xhi
  0.000000  200.000000 ylo yhi
  0.000000  11.071269 zlo zhi
  0.000000  0.000000  0.000000 xy xz yz
  Atoms

*/
void EMHome::gbCnf::writeLmpData(Config& c, string fnm = "out.lmp.init") {
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "# lmp data config " << endl;
  ofs << (int)c.natoms << " atoms" << endl;
  //ofs << (int)c.ntypes << " atom types " << endl;
  ofs << 4 << " atom types " << endl;
  ofs << c.cell[0] << " " << c.cell[3] << " xlo xhi" << endl;
  ofs << c.cell[1] << " " << c.cell[4] << " ylo yhi" << endl;
  ofs << c.cell[2] << " " << c.cell[5] << " zlo zhi" << endl;
  ofs << c.cell[6] << " " << c.cell[7] << " " << c.cell[8] << " xy xz yz"
      << endl;
  ofs << "Atoms" << endl << endl;
  for (int i = 0; i < c.natoms; ++i) {
    auto&& a = c.atoms[i];
    ofs << a.id + 1 << " " << a.tp << " " << a.pst[0] << " " << a.pst[1]
        << " " << a.pst[2] << endl;
  }
}
void EMHome::gbCnf::writeLmpDataDebug(Config& c, string fnm = "out.lmp.init") {
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "# " << setprecision(14) << c.oldEngy << " " << c.engy << endl;
  ofs << (int)c.natoms << " atoms" << endl;
  //ofs << (int)c.ntypes << " atom types " << endl;
  ofs << 4 << " atom types " << endl;
  ofs << c.cell[0] << " " << c.cell[3] << " xlo xhi" << endl;
  ofs << c.cell[1] << " " << c.cell[4] << " ylo yhi" << endl;
  ofs << c.cell[2] << " " << c.cell[5] << " zlo zhi" << endl;
  ofs << c.cell[6] << " " << c.cell[7] << " " << c.cell[8] << " xy xz yz"
      << endl;
  ofs << "Atoms" << endl << endl;
  for (int i = 0; i < c.natoms; ++i) {
    auto&& a = c.atoms[i];
    ofs << a.id + 1 << " " << a.tp << " " << a.pst[0] << " " << a.pst[1]
        << " " << a.pst[2] << endl;
  }
}

/* write cfg data files*/
/*
Number of particles = 16200
A = 4.37576470588235 Angstrom (basic length-scale)
H0(1,1) = 127.5 A
H0(1,2) = 0 A
H0(1,3) = 0 A
H0(2,1) = 0 A
H0(2,2) = 119.501132067411 A
H0(2,3) = 0 A
H0(3,1) = 0 A
H0(3,2) = 0 A
H0(3,3) = 3 A
.NO_VELOCITY.
entry_count = 9
auxiliary[0] = kine [reduced unit]
auxiliary[1] = pote [reduced unit]
auxiliary[2] = s11 [reduced unit]
auxiliary[3] = s22 [reduced unit]
auxiliary[4] = s12 [reduced unit]
auxiliary[5] = hydro [reduced unit]
1.000000
Ar
0.0016667 0.00616 0.5 0 -2 -1.9431e-13 -2.9917e-13 1.6811e-13 -2.4674e-13 

      buffData[0] = meanX;
      buffData[1] = stdX;
      buffData[2] = meanY;
      buffData[3] = stdY;
      buffData[4] = meanZ;
      buffData[5] = stdZ;
      buffData[6] = meanDist;
      buffData[7] = stdDist;
      buffData[8] = meanType;
      buffData[9] = stdTp;

 */
void EMHome::gbCnf::writeCfgData(const Config& c, 
                                 const vector<vector<double>>& data,
                                 const vector<string>& elems, 
                                 string fnm = "out.lmp.init")
{
  ofstream ofs(fnm, std::ofstream::out);
  ofs << "Number of particles =  " << (int)c.atoms.size() << endl;
  ofs << "A = 1.0 Angstrom (basic length-scale)" << endl;
  ofs << "H0(1,1) = " << c.bvx[X] << " A" << endl;
  ofs << "H0(1,2) = " << c.bvx[Y] << " A" << endl;
  ofs << "H0(1,3) = " << c.bvx[Z] << " A" << endl;
  ofs << "H0(2,1) = " << c.bvy[X] << " A" << endl;
  ofs << "H0(2,2) = " << c.bvy[Y] << " A" << endl;
  ofs << "H0(2,3) = " << c.bvy[Z] << " A" << endl;
  ofs << "H0(3,1) = " << c.bvz[X] << " A" << endl;
  ofs << "H0(3,2) = " << c.bvz[Y] << " A" << endl;
  ofs << "H0(3,3) = " << c.bvz[Z] << " A" << endl;
  ofs << ".NO_VELOCITY." << endl;
  ofs << "entry_count = 10" << endl;
  ofs << "auxiliary[0] = stdX" << endl;
  ofs << "auxiliary[1] = stdY" << endl;
  ofs << "auxiliary[2] = stdZ" << endl;
  ofs << "auxiliary[3] = mean_distance" << endl;
  ofs << "auxiliary[4] = std_distance" << endl;
  ofs << "auxiliary[5] = mean_type" << endl;
  ofs << "auxiliary[6] = std_type" << endl;
  double mass0 = findMass(elems[0]);
  double mass1 = findMass(elems[1]);
  for (int i = 0; i < c.atoms.size(); ++i) 
  {
    auto&& a = c.atoms[i];
    if (a.tp == 1)
      ofs << mass0 << endl << elems[0] << endl;
    else if (a.tp == 4)
      ofs << mass1 << endl << elems[1] << endl;
    ofs << a.prl[X] << " " << a.prl[Y] << " " << a.prl[Z] << " ";
    if (data.size())
      ofs << data[i][1] << " " << data[i][3] << " " << data[i][5] << " "
          << data[i][6] << " " << data[i][7] << " " << data[i][8] << " "
          << data[i][9] << endl; 
  }
}

double EMHome::gbCnf::findMass(string x)
{
  /*this "it" is indeed a iterator*/
  auto it = std::find(element.begin(), element.end(), x);
  int index = std::distance(element.begin(), it);
  return mass[index];
}
