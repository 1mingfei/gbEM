/*
 * @Author: chaomy
 * @Date:   2018-07-07 16:58:27
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-10-25 12:36:34
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
