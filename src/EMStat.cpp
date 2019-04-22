/*
 * Author: 1mingfei 
 * Date:   2019-04-17
 * Purpose: functions for EMHome class 
 * estimate mean color
 * To-do: use a k-d Tree to better searching 1NN
 */

#include "gbCnf.h"

/*gather energies from all the structures
 *and calculate probability */
void EMHome::getProb(gbCnf& cnfModifier, double T, double N)
{
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "processing Boltzmann probability\n";
  
  int quotient = NI / nProcs;
  int remainder = NI % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;
  /*calculate total boltzmann energy Z*/
  for (int j = 0; j < nCycle; ++j)
  {
    double* subEngys = new double[nProcs];
    for (int i = (j * nProcs); i < ((j + 1) * nProcs); ++i)
    {
      if (i % nProcs != me) continue;
      boltzEngy = 0.0; //important to initialize here
      if (i >= NI) continue;

      c1 = move(cnfModifier.readLmpData(to_string(i) + ".txt"));
      engy = c1.engy / N;
      boltzEngy = exp(-engy/KB/T);
      cout << "read in config " << i << " on processor " << me 
           << " out of total " << nProcs  << " processors, energy: "
           << engy << " BE " << boltzEngy << endl;
      /*cnf:  oldEngy -- normalized energy
       *      engy    -- Boltzmann energy */
      c1.oldEngy = engy;
      c1.engy = boltzEngy;
      cnfs[i] = move(c1);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&boltzEngy, 1, MPI_DOUBLE, subEngys, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);

    if (me == 0)
    {
      cout << "subset energies: ";
      for (int i = 0; i < nProcs; ++i)
        cout << subEngys[i] << "   ";
      cout << endl;
      for (int i = 0; i < nProcs; ++i)
        sumEngy += subEngys[i];
      cout << "running total: " << sumEngy << endl;
    }
    delete [] subEngys;
  }
  MPI_Bcast(&sumEngy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  cout << "rank " << me << " sum energy: " << sumEngy << endl;
  /*calculate total boltzmann energy Z*/
  for (int i = 0; i < NI; ++i)
  {
    if (i % nProcs != me) continue;
    cnfs[i].oldEngy = cnfs[i].engy/sumEngy;
    cnfModifier.writeLmpDataDebug(cnfs[i], to_string(i) + ".txt");
  }
  /*Note: by now all the structures are stored in cnfs, with oldEngy represent
   *probability. */
}

/*calculate which is the nearest point in reference system and average pos 
 *then calculate type mean and std. */
void EMHome::estimateMean(gbCnf& cnfModifier)
{
  MPI_Barrier(MPI_COMM_WORLD);
  if (me == 0)
    std::cout << "calculating ave configuration based on energies...\n";
 
  c0 = move(cnfModifier.readLmpData(to_string(NI) + ".txt"));
  for (int i = 0; i < NI; ++i)
    cnfs[i] = move(cnfModifier.readLmpData(to_string(i) + ".txt"));

  MPI_Barrier(MPI_COMM_WORLD);
  int nAtom = c0.atoms.size();


  /*here based on atoms*/
  int quotient = nAtom / nProcs;
  int remainder = nAtom % nProcs;
  int nCycle = remainder ? (quotient + 1) : quotient;


  /*allocate largest memory for all nAtom (this is slightly larger than needed) 
   *Note that "nCycle * nProcs >= nAtoms"
   *all data store here after this*/
  double** data = new double* [nCycle * nProcs];
  for (int i = 0; i < (nCycle * nProcs); ++i)
    data[i] = new double [10];

  /*smallest buff for gathering in each cycle*/
  double* buffData = new double [10];

  for (int k = 0; k < nCycle; ++k)
  {
    for (int i = (k * nProcs); i < ((k + 1) * nProcs); ++i)
    {
      if ((me == 0) && (i % nProcs != 0))
      {
        MPI_Recv(&data[i][0], 10, MPI_DOUBLE, (i % nProcs), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      if (i % nProcs != me) continue;

      vector<double> dists(NI, 0.0);
      vector<double> posX(NI, 0.0);
      vector<double> posY(NI, 0.0);
      vector<double> posZ(NI, 0.0);
      vector<double> IDtyp(NI, 0.0);
      vector<double> prob(NI, 0.0);
      
      /*!!!!initialize buffer data here 
       *very important for unequal size of nProcs and leftover iterations*/
      for (int j = 0; j < 10; ++j)
        buffData[j] = 0.0;

      if (i >= nAtom) continue;

      for (int j = 0; j < NI; ++j)
      {
        /*get which atom in cnfs[j] is closest to atom[i] in c0. */
        double dist;
        int closest = cnfModifier.getNNID(c0.atoms[i], cnfs[j], dist);
        Atom tmpAtom = cnfs[j].atoms[closest];
        
        dists[j] = dist;

        /*correct relative position*/
        double tmp = tmpAtom.pst[X] - c0.atoms[i].pst[X];
        if (tmp > c0.length[X]/2.0) tmpAtom.pst[X] -= c0.length[X];
        if (tmp < -c0.length[X]/2.0) tmpAtom.pst[X] += c0.length[X];
        tmp = tmpAtom.pst[Y] - c0.atoms[i].pst[Y];
        if (tmp > c0.length[Y]/2.0) tmpAtom.pst[Y] -= c0.length[Y];
        if (tmp < -c0.length[Y]/2.0) tmpAtom.pst[Y] += c0.length[Y];
        tmp = tmpAtom.pst[Z] - c0.atoms[i].pst[Z];
        if (tmp > c0.length[Z]/2.0) tmpAtom.pst[Z] -= c0.length[Z];
        if (tmp < -c0.length[Z]/2.0) tmpAtom.pst[Z] += c0.length[Z];

        posX[j] = tmpAtom.pst[X];
        posY[j] = tmpAtom.pst[Y];
        posZ[j] = tmpAtom.pst[Z];
        double tmpType = (tmpAtom.tp == 1) ? 0.0 : 1.0;
        IDtyp[j] = tmpType;
        prob[j] = cnfs[j].oldEngy;
      }

      /*get some stats*/
      double meanDist = getProbMean(dists, prob);
      double meanX = getProbMean(posX, prob);
      double meanY = getProbMean(posY, prob);
      double meanZ = getProbMean(posZ, prob);
      double meanType = getProbMean(IDtyp, prob);
      double stdDist = getStdDevProb(dists, prob, meanDist);
      double stdX = getStdDevProb(posX, prob, meanX);
      double stdY = getStdDevProb(posY, prob, meanY);
      double stdZ = getStdDevProb(posZ, prob, meanZ);
      double stdTp = getStdDevProb(IDtyp, prob, meanType);

      cout << std::setprecision(6) << "on rank " << me 
           << " atomID " << c0.atoms[i].id
           << " meanX: " << meanX << " meanY: " << meanY << " meanZ: " << meanZ
           << " stdX: " << stdX << " stdY: " << stdY << " stdZ: " << stdZ
           << " posStd: " << stdDist
           << " typeMean: " << meanType 
           << " typeStd: " << stdTp << endl;

      /*pack data and send*/
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

      if (me != 0)
        MPI_Send(&buffData[0], 10, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      else
      {
        for (int ii = 0; ii < 10; ++ii)
          data[k * nProcs + i % nProcs][ii] = buffData[ii];
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (me == 0)
  {
    //convert data to vector dtype
    vector<vector<double>> VVData;
    for (int i = 0; i < (nCycle * nProcs); ++i)
    {
      vector<double> tmpVec;
#ifdef DEBUG
      cout << i << " ";
#endif
      for (int j = 0; j < 10; ++j)
      {
#ifdef DEBUG
        cout  << std::setprecision(5) << data[i][j] << " ";
#endif
        tmpVec.push_back(data[i][j]);
      }
#ifdef DEBUG
      cout << endl;
#endif
      VVData.push_back(tmpVec);
    }
    /*update c0 position informations*/
    for (int i = 0; i < c0.atoms.size(); ++i)
    {
      c0.atoms[i].pst[X] = data[i][0];
      c0.atoms[i].pst[Y] = data[i][2];
      c0.atoms[i].pst[Z] = data[i][4];
    }
    cnfModifier.cnvpst2prl(c0);
    cnfModifier.writeCfgData(c0, VVData, "final.cfg");
    cout << "final averaged structure: final.cfg" << endl;
  }


  /*free smallest buffer*/
  delete [] buffData;
 
  /*free largest data set*/
  for (int i = 0; i < (nCycle * nProcs); ++i)
    delete [] data[i];
  delete [] data;

}

/*get which atom in c1 is closest to atom[i] in c0.
 *atm is our reference atom, cnf is the sturcture to search
 *calculated dist is returned afterwards. */
int EMHome::gbCnf::getNNID(Atom& atm, Config& cnf, double& dist)
{
  dist = std::numeric_limits<double>::max();
  int res;
  for (int i = 0; i < cnf.atoms.size(); ++i)
  {
    double tmpDist = calDist(cnf.length, cnf.atoms[i], atm);
    if (tmpDist < dist)
    {
      dist = tmpDist;
      res = i;
    }
  }
  return res;
}
