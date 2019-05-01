/*
 * Author: 1mingfei 
 * Date:   2019-04-14
 * Purpose:user can define their own cost function to 
 * evaluate gb location in y
 */

#include "gbCnf.h"

int meanScore(const vector<vector<double>>& Score, \
                     const vector<double>& Loc)
{
  int bestIndex = Score.size();
  double bestScore = std::numeric_limits<double>::max(); //FCC100
  for (int i = 0; i < Score.size(); ++i)
  {
    if (Loc[i] > 0.01)
    {
      if (meanV(Score[i]) < bestScore) //FCC100
      {
        bestScore = meanV(Score[i]); //FCC100
        bestIndex = i;
      }
    }
#ifdef DEBUG
    cout << i << " " << meanV(Score[i]) << " " << Loc[i] << endl;
#endif
  }
  return bestIndex;
}
int stdScore(const vector<vector<double>>& Score, \
                     const vector<double>& Loc)
{
  int bestIndex = Score.size();
  double bestScore = 0.0; //1210TB
  for (int i = 0; i < Score.size(); ++i)
  {
    if (Loc[i] > 0.01)
    {
      if (stddev(Score[i]) > bestScore) //1210TB
      {
        bestScore = stddev(Score[i]); //1210TB
        bestIndex = i;
      }
    }
#ifdef DEBUG
    cout << i << " " << stddev(Score[i]) << " " << Loc[i] << endl;
#endif
  }
  return bestIndex;
}
