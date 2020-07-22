#ifndef __CardioXcomp__drugcompound__
#define __CardioXcomp__drugcompound__


#include "luafile.hpp"
#include <map>
namespace cardioxcomp {
  class DrugCompound
  {
  protected:
    string mName;
    double mDose;
    std::map<std::string,double> mChannelsIC50values; // map from channels name to its corresponding IC50 value
    double mHillCoef;
  public:
    DrugCompound(std::string name, double dose, double hillcoef=1);
    DrugCompound();
    ~DrugCompound();
    std::string getName();
    double getDose();
    double getHillCoef();
    void setChannelIC50(std::string channelsName,double IC50value);
  };
}
#endif

