#include "drugCompound.hpp"



namespace cardioxcomp{
  DrugCompound::DrugCompound()
  {
  }

  DrugCompound::DrugCompound(std::string name, double dose, double hillcoef){
    mName = name;
    mDose= dose;
    mHillCoef=hillcoef;
  }

  DrugCompound::~DrugCompound(){

  }

  std::string DrugCompound::getName()
  {
    return mName;
  }

  double DrugCompound::getDose()
  {
    return mDose;
  }
  double DrugCompound::getHillCoef()
  {
    return mHillCoef;
  }
  void DrugCompound::setChannelIC50(std::string channelsName,double IC50value)
  {
    mChannelsIC50values.insert(std::pair<std::string, double>(channelsName,IC50value));
  }
}
