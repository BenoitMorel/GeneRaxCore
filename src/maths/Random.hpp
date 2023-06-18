#pragma once
#include <random>

class Random {
public:
  Random() = delete;
  static void setSeed(unsigned int seed);
  static int getInt(); 
  /**
   *  Return a uniform random value between min and max (included)
   */
  static int getInt(unsigned int min, unsigned int max); 
  static double getProba();
  static bool getBool();
private:
  static std::mt19937_64 _rng;
  static std::uniform_int_distribution<int> _unii;
  static std::uniform_real_distribution<double> _uniproba;
};
