#include <iostream>
#include <iomanip>
#include <string>

template <typename T, typename T2>
void progressBar(std::ostream& out, T current, T2 max,
                 unsigned short width)
{
  std::string V = std::to_string(current);
  std::string M = std::to_string(max);
  float v = std::stof(V);
  float m = std::stof(M);

  float ratio = v/m;

  unsigned short bars = static_cast<unsigned short>(width * ratio);
  bars = bars > width ? 0 : bars;
  std::cout << "\r\t[";
  for(unsigned short i = 0; i < bars; i++)
  {
    std::cout << "#";
  }
  for(unsigned short i = 0; i < (width - bars); i++)
  {
    std::cout << " ";
  }
  std::cout << std::setprecision(2) <<
    std::setiosflags (std::ios::fixed) << "]\t" <<
      (ratio * 100) << "%\t  ";

}
