#pragma once

#include <iostream>
#include <string>
#include <Rcpp.h>
#include <cmath> //floor

class progressBar {
public:
  progressBar(unsigned short width, unsigned long max, bool init = true) :
  _printed(0),
  _width(width),
  _max(max),
  _done(false),
  _printed_end(false)
  {
    if(init)
      Rcpp::Rcout << "[";
    else
      _done = true;
  }
  void update(unsigned long current);
  bool printed_end() {return _printed_end;}
  void print_end();
private:
  unsigned short _printed;
  unsigned short _width;
  unsigned long _max;
  bool _done;
  bool _printed_end;
};

inline void progressBar::update(unsigned long current)
{
  float a = current;
  float b = _max;
  float r = a/b;
  float w = _width;
  unsigned short x = floor(w * r);
  while(_printed < x && ! _done)
  {
    unsigned int rperc = (r * 100);
    Rcpp::Rcout << "." << rperc << "%.";
    _printed++;
  }
  if(_printed >= _width && ! _done)
  {
    Rcpp::Rcout << "]\n";
    _done = true;
    _printed_end = true;
  }
}

inline void progressBar::print_end()
{
  Rcpp::Rcout << "]\n";
}
