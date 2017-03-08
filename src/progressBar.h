#pragma once

#include <iostream>
#include <string>
#include <Rcpp.h>
#include <cmath> //floor

class progressBar {
public:
  progressBar(unsigned short width, unsigned long max, bool init = true) :
  _printed(0),
  _p_chars(0),
  _width(width),
  _max(max),
  _done(false),
  _printed_end(false),
  _type(1)
  {
    if(init)
    {
      if(_type == 1)
        Rcpp::Rcout
        << "0%                                    100%\n"
        << "[----------------------------------------]\n";
        Rcpp::Rcout << "[";
        _p_chars++;
    }
    else
      _done = true;
  }
  void update(unsigned long current);
  bool printed_end() {return _printed_end;}
  void print_end();
private:
  unsigned short _printed;
  unsigned short _p_chars;
  unsigned short _width;
  unsigned long _max;
  bool _done;
  bool _printed_end;
  char _type;
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
    if(_type == 2)
    {
      unsigned int rperc = (r * 100);
      Rcpp::Rcout << "." << rperc << "%.";
      _printed++;
      _p_chars += 3;
      _p_chars += (rperc > 99 ? 3 : rperc > 9 ? 2 : 1);
      if(_p_chars > 40)
      {
        Rcpp::Rcout << "\n";
        _p_chars = 0;
      }
    }
    else if(_type == 1)
    {
      Rcpp::Rcout << '#';
      _printed++;
    }
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
  while(_printed < 40)
  {
    Rcpp::Rcout << "#";
    _printed++;
  }
  Rcpp::Rcout << "]\n";
}
