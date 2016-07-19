
#ifndef FSTREAM_NAMD_H
#define FSTREAM_NAMD_H

#include <iostream>
#include <sstream>
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#else
#include <io.h>
#endif

class ofstream_namd : public std::ostringstream {

private:
  int fd;
  std::string fname;

public:

  ofstream_namd() : fd(0) { ; }

  explicit ofstream_namd(const char *_fname, std::ios_base::openmode _mode = std::ios_base::out) : fd(0) {
    open(_fname, _mode);
  }

  void open(const char *_fname, std::ios_base::openmode _mode = std::ios_base::out);

  bool is_open() const { return ! ! fd; }

  ofstream_namd& flush();

  void close();

  ~ofstream_namd() {
    if ( fd ) close();
  }

  void seekbegin() {
#if !defined(WIN32) || defined(__CYGWIN)
    lseek(fd, 0, SEEK_SET);
#else
    _lseek(fd, 0, SEEK_SET);
#endif
  }

  bool good() const { return true; }
  bool fail() const { return false; }
  bool bad() const { return false; }
  bool operator!() const { return false; }
  operator bool() const { return true; }
  void clear() { ; }

};

#endif  // FSTREAM_NAMD_H

