#ifndef SIMPLE_TIMINGS_H
#define SIMPLE_TIMINGS_H

#include <string>
#include <vector>

#include <Partition.h>
#include <rru_mpi.h>

#define NAMEWIDTH 5

using namespace std;

class SimpleTimings {

 public:

  typedef int TimerRef;

  SimpleTimings();
  ~SimpleTimings();

  static TimerRef getTimer(const char *nm);
  static void startTimer(TimerRef tr);
  static void stopTimer(TimerRef tr);
  static void timerStats(TimerRef, int nameWidth=NAMEWIDTH);
  static void stopTimerStats(TimerRef, int nameWidth=NAMEWIDTH);
  static void fakeTimer(TimerRef tr, double time);
  static void fakeTimerStats(TimerRef tr,double time,int nameWidth=NAMEWIDTH);
  static void accumStats(int nameWidth=NAMEWIDTH);
  static void timingStats(double t, string name, int nameWidth=NAMEWIDTH);

 private:

  static vector <string> names;
  static vector <double> start;
  static vector <double> stop;
  static vector <double> accum;
};

#endif
