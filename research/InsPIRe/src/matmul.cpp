/*
See https://github.com/ahenzinger/simplepir for original.

MIT License

Copyright (c) 2022, Alexandra Henzinger

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// #define PROFILE

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

#ifdef PROFILE
  #include <sstream>
  #include <iostream>
  #include <functional>
  #include <fcntl.h>
  #include <sys/stat.h>
  #include <sys/wait.h>
  #include <sys/types.h>
  #include <unistd.h>
  #include <signal.h>

  struct System
  {
      static void profile(const std::string& name,std::function<void()> body) {
          std::string filename = name.find(".data") == std::string::npos ? (name + ".data") : name;

          // Launch profiler
          pid_t pid;
          std::stringstream s;
          s << getpid();
          pid = fork();
          if (pid == 0) {
              exit(execl("/usr/bin/perf","perf","stat","-e","cache-references,cache-misses,cycles,instructions,branches,faults,migrations,l1d.replacement,l2_rqsts.all_demand_miss,cycle_activity.stalls_l3_miss","-o",filename.c_str(),"-p",s.str().c_str(),nullptr));
          }

          // Run body
          body();

          // Kill profiler  
          kill(pid,SIGINT);
          waitpid(pid,nullptr,0);
      }

      static void profile(std::function<void()> body) {
          profile("perf.data",body);
      }
  };
#endif

typedef uint32_t Elem;

extern "C" {
  void matMulVecPacked(uint32_t *out, const uint32_t *a, const uint32_t *b,
      size_t aRows, size_t aCols);

  void matMulVecPacked2(Elem *out, const Elem *a, const Elem *b_full,
      size_t aRows, size_t aCols);

  void matMulVecPacked4(Elem *out, const Elem *a, const Elem *b_full,
      size_t aRows, size_t aCols);

  void matMulVecPacked6(Elem *out, const Elem *a, const Elem *b_full,
      size_t aRows, size_t aCols);

  void matMulVecPacked8(Elem *out, const Elem *a, const Elem *b_full,
      size_t aRows, size_t aCols);

  void matMulVecPacked8Alt(Elem *out, const Elem *a, const Elem *b_full,
      size_t aRows, size_t aCols);

  void matMulVecPacked8Orig(Elem *out, const Elem *a, const Elem *b_full,
      size_t aRows, size_t aCols);
}

// Hard-coded, to allow for compiler optimizations:
#define COMPRESSION 4
#define BASIS       8
#define BASIS2      16
#define BASIS3      24
#define MASK        0xff

void matMulVecPacked(uint32_t *out, const uint32_t *a, const uint32_t *b,
    size_t aRows, size_t aCols)
{
  uint32_t db, db2, db3, db4, db5, db6, db7, db8;
  uint32_t val, val2, val3, val4, val5, val6, val7, val8;
  uint32_t tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  size_t index = 0;
  size_t index2;

  for (size_t i = 0; i < aRows; i += 8) {
    tmp  = 0;
    tmp2 = 0;
    tmp3 = 0;
    tmp4 = 0;
    tmp5 = 0;
    tmp6 = 0;
    tmp7 = 0;
    tmp8 = 0;

    index2 = 0;
    for (size_t j = 0; j < aCols; j++) {
      db  = a[index];
      db2 = a[index+1*aCols];
      db3 = a[index+2*aCols];
      db4 = a[index+3*aCols];
      db5 = a[index+4*aCols];
      db6 = a[index+5*aCols];
      db7 = a[index+6*aCols];
      db8 = a[index+7*aCols];

      val  = db & MASK;
      val2 = db2 & MASK;
      val3 = db3 & MASK;
      val4 = db4 & MASK;
      val5 = db5 & MASK;
      val6 = db6 & MASK;
      val7 = db7 & MASK;
      val8 = db8 & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      index2 += 1;

      val  = (db >> BASIS) & MASK;
      val2 = (db2 >> BASIS) & MASK;
      val3 = (db3 >> BASIS) & MASK;
      val4 = (db4 >> BASIS) & MASK;
      val5 = (db5 >> BASIS) & MASK;
      val6 = (db6 >> BASIS) & MASK;
      val7 = (db7 >> BASIS) & MASK;
      val8 = (db8 >> BASIS) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      index2 += 1;

      val  = (db >> BASIS2) & MASK;
      val2 = (db2 >> BASIS2) & MASK;
      val3 = (db3 >> BASIS2) & MASK;
      val4 = (db4 >> BASIS2) & MASK;
      val5 = (db5 >> BASIS2) & MASK;
      val6 = (db6 >> BASIS2) & MASK;
      val7 = (db7 >> BASIS2) & MASK;
      val8 = (db8 >> BASIS2) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      index2 += 1;

      val  = (db >> BASIS3) & MASK;
      val2 = (db2 >> BASIS3) & MASK;
      val3 = (db3 >> BASIS3) & MASK;
      val4 = (db4 >> BASIS3) & MASK;
      val5 = (db5 >> BASIS3) & MASK;
      val6 = (db6 >> BASIS3) & MASK;
      val7 = (db7 >> BASIS3) & MASK;
      val8 = (db8 >> BASIS3) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      index2 += 1;
      index += 1;
    }
    out[i]   += tmp;
    out[i+1] += tmp2;
    out[i+2] += tmp3;
    out[i+3] += tmp4;
    out[i+4] += tmp5;
    out[i+5] += tmp6;
    out[i+6] += tmp7;
    out[i+7] += tmp8;
    index += aCols*7;
  }
}

void matMulVecPacked2(Elem *out, const Elem *a, const Elem *b_full,
    size_t aRows, size_t aCols)
{
  // in this variant, b has 2 columns

  Elem db, db2, db3, db4, db5, db6, db7, db8;
  Elem val, val2, val3, val4, val5, val6, val7, val8;
  Elem tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  Elem x_tmp, x_tmp2, x_tmp3, x_tmp4, x_tmp5, x_tmp6, x_tmp7, x_tmp8;
  size_t index = 0;
  size_t index2;

  const Elem *b = b_full;
  const Elem *b2 = b_full + aCols*COMPRESSION;

  for (size_t i = 0; i < aRows; i += 8) {
    tmp  = 0;
    tmp2 = 0;
    tmp3 = 0;
    tmp4 = 0;
    tmp5 = 0;
    tmp6 = 0;
    tmp7 = 0;
    tmp8 = 0;

    x_tmp  = 0;
    x_tmp2 = 0;
    x_tmp3 = 0;
    x_tmp4 = 0;
    x_tmp5 = 0;
    x_tmp6 = 0;
    x_tmp7 = 0;
    x_tmp8 = 0;

    index2 = 0;
    for (size_t j = 0; j < aCols; j++) {
      db  = a[index];
      db2 = a[index+1*aCols];
      db3 = a[index+2*aCols];
      db4 = a[index+3*aCols];
      db5 = a[index+4*aCols];
      db6 = a[index+5*aCols];
      db7 = a[index+6*aCols];
      db8 = a[index+7*aCols];

      val  = db & MASK;
      val2 = db2 & MASK;
      val3 = db3 & MASK;
      val4 = db4 & MASK;
      val5 = db5 & MASK;
      val6 = db6 & MASK;
      val7 = db7 & MASK;
      val8 = db8 & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      index2 += 1;

      val  = (db >> BASIS) & MASK;
      val2 = (db2 >> BASIS) & MASK;
      val3 = (db3 >> BASIS) & MASK;
      val4 = (db4 >> BASIS) & MASK;
      val5 = (db5 >> BASIS) & MASK;
      val6 = (db6 >> BASIS) & MASK;
      val7 = (db7 >> BASIS) & MASK;
      val8 = (db8 >> BASIS) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      index2 += 1;

      val  = (db >> BASIS2) & MASK;
      val2 = (db2 >> BASIS2) & MASK;
      val3 = (db3 >> BASIS2) & MASK;
      val4 = (db4 >> BASIS2) & MASK;
      val5 = (db5 >> BASIS2) & MASK;
      val6 = (db6 >> BASIS2) & MASK;
      val7 = (db7 >> BASIS2) & MASK;
      val8 = (db8 >> BASIS2) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      index2 += 1;

      val  = (db >> BASIS3) & MASK;
      val2 = (db2 >> BASIS3) & MASK;
      val3 = (db3 >> BASIS3) & MASK;
      val4 = (db4 >> BASIS3) & MASK;
      val5 = (db5 >> BASIS3) & MASK;
      val6 = (db6 >> BASIS3) & MASK;
      val7 = (db7 >> BASIS3) & MASK;
      val8 = (db8 >> BASIS3) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      index2 += 1;

      index += 1;
    }
    out[2*i]    += tmp;
    out[2*i+1]  += x_tmp;
    out[2*i+2]  += tmp2;
    out[2*i+3]  += x_tmp2;
    out[2*i+4]  += tmp3;
    out[2*i+5]  += x_tmp3;
    out[2*i+6]  += tmp4;
    out[2*i+7]  += x_tmp4;
    out[2*i+8]  += tmp5;
    out[2*i+9]  += x_tmp5;
    out[2*i+10] += tmp6;
    out[2*i+11] += x_tmp6;
    out[2*i+12] += tmp7;
    out[2*i+13] += x_tmp7;
    out[2*i+14] += tmp8;
    out[2*i+15] += x_tmp8;
    index += aCols*7;
  }
}

void trueMatMulVecPacked4(Elem *out, const Elem *a, const Elem *b_full,
    size_t aRows, size_t aCols)
{
  // in this variant, b has 4 columns

  Elem db, db2, db3, db4, db5, db6, db7, db8;
  Elem val, val2, val3, val4, val5, val6, val7, val8;
  Elem tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  Elem x_tmp, x_tmp2, x_tmp3, x_tmp4, x_tmp5, x_tmp6, x_tmp7, x_tmp8;
  Elem y_tmp, y_tmp2, y_tmp3, y_tmp4, y_tmp5, y_tmp6, y_tmp7, y_tmp8;
  Elem z_tmp, z_tmp2, z_tmp3, z_tmp4, z_tmp5, z_tmp6, z_tmp7, z_tmp8;
  size_t index = 0;
  size_t index2;

  const Elem *b = b_full;
  const Elem *b2 = b_full + aCols*COMPRESSION;
  const Elem *b3 = b_full + 2*aCols*COMPRESSION;
  const Elem *b4 = b_full + 3*aCols*COMPRESSION;

  for (size_t i = 0; i < aRows; i += 8) {
    tmp  = 0;
    tmp2 = 0;
    tmp3 = 0;
    tmp4 = 0;
    tmp5 = 0;
    tmp6 = 0;
    tmp7 = 0;
    tmp8 = 0;

    x_tmp  = 0;
    x_tmp2 = 0;
    x_tmp3 = 0;
    x_tmp4 = 0;
    x_tmp5 = 0;
    x_tmp6 = 0;
    x_tmp7 = 0;
    x_tmp8 = 0;

    y_tmp  = 0;
    y_tmp2 = 0;
    y_tmp3 = 0;
    y_tmp4 = 0;
    y_tmp5 = 0;
    y_tmp6 = 0;
    y_tmp7 = 0;
    y_tmp8 = 0;

    z_tmp  = 0;
    z_tmp2 = 0;
    z_tmp3 = 0;
    z_tmp4 = 0;
    z_tmp5 = 0;
    z_tmp6 = 0;
    z_tmp7 = 0;
    z_tmp8 = 0;

    index2 = 0;
    for (size_t j = 0; j < aCols; j++) {
      db  = a[index];
      db2 = a[index+1*aCols];
      db3 = a[index+2*aCols];
      db4 = a[index+3*aCols];
      db5 = a[index+4*aCols];
      db6 = a[index+5*aCols];
      db7 = a[index+6*aCols];
      db8 = a[index+7*aCols];

      val  = db & MASK;
      val2 = db2 & MASK;
      val3 = db3 & MASK;
      val4 = db4 & MASK;
      val5 = db5 & MASK;
      val6 = db6 & MASK;
      val7 = db7 & MASK;
      val8 = db8 & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];

      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];

      index2 += 1;

      val  = (db >> BASIS) & MASK;
      val2 = (db2 >> BASIS) & MASK;
      val3 = (db3 >> BASIS) & MASK;
      val4 = (db4 >> BASIS) & MASK;
      val5 = (db5 >> BASIS) & MASK;
      val6 = (db6 >> BASIS) & MASK;
      val7 = (db7 >> BASIS) & MASK;
      val8 = (db8 >> BASIS) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];

      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];

      index2 += 1;

      val  = (db >> BASIS2) & MASK;
      val2 = (db2 >> BASIS2) & MASK;
      val3 = (db3 >> BASIS2) & MASK;
      val4 = (db4 >> BASIS2) & MASK;
      val5 = (db5 >> BASIS2) & MASK;
      val6 = (db6 >> BASIS2) & MASK;
      val7 = (db7 >> BASIS2) & MASK;
      val8 = (db8 >> BASIS2) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];

      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];

      index2 += 1;

      val  = (db >> BASIS3) & MASK;
      val2 = (db2 >> BASIS3) & MASK;
      val3 = (db3 >> BASIS3) & MASK;
      val4 = (db4 >> BASIS3) & MASK;
      val5 = (db5 >> BASIS3) & MASK;
      val6 = (db6 >> BASIS3) & MASK;
      val7 = (db7 >> BASIS3) & MASK;
      val8 = (db8 >> BASIS3) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];

      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];

      index2 += 1;
      
      index += 1;
    }
    out[4*i]    += tmp;
    out[4*i+1]  += x_tmp;
    out[4*i+2]  += y_tmp;
    out[4*i+3]  += z_tmp;
    out[4*i+4]  += tmp2;
    out[4*i+5]  += x_tmp2;
    out[4*i+6]  += y_tmp2;
    out[4*i+7]  += z_tmp2;
    out[4*i+8]  += tmp3;
    out[4*i+9]  += x_tmp3;
    out[4*i+10] += y_tmp3;
    out[4*i+11] += z_tmp3;
    out[4*i+12] += tmp4;
    out[4*i+13] += x_tmp4;
    out[4*i+14] += y_tmp4;
    out[4*i+15] += z_tmp4;
    out[4*i+16] += tmp5;
    out[4*i+17] += x_tmp5;
    out[4*i+18] += y_tmp5;
    out[4*i+19] += z_tmp5;
    out[4*i+20] += tmp6;
    out[4*i+21] += x_tmp6;
    out[4*i+22] += y_tmp6;
    out[4*i+23] += z_tmp6;
    out[4*i+24] += tmp7;
    out[4*i+25] += x_tmp7;
    out[4*i+26] += y_tmp7;
    out[4*i+27] += z_tmp7;
    out[4*i+28] += tmp8;
    out[4*i+29] += x_tmp8;
    out[4*i+30] += y_tmp8;
    out[4*i+31] += z_tmp8;
    index += aCols*7;
  }
}

void matMulVecPacked4(Elem *out, const Elem *a, const Elem *b_full,
    size_t aRows, size_t aCols)
{
  #ifdef PROFILE
    System::profile("matmul4", [&]() {
  #endif

  trueMatMulVecPacked4(out, a, b_full, aRows, aCols);
  
  #ifdef PROFILE
    });
  #endif
}

void matMulVecPacked6(Elem *__restrict__ out, const Elem *__restrict__ a, const Elem *__restrict__ b_full,
    size_t aRows, size_t aCols)
{
  // in this variant, b has 6 columns

  Elem db, db2, db3, db4, db5, db6, db7, db8;
  Elem val, val2, val3, val4, val5, val6, val7, val8;
  Elem tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  Elem x_tmp, x_tmp2, x_tmp3, x_tmp4, x_tmp5, x_tmp6, x_tmp7, x_tmp8;
  Elem y_tmp, y_tmp2, y_tmp3, y_tmp4, y_tmp5, y_tmp6, y_tmp7, y_tmp8;
  Elem z_tmp, z_tmp2, z_tmp3, z_tmp4, z_tmp5, z_tmp6, z_tmp7, z_tmp8;
  Elem r_tmp, r_tmp2, r_tmp3, r_tmp4, r_tmp5, r_tmp6, r_tmp7, r_tmp8;
  Elem s_tmp, s_tmp2, s_tmp3, s_tmp4, s_tmp5, s_tmp6, s_tmp7, s_tmp8;
  size_t index = 0;
  size_t index2;

  const Elem *b = b_full;
  const Elem *b2 = b_full + aCols*COMPRESSION;
  const Elem *b3 = b_full + 2*aCols*COMPRESSION;
  const Elem *b4 = b_full + 3*aCols*COMPRESSION;
  const Elem *b5 = b_full + 4*aCols*COMPRESSION;
  const Elem *b6 = b_full + 5*aCols*COMPRESSION;

  for (size_t i = 0; i < aRows; i += 8) {
    tmp  = 0;
    tmp2 = 0;
    tmp3 = 0;
    tmp4 = 0;
    tmp5 = 0;
    tmp6 = 0;
    tmp7 = 0;
    tmp8 = 0;

    x_tmp  = 0;
    x_tmp2 = 0;
    x_tmp3 = 0;
    x_tmp4 = 0;
    x_tmp5 = 0;
    x_tmp6 = 0;
    x_tmp7 = 0;
    x_tmp8 = 0;

    y_tmp  = 0;
    y_tmp2 = 0;
    y_tmp3 = 0;
    y_tmp4 = 0;
    y_tmp5 = 0;
    y_tmp6 = 0;
    y_tmp7 = 0;
    y_tmp8 = 0;

    z_tmp  = 0;
    z_tmp2 = 0;
    z_tmp3 = 0;
    z_tmp4 = 0;
    z_tmp5 = 0;
    z_tmp6 = 0;
    z_tmp7 = 0;
    z_tmp8 = 0;

    r_tmp  = 0;
    r_tmp2 = 0;
    r_tmp3 = 0;
    r_tmp4 = 0;
    r_tmp5 = 0;
    r_tmp6 = 0;
    r_tmp7 = 0;
    r_tmp8 = 0;

    s_tmp  = 0;
    s_tmp2 = 0;
    s_tmp3 = 0;
    s_tmp4 = 0;
    s_tmp5 = 0;
    s_tmp6 = 0;
    s_tmp7 = 0;
    s_tmp8 = 0;


    index2 = 0;
    for (size_t j = 0; j < aCols; j++) {
      db  = a[index];
      db2 = a[index+1*aCols];
      db3 = a[index+2*aCols];
      db4 = a[index+3*aCols];
      db5 = a[index+4*aCols];
      db6 = a[index+5*aCols];
      db7 = a[index+6*aCols];
      db8 = a[index+7*aCols];

      val  = db & MASK;
      val2 = db2 & MASK;
      val3 = db3 & MASK;
      val4 = db4 & MASK;
      val5 = db5 & MASK;
      val6 = db6 & MASK;
      val7 = db7 & MASK;
      val8 = db8 & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];
      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];
      r_tmp  += val*b5[index2];
      r_tmp2 += val2*b5[index2];
      r_tmp3 += val3*b5[index2];
      r_tmp4 += val4*b5[index2];
      r_tmp5 += val5*b5[index2];
      r_tmp6 += val6*b5[index2];
      r_tmp7 += val7*b5[index2];
      r_tmp8 += val8*b5[index2];
      s_tmp  += val*b6[index2];
      s_tmp2 += val2*b6[index2];
      s_tmp3 += val3*b6[index2];
      s_tmp4 += val4*b6[index2];
      s_tmp5 += val5*b6[index2];
      s_tmp6 += val6*b6[index2];
      s_tmp7 += val7*b6[index2];
      s_tmp8 += val8*b6[index2];
      index2 += 1;

      val  = (db >> BASIS) & MASK;
      val2 = (db2 >> BASIS) & MASK;
      val3 = (db3 >> BASIS) & MASK;
      val4 = (db4 >> BASIS) & MASK;
      val5 = (db5 >> BASIS) & MASK;
      val6 = (db6 >> BASIS) & MASK;
      val7 = (db7 >> BASIS) & MASK;
      val8 = (db8 >> BASIS) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];
      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];
      r_tmp  += val*b5[index2];
      r_tmp2 += val2*b5[index2];
      r_tmp3 += val3*b5[index2];
      r_tmp4 += val4*b5[index2];
      r_tmp5 += val5*b5[index2];
      r_tmp6 += val6*b5[index2];
      r_tmp7 += val7*b5[index2];
      r_tmp8 += val8*b5[index2];
      s_tmp  += val*b6[index2];
      s_tmp2 += val2*b6[index2];
      s_tmp3 += val3*b6[index2];
      s_tmp4 += val4*b6[index2];
      s_tmp5 += val5*b6[index2];
      s_tmp6 += val6*b6[index2];
      s_tmp7 += val7*b6[index2];
      s_tmp8 += val8*b6[index2];
      index2 += 1;

      val  = (db >> BASIS2) & MASK;
      val2 = (db2 >> BASIS2) & MASK;
      val3 = (db3 >> BASIS2) & MASK;
      val4 = (db4 >> BASIS2) & MASK;
      val5 = (db5 >> BASIS2) & MASK;
      val6 = (db6 >> BASIS2) & MASK;
      val7 = (db7 >> BASIS2) & MASK;
      val8 = (db8 >> BASIS2) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];
      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];
      r_tmp  += val*b5[index2];
      r_tmp2 += val2*b5[index2];
      r_tmp3 += val3*b5[index2];
      r_tmp4 += val4*b5[index2];
      r_tmp5 += val5*b5[index2];
      r_tmp6 += val6*b5[index2];
      r_tmp7 += val7*b5[index2];
      r_tmp8 += val8*b5[index2];
      s_tmp  += val*b6[index2];
      s_tmp2 += val2*b6[index2];
      s_tmp3 += val3*b6[index2];
      s_tmp4 += val4*b6[index2];
      s_tmp5 += val5*b6[index2];
      s_tmp6 += val6*b6[index2];
      s_tmp7 += val7*b6[index2];
      s_tmp8 += val8*b6[index2];
      index2 += 1;

      val  = (db >> BASIS3) & MASK;
      val2 = (db2 >> BASIS3) & MASK;
      val3 = (db3 >> BASIS3) & MASK;
      val4 = (db4 >> BASIS3) & MASK;
      val5 = (db5 >> BASIS3) & MASK;
      val6 = (db6 >> BASIS3) & MASK;
      val7 = (db7 >> BASIS3) & MASK;
      val8 = (db8 >> BASIS3) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];
      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];
      r_tmp  += val*b5[index2];
      r_tmp2 += val2*b5[index2];
      r_tmp3 += val3*b5[index2];
      r_tmp4 += val4*b5[index2];
      r_tmp5 += val5*b5[index2];
      r_tmp6 += val6*b5[index2];
      r_tmp7 += val7*b5[index2];
      r_tmp8 += val8*b5[index2];
      s_tmp  += val*b6[index2];
      s_tmp2 += val2*b6[index2];
      s_tmp3 += val3*b6[index2];
      s_tmp4 += val4*b6[index2];
      s_tmp5 += val5*b6[index2];
      s_tmp6 += val6*b6[index2];
      s_tmp7 += val7*b6[index2];
      s_tmp8 += val8*b6[index2];
      index2 += 1;

      index += 1;
    }
    out[6*i+0]  += tmp;
    out[6*i+1]  += x_tmp;
    out[6*i+2]  += y_tmp;
    out[6*i+3]  += z_tmp;
    out[6*i+4]  += r_tmp;
    out[6*i+5]  += s_tmp;
    out[6*i+6]  += tmp2;
    out[6*i+7]  += x_tmp2;
    out[6*i+8]  += y_tmp2;
    out[6*i+9]  += z_tmp2;
    out[6*i+10] += r_tmp2;
    out[6*i+11] += s_tmp2;
    out[6*i+12] += tmp3;
    out[6*i+13] += x_tmp3;
    out[6*i+14] += y_tmp3;
    out[6*i+15] += z_tmp3;
    out[6*i+16] += r_tmp3;
    out[6*i+17] += s_tmp3;    
    out[6*i+18] += tmp4;
    out[6*i+19] += x_tmp4;
    out[6*i+20] += y_tmp4;
    out[6*i+21] += z_tmp4;
    out[6*i+22] += r_tmp4;
    out[6*i+23] += s_tmp4;
    out[6*i+24] += tmp5;
    out[6*i+25] += x_tmp5;
    out[6*i+26] += y_tmp5;
    out[6*i+27] += z_tmp5;
    out[6*i+28] += r_tmp5;
    out[6*i+29] += s_tmp5;
    out[6*i+30] += tmp6;
    out[6*i+31] += x_tmp6;
    out[6*i+32] += y_tmp6;
    out[6*i+33] += z_tmp6;
    out[6*i+34] += r_tmp6;
    out[6*i+35] += s_tmp6;
    out[6*i+36] += tmp7;
    out[6*i+37] += x_tmp7;
    out[6*i+38] += y_tmp7;
    out[6*i+39] += z_tmp7;
    out[6*i+40] += r_tmp7;
    out[6*i+41] += s_tmp7;
    out[6*i+42] += tmp8;
    out[6*i+43] += x_tmp8;
    out[6*i+44] += y_tmp8;
    out[6*i+45] += z_tmp8;
    out[6*i+46] += r_tmp8;
    out[6*i+47] += s_tmp8;
    index += aCols*7;
  }
}

void trueMatMulVecPacked8(Elem *__restrict__ out, const Elem *__restrict__ a, 
    const Elem *__restrict__ b,
    const Elem *__restrict__ b2,
    const Elem *__restrict__ b3,
    const Elem *__restrict__ b4,
    const Elem *__restrict__ b5,
    const Elem *__restrict__ b6,
    const Elem *__restrict__ b7,
    const Elem *__restrict__ b8,
    size_t aRows, size_t aCols)
{
  // in this variant, b has 8 columns
  
  Elem db, db2, db3, db4, db5, db6, db7, db8;
  Elem val, val2, val3, val4, val5, val6, val7, val8;
  Elem tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
  Elem x_tmp, x_tmp2, x_tmp3, x_tmp4, x_tmp5, x_tmp6, x_tmp7, x_tmp8;
  Elem y_tmp, y_tmp2, y_tmp3, y_tmp4, y_tmp5, y_tmp6, y_tmp7, y_tmp8;
  Elem z_tmp, z_tmp2, z_tmp3, z_tmp4, z_tmp5, z_tmp6, z_tmp7, z_tmp8;
  Elem r_tmp, r_tmp2, r_tmp3, r_tmp4, r_tmp5, r_tmp6, r_tmp7, r_tmp8;
  Elem s_tmp, s_tmp2, s_tmp3, s_tmp4, s_tmp5, s_tmp6, s_tmp7, s_tmp8;
  Elem t_tmp, t_tmp2, t_tmp3, t_tmp4, t_tmp5, t_tmp6, t_tmp7, t_tmp8;
  Elem u_tmp, u_tmp2, u_tmp3, u_tmp4, u_tmp5, u_tmp6, u_tmp7, u_tmp8;
  size_t index = 0;
  size_t index2;

  for (size_t i = 0; i < aRows; i += 8) {
    tmp  = 0;
    tmp2 = 0;
    tmp3 = 0;
    tmp4 = 0;
    tmp5 = 0;
    tmp6 = 0;
    tmp7 = 0;
    tmp8 = 0;

    x_tmp  = 0;
    x_tmp2 = 0;
    x_tmp3 = 0;
    x_tmp4 = 0;
    x_tmp5 = 0;
    x_tmp6 = 0;
    x_tmp7 = 0;
    x_tmp8 = 0;

    y_tmp  = 0;
    y_tmp2 = 0;
    y_tmp3 = 0;
    y_tmp4 = 0;
    y_tmp5 = 0;
    y_tmp6 = 0;
    y_tmp7 = 0;
    y_tmp8 = 0;

    z_tmp  = 0;
    z_tmp2 = 0;
    z_tmp3 = 0;
    z_tmp4 = 0;
    z_tmp5 = 0;
    z_tmp6 = 0;
    z_tmp7 = 0;
    z_tmp8 = 0;

    r_tmp  = 0;
    r_tmp2 = 0;
    r_tmp3 = 0;
    r_tmp4 = 0;
    r_tmp5 = 0;
    r_tmp6 = 0;
    r_tmp7 = 0;
    r_tmp8 = 0;

    s_tmp  = 0;
    s_tmp2 = 0;
    s_tmp3 = 0;
    s_tmp4 = 0;
    s_tmp5 = 0;
    s_tmp6 = 0;
    s_tmp7 = 0;
    s_tmp8 = 0;

    t_tmp  = 0;
    t_tmp2 = 0;
    t_tmp3 = 0;
    t_tmp4 = 0;
    t_tmp5 = 0;
    t_tmp6 = 0;
    t_tmp7 = 0;
    t_tmp8 = 0;

    u_tmp  = 0;
    u_tmp2 = 0;
    u_tmp3 = 0;
    u_tmp4 = 0;
    u_tmp5 = 0;
    u_tmp6 = 0;
    u_tmp7 = 0;
    u_tmp8 = 0;


    index2 = 0;
    for (size_t j = 0; j < aCols; j++) {
      db  = a[index];
      db2 = a[index+1*aCols];
      db3 = a[index+2*aCols];
      db4 = a[index+3*aCols];
      db5 = a[index+4*aCols];
      db6 = a[index+5*aCols];
      db7 = a[index+6*aCols];
      db8 = a[index+7*aCols];

      val  = db & MASK;
      val2 = db2 & MASK;
      val3 = db3 & MASK;
      val4 = db4 & MASK;
      val5 = db5 & MASK;
      val6 = db6 & MASK;
      val7 = db7 & MASK;
      val8 = db8 & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];
      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];
      r_tmp  += val*b5[index2];
      r_tmp2 += val2*b5[index2];
      r_tmp3 += val3*b5[index2];
      r_tmp4 += val4*b5[index2];
      r_tmp5 += val5*b5[index2];
      r_tmp6 += val6*b5[index2];
      r_tmp7 += val7*b5[index2];
      r_tmp8 += val8*b5[index2];
      s_tmp  += val*b6[index2];
      s_tmp2 += val2*b6[index2];
      s_tmp3 += val3*b6[index2];
      s_tmp4 += val4*b6[index2];
      s_tmp5 += val5*b6[index2];
      s_tmp6 += val6*b6[index2];
      s_tmp7 += val7*b6[index2];
      s_tmp8 += val8*b6[index2];
      t_tmp  += val*b7[index2];
      t_tmp2 += val2*b7[index2];
      t_tmp3 += val3*b7[index2];
      t_tmp4 += val4*b7[index2];
      t_tmp5 += val5*b7[index2];
      t_tmp6 += val6*b7[index2];
      t_tmp7 += val7*b7[index2];
      t_tmp8 += val8*b7[index2];
      u_tmp  += val*b8[index2];
      u_tmp2 += val2*b8[index2];
      u_tmp3 += val3*b8[index2];
      u_tmp4 += val4*b8[index2];
      u_tmp5 += val5*b8[index2];
      u_tmp6 += val6*b8[index2];
      u_tmp7 += val7*b8[index2];
      u_tmp8 += val8*b8[index2];


      index2 += 1;

      val  = (db >> BASIS) & MASK;
      val2 = (db2 >> BASIS) & MASK;
      val3 = (db3 >> BASIS) & MASK;
      val4 = (db4 >> BASIS) & MASK;
      val5 = (db5 >> BASIS) & MASK;
      val6 = (db6 >> BASIS) & MASK;
      val7 = (db7 >> BASIS) & MASK;
      val8 = (db8 >> BASIS) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];
      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];
      r_tmp  += val*b5[index2];
      r_tmp2 += val2*b5[index2];
      r_tmp3 += val3*b5[index2];
      r_tmp4 += val4*b5[index2];
      r_tmp5 += val5*b5[index2];
      r_tmp6 += val6*b5[index2];
      r_tmp7 += val7*b5[index2];
      r_tmp8 += val8*b5[index2];
      s_tmp  += val*b6[index2];
      s_tmp2 += val2*b6[index2];
      s_tmp3 += val3*b6[index2];
      s_tmp4 += val4*b6[index2];
      s_tmp5 += val5*b6[index2];
      s_tmp6 += val6*b6[index2];
      s_tmp7 += val7*b6[index2];
      s_tmp8 += val8*b6[index2];
      t_tmp  += val*b7[index2];
      t_tmp2 += val2*b7[index2];
      t_tmp3 += val3*b7[index2];
      t_tmp4 += val4*b7[index2];
      t_tmp5 += val5*b7[index2];
      t_tmp6 += val6*b7[index2];
      t_tmp7 += val7*b7[index2];
      t_tmp8 += val8*b7[index2];
      u_tmp  += val*b8[index2];
      u_tmp2 += val2*b8[index2];
      u_tmp3 += val3*b8[index2];
      u_tmp4 += val4*b8[index2];
      u_tmp5 += val5*b8[index2];
      u_tmp6 += val6*b8[index2];
      u_tmp7 += val7*b8[index2];
      u_tmp8 += val8*b8[index2];

      index2 += 1;

      val  = (db >> BASIS2) & MASK;
      val2 = (db2 >> BASIS2) & MASK;
      val3 = (db3 >> BASIS2) & MASK;
      val4 = (db4 >> BASIS2) & MASK;
      val5 = (db5 >> BASIS2) & MASK;
      val6 = (db6 >> BASIS2) & MASK;
      val7 = (db7 >> BASIS2) & MASK;
      val8 = (db8 >> BASIS2) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];
      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];
      r_tmp  += val*b5[index2];
      r_tmp2 += val2*b5[index2];
      r_tmp3 += val3*b5[index2];
      r_tmp4 += val4*b5[index2];
      r_tmp5 += val5*b5[index2];
      r_tmp6 += val6*b5[index2];
      r_tmp7 += val7*b5[index2];
      r_tmp8 += val8*b5[index2];
      s_tmp  += val*b6[index2];
      s_tmp2 += val2*b6[index2];
      s_tmp3 += val3*b6[index2];
      s_tmp4 += val4*b6[index2];
      s_tmp5 += val5*b6[index2];
      s_tmp6 += val6*b6[index2];
      s_tmp7 += val7*b6[index2];
      s_tmp8 += val8*b6[index2];
      t_tmp  += val*b7[index2];
      t_tmp2 += val2*b7[index2];
      t_tmp3 += val3*b7[index2];
      t_tmp4 += val4*b7[index2];
      t_tmp5 += val5*b7[index2];
      t_tmp6 += val6*b7[index2];
      t_tmp7 += val7*b7[index2];
      t_tmp8 += val8*b7[index2];
      u_tmp  += val*b8[index2];
      u_tmp2 += val2*b8[index2];
      u_tmp3 += val3*b8[index2];
      u_tmp4 += val4*b8[index2];
      u_tmp5 += val5*b8[index2];
      u_tmp6 += val6*b8[index2];
      u_tmp7 += val7*b8[index2];
      u_tmp8 += val8*b8[index2];

      index2 += 1;

      val  = (db >> BASIS3) & MASK;
      val2 = (db2 >> BASIS3) & MASK;
      val3 = (db3 >> BASIS3) & MASK;
      val4 = (db4 >> BASIS3) & MASK;
      val5 = (db5 >> BASIS3) & MASK;
      val6 = (db6 >> BASIS3) & MASK;
      val7 = (db7 >> BASIS3) & MASK;
      val8 = (db8 >> BASIS3) & MASK;
      tmp  += val*b[index2];
      tmp2 += val2*b[index2];
      tmp3 += val3*b[index2];
      tmp4 += val4*b[index2];
      tmp5 += val5*b[index2];
      tmp6 += val6*b[index2];
      tmp7 += val7*b[index2];
      tmp8 += val8*b[index2];
      x_tmp  += val*b2[index2];
      x_tmp2 += val2*b2[index2];
      x_tmp3 += val3*b2[index2];
      x_tmp4 += val4*b2[index2];
      x_tmp5 += val5*b2[index2];
      x_tmp6 += val6*b2[index2];
      x_tmp7 += val7*b2[index2];
      x_tmp8 += val8*b2[index2];
      y_tmp  += val*b3[index2];
      y_tmp2 += val2*b3[index2];
      y_tmp3 += val3*b3[index2];
      y_tmp4 += val4*b3[index2];
      y_tmp5 += val5*b3[index2];
      y_tmp6 += val6*b3[index2];
      y_tmp7 += val7*b3[index2];
      y_tmp8 += val8*b3[index2];
      z_tmp  += val*b4[index2];
      z_tmp2 += val2*b4[index2];
      z_tmp3 += val3*b4[index2];
      z_tmp4 += val4*b4[index2];
      z_tmp5 += val5*b4[index2];
      z_tmp6 += val6*b4[index2];
      z_tmp7 += val7*b4[index2];
      z_tmp8 += val8*b4[index2];
      r_tmp  += val*b5[index2];
      r_tmp2 += val2*b5[index2];
      r_tmp3 += val3*b5[index2];
      r_tmp4 += val4*b5[index2];
      r_tmp5 += val5*b5[index2];
      r_tmp6 += val6*b5[index2];
      r_tmp7 += val7*b5[index2];
      r_tmp8 += val8*b5[index2];
      s_tmp  += val*b6[index2];
      s_tmp2 += val2*b6[index2];
      s_tmp3 += val3*b6[index2];
      s_tmp4 += val4*b6[index2];
      s_tmp5 += val5*b6[index2];
      s_tmp6 += val6*b6[index2];
      s_tmp7 += val7*b6[index2];
      s_tmp8 += val8*b6[index2];
      t_tmp  += val*b7[index2];
      t_tmp2 += val2*b7[index2];
      t_tmp3 += val3*b7[index2];
      t_tmp4 += val4*b7[index2];
      t_tmp5 += val5*b7[index2];
      t_tmp6 += val6*b7[index2];
      t_tmp7 += val7*b7[index2];
      t_tmp8 += val8*b7[index2];
      u_tmp  += val*b8[index2];
      u_tmp2 += val2*b8[index2];
      u_tmp3 += val3*b8[index2];
      u_tmp4 += val4*b8[index2];
      u_tmp5 += val5*b8[index2];
      u_tmp6 += val6*b8[index2];
      u_tmp7 += val7*b8[index2];
      u_tmp8 += val8*b8[index2];

      index2 += 1;
      
      index += 1;
    }
    out[8*i+0]    += tmp;
    out[8*i+1]  += x_tmp;
    out[8*i+2]  += y_tmp;
    out[8*i+3]  += z_tmp;
    out[8*i+4]  += r_tmp;
    out[8*i+5]  += s_tmp;
    out[8*i+6]  += t_tmp;
    out[8*i+7]  += u_tmp;
    out[8*i+8]  += tmp2;
    out[8*i+9]  += x_tmp2;
    out[8*i+10]  += y_tmp2;
    out[8*i+11]  += z_tmp2;
    out[8*i+12]  += r_tmp2;
    out[8*i+13]  += s_tmp2;
    out[8*i+14]  += t_tmp2;
    out[8*i+15]  += u_tmp2;
    out[8*i+16]  += tmp3;
    out[8*i+17]  += x_tmp3;
    out[8*i+18]  += y_tmp3;
    out[8*i+19]  += z_tmp3;
    out[8*i+20]  += r_tmp3;
    out[8*i+21]  += s_tmp3;
    out[8*i+22]  += t_tmp3;
    out[8*i+23]  += u_tmp3;
    out[8*i+24]  += tmp4;
    out[8*i+25]  += x_tmp4;
    out[8*i+26]  += y_tmp4;
    out[8*i+27]  += z_tmp4;
    out[8*i+28]  += r_tmp4;
    out[8*i+29]  += s_tmp4;
    out[8*i+30]  += t_tmp4;
    out[8*i+31]  += u_tmp4;
    out[8*i+32]  += tmp5;
    out[8*i+33]  += x_tmp5;
    out[8*i+34]  += y_tmp5;
    out[8*i+35]  += z_tmp5;
    out[8*i+36]  += r_tmp5;
    out[8*i+37]  += s_tmp5;
    out[8*i+38]  += t_tmp5;
    out[8*i+39]  += u_tmp5;
    out[8*i+40]  += tmp6;
    out[8*i+41]  += x_tmp6;
    out[8*i+42]  += y_tmp6;
    out[8*i+43]  += z_tmp6;
    out[8*i+44]  += r_tmp6;
    out[8*i+45]  += s_tmp6;
    out[8*i+46]  += t_tmp6;
    out[8*i+47]  += u_tmp6;
    out[8*i+48]  += tmp7;
    out[8*i+49]  += x_tmp7;
    out[8*i+50]  += y_tmp7;
    out[8*i+51]  += z_tmp7;
    out[8*i+52]  += r_tmp7;
    out[8*i+53]  += s_tmp7;
    out[8*i+54]  += t_tmp7;
    out[8*i+55]  += u_tmp7;
    out[8*i+56]  += tmp8;
    out[8*i+57]  += x_tmp8;
    out[8*i+58]  += y_tmp8;
    out[8*i+59]  += z_tmp8;
    out[8*i+60]  += r_tmp8;
    out[8*i+61]  += s_tmp8;
    out[8*i+62]  += t_tmp8;
    out[8*i+63]  += u_tmp8;
    index += aCols*7;
  }
}

void matMulVecPacked8(Elem *__restrict__ out, const Elem *__restrict__ a, const Elem *__restrict__ b_full,
    size_t aRows, size_t aCols)
{
  #ifdef PROFILE
    System::profile("matmul4", [&]() {
  #endif

  trueMatMulVecPacked8(out, a, b_full, b_full + (1*aCols*COMPRESSION), b_full + (2*aCols*COMPRESSION), b_full + (3*aCols*COMPRESSION), b_full + (4*aCols*COMPRESSION), b_full + (5*aCols*COMPRESSION), b_full + (6*aCols*COMPRESSION), b_full + (7*aCols*COMPRESSION), aRows, aCols);
  
  #ifdef PROFILE
    });
  #endif
}