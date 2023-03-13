
#ifndef MULTIARRAY_H
#define MULTIARRAY_H

#include <cassert>
#include <cstdio>   // for printf
#include <cstdlib>  // for malloc

template <class T>
class MultiArray
{
  //class to hold an up-to-six dimensional array of whatever T is.  Any indices not used are flattened.  This should probably be replaced with sets of TH3s... but the intention was to take advantage of indices for all elements being the same, to avoid unpacking and re-packing TVectors, etc, and to get rid of any other overhead that might be showing up in the TH3 implementation.
  //it does make it more annoying to interpolate, though.
 public:
  static const int MAX_DIM = 6;
  int dim;
  int n[6];
  long int length;
  T *field;

  MultiArray(int a = 0, int b = 0, int c = 0, int d = 0, int e = 0, int f = 0)
  {
    int n_[6];
    for (int i = 0; i < MAX_DIM; i++)
      n[i] = 0;
    n_[0] = a;
    n_[1] = b;
    n_[2] = c;
    n_[3] = d;
    n_[4] = e;
    n_[5] = f;
    length = 1;
    dim = MAX_DIM;
    for (int i = 0; i < dim; i++)
    {
      if (n_[i] < 1)
      {
        dim = i;
        break;
      }
      n[i] = n_[i];
      length *= n[i];
    }
    field = static_cast<T *>(malloc(length * sizeof(T)));
    //note that since we don't know what T is, we can't safely zero it.  Someone else will have to do that.
  }
  //! delete copy ctor and assignment opertor (cppcheck)
  explicit MultiArray(const MultiArray &) = delete;
  MultiArray &operator=(const MultiArray &) = delete;

  ~MultiArray()
  {
    free(field);
  }

  void Add(int a, int b, int c, T in)
  {
    Add(a, b, c, 0, 0, 0, in);
    return;
  };

  void Add(int a, int b, int c, int d, int e, int f, T in)
  {
    int n_[6];
    n_[0] = a;
    n_[1] = b;
    n_[2] = c;
    n_[3] = d;
    n_[4] = e;
    n_[5] = f;
    long int index = n_[0];
    for (int i = 1; i < dim; i++)
    {
      index = (index * n[i]) + n_[i];
    }
    field[index] = field[index] + in;
    return;
  }

  T Get(int a = 0, int b = 0, int c = 0, int d = 0, int e = 0, int f = 0)
  {
    int n_[6];
    n_[0] = a;
    n_[1] = b;
    n_[2] = c;
    n_[3] = d;
    n_[4] = e;
    n_[5] = f;
    long int index = 0;
    for (int i = 0; i < dim; i++)
    {
      if (n[i] <= n_[i] || n_[i] < 0)
      {  //check bounds
        printf("asking for el %d %d %d %d %d %d.  %dth element is outside of bounds 0<x<%d\n", n_[0], n_[1], n_[2], n_[3], n_[4], n_[5], n_[i], n[i]);
        assert(false);
      }
      index = (index * n[i]) + n_[i];
    }
    return field[index];
  }

  T *GetPtr(int a = 0, int b = 0, int c = 0, int d = 0, int e = 0, int f = 0)
  {  //faster for repeated access.
    int n_[6];
    n_[0] = a;
    n_[1] = b;
    n_[2] = c;
    n_[3] = d;
    n_[4] = e;
    n_[5] = f;
    long int index = n_[0];
    for (int i = 1; i < dim; i++)
    {
      index = (index * n[i]) + n_[i];
    }
    return &(field[index]);
  }

  T *GetFlat(int a = 0)
  {  //get the value at position a in the 1D equivalent, assuming the math is done elsewhere, or we're just going straight through the thing.
    if (a < 0 || a >= length)
    {
      printf("tried to seek element %d of multiarray, but bounds are 0<a<%ld\n", a, length);
      assert(a < 0 || a >= length);  //check bounds
    }
    return &(field[a]);
  }

  int Length()
  {
    return (int) length;
  }

  void Set(int a, int b, int c, T in)
  {
    Set(a, b, c, 0, 0, 0, in);
    return;
  };

  void Set(int a, int b, int c, int d, int e, int f, T in)
  {
    int n_[6];
    n_[0] = a;
    n_[1] = b;
    n_[2] = c;
    n_[3] = d;
    n_[4] = e;
    n_[5] = f;
    long int index = n_[0];
    for (int i = 1; i < dim; i++)
    {
      index = (index * n[i]) + n_[i];
    }
    field[index] = in;
    return;
  }

  void SetAll(T in)
  {
    //this assumes there's an '=' operator for T, but that's generally true.
    for (long int i = 0; i < length; i++)
    {
      field[i] = in;
    }
    return;
  }
};
#endif  //MULTIARRAY_H
