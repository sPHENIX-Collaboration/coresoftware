#ifndef __fastvec_h__
#define __fastvec_h__

#include <vector>
#include <cstring>
#include <iostream>


// a vector which resides on the stack when small enough,
// otherwise it resides on the heap
class fastvec
{
  public:
    fastvec() : size(0){}
    ~fastvec(){}
    
    void push_back(unsigned int a)
    {
      if(size > 16383)
      {
        vec.push_back(a);
        size += 1;
      }
      else
      {
        arr[size] = a;
        size += 1;
      }
    }
    
    void push_back(unsigned int* a, unsigned int num)
    {
      if( (size+num) <  16384)
      {
        memcpy(&(arr[size]), a, (num<<2));
        size += num;
      }
      else
      {
        for(unsigned int i=0;i<num;++i)
        {
          push_back(a[i]);
        }
      }
    }
    
    unsigned int& operator[](unsigned int idx)
    {
//       unsigned long int address = (( (unsigned long int)(&(arr[idx])) & ((unsigned long int)((idx > 16383) - 1))) ^ ( (unsigned long int)(&(vec[idx-16384])) & ((unsigned long int)((idx <= 16383) - 1))));
//       
//       return ( *( (unsigned int*)(address) ) );
      
      if(idx <= 16383)
      {
        return arr[idx];
      }
      else
      {
        return vec[idx-16384];
      }
    }
    
    void operator=(const fastvec& other)
    {
      size = other.size;
      if(size < 16384)
      {
        for(unsigned int i=0;i<size;++i)
        {
          arr[i] = other.arr[i];
        }
      }
      else
      {
        for(unsigned int i=0;i<16384;++i)
        {
          arr[i] = other.arr[i];
        }
        vec = other.vec;
      }
    }
    
    void clear()
    {
      size = 0;
      vec.clear();
    }
    
    std::vector<unsigned int> vec;
    unsigned int arr[16384];
    unsigned int size;
};


class fastvec2d
{
  public:
    fastvec2d(unsigned int second_dimension_size) : size(0), size2(second_dimension_size)
    {
      nstack = (16384/size2);
    }
    ~fastvec2d(){}
    
    void fill(unsigned int* buf, unsigned int bufsize)
    {
      if(size < nstack)
      {
        arr_size[size] = bufsize;
        memcpy(&(arr[size*size2]), &(buf[0]), (bufsize<<2));
        size += 1;
      }
      else
      {
        std::vector<unsigned int> tempvec;
        for(unsigned int i=0;i<bufsize;++i)
        {
          tempvec.push_back(buf[i]);
        }
        vec.push_back(tempvec);
        size += 1;
      }
    }
    
    unsigned int n_entries(unsigned int idx)
    {
      if(idx < nstack)
      {
        return arr_size[idx];
      }
      else
      {
        return vec[idx - nstack].size();
      }
    }
    
    unsigned int& operator()(unsigned int idx1, unsigned int idx2)
    {
      if(idx1 < nstack)
      {
        return arr[idx1*size2 + idx2];
      }
      else
      {
        return vec[idx1 - nstack][idx2];
      }
    }
    
    void fetch(unsigned int begin, unsigned int end, unsigned int* result_arr, unsigned int* result_size)
    {
      if(end < nstack)
      {
        memcpy(result_arr, &(arr[begin*size2]), ((1+end-begin)*size2)*4);
        memcpy(result_size, &(arr_size[begin]), (1+end-begin)*4);
      }
      else
      {
        for(unsigned int i=begin;i<=end;++i)
        {
          unsigned int entries = n_entries(i);
          result_size[i-begin] = entries;
          for(unsigned int j=0;j<entries;++j)
          {
            result_arr[(i-begin)*size2 + j] = (*this)(i,j);
          }
        }
      }
    }
    
    
    unsigned int arr[16384];
    unsigned int arr_size[16384];
    unsigned int size;
    unsigned int size2;
    unsigned int nstack;
    std::vector<std::vector<unsigned int> > vec;
};



#endif
