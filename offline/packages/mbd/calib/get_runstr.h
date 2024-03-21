#ifndef __GET_RUNSTR__
#define __GET_RUNSTR__

//namespace RUNSTR {

TString get_runstr(const char *fname)
{
  TString name = fname;
  name.ReplaceAll(".prdf","");
  name.ReplaceAll(".root","");
  int index = name.Last('/');
  if ( index > 0 )
  {
    name.Remove(0,index+1);
  }
  index = name.First('-');
  if ( index > 0 )
  {
    name.Remove(0,index+1);
  }
  //cout << "aaa " << name << endl;

  /*
  index = name.Last('-');
  if ( index > 0 )
  {
    name.Remove(index,name.Length());
  }
  //cout << "bbb " << name << endl;
  */

  return name;

}

int get_runnumber(const char *fname)
{
  TString str = get_runstr(fname);
  cout << str << endl;
  int index = str.Last('-');
  if ( index > 0 )
  {
    str.Remove(index,str.Length());
  }
  cout << " get_runnumber " << str << "\t" << str.Atoi() << endl;
  return str.Atoi();
}

//} // namespace RUNSTR

#endif  // __get_runstr__

