#ifndef INTT_FELIX_MAP_H
#define INTT_FELIX_MAP_H

namespace InttNameSpace
{
  struct Online_s;
}
namespace InttNameSpace
{
  struct RawData_s;
}

namespace InttFelix
{
  int RawDataToOnline(struct InttNameSpace::RawData_s const&, struct InttNameSpace::Online_s&);
  int OnlineToRawData(struct InttNameSpace::Online_s const&, struct InttNameSpace::RawData_s&);
};  // namespace InttFelix

#endif  // INTT_FELIX_MAP_H
