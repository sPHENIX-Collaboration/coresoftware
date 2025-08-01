#include "SingleTriggeredInput.h"
#include "SingleGl1TriggeredInput.h"

#include <frog/FROG.h>

#include <ffarawobjects/CaloPacketContainerv1.h>
#include <ffarawobjects/CaloPacketv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TSystem.h>

#include <cstdint>   // for uint64_t
#include <iostream>  // for operator<<, basic_ostream, endl
#include <ranges>
#include <set>
#include <unordered_set>
#include <utility>  // for pair
#include <vector>

SingleTriggeredInput::SingleTriggeredInput(const std::string &name)
  : Fun4AllBase(name)
{
  m_bclkarray.fill(std::numeric_limits<uint64_t>::max());
  m_bclkdiffarray.fill(std::numeric_limits<uint64_t>::max());
}

SingleTriggeredInput::~SingleTriggeredInput()
{
  for (auto& [pid, dq] : m_PacketEventDeque)
  {
    while (!dq.empty())
    {
      delete dq.front();
      dq.pop_front();
    }
  }

  for (auto& [pid, evt] : m_PacketEventBackup)
  {
    delete evt;
  }

  delete m_EventIterator;
}

bool SingleTriggeredInput::CheckFemDiffIdx(int pid, size_t index, const std::deque<Event*>& events, uint64_t gl1diffidx)
{
  if (index >= events.size())
  {
    return false;
  }

  Packet* pkt_prev = events[index-1]->getPacket(pid);
  Packet* pkt_curr = events[index]->getPacket(pid);
  if (!pkt_prev || !pkt_curr)
  {
    delete pkt_prev;
    delete pkt_curr;
    return false;
  }

  auto get_majority_femclk = [](Packet* pkt) -> uint16_t {
    int nmod = pkt->iValue(0, "NRMODULES");
    std::map<uint16_t, int> counts;
    for (int j = 0; j < nmod; ++j)
    {
      uint16_t clk = static_cast<uint16_t>(pkt->iValue(j, "FEMCLOCK"));
      counts[clk]++;
    }
    if (counts.empty())
    {
      return std::numeric_limits<uint16_t>::max();
    }
    return std::max_element(counts.begin(), counts.end(), [](const auto& a, const auto& b) { return a.second < b.second; })->first;
  };

  uint16_t clk_prev = get_majority_femclk(pkt_prev);
  uint16_t clk_curr = get_majority_femclk(pkt_curr);

  delete pkt_prev;
  delete pkt_curr;

  if (clk_prev == std::numeric_limits<uint16_t>::max() || clk_curr == std::numeric_limits<uint16_t>::max())
  {
    return false;
  }

  uint16_t femdiff = static_cast<uint16_t>(clk_curr - clk_prev);
  gl1diffidx = static_cast<uint16_t>(gl1diffidx & 0xFFFFU);
  return (femdiff == gl1diffidx);
}

bool SingleTriggeredInput::CheckPoolAlignment(int pid, const std::array<uint64_t, pooldepth>& sebdiff, const std::array<uint64_t, pooldepth>& gl1diff, std::vector<int>& bad_indices, int& shift, bool& CurrentPoolLastDiffBad, bool PrevPoolLastDiffBad)
{
  bad_indices.clear();
  shift = 0;
  CurrentPoolLastDiffBad=false;

  if (std::equal(sebdiff.begin(), sebdiff.end(), gl1diff.begin()))
  {
    return true;
  }

  //Finding intermittent corrupted data
  size_t n = sebdiff.size();
  std::vector<int> bad_diff_indices;
  for (size_t i = 0; i < n; ++i)
  {
    if ( sebdiff[i] != gl1diff[i] )
    {
      if ( !m_packetclk_copy_runs )
      {
        //backup procedure to recover stuck 16bit XMIT clock
        size_t idxcheck =  i == 0  ? i+1 : i;
        bool passFemDiffCheckIdx = CheckFemDiffIdx(pid, idxcheck, m_PacketEventDeque[pid], gl1diff[idxcheck]);
        if ( passFemDiffCheckIdx )
        {
          m_OverrideWithRepClock.insert(pid);
          continue;
        }
      } 
      bad_diff_indices.push_back(i);
    }
  }

  if (bad_diff_indices.empty())
  {
    if ( Verbosity() > 0 )
    {
      std::cout << Name() << " recovered from bad XMIT clocks. Merging pool" << std::endl;
    }
    return true;
  }

  bool move_to_shift_algo = false;
  if(bad_diff_indices.size() >=5)
  {
    std::cout << std::endl;
    std::cout << "----------------- " << Name() << " -----------------" << std::endl;
    std::cout << "More than 5 diffs are bad.. try shifting algorithm" << std::endl;
    move_to_shift_algo = true;
  }
  if(!move_to_shift_algo)
  {
    std::cout << std::endl;
    std::cout << "----------------- " << Name() << " -----------------" << std::endl;
  }

  size_t idx = 0;
  while (idx < bad_diff_indices.size() && !move_to_shift_algo)
  {
    int start = bad_diff_indices[idx];
    int end = start;
    while ((idx+1) < bad_diff_indices.size() && bad_diff_indices[idx+1] == end + 1)
    {
      ++idx;
      ++end;
    }

    int length = end - start + 1;
    if(length<=0)
    {
      std::cout << Name() << ": length of bad diffs is <=0. This should not happen... something very wrong. rejecting the pool" << std::endl;
      return false;
    }
    if(length>=5)
    {
      std::cout << Name() << ": length of bad diffs >=5 with bad_diff_indices.size() " << bad_diff_indices.size() << ". This should not have happened.. rejecting pool" << std::endl;
      return false;
    }

    if(start==static_cast<int>(pooldepth - 1))
    {
      bad_indices.push_back(start);
      CurrentPoolLastDiffBad= true;
    }
    else if (start==0)
    {
      if (PrevPoolLastDiffBad)
      {
        for (int j = start; j < end; ++j)
        {
          bad_indices.push_back(j);
        }
      }
      else
      {
        if (length == 1)
        {
          std::cout << Name() << ": diff[0] alone bad. isolated bad diff which should not happen.. rejecting pool" << std::endl;
          return false;
        }
        for (int j = start; j < end; ++j)
        {
          bad_indices.push_back(j);
        }
      }
    }
    else if (start < static_cast<int>(pooldepth - 1) && start >0)
    {
      if(length==1)
      {
        std::cout << Name() << ": Isolated bad diff[" << start << "] - rejecting pool" << std::endl;
        return false;
      }
      if(length>=2)
      {
        for (int j = start; j < end; ++j)
        {
          bad_indices.push_back(j);
        }
      }
    }
    else
    {
      std::cout << Name() << ": no categories assigned for length " << length << " and start / end " << start << " / " << end << " rejecting pool" << std::endl; 
      return false;
    }
    ++idx;
  }

  if (!move_to_shift_algo)
  {
    size_t nbads = bad_indices.size();

    if (nbads == 0 || nbads >= 4)
    {
      std::cout << Name() << ": unexpected number of bad events = " << nbads << " – rejecting pool" << std::endl;
      return false;
    }

    std::cout << Name() << ": intermittent bad events = " << nbads << " – do not try shifting algorithm" << std::endl;
    return true;
  }


  //Try shift 
  if(!move_to_shift_algo)
  {
    std::cout << Name() << ": Unexpected shift flag = " << move_to_shift_algo << ". Something went wrong - rejecting pool" << std::endl;
    return false;
  }
  
  if(move_to_shift_algo)
  {
    std::cout << Name() << ": Inconsistent diffs of " << bad_diff_indices.size() << ". Trying now shifting events to resynchronize" << std::endl;
  }

  bool match = true;
  bool first_pool = (gl1diff[0] == std::numeric_limits<uint64_t>::max());
  size_t start = first_pool ? 2 : 1;

  for (size_t i = start; i < pooldepth; ++i)
  {
    if (sebdiff[i] != gl1diff[i - 1])
    {
      match= false;
      break;
    }
  }
  if (match)
  {
    shift = -1;
    return true;
  }

  match = true;
  start = first_pool ? 1 : 0;
  for (size_t i = start; i < pooldepth - 1; ++i)
  {
    if (sebdiff[i] != gl1diff[i + 1])
    {
      match= false;
      break;
    }
  }
  if (match)
  {
    shift = 1;
    return true;
  }

  return false;
}

int SingleTriggeredInput::fileopen(const std::string &filenam)
{
  std::cout << PHWHERE << "trying to open " << filenam << std::endl;
  if (IsOpen())
  {
    std::cout << "Closing currently open file "
              << FileName()
              << " and opening " << filenam << std::endl;
    fileclose();
  }
  FileName(filenam);
  FROG frog;
  std::string fname = frog.location(FileName());
  if (Verbosity() > 0)
  {
    std::cout << Name() << ": opening file " << FileName() << std::endl;
  }
  int status = 0;
  m_EventIterator = new fileEventiterator(fname.c_str(), status);
  if (status)
  {
    delete m_EventIterator;
    m_EventIterator = nullptr;
    std::cout << PHWHERE << Name() << ": could not open file " << fname << std::endl;
    return -1;
  }
  IsOpen(1);
  AddToFileOpened(fname);  // add file to the list of files which were opened
  return 0;
}

int SingleTriggeredInput::fileclose()
{
  if (!IsOpen())
  {
    std::cout << Name() << ": fileclose: No Input file open" << std::endl;
    return -1;
  }
  delete m_EventIterator;
  m_EventIterator = nullptr;
  IsOpen(0);
  UpdateFileList();
  return 0;
}

int SingleTriggeredInput::FillEventVector()
{
  while (GetEventIterator() == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return -1;
    }
  }

  bool allPacketEventDequeEmpty = true;
  int representative_pid = -1;
  for (int pid : m_PacketSet)
  {
    if ( !m_PacketEventDeque[pid].empty() )
    {
      allPacketEventDequeEmpty = false;
      break;
    }
    uint64_t tmp = m_bclkarray_map[pid][pooldepth];
    m_bclkarray_map[pid].fill(std::numeric_limits<uint64_t>::max());
    m_bclkarray_map[pid][0] = tmp;
    m_bclkdiffarray_map[pid].fill(std::numeric_limits<uint64_t>::max());

    if ( representative_pid == -1 ) 
    {
      representative_pid = pid;
    }
  }
  if ( !allPacketEventDequeEmpty )
  {
    return 0;
  }

  size_t i{0};
  std::map<int, Event*> m_ShiftedEvents;

  while (i < pooldepth)
  {
    Event* evt{nullptr};
    if (this != Gl1Input())
    {
      auto* gl1 = dynamic_cast<SingleGl1TriggeredInput*>(Gl1Input());
      if (gl1)
      {
        int nskip = gl1->GetGl1SkipArray()[i];
        if (m_Gl1PacketOneSkipActiveTrace)
        {
          int gl1packetcountdiff = static_cast<int>(gl1->GetPacketNumbers()[i] - m_Gl1PacketNumberOneSkip);
          if ( nskip == -1 && gl1packetcountdiff < m_Gl1PacketOneSkipCount )
          {
            m_Gl1PacketOneSkipActiveTrace = false;
            if (Verbosity() > 0)
            {
              std::cout << Name() << ": GL1 stuck detected at " << gl1->GetPacketNumbers()[i] << " after skipping " << gl1packetcountdiff << ". Clearing trace." << std::endl;
            }
          }
          else if ( gl1packetcountdiff >= m_Gl1PacketOneSkipCount )
          {
            m_Gl1PacketOneSkipActiveTrace = false;
            if (Verbosity() > 0)
            {
              std::cout <<  Name() << ": No stuck found before " << m_Gl1PacketOneSkipCount << " events. Skipping one SEB event now." << std::endl;
            }
            Event* skip_evt = GetEventIterator()->getNextEvent();
            while (!skip_evt)
            {
              fileclose();
              if (!OpenNextFile())
              {
                FilesDone(1);
                return -1;
              }
              skip_evt = GetEventIterator()->getNextEvent();
            }
          }
        }
        if (nskip == 1 && m_packetclk_copy_runs == true)
        {
          m_Gl1PacketOneSkipActiveTrace = true;
          m_Gl1PacketNumberOneSkip = gl1->GetPacketNumbers()[i];
          if (Verbosity() > 0)
          {
            std::cout << Name() << ": GL1 packet skip for " << m_Gl1PacketNumberOneSkip << ". Start tracing Gl1 packet number." << std::endl;
          }
          nskip = 0;
        }

        while (nskip > 0)
        {
          Event* skip_evt = GetEventIterator()->getNextEvent();
          while (!skip_evt)
          {
            fileclose();
            if (!OpenNextFile())
            {
              FilesDone(1);
              return -1;
            }
            skip_evt = GetEventIterator()->getNextEvent();
          }
          if (Verbosity() > 0)
          {
            std::cout << Name() << ": Skipping SEB events because of GL1 packet number diff : " << nskip
              << ", with event sequence number " << skip_evt->getEvtSequence() << std::endl;
          }

          Packet* pkt = skip_evt->getPacket(representative_pid);
          if (!pkt)
          {
            delete skip_evt;
            continue;
          }

          FillPacketClock(skip_evt, pkt, i);
          delete pkt;

          uint64_t seb_diff = m_bclkdiffarray_map[representative_pid][i];
          int gl1pid = Gl1Input()->m_bclkdiffarray_map.begin()->first;
          uint64_t gl1_diff = Gl1Input()->m_bclkdiffarray_map[gl1pid][i];

          if (seb_diff == gl1_diff)
          {
            if (Verbosity() > 0)
            {
              std::cout << Name() << ": Early stop of SEB skip after " << (gl1->GetGl1SkipArray()[i] - nskip)
                << " from intial " << gl1->GetGl1SkipArray()[i] << " events." << std::endl;
            }
            evt = skip_evt;
            break;
          }
          delete skip_evt;
          nskip--;
        }
      }
    }

    if (!evt)
    {
      evt = GetEventIterator()->getNextEvent();
      while (!evt)
      {
        fileclose();
        if (!OpenNextFile())
        {
          FilesDone(1);
          return -1;
        }
        evt = GetEventIterator()->getNextEvent();
      }
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      if (Verbosity() > 0)
      {
        std::cout << Name() << " dropping non data event: " << evt->getEvtSequence() << std::endl;
      }
      delete evt;
      continue;
    }
    evt->convert();
    
    if (firstcall)
    {
      std::cout << "Creating DSTs first call" << std::endl;
      CreateDSTNodes(evt);
      int run = evt->getRunNumber();
      m_packetclk_copy_runs = (run >= 44000 && run < 56079); 
      firstcall = false;
    }

    for (int pid : m_PacketSet)
    {
      Event *thisevt = evt;
      if (m_PacketShiftOffset[pid] == 1)
      {
        if (i==0)
        {
          thisevt = m_PacketEventBackup[pid];
          m_ShiftedEvents[pid] = evt;
        }
        else if (i > 0)
        {
          thisevt = m_ShiftedEvents[pid];
          m_ShiftedEvents[pid] = evt;
          if (i == pooldepth -1)
          {
            m_PacketEventBackup[pid] = evt;
          }
        }
      }
      
      Packet* pkt = thisevt->getPacket(pid);
      if (!pkt)
      {
        continue;
      }
      FillPacketClock(thisevt, pkt, i);
      m_PacketEventDeque[pid].push_back(thisevt);
      delete pkt;
    
      if (representative_pid == -1 && m_PacketShiftOffset[pid] == 0)
      {
        representative_pid = pid;
      }
    }
    i++;
  }

  size_t minSize = pooldepth;
  for (const auto& [pid, dq] : m_PacketEventDeque)
  {
    minSize = std::min(dq.size(), minSize);
  }
  return minSize;
}

uint64_t SingleTriggeredInput::GetClock(Event *evt, int pid)
{
  Packet* packet = evt->getPacket(pid);
  if (!packet)
  {
    std::cout << Name() << ": Missing packet " << pid << " in event " << evt->getEvtSequence() << std::endl;
    return std::numeric_limits<uint64_t>::max();
  }
  uint64_t clkval = static_cast<uint64_t>(packet->lValue(0, "CLOCK"));
  uint64_t clk = clkval & 0xFFFFFFFFU;
  delete packet;
  return clk;
}

void SingleTriggeredInput::FillPacketClock(Event* evt, Packet* pkt, size_t event_index)
{
  if (!pkt)
  {
    return;
  }
  int pid = pkt->getIdentifier();

  if (m_bclkarray_map.find(pid) == m_bclkarray_map.end())
  {
    m_bclkarray_map[pid].fill(std::numeric_limits<uint64_t>::max());
    m_bclkdiffarray_map[pid].fill(std::numeric_limits<uint64_t>::max());
  }

  auto& clkarray = m_bclkarray_map[pid];
  auto& diffarray = m_bclkdiffarray_map[pid];


  // Special handling for FEM-copied clocks
  if (m_packetclk_copy_runs && m_CorrectCopiedClockPackets.count(pid))
  {
    if (event_index == 0)
    {
      clkarray[event_index+1] = m_PreviousValidBCOMap[pid];
    }
    else if (event_index >=1) 
    {
      Event* shifted_evt = m_PacketEventDeque[pid][event_index - 1];
      clkarray[event_index+1] = GetClock(shifted_evt, pid);
    }
    
    uint64_t prev = clkarray[event_index];
    uint64_t curr = clkarray[event_index + 1];

    if (prev == std::numeric_limits<uint64_t>::max() || curr == std::numeric_limits<uint64_t>::max())
    {
      diffarray[event_index] = std::numeric_limits<uint64_t>::max();
    }
    else
    {
      diffarray[event_index] = ComputeClockDiff(curr, prev);
    }

    return;
  }


  uint64_t clk = GetClock(evt, pid);
  if (clk == std::numeric_limits<uint64_t>::max())
  {
    std::cout << Name() << ": Bad clock for packet " << pid << " at event index " << event_index << std::endl;
    return;
  }

  clkarray[event_index + 1] = clk;

  uint64_t prev = clkarray[event_index];
  if(prev == std::numeric_limits<uint64_t>::max())
  {
    static std::unordered_set<int> warned;

    if (warned.find(pid) == warned.end())
    {
      std::cout << Name() << ": First pool for pacekt " << pid << " – skipping first diff because of no previous clock" << std::endl;
      warned.insert(pid);
    }
    else
    {
      std::cout << "prev clock is max something is wrong... : " << event_index << std::endl;
    }
    diffarray[event_index] = std::numeric_limits<uint64_t>::max();
  }
  else
  {
    diffarray[event_index] = ComputeClockDiff(clk, prev);
  }

  if (auto* gl1 = dynamic_cast<SingleGl1TriggeredInput*>(this))
  {
    int packet_number = pkt->iValue(0);
    gl1->SetPacketNumbers(gl1->GetCurrentPacketNumber(), packet_number);
    if ( event_index < pooldepth )
    {
      gl1->SetGl1PacketNumber(event_index, packet_number);
    }

    int skip_count = 0;
    if (gl1->GetLastPacketNumber() != 0)
    {
      int diff = gl1->GetCurrentPacketNumber() - gl1->GetLastPacketNumber() ;
      skip_count = diff - 1;
    }

    if (event_index < pooldepth)
    {
      gl1->SetGl1SkipAtIndex(event_index, skip_count);
    }
  }
}

void SingleTriggeredInput::FillPool()
{
  if (AllDone() || EventAlignmentProblem())
  {
    return;
  }

  bool all_packets_bad = !m_PacketAlignmentProblem.empty() && std::all_of(m_PacketAlignmentProblem.begin(), m_PacketAlignmentProblem.end(), [](const std::pair<const int, bool> &entry) -> bool { return entry.second;});
  if (all_packets_bad)
  {
    std::cout << Name() << ": ALL packets are marked as bad. Stop combining for this SEB." << std::endl;
    EventAlignmentProblem(1);
    return;
  }

  if (!FilesDone())
  {
    int eventvectorsize = FillEventVector();
    if (eventvectorsize != 0)
    {
      if (Gl1Input()->m_bclkdiffarray_map.empty())
      {
        std::cout << Name() << " : GL1 clock map is empty!" << std::endl;
        return;
      }
      m_OverrideWithRepClock.clear();

      int gl1pid = Gl1Input()->m_bclkdiffarray_map.begin()->first;
      const auto& gl1diff = Gl1Input()->m_bclkdiffarray_map.at(gl1pid);

      bool allgl1max = std::all_of(gl1diff.begin(), gl1diff.end(), [](uint64_t val) {
          return val == std::numeric_limits<uint64_t>::max();
        });
      if (allgl1max)
      {
        std::cout << Name() << " : GL1 clock diffs all filled with max 64 bit values for PID " << gl1pid << " return and try next pool" << std::endl;
        return;
      }
      m_DitchPackets.clear();

      for (auto& [pid, _] : m_PrevPoolLastDiffBad)
      {
        m_PrevPoolLastDiffBad[pid] = false;
      }

      for (const auto& [pid, sebdiff] : m_bclkdiffarray_map)
      {
        size_t packetpoolsize = m_PacketEventDeque[pid].size();
        if(packetpoolsize==0)
        {
          std::cout << Name() << ": packet pool size is zero.... something is wrong" << std::endl;
          return;
        }

        if(m_PacketAlignmentProblem[pid])
        {
          continue;
        }
        std::vector<int> bad_indices;
        int shift = 0;

        bool CurrentPoolLastDiffBad = false;
        bool PrevPoolLastDiffBad = m_PrevPoolLastDiffBad[pid];

        bool aligned = CheckPoolAlignment(pid, sebdiff, gl1diff, bad_indices, shift, CurrentPoolLastDiffBad, PrevPoolLastDiffBad);
        
        if (aligned)
        {
          m_PrevPoolLastDiffBad[pid] = CurrentPoolLastDiffBad;
          if (!bad_indices.empty())
          {
            std::cout << Name() << ": Packet " << pid << " has bad indices: ";
            for (int bi : bad_indices){
              std::cout << bi << " ";
              m_DitchPackets[pid].insert(bi);
            }
            std::cout << std::endl;

            std::cout << "full print out of gl1 vs seb clocks " << std::endl;
            for (size_t i = 0; i <= pooldepth; ++i)
            {
              uint64_t gl1_clk = Gl1Input()->m_bclkarray_map[gl1pid][i];
              uint64_t seb_clk = m_bclkarray_map[pid][i];
              std::cout << "pool index i " << i << ", gl1 / seb : " << gl1_clk << " / " << seb_clk;
              if(i<pooldepth){
                uint64_t gl1_diff = Gl1Input()->m_bclkdiffarray_map[gl1pid][i];
                uint64_t seb_diff = m_bclkdiffarray_map[pid][i];
                std::cout << " -> diff of gl1 vs seb : " << gl1_diff << " " << seb_diff << std::endl;
              }
              else if(i==pooldepth)
              {
                std::cout << std::endl;
              }
            }
          }

          if (shift == -1)
          {
            std::cout << Name() << ": Packet " << pid << " shifted by -1 with dropping the first seb event" << std::endl;
            if(m_PacketShiftOffset[pid] == -1)
            {
              std::cout << "Packet " << pid << " requires an additional shift -1. Lets not handle this for the moment.. stop combining" << std::endl;
              m_PacketAlignmentProblem[pid] = true;
            } 

            if (!m_PacketEventDeque[pid].empty())
            {
              m_PacketEventDeque[pid].pop_front();
            }
            else
            {
              std::cout << Name() << ": ERROR — shift -1 requested but packet deque is empty!" << std::endl;
              continue;
            }

            for (size_t i = 0; i < packetpoolsize - 1; ++i)
            {
              m_bclkarray_map[pid][i] = m_bclkarray_map[pid][i+1];
            }

            for (size_t i = 0; i < packetpoolsize; ++i)
            {
              m_bclkdiffarray_map[pid][i] = ComputeClockDiff(m_bclkarray_map[pid][i+1], m_bclkarray_map[pid][i]);
            }
            Event* evt = GetEventIterator()->getNextEvent();
            if (evt)
            {
              evt->convert();
              std::vector<Packet*> pktvec = evt->getPacketVector();
              for (Packet* pkt : pktvec)
              {
                if (pkt->getIdentifier() == pid)
                {
                  FillPacketClock(evt, pkt, packetpoolsize - 1);
                  m_PacketEventDeque[pid].push_back(evt);
                }
                delete pkt;
              }
            }
            else
            {
              std::cout << Name() << ": Cannot refill after shift -1" << std::endl;
              FilesDone(1);
              return;
            }
            m_PacketShiftOffset[pid] -= 1;
          }
          else if (shift == 1)
          {
            std::cout << Name() << ": Packet " << pid << " requires shift +1 (insert dummy at front)" << std::endl;
            
            if (m_packetclk_copy_runs)
            {
              std::cout << Name() << " : runs where clocks are copied from the first XMIT. Checking FEM clock diff" << std::endl;
              if (FemClockAlignment(pid, m_PacketEventDeque[pid], gl1diff))
              {
                std::cout << Name() << " : Packet identified as aligned with FEM clocks. Apply shift only on packet clocks." << std::endl;
                m_CorrectCopiedClockPackets.insert(pid);
                m_DitchPackets[pid].insert(0);

                Event* evt0 = m_PacketEventDeque[pid].front();
                m_PreviousValidBCOMap[pid] = GetClock(evt0, pid);
                m_bclkarray_map[pid][pooldepth] = m_bclkarray_map[pid][pooldepth - 1];
                continue;
              }
              std::cout << Name() << " : Packet identified as misaligned also with FEMs. Do normal recovery process" << std::endl;
            }

            if(m_PacketShiftOffset[pid] == 1)
            {
              std::cout << "Packet " << pid << " requires an additional shift +1. Lets not handle this for the moment.. stop combining" << std::endl;
              m_PacketAlignmentProblem[pid] = true;
            }

            for (size_t i = pooldepth; i > 0; --i)
            {
              m_bclkarray_map[pid][i] = m_bclkarray_map[pid][i-1];
            }
            for (size_t i = 1 ; i < pooldepth; ++i)
            {
              m_bclkdiffarray_map[pid][i] = ComputeClockDiff(m_bclkarray_map[pid][i+1], m_bclkarray_map[pid][i]);
            }

            m_bclkarray_map[pid][0] = 0;
            m_bclkdiffarray_map[pid][0] = 0;
            m_DitchPackets[pid].insert(0);

            if (!m_PacketEventDeque[pid].empty())
            {
              m_PacketEventBackup[pid] = m_PacketEventDeque[pid].back();
              Event* dummy_event = m_PacketEventDeque[pid][0]; 
              m_PacketEventDeque[pid].push_front(dummy_event);
              m_PacketEventDeque[pid].pop_back();
            }
            else
            {
              std::cout << Name() << ": m_PacketEventDeque is empty, cannot insert dummy event!" << std::endl;
              return;
            }

            m_PacketShiftOffset[pid] += 1;
            std::cout << std::endl;
          }
        }
        else
        {
          std::cout << Name() << ": Alignment failed for packet " << pid
            << " (retry count = " << m_PacketAlignmentFailCount[pid] << ")" << std::endl;
          std::cout << "full print out of gl1 vs seb clocks " << std::endl;
          for (size_t i = 0; i <= pooldepth; ++i)
          {
            uint64_t gl1_clk = Gl1Input()->m_bclkarray_map[gl1pid][i];
            uint64_t seb_clk = m_bclkarray_map[pid][i];
            std::cout << "pool index i " << i << ", gl1 / seb : " << gl1_clk << " / " << seb_clk;
            if(i<pooldepth){
              uint64_t gl1_diff = Gl1Input()->m_bclkdiffarray_map[gl1pid][i];
              uint64_t seb_diff = m_bclkdiffarray_map[pid][i];
              std::cout << " -- diff of gl1 vs seb : " << gl1_diff << " " << seb_diff << std::endl;
            }
            else if(i==pooldepth)
            {
              std::cout << std::endl;
            }
          }

          m_PacketAlignmentFailCount[pid]++;
          for (size_t i = 0; i < pooldepth; ++i)
          {
            m_DitchPackets[pid].insert(i);
          }

          if (m_PacketAlignmentFailCount[pid] >= m_max_alignment_retries)
          {
            std::cout << Name() << ": Max retries reached — permanently ditching packet " << pid << std::endl;
            m_PacketAlignmentFailCount[pid] = 0; 
            m_PacketAlignmentProblem[pid] = true;
          }
        }
      }
    }
  }
  return;
}

void SingleTriggeredInput::CreateDSTNodes(Event *evt)
{
  std::string CompositeNodeName = "Packets";
  if (KeepMyPackets())
  {
    CompositeNodeName = "PacketsKeep";
  }
  PHNodeIterator iter(m_topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    m_topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", CompositeNodeName));
  if (!detNode)
  {
    detNode = new PHCompositeNode(CompositeNodeName);
    dstNode->addNode(detNode);
  }
  std::vector<Packet *> pktvec = evt->getPacketVector();
  for (auto *piter : pktvec)
  {
    int packet_id = piter->getIdentifier();
    m_PacketSet.insert(packet_id);
    std::string PacketNodeName = std::to_string(packet_id);
    CaloPacket *calopacket = findNode::getClass<CaloPacket>(detNode, PacketNodeName);
    if (!calopacket)
    {
      calopacket = new CaloPacketv1();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(calopacket, PacketNodeName, "PHObject");
      detNode->addNode(newNode);
    }
    m_PacketShiftOffset.try_emplace(packet_id, 0);
    delete piter;
  }
}

bool SingleTriggeredInput::FemClockAlignment(int pid, const std::deque<Event*>& events, const std::array<uint64_t, pooldepth>& gl1diff)
{
  if (events.size() < pooldepth)
  {
    std::cout << Name() << ": Not enough events for FEMClockAlignment check for packet " << pid << std::endl;
    return false;
  }

  uint64_t prev_clk = std::numeric_limits<uint64_t>::max();

  for (size_t i = 0; i < pooldepth; ++i)
  {
    Event* evt = events[i];
    Packet* pkt = evt->getPacket(pid);
    if (!pkt)
    {
      continue;
    }

    int nmod = pkt->iValue(0, "NRMODULES");
    std::map<int, int> clk_count;

    for (int j = 0; j < nmod; ++j)
    {
      int femclk = static_cast<uint16_t>(pkt->iValue(j, "FEMCLOCK"));
      clk_count[femclk]++;
    }

    delete pkt;

    if (clk_count.empty())
    {
      continue;
    }

    int majority_clk = std::max_element(
        clk_count.begin(), clk_count.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; })->first;

    if (clk_count[majority_clk] < 2)
    {
      std::cout << Name() << ": FemClockAlignment — no majority FEM clocks for packet " << pid << " at pool index " << i << std::endl;
      return false;
    }


    if (i >= 1 && prev_clk != std::numeric_limits<uint64_t>::max() && gl1diff[i] != std::numeric_limits<uint64_t>::max())
    {
      uint16_t fem_diff = static_cast<uint16_t>(ComputeClockDiff(majority_clk, prev_clk) & 0xFFFFU);
      uint16_t gl1_diff = static_cast<uint16_t>(gl1diff[i] & 0xFFFFU);

      std::cout << "i " << i << " curr_fem_clk / prev_fem_clk " << majority_clk << " " << prev_clk << " , fem_diff / gl1_diff " << fem_diff << " " << gl1_diff << std::endl;
      if (fem_diff != gl1_diff)
      {
        return false;
      }
    }

    prev_clk = majority_clk;
  }

  return true;
}

int SingleTriggeredInput::FemEventNrClockCheck(OfflinePacket *pkt)
{
  CaloPacket *calopkt = dynamic_cast<CaloPacket *>(pkt);
  if (!calopkt)
  {
    return 0;
  }
  // make sure all clocks of the FEM are fine,
  int nrModules = calopkt->iValue(0, "NRMODULES");
  std::set<int> EventNoSet;
  for (int j = 0; j < nrModules; j++)
  {
    EventNoSet.insert(calopkt->iValue(j, "FEMEVTNR"));
  }
  size_t femeventnumbers = EventNoSet.size();
  if (femeventnumbers > 1)
  {
    int goodfemevent = 0;  // store the good event number so we can insert it in the set, but only if the clock counters agree
    if (femeventnumbers == 2)
    {
      // find the outlier if we have a 2:1 decision
      std::map<int, int> EventMap;
      std::map<int, int> BadModuleMap;
      for (int j = 0; j < nrModules; j++)
      {
        EventMap[calopkt->iValue(j, "FEMEVTNR")]++;
        BadModuleMap[calopkt->iValue(j, "FEMEVTNR")] = j;
      }
      for (const auto iter : EventMap)
      {
        if (iter.second == 1)
        {
          calopkt->setFemStatus(BadModuleMap[iter.first], CaloPacket::BAD_EVENTNR);
        }
        else
        {
          goodfemevent = iter.first;
        }
      }
    }
    else 
    {
      for (int j = 0; j < nrModules; j++)
      {
        calopkt->setFemStatus(j, CaloPacket::BAD_EVENTNR);
      }
    }
    std::set<int> FemClockSet;
    for (int j = 0; j < nrModules; j++)
    {
      FemClockSet.insert(calopkt->iValue(j, "FEMCLOCK"));
    }
    if (FemClockSet.size() == 1)
    {
      static int icnt = 0;
      if (icnt < 10)
      {
        icnt++;
        std::cout << "Packet " << calopkt->getIdentifier() << " has not unique event numbers"
                  << " but FEM Clock counters are identical" << std::endl;
      }
      if (goodfemevent > 0)
      {
        m_FEMEventNrSet.insert(goodfemevent);
      }
      return 1;
    }
    static int icnt = 0;
    if (icnt < 1000)
    {
      icnt++;
      std::cout << "resetting packet " << calopkt->getIdentifier()
                << " with fem event and clock mismatch" << std::endl;
      std::map<int, int> ClockMap;
      std::map<int, int> EventMap;
      for (int j = 0; j < nrModules; j++)
      {
        EventMap[calopkt->iValue(j, "FEMEVTNR")]++;
        ClockMap[calopkt->iValue(j, "FEMCLOCK")]++;
      }
      for (const auto iterA : EventMap)
      {
        std::cout << "Event Nr : " << iterA.first << " shows up " << iterA.second << " times"
                  << std::hex << ", Event Nr 0x" << iterA.first << std::dec << std::endl;
      }
      for (const auto iterA : ClockMap)
      {
        std::cout << "Clock : 0x" << std::hex << iterA.first << std::dec
                 << " shows up " << iterA.second << " times" << std::endl;
      }
    }
    return -1;
  }
  m_FEMEventNrSet.insert(*(EventNoSet.begin()));
  return 0;
}

void SingleTriggeredInput::dumpdeque()
{
  const auto *iter1 = clkdiffbegin();
  const auto *iter2 = Gl1Input()->clkdiffbegin();
  while (iter1 != clkdiffend())
  {
    std::cout << Name() << " clk: 0x" << std::hex << *iter1
              << " Gl1 clk: 0x" << *iter2 << std::dec << std::endl;
    iter1++;
    iter2++;
  }
  return;
}

int SingleTriggeredInput::ReadEvent()
{
  for (const auto& [pid, dq] : m_PacketEventDeque)
  {
    if (dq.empty())
    {
      if (!EventAlignmentProblem())
      {
        std::cout << Name() << ": Packet " << pid << " has empty deque — all events done" << std::endl;
        AllDone(1);
      }
      return -1;
    }
  }

  if (Verbosity() > 1)
  {
    size_t size = m_PacketEventDeque.begin()->second.size();
    std::cout << "deque size: " << size << std::endl;
  }

  auto *ref_evt = m_PacketEventDeque.begin()->second.front();
  RunNumber(ref_evt->getRunNumber());

  uint64_t event_number = ref_evt->getEvtSequence();
  if(event_number % 10000==0)
  {
    std::cout << "processed events : " << event_number << std::endl;
  }

  m_FEMEventNrSet.clear();

  bool all_packets_unshifted = std::all_of(
      m_PacketShiftOffset.begin(), m_PacketShiftOffset.end(),
      [](const std::pair<int, int>& p) { return p.second == 0; });

  std::set<Event*> events_to_delete;

  for (auto& [pid, dq] : m_PacketEventDeque)
  {
    if(m_PacketAlignmentProblem[pid]) 
    {
      continue;
    }
    Event* evt = dq.front();
    Packet* packet = evt->getPacket(pid);

    int packet_id = packet->getIdentifier();
    if (packet_id != pid)
    {
      std::cout << Name() << ": packet id mismatch... Should never happen. Abort combining" << std::endl;
      EventAlignmentProblem(1);
      delete packet;
      return -1;
    }

    CaloPacket *newhit = findNode::getClass<CaloPacket>(m_topNode, packet_id);
    newhit->Reset();

    if (m_DitchPackets.count(packet_id) && m_DitchPackets[packet_id].count(0))
    {
      newhit->setStatus(OfflinePacket::PACKET_DROPPED);
      newhit->setIdentifier(packet_id);
      std::cout << "ditching packet " << packet_id << " from prdf event " << evt->getEvtSequence() << std::endl;
      delete packet;
      continue;
    }

    newhit->setStatus(OfflinePacket::PACKET_OK);
    if (m_OverrideWithRepClock.count(packet_id))
    {
      newhit->setStatus(OfflinePacket::PACKET_CORRUPT);
    }
    newhit->setPacketEvtSequence(packet->iValue(0, "EVTNR"));
    int nr_modules = packet->iValue(0, "NRMODULES");
    int nr_channels = packet->iValue(0, "CHANNELS");
    int nr_samples = packet->iValue(0, "SAMPLES");
    newhit->setNrModules(nr_modules);
    newhit->setNrChannels(nr_channels);
    newhit->setNrSamples(nr_samples);
    newhit->setIdentifier(packet_id);
    if (m_packetclk_copy_runs && m_CorrectCopiedClockPackets.count(packet_id))
    {
      uint64_t prev_packet_clock = m_PreviousValidBCOMap[packet_id];
      newhit->setBCO(prev_packet_clock);
      m_PreviousValidBCOMap[packet_id] = GetClock(evt,packet_id);
    }
    else
    {
      newhit->setBCO(packet->lValue(0, "CLOCK"));
    }

    for (int ifem = 0; ifem < nr_modules; ifem++)
    {
      newhit->setFemClock(ifem, packet->iValue(ifem, "FEMCLOCK"));
      newhit->setFemEvtSequence(ifem, packet->iValue(ifem, "FEMEVTNR"));
      newhit->setFemSlot(ifem, packet->iValue(ifem, "FEMSLOT"));
      newhit->setChecksumLsb(ifem, packet->iValue(ifem, "CHECKSUMLSB"));
      newhit->setChecksumMsb(ifem, packet->iValue(ifem, "CHECKSUMMSB"));
      newhit->setCalcChecksumLsb(ifem, packet->iValue(ifem, "CALCCHECKSUMLSB"));
      newhit->setCalcChecksumMsb(ifem, packet->iValue(ifem, "CALCCHECKSUMMSB"));
      newhit->setFemStatus(ifem, CaloPacket::FEM_OK);
    }
    for (int ipmt = 0; ipmt < nr_channels; ipmt++)
    {
      bool isSuppressed = packet->iValue(ipmt, "SUPPRESSED");
      newhit->setSuppressed(ipmt, isSuppressed);
      if (isSuppressed)
      {
        newhit->setPre(ipmt, packet->iValue(ipmt, "PRE"));
        newhit->setPost(ipmt, packet->iValue(ipmt, "POST"));
      }
      else
      {
        for (int isamp = 0; isamp < nr_samples; isamp++)
        {
          newhit->setSample(ipmt, isamp, packet->iValue(isamp, ipmt));
        }
      }
    }
    delete packet;
    int iret = FemEventNrClockCheck(newhit);
    if (iret < 0)
    {
      std::cout << Name() <<" : failed on FemEventNrClockCheck reset calo packet " << std::endl;
      newhit->Reset();
    }

    if (all_packets_unshifted || m_PacketShiftOffset[pid] == 1)
    {
      events_to_delete.insert(evt);
    }
  }

  for(Event *evtdelete : events_to_delete)
  {
    delete evtdelete;
  }

  for (auto& [pid, idx_set] : m_DitchPackets)
  {
    std::set<int> new_set;
    for (int idx : idx_set)
    {
      if (idx > 0)
      {
        new_set.insert(idx - 1);
      }
    }
    idx_set = std::move(new_set);
  }

  for (auto& [pid, dq] : m_PacketEventDeque)
  {
    if (!dq.empty())
    {
      dq.pop_front();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

