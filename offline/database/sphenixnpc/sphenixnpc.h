#ifndef DBTOOLS_SPHENIXNPC_H
#define DBTOOLS_SPHENIXNPC_H


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <nopayloadclient/nopayloadclient.hpp>
#pragma GCC diagnostic pop

#include <nlohmann/json.hpp>

#include <iostream>

using json = nlohmann::json;

//namespace sphenixnpc {

class sphenixnpc :  public nopayloadclient::Client
{

public: 
  static sphenixnpc *instance(const std::string &globaltag = "NONE");
   ~sphenixnpc(); 
    json getUrlDict(long long iov);
    json createGlobalTag(const std::string &tagname);
    
    json deleteGlobalTag(const std::string&);
    json get(std::string pl_type, long long iov);
    json insertPayload(std::string pl_type, std::string file_url,
                       long long iov_start);
    json insertPayload(std::string pl_type, std::string file_url,
                       long long iov_start, long long iov_end);
    json setGlobalTag(std::string name);
    json clearCache();
    std::string getCalibrationFile(const std::string &type, uint64_t iov);
    int insertcalib(const std::string &fname, const std::string &payloadtype, uint64_t iov_start);
    int insertcalib(const std::string &fname, const std::string &payloadtype, uint64_t iov_start, uint64_t iov_end);
    void Verbosity(int i) {m_Verbosity = i;}
      int Verbosity() const {return m_Verbosity;}

private:
    static sphenixnpc *__instance;
    int m_Verbosity = 0;
    json url_dict_; // valid until global tag is switched
    std::string m_CachedGlobalTag;
};

//}
#endif // DBTOOLS_SPHENIXNPC_H
