#include "sphenixnpc.h"

sphenixnpc *sphenixnpc::__instance = nullptr;

sphenixnpc *sphenixnpc::instance(const std::string &globaltag)
{
  if (! __instance) 
  {
    __instance = new sphenixnpc();
  }
  __instance->setGlobalTag(globaltag);
  return __instance;
}

sphenixnpc::~sphenixnpc()
{
  clearCache();
  __instance = nullptr;
}

json sphenixnpc::createGlobalTag(const std::string &tagname)
{
  setGlobalTag(tagname);
  return nopayloadclient::Client::createGlobalTag();
}

json sphenixnpc::deleteGlobalTag(const std::string &tagname)
{
  setGlobalTag(tagname);
  return nopayloadclient::Client::deleteGlobalTag();
}

json sphenixnpc::getUrlDict(long long iov) {
    return nopayloadclient::Client::getUrlDict(0, iov);
}

json sphenixnpc::get(std::string pl_type, long long iov) {
    if (url_dict_.is_null()) {
        json resp = getUrlDict(iov);
        if (resp["code"] != 0) return resp;
        url_dict_ = resp["msg"];
    }
    if (not url_dict_.contains(pl_type)) {
        return DataBaseException("No payload with type " + pl_type + " exists.").jsonify();
    }
    return makeResp(url_dict_[pl_type]);
}

json sphenixnpc::insertPayload(std::string pl_type, std::string file_url,
                           long long iov_start) {
    return nopayloadclient::Client::insertPayload(pl_type, file_url, 0, iov_start);
}

json sphenixnpc::insertPayload(std::string pl_type, std::string file_url,
                           long long iov_start, long long iov_end) {
    return nopayloadclient::Client::insertPayload(pl_type, file_url, 0, iov_start, 0, iov_end);
}

json sphenixnpc::setGlobalTag(std::string name) {
  if (name != m_CachedGlobalTag)
  {
    url_dict_ = json {};
    m_CachedGlobalTag = name;
  }
    return nopayloadclient::Client::setGlobalTag(name);
}

json sphenixnpc::clearCache() {
    url_dict_ = json {};
    return nopayloadclient::Client::clearCache();
}

std::string sphenixnpc::getCalibrationFile(const std::string &type, uint64_t iov)
{
  nlohmann::json result = get(type, iov);
  return result.at("msg");
}

int sphenixnpc::insertcalib(const std::string &payloadtype, const std::string &fname, uint64_t iov_start)
{
  json ret = insertPayload(payloadtype, fname, iov_start);
  if (Verbosity())
  {
    std::cout << ret << std::endl;
  }
  return 0;
}

int sphenixnpc::insertcalib(const std::string &payloadtype, const std::string &fname, uint64_t iov_start, uint64_t iov_end)
{
  json ret = insertPayload(payloadtype, fname, iov_start, iov_end);
  if (Verbosity())
  {
    std::cout << ret << std::endl;
  }
  return 0;
}
