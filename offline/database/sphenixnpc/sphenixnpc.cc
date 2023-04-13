#include "sphenixnpc.h"

sphenixnpc *sphenixnpc::__instance = nullptr;

sphenixnpc *sphenixnpc::instance(const std::string &globaltag)
{
  if (!__instance)
  {
    __instance = new sphenixnpc();
  }
  __instance->cache_set_GlobalTag(globaltag);
  return __instance;
}

sphenixnpc::~sphenixnpc()
{
  __instance = nullptr;
}

int sphenixnpc::createGlobalTag(const std::string &tagname)
{
  setGlobalTag(tagname);
  nlohmann::json result = nopayloadclient::Client::createGlobalTag();
  if (Verbosity())
  {
    std::cout << result << std::endl;
  }
  return 0;
}

int sphenixnpc::deleteGlobalTag(const std::string &tagname)
{
  setGlobalTag(tagname);
  nlohmann::json result = nopayloadclient::Client::deleteGlobalTag();
  if (Verbosity())
  {
    std::cout << result << std::endl;
  }
  return 0;
}

nlohmann::json sphenixnpc::getUrlDict(long long iov)
{
  return nopayloadclient::Client::getUrlDict(0, iov);
}

nlohmann::json sphenixnpc::get(const std::string &pl_type, long long iov)
{
  if (url_dict_.is_null())
  {
    nlohmann::json resp = getUrlDict(iov);
    if (resp["code"] != 0)
    {
      return resp;
    }
    url_dict_ = resp["msg"];
  }
  if (not url_dict_.contains(pl_type))
  {
    return nopayloadclient::DataBaseException("No payload with type " + pl_type + " exists.").jsonify();
  }
  return makeResp(url_dict_[pl_type]);
}

nlohmann::json sphenixnpc::cache_set_GlobalTag(const std::string &name)
{
  if (name != m_CachedGlobalTag)
  {
    url_dict_ = nlohmann::json{};
    m_CachedGlobalTag = name;
  }
  return nopayloadclient::Client::setGlobalTag(name);
}

nlohmann::json sphenixnpc::clearCache()
{
  url_dict_ = nlohmann::json{};
  return nopayloadclient::Client::clearCache();
}

std::string sphenixnpc::getCalibrationFile(const std::string &type, uint64_t iov)
{
  nlohmann::json result = get(type, iov);
  return result.at("msg");
}

int sphenixnpc::insertcalib(const std::string &pl_type, const std::string &file_url, uint64_t iov_start)
{
  nlohmann::json ret = nopayloadclient::Client::insertPayload(pl_type, file_url, 0, iov_start);
  if (Verbosity())
  {
    std::cout << ret << std::endl;
  }
  return 0;
}

int sphenixnpc::insertcalib(const std::string &pl_type, const std::string &file_url, uint64_t iov_start, uint64_t iov_end)
{
  nlohmann::json ret = nopayloadclient::Client::insertPayload(pl_type, file_url, 0, iov_start, 0, iov_end);
  if (Verbosity())
  {
    std::cout << ret << std::endl;
  }
  return 0;
}
