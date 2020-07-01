#ifndef PHOOL_PHLOG_H
#define PHOOL_PHLOG_H

#include <log4cpp/Appender.hh>
#include <log4cpp/Category.hh>
#include <log4cpp/FileAppender.hh>
#include <log4cpp/OstreamAppender.hh>
#include <log4cpp/PatternLayout.hh>
#include <log4cpp/Priority.hh>
#include <log4cpp/PropertyConfigurator.hh>

#include <pwd.h>
#include <sys/types.h>
#include <unistd.h>

class PHLog
{
 public:
  static void Init()
  {
    const std::string initFileName = "phlog.conf";
    const char* env_file = std::getenv("PHLOGCONF");
    if (_file_exists(initFileName))
    {
      log4cpp::PropertyConfigurator::configure(initFileName);
    }
    else if (_file_exists(std::string(_get_homedir()) + "/" + initFileName))
    {
      log4cpp::PropertyConfigurator::configure(std::string(_get_homedir()) + "/" + initFileName);
    }
    else if (env_file && _file_exists(env_file))
    {
      log4cpp::PropertyConfigurator::configure(env_file);
    }
    else
    {
      log4cpp::PatternLayout* layout1 = new log4cpp::PatternLayout();
      layout1->setConversionPattern("%d [%p] %c: %m%n");
      log4cpp::Appender* appender1 = new log4cpp::OstreamAppender("console", &std::cout);
      appender1->setLayout(layout1);
      log4cpp::Category& root = log4cpp::Category::getRoot();
      root.setPriority(log4cpp::Priority::FATAL);
      root.addAppender(appender1);
    }
  };

  static void Deinit()
  {
    log4cpp::Category::shutdown();
  };

  static log4cpp::Category& get(const std::string& category = "")
  {
    static bool phlog_is_initialized = false;
    if (!phlog_is_initialized)
    {
      Init();
      phlog_is_initialized = true;
    }
    return log4cpp::Category::getInstance(category);
  }

 private:
  static const char* _get_homedir()
  {
    struct passwd* pw = getpwuid(getuid());
    return pw->pw_dir;
  };

  static inline bool _file_exists(const std::string& filename)
  {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
  };
};

// conditional logging

#define LOG_DEBUG_IF(category, cond)                                              \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::DEBUG) && (cond)) \
  PHLog::get(category) << log4cpp::Priority::DEBUG

#define LOG_INFO_IF(category, cond)                                              \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::INFO) && (cond)) \
  PHLog::get(category) << log4cpp::Priority::INFO

#define LOG_NOTICE_IF(category, cond)                                              \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::NOTICE) && (cond)) \
  PHLog::get(category) << log4cpp::Priority::NOTICE

#define LOG_WARN_IF(category, cond)                                              \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::WARN) && (cond)) \
  PHLog::get(category) << log4cpp::Priority::WARN

#define LOG_ERROR_IF(category, cond)                                              \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::ERROR) && (cond)) \
  PHLog::get(category) << log4cpp::Priority::ERROR

#define LOG_CRIT_IF(category, cond)                                              \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::CRIT) && (cond)) \
  PHLog::get(category) << log4cpp::Priority::CRIT

#define LOG_ALERT_IF(category, cond)                                              \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::ALERT) && (cond)) \
  PHLog::get(category) << log4cpp::Priority::ALERT

#define LOG_FATAL_IF(category, cond)                                              \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::FATAL) && (cond)) \
  PHLog::get(category) << log4cpp::Priority::FATAL

// non-conditional logging

#define LOG_DEBUG(category)                                             \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::DEBUG)) \
  PHLog::get(category) << log4cpp::Priority::DEBUG

#define LOG_INFO(category)                                             \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::INFO)) \
  PHLog::get(category) << log4cpp::Priority::INFO

#define LOG_NOTICE(category)                                             \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::NOTICE)) \
  PHLog::get(category) << log4cpp::Priority::NOTICE

#define LOG_WARN(category)                                             \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::WARN)) \
  PHLog::get(category) << log4cpp::Priority::WARN

#define LOG_ERROR(category)                                             \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::ERROR)) \
  PHLog::get(category) << log4cpp::Priority::ERROR

#define LOG_CRIT(category)                                             \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::CRIT)) \
  PHLog::get(category) << log4cpp::Priority::CRIT

#define LOG_ALERT(category)                                             \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::ALERT)) \
  PHLog::get(category) << log4cpp::Priority::ALERT

#define LOG_FATAL(category)                                             \
  if (PHLog::get(category).isPriorityEnabled(log4cpp::Priority::FATAL)) \
  PHLog::get(category) << log4cpp::Priority::FATAL

#endif  // PHOOL_PHLOG_H
