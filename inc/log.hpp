#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>

#ifndef HPP_LOG
#define HPP_LOG
namespace logging = boost::log;

/**
 * Sets up default logging at info and above
 * and allows for logging at debug and above
 * @param log_level must be debug to be useful, or left blank
 */
void init_logging(const std::string& log_level = "");

#endif
