#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem.hpp>
#include <sdsl/suffix_arrays.hpp>

#ifndef HPP_UTILS
#define HPP_UTILS
namespace logging = boost::log;
namespace filesystem = boost::filesystem;

/**
 * Sets up default logging at info and above
 * and allows for logging at debug and above
 * @param log_level must be debug to be useful, or left blank
 */
void init_logging(const std::string& log_level = "");

/**
 * Convert a nucleotide character into its integer representation.
 */
int encode_char(const char &c);

#endif
