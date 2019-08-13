#include "log.hpp"

void init_logging(const std::string& log_level){
    auto g_log_level{logging::trivial::info};
    if (log_level == "debug")
        g_log_level = logging::trivial::debug;
    boost::log::core::get()->set_filter
            (
                    logging::trivial::severity >= g_log_level
            );
}
