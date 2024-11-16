#include "../include/client_logger_builder.h"
#include "C:/OpSys_cw/allocator/common/include/not_implemented.h"
#include "../include/client_logger.h"
#include "../include/client_logger_builder.h"
#include <fstream>
#include <vector>


client_logger_builder::client_logger_builder()
{
    _format = "%s %m %t %d";
}

client_logger_builder::client_logger_builder(std::string const &format) : _format(format)
{

}

client_logger_builder::client_logger_builder(
    client_logger_builder const &other) : streams(other.streams), _format(other._format), _console_sev(other._console_sev)

{

}

client_logger_builder &client_logger_builder::operator=(
    client_logger_builder const &other)
{
    if (this != &other) {
        streams = other.streams;
        _format = other._format;
        _console_sev = other._console_sev;
    }

    return *this;
}

client_logger_builder::client_logger_builder(
    client_logger_builder &&other) noexcept : streams(std::move(other.streams)), _format(std::move(other._format)), _console_sev(std::move(other._console_sev))
{

}

client_logger_builder &client_logger_builder::operator=(
    client_logger_builder &&other) noexcept
{
    if (this != &other) {
        streams = std::move(other.streams);
        _format = std::move(other._format);
        _console_sev = std::move(other._console_sev);
    }

    return *this;
}

client_logger_builder::~client_logger_builder() noexcept
{

}

logger_builder *client_logger_builder::add_file_stream(
    std::string const &stream_file_path,
    logger::severity severity)
{
    auto absolute_path = std::filesystem::absolute(stream_file_path);
    streams[absolute_path.string()].emplace(severity);

    return this;
}

logger_builder *client_logger_builder::add_console_stream(
    logger::severity severity)
{
    _console_sev.emplace(severity);

    return this;
}

logger_builder *client_logger_builder::clear()
{
    streams.clear();
    _console_sev.clear();
    return this;
}

logger *client_logger_builder::build() const
{
    client_logger* client_log = new client_logger(_format);
    
    auto map_iterator = streams.begin();
    while (map_iterator != streams.end())
    {
        auto set_iterator = (*map_iterator).second.begin();
        while (set_iterator != (*map_iterator).second.end())
        {
            (*client_log).add_file_stream_to_log_map_streams((*map_iterator).first, (*set_iterator));
            ++set_iterator;
        }
        ++map_iterator;
    }

    for (auto &sev: _console_sev)
    {
        (*client_log).add_console_stream_to_set(sev);
    }

    return client_log;
}