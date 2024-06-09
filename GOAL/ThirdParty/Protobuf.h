////////////////////////////////
/// usage : 1.	data format converters.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_THIRD_PARTY_PROTOBUF_H
#define CN_HUST_GOAL_THIRD_PARTY_PROTOBUF_H


#include "GOAL/Flag.h"

#if _PLUGIN_PROTOBUF

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#pragma warning(push)
#pragma warning(disable: 26495) // Warning C26495 Variable is uninitialized.Always initialize a member variable(type.6).
#pragma warning(disable: 26451) // Warning C26451 Arithmetic overflow : Using operator '*' on a 4 byte value and then casting the result to a 8 byte value.Cast the value to the wider type before calling operator '*' to avoid overflow(io.2).
#include <google/protobuf/util/json_util.h>
#include <google/protobuf/io/gzip_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#pragma warning(pop)


namespace goal {

template<typename T>
std::string protobufToJson(const T &obj, bool pretty = false) {
    std::string data;

    google::protobuf::util::JsonPrintOptions options;
    options.add_whitespace = pretty;
    options.always_print_primitive_fields = true;
    options.preserve_proto_field_names = true;

    google::protobuf::util::MessageToJsonString(obj, &data, options);

    return data;
}

template<typename T>
void jsonToProtobuf(const std::string &data, T &obj) {
    google::protobuf::util::JsonParseOptions options;
    google::protobuf::util::JsonStringToMessage(data, &obj, options);
}

template<typename T>
bool loadJson(const std::string &path, T &obj) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) { return false; }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    jsonToProtobuf(oss.str(), obj);
    return true;
}

template<typename T>
bool saveJson(const std::string &path, const T &obj) {
    std::ofstream ofs(path);
    if (!ofs.is_open()) { return false; }
    ofs << protobufToJson(obj);
    return true;
}

template<typename T>
bool loadProtobufZip(const std::string &path, T &obj) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open()) { return false; }
    google::protobuf::io::IstreamInputStream iis(&ifs);

    google::protobuf::io::GzipInputStream gis(&iis, google::protobuf::io::GzipInputStream::Format::GZIP);

    return obj.ParseFromZeroCopyStream(&gis);
    return true;
}

template<typename T>
bool saveProtobufZip(const std::string &path, const T &obj) {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.is_open()) { return false; }
    google::protobuf::io::OstreamOutputStream fos(&ofs);
    
    google::protobuf::io::GzipOutputStream::Options options;
    options.compression_level = 9;
    options.format = google::protobuf::io::GzipOutputStream::Format::GZIP;
    google::protobuf::io::GzipOutputStream gos(&fos, options);
    
    return obj.SerializeToZeroCopyStream(&gos);
    return true;
}

template<typename T>
bool loadProtobuf(const std::string& path, T& obj) {
    #if _PLUGIN_ZLIB
    return loadProtobufZip(path, obj);
    #else
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.is_open()) { return false; }
    return obj.ParseFromIstream(&ifs);
    #endif // _PLUGIN_ZLIB
}

template<typename T>
bool saveProtobuf(const std::string& path, const T& obj) {
    #if _PLUGIN_ZLIB
    return saveProtobufZip(path, obj);
    #else
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.is_open()) { return false; }
    return obj.SerializeToOstream(&ofs);
    #endif // _PLUGIN_ZLIB
}

}

#endif // _PLUGIN_PB


#endif // CN_HUST_GOAL_THIRD_PARTY_PROTOBUF_H
