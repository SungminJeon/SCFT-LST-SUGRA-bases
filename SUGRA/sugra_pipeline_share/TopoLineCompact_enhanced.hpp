#pragma once
#include "Topology_enhanced.h"
#include <string>
#include <vector>
#include <sstream>
#include <cctype>
#include <tuple>

// Enhanced line-compact format with External curve support
// Format: kinds | bparams | S=(u,v);... | I=(u,v);... | sp=... | ip=... | E=(parent,type,port,eid);... | ep=...
//
// E connections: (parent_id, parent_type, port_idx, external_id)
//   - parent_id: index of parent object
//   - parent_type: 0=block, 1=sidelink, 2=instanton
//   - port_idx: port index on parent object
//   - external_id: index of external curve
// ep: external parameters CSV

class TopoLineCompact_enhanced {
public:
    static std::string serialize(const Topology_enhanced& T);
    static bool deserialize(const std::string& line, Topology_enhanced& out);
    
    // Validation
    static bool validate(const Topology_enhanced& T);
    static std::string getValidationErrors(const Topology_enhanced& T);
    
private:
    // Helper functions
    static int kindToInt(LKind k);
    static LKind intToKind(int k);
    static void trim(std::string& s);
    static std::vector<std::string> split(const std::string& s, char sep);
    static std::vector<int> parse_csv_ints(const std::string& s);
    static std::vector<std::pair<int,int>> parse_pairs(const std::string& s);
    static std::vector<std::tuple<int,int,int,int>> parse_quads(const std::string& s);
};
