#include "TopoLineCompact_enhanced.hpp"
#include <algorithm>
#include <stdexcept>

// ===== Helper functions =====
int TopoLineCompact_enhanced::kindToInt(LKind k) {
    switch (k) {
        case LKind::g: return 0;
        case LKind::L: return 1;
        case LKind::S: return 2;
        case LKind::I: return 3;
        case LKind::E: return 4;  // Enhanced: External
    }
    return -1;
}

LKind TopoLineCompact_enhanced::intToKind(int k) {
    switch (k) {
        case 0: return LKind::g;
        case 1: return LKind::L;
        case 2: return LKind::S;
        case 3: return LKind::I;
        case 4: return LKind::E;  // Enhanced: External
    }
    return LKind::g;
}

void TopoLineCompact_enhanced::trim(std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) ++a;
    while (b > a && std::isspace((unsigned char)s[b-1])) --b;
    s = s.substr(a, b - a);
}

std::vector<std::string> TopoLineCompact_enhanced::split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

std::vector<int> TopoLineCompact_enhanced::parse_csv_ints(const std::string& s) {
    std::vector<int> v;
    if (s.empty()) return v;
    for (auto& t : split(s, ',')) {
        std::string u = t;
        trim(u);
        if (!u.empty()) v.push_back(std::stoi(u));
    }
    return v;
}

std::vector<std::pair<int,int>> TopoLineCompact_enhanced::parse_pairs(const std::string& s) {
    // "(u,v);(u,v)" format
    std::vector<std::pair<int,int>> v;
    if (s.empty()) return v;
    
    for (auto& t : split(s, ';')) {
        std::string u = t;
        trim(u);
        if (u.empty()) continue;
        
        // Remove parentheses
        if (u.front() == '(' && u.back() == ')') {
            u = u.substr(1, u.size() - 2);
        }
        
        auto xy = split(u, ',');
        if (xy.size() != 2) continue;
        
        v.emplace_back(std::stoi(xy[0]), std::stoi(xy[1]));
    }
    return v;
}

std::vector<std::tuple<int,int,int,int>> TopoLineCompact_enhanced::parse_quads(const std::string& s) {
    // "(a,b,c,d);(a,b,c,d)" format for External connections
    std::vector<std::tuple<int,int,int,int>> v;
    if (s.empty()) return v;
    
    for (auto& t : split(s, ';')) {
        std::string u = t;
        trim(u);
        if (u.empty()) continue;
        
        // Remove parentheses
        if (u.front() == '(' && u.back() == ')') {
            u = u.substr(1, u.size() - 2);
        }
        
        auto parts = split(u, ',');
        if (parts.size() != 4) continue;
        
        v.emplace_back(
            std::stoi(parts[0]),  // parent_id
            std::stoi(parts[1]),  // parent_type
            std::stoi(parts[2]),  // port_idx
            std::stoi(parts[3])   // external_id
        );
    }
    return v;
}

// ===== Serialization =====
std::string TopoLineCompact_enhanced::serialize(const Topology_enhanced& T) {
    // Build kinds and bparams
    std::ostringstream kinds, bparams;
    for (size_t i = 0; i < T.block.size(); ++i) {
        if (i) { 
            kinds << ","; 
            bparams << ","; 
        }
        kinds << kindToInt(T.block[i].kind);
        bparams << T.block[i].param;
    }
    
    // S/I connections
    std::ostringstream Sc, Ic;
    for (size_t i = 0; i < T.s_connection.size(); ++i) {
        if (i) Sc << ";";
        Sc << "(" << T.s_connection[i].u << "," << T.s_connection[i].v << ")";
    }
    for (size_t i = 0; i < T.i_connection.size(); ++i) {
        if (i) Ic << ";";
        Ic << "(" << T.i_connection[i].u << "," << T.i_connection[i].v << ")";
    }
    
    // S/I parameters
    std::ostringstream Sp, Ip;
    for (size_t i = 0; i < T.side_links.size(); ++i) {
        if (i) Sp << ",";
        Sp << T.side_links[i].param;
    }
    for (size_t i = 0; i < T.instantons.size(); ++i) {
        if (i) Ip << ",";
        Ip << T.instantons[i].param;
    }
    
    // Enhanced: External connections and parameters
    std::ostringstream Ec, Ep;
    for (size_t i = 0; i < T.e_connection.size(); ++i) {
        if (i) Ec << ";";
        const auto& e = T.e_connection[i];
        Ec << "(" << e.parent_id << "," << e.parent_type << "," 
           << e.port_idx << "," << e.external_id << ")";
    }
    for (size_t i = 0; i < T.externals.size(); ++i) {
        if (i) Ep << ",";
        Ep << T.externals[i].param;
    }
    
    // Build final output
    std::ostringstream out;
    out << kinds.str() << " | "
        << bparams.str() << " | "
        << "S=" << Sc.str() << " | "
        << "I=" << Ic.str() << " | "
        << "sp=" << Sp.str() << " | "
        << "ip=" << Ip.str();
    
    // Only add External fields if there are externals
    if (!T.externals.empty() || !T.e_connection.empty()) {
        out << " | E=" << Ec.str() << " | ep=" << Ep.str();
    }
    
    return out.str();
}

// ===== Deserialization =====
bool TopoLineCompact_enhanced::deserialize(const std::string& line, Topology_enhanced& out) {
    auto parts = split(line, '|');
    if (parts.size() < 6) return false;  // Minimum required fields
    
    for (auto& p : parts) trim(p);

    // Parse kinds and params
    auto kinds_csv = parts[0];
    auto bparams_csv = parts[1];
    auto kinds_vec = parse_csv_ints(kinds_csv);
    auto bparams_vec = parse_csv_ints(bparams_csv);
    
    if (kinds_vec.size() != bparams_vec.size()) return false;

    // Initialize topology
    out.Initialize();
    for (size_t i = 0; i < kinds_vec.size(); ++i) {
        const LKind k = intToKind(kinds_vec[i]);
        // Chain assumption: connect to the right
        if (i == 0) {
            out.addBlock(k, bparams_vec[i]);
        } else {
            out.addBlockRight(k, bparams_vec[i]);
        }
    }

    // Helper lambda to read fields
    auto rd = [&](const std::string& tag) -> std::string {
        for (size_t i = 2; i < parts.size(); ++i) {
            auto s = parts[i];
            if (s.rfind(tag, 0) == 0) {
                auto t = s.substr(tag.size());
                trim(t);
                return t;
            }
        }
        return std::string();
    };

    // Parse connections and parameters
    const std::string Sraw = rd("S=");
    const std::string Iraw = rd("I=");
    const std::string SPar = rd("sp=");
    const std::string IPar = rd("ip=");
    const std::string Eraw = rd("E=");   // Enhanced: External connections
    const std::string EPar = rd("ep=");   // Enhanced: External parameters

    // Parse parameters first
    auto sp = parse_csv_ints(SPar);
    auto ip = parse_csv_ints(IPar);
    auto ep = parse_csv_ints(EPar);  // Enhanced
    
    // Add side links and instantons
    for (int v : sp) out.addSideLink(v);
    for (int v : ip) out.addInstanton(v);
    for (int v : ep) out.addExternal(v);  // Enhanced
    
    // Parse connections
    std::vector<SideLinkStructure> s_conn;
    for (auto [u, sid] : parse_pairs(Sraw)) {
        s_conn.push_back({u, sid});
    }
    out.LinkingSideLinks(std::move(s_conn));
    
    std::vector<InstantonStructure> i_conn;
    for (auto [u, iid] : parse_pairs(Iraw)) {
        i_conn.push_back({u, iid});
    }
    out.LinkingInstantonStructure(std::move(i_conn));
    
    // Enhanced: Parse external connections
    if (!Eraw.empty()) {
        std::vector<ExternalStructure> e_conn;
        for (auto [parent_id, parent_type, port_idx, eid] : parse_quads(Eraw)) {
            e_conn.push_back({parent_id, parent_type, port_idx, eid});
        }
        out.LinkingExternals(std::move(e_conn));
    }
    
    return true;
}

// ===== Validation =====
bool TopoLineCompact_enhanced::validate(const Topology_enhanced& T) {
    // Validate all connections reference valid objects
    
    // Check S-connections
    for (const auto& sc : T.s_connection) {
        if (sc.u < 0 || sc.u >= static_cast<int>(T.block.size())) return false;
        if (sc.v < 0 || sc.v >= static_cast<int>(T.side_links.size())) return false;
    }
    
    // Check I-connections
    for (const auto& ic : T.i_connection) {
        if (ic.u < 0 || ic.u >= static_cast<int>(T.block.size())) return false;
        if (ic.v < 0 || ic.v >= static_cast<int>(T.instantons.size())) return false;
    }
    
    // Check E-connections (Enhanced)
    for (const auto& ec : T.e_connection) {
        // Check external exists
        if (ec.external_id < 0 || ec.external_id >= static_cast<int>(T.externals.size())) {
            return false;
        }
        
        // Check parent exists based on type
        switch (ec.parent_type) {
            case 0:  // Block
                if (ec.parent_id < 0 || ec.parent_id >= static_cast<int>(T.block.size())) {
                    return false;
                }
                break;
            case 1:  // SideLink
                if (ec.parent_id < 0 || ec.parent_id >= static_cast<int>(T.side_links.size())) {
                    return false;
                }
                break;
            case 2:  // Instanton
                if (ec.parent_id < 0 || ec.parent_id >= static_cast<int>(T.instantons.size())) {
                    return false;
                }
                break;
            default:
                return false;
        }
    }
    
    return true;
}

std::string TopoLineCompact_enhanced::getValidationErrors(const Topology_enhanced& T) {
    std::ostringstream errors;
    
    // Check S-connections
    for (const auto& sc : T.s_connection) {
        if (sc.u < 0 || sc.u >= static_cast<int>(T.block.size())) {
            errors << "S-connection: Invalid block index " << sc.u << "\n";
        }
        if (sc.v < 0 || sc.v >= static_cast<int>(T.side_links.size())) {
            errors << "S-connection: Invalid sidelink index " << sc.v << "\n";
        }
    }
    
    // Check I-connections
    for (const auto& ic : T.i_connection) {
        if (ic.u < 0 || ic.u >= static_cast<int>(T.block.size())) {
            errors << "I-connection: Invalid block index " << ic.u << "\n";
        }
        if (ic.v < 0 || ic.v >= static_cast<int>(T.instantons.size())) {
            errors << "I-connection: Invalid instanton index " << ic.v << "\n";
        }
    }
    
    // Check E-connections (Enhanced)
    for (const auto& ec : T.e_connection) {
        if (ec.external_id < 0 || ec.external_id >= static_cast<int>(T.externals.size())) {
            errors << "E-connection: Invalid external_id " << ec.external_id << "\n";
        }
        
        switch (ec.parent_type) {
            case 0:
                if (ec.parent_id < 0 || ec.parent_id >= static_cast<int>(T.block.size())) {
                    errors << "E-connection: Invalid block parent_id " << ec.parent_id << "\n";
                }
                break;
            case 1:
                if (ec.parent_id < 0 || ec.parent_id >= static_cast<int>(T.side_links.size())) {
                    errors << "E-connection: Invalid sidelink parent_id " << ec.parent_id << "\n";
                }
                break;
            case 2:
                if (ec.parent_id < 0 || ec.parent_id >= static_cast<int>(T.instantons.size())) {
                    errors << "E-connection: Invalid instanton parent_id " << ec.parent_id << "\n";
                }
                break;
            default:
                errors << "E-connection: Invalid parent_type " << ec.parent_type << "\n";
        }
    }
    
    return errors.str();
}
