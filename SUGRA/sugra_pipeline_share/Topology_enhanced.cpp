#include "Topology_enhanced.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>

// ========== Helper functions ==========
const char* LKindToString(LKind k) {
    switch (k) {
        case LKind::g: return "g";
        case LKind::L: return "L";
        case LKind::S: return "S";
        case LKind::I: return "I";
        case LKind::E: return "E";  // ✨ NEW
    }
    return "?";
}

// ========== Basic operations ==========
int Topology_enhanced::addBlock(LKind kind, int param) {
    block.push_back(Block{kind, param});
    return static_cast<int>(block.size()) - 1;
}

int Topology_enhanced::addBlockRight(LKind kind, int param) {
    const int newId = addBlock(kind, param);
    if (newId > 0) {
        l_connection.push_back(InteriorStructure{newId - 1, newId});
    }
    return newId;
}

int Topology_enhanced::addSideLink(int param) {
    side_links.push_back(SideLinks{param});
    return static_cast<int>(side_links.size()) - 1;
}

int Topology_enhanced::addInstanton(int param) {
    instantons.push_back(Instantons{param});
    return static_cast<int>(instantons.size()) - 1;
}

// ✨ NEW: Add external curve
int Topology_enhanced::addExternal(int param) {
    externals.push_back(External{param});
    return static_cast<int>(externals.size()) - 1;
}

int Topology_enhanced::addDecoration(LKind kind, int param, int blockidx) {
    switch (kind) {
        case LKind::S: {
            const int sid = addSideLink(param);
            s_connection.push_back({blockidx, sid});
            return sid;
        }
        case LKind::I: {
            const int iid = addInstanton(param);
            i_connection.push_back({blockidx, iid});
            return iid;
        }
        case LKind::E: {  // ✨ NEW: External as decoration
            const int eid = addExternal(param);
            e_connection.push_back({blockidx, 0, 0, eid});  // Default to port 0
            return eid;
        }
        default:
            return -1;
    }
}

// ✨ NEW: External attachment methods
bool Topology_enhanced::attachExternalToBlock(int external_id, int block_id, int port_idx) {
    if (external_id < 0 || external_id >= static_cast<int>(externals.size())) return false;
    if (block_id < 0 || block_id >= static_cast<int>(block.size())) return false;
    
    e_connection.push_back({block_id, 0, port_idx, external_id});
    return true;
}

bool Topology_enhanced::attachExternalToSideLink(int external_id, int sidelink_id, int port_idx) {
    if (external_id < 0 || external_id >= static_cast<int>(externals.size())) return false;
    if (sidelink_id < 0 || sidelink_id >= static_cast<int>(side_links.size())) return false;
    
    e_connection.push_back({sidelink_id, 1, port_idx, external_id});
    return true;
}

bool Topology_enhanced::attachExternalToInstanton(int external_id, int instanton_id, int port_idx) {
    if (external_id < 0 || external_id >= static_cast<int>(externals.size())) return false;
    if (instanton_id < 0 || instanton_id >= static_cast<int>(instantons.size())) return false;
    
    e_connection.push_back({instanton_id, 2, port_idx, external_id});
    return true;
}

bool Topology_enhanced::attachExternal(int external_id, int parent_id, int parent_type, int port_idx) {
    if (external_id < 0 || external_id >= static_cast<int>(externals.size())) return false;
    
    switch (parent_type) {
        case 0: return attachExternalToBlock(external_id, parent_id, port_idx);
        case 1: return attachExternalToSideLink(external_id, parent_id, port_idx);
        case 2: return attachExternalToInstanton(external_id, parent_id, port_idx);
        default: return false;
    }
}

// ========== Connection management ==========
int Topology_enhanced::LinkingBlocks(std::vector<InteriorStructure> lcon) {
    l_connection = std::move(lcon);
    return static_cast<int>(l_connection.size());
}

int Topology_enhanced::LinkingSideLinks(std::vector<SideLinkStructure> scon) {
    s_connection = std::move(scon);
    return static_cast<int>(s_connection.size());
}

int Topology_enhanced::LinkingInstantonStructure(std::vector<InstantonStructure> icon) {
    i_connection = std::move(icon);
    return static_cast<int>(i_connection.size());
}

int Topology_enhanced::LinkingExternals(std::vector<ExternalStructure> econ) {
    e_connection = std::move(econ);
    return static_cast<int>(e_connection.size());
}

// ========== Initialization ==========
void Topology_enhanced::Initialize() {
    name.clear();
    block.clear();
    side_links.clear();
    instantons.clear();
    externals.clear();  // ✨ NEW
    l_connection.clear();
    s_connection.clear();
    i_connection.clear();
    e_connection.clear();  // ✨ NEW
}

void Topology_enhanced::InitializeBlocks() {
    block.clear();
}

void Topology_enhanced::InitializeDecorations() {
    side_links.clear();
    instantons.clear();
    externals.clear();  // ✨ NEW
}

void Topology_enhanced::InitializeGluings() {
    l_connection.clear();
    s_connection.clear();
    i_connection.clear();
    e_connection.clear();  // ✨ NEW
}

void Topology_enhanced::InitializeExternals() {
    externals.clear();
    e_connection.clear();
}

// ========== Query functions ==========
bool Topology_enhanced::hasNode(int id) const {
    return (id >= 0 && id < static_cast<int>(block.size()));
}

bool Topology_enhanced::hasExternal(int id) const {
    return (id >= 0 && id < static_cast<int>(externals.size()));
}

std::vector<int> Topology_enhanced::getExternalsOnBlock(int block_id) const {
    std::vector<int> result;
    for (const auto& ec : e_connection) {
        if (ec.parent_type == 0 && ec.parent_id == block_id) {
            result.push_back(ec.external_id);
        }
    }
    return result;
}

std::vector<int> Topology_enhanced::getExternalsOnPort(int parent_id, int parent_type, int port) const {
    std::vector<int> result;
    for (const auto& ec : e_connection) {
        if (ec.parent_type == parent_type && 
            ec.parent_id == parent_id && 
            ec.port_idx == port) {
            result.push_back(ec.external_id);
        }
    }
    return result;
}

int Topology_enhanced::getMaxPortIndex(int parent_id, int parent_type) const {
    int max_port = -1;
    for (const auto& ec : e_connection) {
        if (ec.parent_type == parent_type && ec.parent_id == parent_id) {
            max_port = std::max(max_port, ec.port_idx);
        }
    }
    return max_port;
}

int Topology_enhanced::countExternalsOnObject(int parent_id, int parent_type) const {
    int count = 0;
    for (const auto& ec : e_connection) {
        if (ec.parent_type == parent_type && ec.parent_id == parent_id) {
            ++count;
        }
    }
    return count;
}

// ========== Theory checks (stubs) ==========
bool Topology_enhanced::isSCFT() const { return false; }
bool Topology_enhanced::isLST() const { return false; }
bool Topology_enhanced::isSUGRA() const { return false; }

// ========== File I/O ==========
bool Topology_enhanced::saveToFile(const std::string& path) const {
    std::ofstream out(path);
    if (!out) return false;

    out << "name " << std::quoted(name) << "\n";

    // Blocks
    out << "blocks " << block.size() << "\n";
    for (size_t i = 0; i < block.size(); ++i) {
        out << "  [" << i << "] kind " << LKindToString(block[i].kind)
            << " param " << block[i].param << "\n";
    }

    // Side links
    out << "side_links " << side_links.size() << "\n";
    for (size_t i = 0; i < side_links.size(); ++i) {
        out << "  [" << i << "] param " << side_links[i].param << "\n";
    }

    // Instantons
    out << "instantons " << instantons.size() << "\n";
    for (size_t i = 0; i < instantons.size(); ++i) {
        out << "  [" << i << "] param " << instantons[i].param << "\n";
    }

    // ✨ NEW: Externals
    out << "externals " << externals.size() << "\n";
    for (size_t i = 0; i < externals.size(); ++i) {
        out << "  [" << i << "] param " << externals[i].param << "\n";
    }

    // Connections
    out << "l_connection " << l_connection.size() << "\n";
    for (const auto& e : l_connection) {
        out << "  (" << e.u << ", " << e.v << ")\n";
    }

    out << "s_connection " << s_connection.size() << "\n";
    for (const auto& e : s_connection) {
        out << "  (node " << e.u << ", sid " << e.v << ")\n";
    }

    out << "i_connection " << i_connection.size() << "\n";
    for (const auto& e : i_connection) {
        out << "  (node " << e.u << ", iid " << e.v << ")\n";
    }

    // ✨ NEW: External connections with port info
    out << "e_connection " << e_connection.size() << "\n";
    for (const auto& e : e_connection) {
        out << "  (parent " << e.parent_id 
            << ", type " << e.parent_type 
            << ", port " << e.port_idx 
            << ", eid " << e.external_id << ")\n";
    }

    return true;
}

// ========== Debug output ==========
std::ostream& operator<<(std::ostream& os, const Topology_enhanced& T) {
    os << "Topology_enhanced " << std::quoted(T.name) << "\n";

    os << "  Blocks (" << T.block.size() << ")\n";
    for (size_t i = 0; i < T.block.size(); ++i) {
        const auto& b = T.block[i];
        os << "    [" << i << "] kind=" << LKindToString(b.kind)
           << " param=" << b.param << "\n";
    }

    os << "  SideLinks (" << T.side_links.size() << ")\n";
    for (size_t i = 0; i < T.side_links.size(); ++i) {
        os << "    [" << i << "] param=" << T.side_links[i].param << "\n";
    }

    os << "  Instantons (" << T.instantons.size() << ")\n";
    for (size_t i = 0; i < T.instantons.size(); ++i) {
        os << "    [" << i << "] param=" << T.instantons[i].param << "\n";
    }

    // ✨ NEW: Externals
    os << "  Externals (" << T.externals.size() << ")\n";
    for (size_t i = 0; i < T.externals.size(); ++i) {
        os << "    [" << i << "] param=" << T.externals[i].param << "\n";
    }

    os << "  L-Connections (" << T.l_connection.size() << ")\n";
    for (const auto& e : T.l_connection) {
        os << "    " << e.u << " - " << e.v << "\n";
    }

    os << "  S-Connections (" << T.s_connection.size() << ")\n";
    for (const auto& e : T.s_connection) {
        os << "    node " << e.u << " ~ sid " << e.v << "\n";
    }

    os << "  I-Connections (" << T.i_connection.size() << ")\n";
    for (const auto& e : T.i_connection) {
        os << "    node " << e.u << " ~ iid " << e.v << "\n";
    }

    // ✨ NEW: External connections
    os << "  E-Connections (" << T.e_connection.size() << ")\n";
    for (const auto& e : T.e_connection) {
        os << "    External[" << e.external_id << "] on ";
        switch (e.parent_type) {
            case 0: os << "Block[" << e.parent_id << "]"; break;
            case 1: os << "SideLink[" << e.parent_id << "]"; break;
            case 2: os << "Instanton[" << e.parent_id << "]"; break;
        }
        os << " port " << e.port_idx << "\n";
    }

    return os;
}
