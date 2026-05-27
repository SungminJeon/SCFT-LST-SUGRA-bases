#pragma once
#include <string>
#include <vector>
#include <utility>
#include <cstdint>
#include <ostream>

// 블록 종류: g, L, S, I, E (External 추가)
enum class LKind : uint8_t { g, L, S, I, E };

// 블록: 종류 + 정수 파라미터
struct Block {
    LKind kind;
    int   param; // e.g., -2, -1, 0, ...
};

// Decoration
struct SideLinks { 
    int param;
};

struct Instantons {
    int param;
};

// ✨ NEW: External curves
struct External {
    int param;  // External curve parameter (typically 0 or specific value)
};

// link structure of nodes and interior links
struct InteriorStructure {
    int u; // node id
    int v; // node id
};

// link structure of nodes and side links
struct SideLinkStructure {
    int u; // node id
    int v; // sidelink id
};

// link structure of instantons and nodes
struct InstantonStructure {
    int u; // node id
    int v; // instanton id
};

// ✨ NEW: External curve connection with port support
struct ExternalStructure {
    int parent_id;      // ID of parent object (block, sidelink, or instanton)
    int parent_type;    // 0=block, 1=sidelink, 2=instanton
    int port_idx;       // Port index on the parent object
    int external_id;    // ID of the external curve
};

// ✨ Enhanced Topology with External support
class Topology_enhanced {
public:
    std::string name; // Topology identifier

    // Core components
    std::vector<Block>         block;        // Set of Blocks (nodes, interior links)
    std::vector<SideLinks>     side_links;   // Set of side links
    std::vector<Instantons>    instantons;   // Set of instantons
    std::vector<External>      externals;    // ✨ NEW: Set of external curves

    // Linking structures
    std::vector<InteriorStructure>  l_connection;  // Interior connections
    std::vector<SideLinkStructure>  s_connection;  // Sidelink connections
    std::vector<InstantonStructure> i_connection;  // Instanton connections
    std::vector<ExternalStructure>  e_connection;  // ✨ NEW: External connections with ports

    // === Utility functions ===
    
    // Basic block operations
    int addBlock(LKind kind, int param);
    int addBlockRight(LKind kind, int param);
    
    // Decoration operations
    int addDecoration(LKind kind, int param, int blockidx);
    int addSideLink(int param);
    int addInstanton(int param);
    int addExternal(int param = 0);  // ✨ NEW
    
    // ✨ NEW: External attachment with port specification
    // Attach external to a block at specific port
    bool attachExternalToBlock(int external_id, int block_id, int port_idx = 0);
    
    // Attach external to a sidelink at specific port
    bool attachExternalToSideLink(int external_id, int sidelink_id, int port_idx = 0);
    
    // Attach external to an instanton at specific port
    bool attachExternalToInstanton(int external_id, int instanton_id, int port_idx = 0);
    
    // ✨ NEW: Generic external attachment
    bool attachExternal(int external_id, int parent_id, int parent_type, int port_idx = 0);
    
    // Connection management
    int LinkingBlocks(std::vector<InteriorStructure> lcon);
    int LinkingSideLinks(std::vector<SideLinkStructure> scon);
    int LinkingInstantonStructure(std::vector<InstantonStructure> icon);
    int LinkingExternals(std::vector<ExternalStructure> econ);  // ✨ NEW

    // Initialization
    void Initialize();
    void InitializeBlocks();
    void InitializeDecorations();
    void InitializeGluings();
    void InitializeExternals();  // ✨ NEW
    
    // Query functions
    bool hasNode(int id) const;
    bool hasExternal(int id) const;  // ✨ NEW
    std::vector<int> getExternalsOnBlock(int block_id) const;  // ✨ NEW
    std::vector<int> getExternalsOnPort(int parent_id, int parent_type, int port) const;  // ✨ NEW
    
    // ✨ NEW: Port information
    int getMaxPortIndex(int parent_id, int parent_type) const;
    int countExternalsOnObject(int parent_id, int parent_type) const;
    
    // Theory checks (stubs for now)
    bool isSCFT() const;
    bool isLST() const;
    bool isSUGRA() const;
    bool hasExternalCurves() const { return !externals.empty(); }  // ✨ NEW
    
    // File I/O
    bool saveToFile(const std::string& path) const;
    bool loadFromFile(const std::string& path);
    
    // Conversion
    static Topology_enhanced fromBasicTopology(const struct Topology& basic);
};

// Debug output operator
std::ostream& operator<<(std::ostream& os, const Topology_enhanced& T);

// Helper function to convert LKind to string
const char* LKindToString(LKind k);
