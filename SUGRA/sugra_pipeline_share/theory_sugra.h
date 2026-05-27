// theory_sugra.h
// SUGRA construction support: catalog reader + IF reconstruction
//
// Reads unified catalog (NN + Nk), reconstructs intersection form from topo string.
//
// Usage:
//   #include "theory_sugra.h"
//   auto catalog = load_catalog("unified.cat");
//   for (auto& entry : catalog) {
//       Eigen::MatrixXi IF = reconstruct_IF(entry);
//       // ... attach external curves, check anomaly, etc.
//   }
//
// Compile with: Tensor.C, TopoLineCompact_enhanced.cpp, Topology_enhanced.cpp

#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <functional>

#include "Tensor.h"
#include "Theory_enhanced.h"
#include "TopoLineCompact_enhanced.hpp"
#include "Topology_enhanced.h"

// ============================================================================
// Catalog Entry
// ============================================================================

struct CatalogEntry {
    int id;
    std::string type;   // "NN" or "Nk" (e.g. "N1", "N8", ...)
    int T;              // tensor multiplets = sig_neg of IF
    std::string topo;   // reconstruction string
    
    bool is_nn() const { return type == "NN"; }
    int node_count() const {
        if (is_nn()) return 0;
        return std::stoi(type.substr(1));
    }
};

// ============================================================================
// Catalog I/O
// ============================================================================

inline std::vector<CatalogEntry> load_catalog(const std::string& filename) {
    std::vector<CatalogEntry> catalog;
    std::ifstream fin(filename);
    if (!fin) {
        std::cerr << "Error: Cannot open catalog " << filename << "\n";
        return catalog;
    }
    
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        // Parse: id | type | T | topo...
        // First 3 fields pipe-delimited, topo = everything after 3rd '|'
        size_t pos = 0;
        int pipe_count = 0;
        for (size_t i = 0; i < line.size(); i++) {
            if (line[i] == '|') {
                pipe_count++;
                if (pipe_count == 3) { pos = i + 2; break; }
            }
        }
        if (pipe_count < 3) continue;
        
        std::istringstream iss(line);
        CatalogEntry entry;
        std::string sep;
        iss >> entry.id >> sep >> entry.type >> sep >> entry.T >> sep;
        entry.topo = line.substr(pos);
        
        // Trim trailing whitespace
        while (!entry.topo.empty() && (entry.topo.back() == ' ' || entry.topo.back() == '\n' || entry.topo.back() == '\r'))
            entry.topo.pop_back();
        
        catalog.push_back(entry);
    }
    
    return catalog;
}

// ============================================================================
// NN Reconstruction: Noble Molecule Builders
// ============================================================================

namespace nn_detail {

struct NobleMolecule {
    std::string name;
    std::function<Tensor()> builder;
};

inline const std::vector<NobleMolecule>& get_noble_molecules() {
    static std::vector<NobleMolecule> molecules = {
        // noble 2-Molecules with no -5 curves
        {"1{1,3}22", []{ Tensor t; t.AT(-1); t.ATS(-1,-3); t.AT(-2); t.AT(-2); return t; }},
        
        // noble 0-Molecules with no -5 curves
        {"2{2,3}13", []{ Tensor t; t.AT(-2); t.ATS(-2,-3); t.AT(-1); t.AT(-3); return t; }},
        {"2{2,3}132", []{ Tensor t; t.AT(-2); t.ATS(-2,-3); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2{2,3}1322", []{ Tensor t; t.AT(-2); t.ATS(-2,-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"23132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"223132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"3132", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2132", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"3123", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-2); t.AT(-3); return t; }},
        {"2231322", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"31322", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"21322", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"313", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); return t; }},
        {"23213", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-3); return t; }},
        {"213", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-3); return t; }},

        // noble 4-Molecules with one -5 curve
        {"1{1,{1,5}}1", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-1); t.not_intersect(4,5); t.intersect(5,3); return t; }},
        {"1{1{1,{1,5}}}13", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-1); t.not_intersect(4,5); t.intersect(5,3); t.AT(-3); t.not_intersect(6,5); t.intersect(6,3); return t; }},
        {"1{1{1,{1,5}}}132", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-1); t.not_intersect(4,5); t.intersect(5,3); t.AT(-3); t.not_intersect(6,5); t.intersect(6,3); t.AT(-2); return t; }},
        {"1{1,5}1", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); return t; }},

        // noble 2-Molecules with one -5 curve
        {"151{1,3}2", []{ Tensor t; t.AT(-1); t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
        {"1{1,5}12", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-2); return t; }},
        {"31{1,5}1", []{ Tensor t; t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); return t; }},
        {"231{1,5}1", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); return t; }},
        {"2231{1,5}1", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); return t; }},
        {"2231{1,5}131", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-1); return t; }},

        // noble 1-Molecules with one -5 curve
        {"51{1,3}2", []{ Tensor t; t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
        {"3151{1,3}2", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
        {"23151{1,3}2", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
        {"223151{1,3}2", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); return t; }},
        {"51{1,3}22", []{ Tensor t; t.AT(-5); t.AT(-1); t.ATS(-1,-3); t.AT(-2); t.AT(-2); return t; }},
        {"512{1,3}2", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-2); t.ATS(-1,-3); t.AT(-2); return t; }},
        {"31{1,5}12", []{ Tensor t; t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-2); return t; }},
        {"231{1,5}12", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-2); return t; }},
        {"2231{1,5}12", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-2); return t; }},
        {"31{1,5}13", []{ Tensor t; t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); return t; }},
        {"231{1,5}13", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); return t; }},
        {"2231{1,5}13", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); return t; }},
        {"231{1,5}132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2231{1,5}132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2231{1,5}1322", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"13215132", []{ Tensor t; t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"223151231", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); return t; }},
        {"3151231", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); return t; }},

        // noble 0-Molecules with one -5 curve
        {"31{3{1,5}13}", []{ Tensor t; t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); return t; }},
        {"231{3{1,5}13}", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); return t; }},
        {"2231{3{1,5}13}", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); return t; }},
        {"231{3{1,5}132}", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2231{3{1,5}132}", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS2(-3,-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"3215", []{ Tensor t; t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"2315", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"32215", []{ Tensor t; t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"22315", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"315", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"23215", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"215", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"2215", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"22215", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"313215", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"231315", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"31315", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},

        // noble molecules with -5 curves (longer chains)
        {"3215132", []{ Tensor t; t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2315132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"22315132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"315132", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"23215132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"215132", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2215132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"22315123", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); return t; }},
        {"315123", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); return t; }},
        {"215123", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); return t; }},
        {"223151322", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"3151322", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"232151322", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"2151322", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"22151322", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"31513", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
        {"2321513", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
        {"21513", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
        {"221513", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
        {"2151232", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); return t; }},
        {"21512", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); return t; }},
        {"31315132", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"313151322", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"3131513", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},

        // noble 1-Molecules with two -5 curves
        {"2231513151321", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); return t; }},
        {"13151315132", []{ Tensor t; t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"131513151322", []{ Tensor t; t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},

        // noble 0-Molecules with two -5 curves
        {"513215", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"5132215", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"51315", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"5123215", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"231513215", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"2231513215", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"31513215", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"21513215", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"31512315", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"2315132215", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"22315132215", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"315132215", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"32151315", []{ Tensor t; t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"23151315", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"223151315", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"3151315", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"2151315", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"23151315132", []{ Tensor t; t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"223151315132", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"3151315132", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"2151315132", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); return t; }},
        {"22315131513222", []{ Tensor t; t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"31513151322", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"21513151322", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); return t; }},
        {"315131513", []{ Tensor t; t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},
        {"215131513", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); return t; }},

        // noble 0-Molecules with three -5 curves
        {"5131513215", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); return t; }},
        {"513151315", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},

        // additional noble molecules
        {"1{1,5}13151", []{ Tensor t; t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
        {"15132151", []{ Tensor t; t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
        {"21513151", []{ Tensor t; t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
        {"1513151315", []{ Tensor t; t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); return t; }},
        {"51232151", []{ Tensor t; t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); return t; }},
    };
    return molecules;
}

// Parse NN topo: "SL:param:port:extra:intnum" or "NM:name:port:extra:intnum"
struct NNTopoInfo {
    bool is_sidelink;
    int sidelink_param;
    std::string noble_name;
    int port;
    int extra_curve;
    int int_num;
    bool valid;
};

inline NNTopoInfo parse_nn_topo(const std::string& topo) {
    NNTopoInfo info;
    info.valid = false;
    
    std::vector<std::string> parts;
    std::istringstream iss(topo);
    std::string token;
    while (std::getline(iss, token, ':')) parts.push_back(token);
    
    if (parts.size() != 5) return info;
    
    if (parts[0] == "SL") {
        info.is_sidelink = true;
        info.sidelink_param = std::stoi(parts[1]);
    } else if (parts[0] == "NM") {
        info.is_sidelink = false;
        info.noble_name = parts[1];
    } else {
        return info;
    }
    
    info.port = std::stoi(parts[2]);
    info.extra_curve = std::stoi(parts[3]);
    info.int_num = std::stoi(parts[4]);
    info.valid = true;
    return info;
}

} // namespace nn_detail

// ============================================================================
// Nk Reconstruction: build IF from topocompactline via TheoryGraph
// ============================================================================

namespace nk_detail {

inline Spec make_spec(LKind kind, int param) {
    switch (kind) {
        case LKind::g: return n(param);
        case LKind::L: return i(param);
        case LKind::S: return s(param);
        case LKind::I: return s(param);
        case LKind::E: return e(param);
        default: return n(param);
    }
}

inline Eigen::MatrixXi build_IF_from_topology(const Topology_enhanced& T) {
    TheoryGraph G;
    std::vector<NodeRef> blockNodes, sideNodes, instNodes, extNodes;

    for (size_t i = 0; i < T.block.size(); ++i)
        blockNodes.push_back(G.add(make_spec(T.block[i].kind, T.block[i].param)));
    for (size_t i = 0; i < T.side_links.size(); ++i)
        sideNodes.push_back(G.add(s(T.side_links[i].param)));
    for (size_t i = 0; i < T.instantons.size(); ++i)
        instNodes.push_back(G.add(s(T.instantons[i].param)));
    for (size_t i = 0; i < T.externals.size(); ++i)
        extNodes.push_back(G.add(e(T.externals[i].param)));

    for (const auto& conn : T.l_connection)
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)blockNodes.size())
            G.connect(blockNodes[conn.u], blockNodes[conn.v]);
    for (const auto& conn : T.s_connection)
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)sideNodes.size())
            G.connect(sideNodes[conn.v], blockNodes[conn.u]);
    for (const auto& conn : T.i_connection)
        if (conn.u >= 0 && conn.u < (int)blockNodes.size() &&
            conn.v >= 0 && conn.v < (int)instNodes.size())
            G.connect(instNodes[conn.v], blockNodes[conn.u]);
    for (const auto& conn : T.e_connection) {
        if (conn.external_id < 0 || conn.external_id >= (int)extNodes.size()) continue;
        NodeRef* p = nullptr;
        if (conn.parent_type == 0 && conn.parent_id >= 0 && conn.parent_id < (int)blockNodes.size())
            p = &blockNodes[conn.parent_id];
        else if (conn.parent_type == 1 && conn.parent_id >= 0 && conn.parent_id < (int)sideNodes.size())
            p = &sideNodes[conn.parent_id];
        else if (conn.parent_type == 2 && conn.parent_id >= 0 && conn.parent_id < (int)instNodes.size())
            p = &instNodes[conn.parent_id];
        if (p) G.connect(extNodes[conn.external_id], AttachmentPoint(-1), *p, AttachmentPoint(conn.port_idx));
    }
    return G.ComposeIF_Gluing();
}

} // namespace nk_detail

// ============================================================================
// Main API: reconstruct IF from any catalog entry
// ============================================================================

inline Eigen::MatrixXi reconstruct_IF(const CatalogEntry& entry) {
    if (entry.is_nn()) {
        // NN path: parse topo → build base Tensor → ATE → IF
        auto info = nn_detail::parse_nn_topo(entry.topo);
        if (!info.valid) return Eigen::MatrixXi();
        
        Tensor base;
        if (info.is_sidelink) {
            base = build_tensor(s(info.sidelink_param));
        } else {
            bool found = false;
            for (const auto& nm : nn_detail::get_noble_molecules()) {
                if (nm.name == info.noble_name) {
                    base = nm.builder();
                    found = true;
                    break;
                }
            }
            if (!found) return Eigen::MatrixXi();
        }
        
        if (base.GetT() == 0) return Eigen::MatrixXi();
        
        Tensor t_ext;
        t_ext.SetIF(base.GetIntersectionForm());
        t_ext.ATE(info.extra_curve, info.port, info.int_num);
        return t_ext.GetIntersectionForm();
        
    } else {
        // Nk path: deserialize topocompactline → TheoryGraph → IF
        Topology_enhanced topo;
        if (!TopoLineCompact_enhanced::deserialize(entry.topo, topo))
            return Eigen::MatrixXi();
        return nk_detail::build_IF_from_topology(topo);
    }
}

// Convenience: get Tensor object (for further operations like blowdown)
inline Tensor reconstruct_tensor(const CatalogEntry& entry) {
    Tensor result;
    Eigen::MatrixXi IF = reconstruct_IF(entry);
    if (IF.rows() > 0) result.SetIF(IF);
    return result;
}

// ============================================================================
// Query helpers
// ============================================================================

inline std::vector<CatalogEntry> filter_by_type(const std::vector<CatalogEntry>& catalog,
                                                 const std::string& type) {
    std::vector<CatalogEntry> result;
    for (const auto& e : catalog)
        if (e.type == type) result.push_back(e);
    return result;
}

inline std::vector<CatalogEntry> filter_by_T(const std::vector<CatalogEntry>& catalog,
                                              int T_min, int T_max) {
    std::vector<CatalogEntry> result;
    for (const auto& e : catalog)
        if (e.T >= T_min && e.T <= T_max) result.push_back(e);
    return result;
}

inline std::vector<CatalogEntry> filter_nn(const std::vector<CatalogEntry>& catalog) {
    return filter_by_type(catalog, "NN");
}

inline std::vector<CatalogEntry> filter_nk(const std::vector<CatalogEntry>& catalog) {
    std::vector<CatalogEntry> result;
    for (const auto& e : catalog)
        if (!e.is_nn()) result.push_back(e);
    return result;
}

inline void print_catalog_summary(const std::vector<CatalogEntry>& catalog) {
    std::map<std::string, int> type_counts;
    for (const auto& e : catalog) type_counts[e.type]++;
    
    std::cout << "Catalog: " << catalog.size() << " entries\n";
    for (const auto& [type, count] : type_counts)
        std::cout << "  " << type << ": " << count << "\n";
}

// end of theory_sugra.h
