// cat2tex.cpp — standalone catalog → LaTeX converter
// Reads any cat_Next/ directory and generates results.tex
//
// Usage: ./cat2tex <cat_dir> [max_entries=500]
//
// Works with Phase 1, Phase 2, or any iteration output.

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <map>
#include <algorithm>
#include <filesystem>

// ── Spec labels (for display) ──
struct SpecLabel { std::string tag, label; };
inline std::vector<SpecLabel> spec_labels() {
    std::vector<SpecLabel> v;
    auto a=[&](const std::string&t,const std::string&l){v.push_back({t,l});};
    a("su2","$\\mathfrak{su}_2$");
    a("su2n3","$\\mathfrak{su}_2{\\to}(-3)$"); a("su8","$\\mathfrak{su}_8{\\to}\\mathfrak{so}_{16}$");
    a("su3","$\\mathfrak{su}_3$");
    a("su3n2","$(-3){\\to}(-2)$");
    a("so8","$\\mathfrak{so}_8$"); a("so16n2","$\\mathfrak{so}_{16}{\\to}\\mathfrak{su}_8$");
    a("f4","$\\mathfrak{f}_4$"); a("e6","$\\mathfrak{e}_6$");
    a("e7p","$\\mathfrak{e}_7'$"); a("e7","$\\mathfrak{e}_7$");
    a("e8","$\\mathfrak{e}_8$");
    a("hat1m1","$\\hat{1}{\\to}(-1)$"); a("hat1m2","$\\hat{1}{\\to}(-2)$");
    a("2","$(-2)$");                                   // no gauge: bare (-2) curve
    a("2mix","$(-2)$");                                // same — bare (-2) with (-1) connections
    return v;
}

// ── Entry parsed from .cat ──
struct Entry {
    int entry_id = 0;
    std::string combo_key;
    struct Ext { int curve_idx=-1, spec_id=-1; std::string tag; int ext_si=0, target_si=0, int_num=1; bool is_hat1=false; };
    std::vector<Ext> externals;
    int catalog_id=0; std::string catalog_type; int base_T=0;
    int T=0, H_charged=0, V=0, H_neutral=0, det=0, sig_pos=0, sig_neg=0, sig_zero=0;
    Eigen::MatrixXi IF;
};

std::vector<Entry> load_cat_file(const std::string& fn) {
    std::vector<Entry> es; std::ifstream f(fn); if(!f) return es;
    std::string line; Entry c; bool in=false;
    int ifl=0,ifr=0,ifn=0,el=0;
    while(std::getline(f,line)){
        if(line.empty()||line[0]=='#') continue;
        if(line.substr(0,6)=="ENTRY "){c=Entry();c.entry_id=std::stoi(line.substr(6));in=true;}
        else if(in&&line.substr(0,5)=="SPEC "){
            std::istringstream ss(line.substr(5));std::string tag;int sid,es2,ts,inum,h1;
            ss>>sid>>tag>>es2>>ts>>inum>>h1;
            if(c.externals.empty()&&c.combo_key.empty()){c.combo_key=tag;
                c.externals.push_back({-1,sid,tag,es2,ts,inum,h1!=0});}}
        else if(in&&line.substr(0,6)=="COMBO "){c.combo_key=line.substr(6);
            while(!c.combo_key.empty()&&c.combo_key.back()<=' ')c.combo_key.pop_back();}
        else if(in&&line.substr(0,10)=="EXTERNALS "){el=std::stoi(line.substr(10));c.externals.clear();}
        else if(in&&el>0){Entry::Ext ei;std::istringstream ss(line);int h1;
            ss>>ei.curve_idx>>ei.spec_id>>ei.tag>>ei.ext_si>>ei.target_si>>ei.int_num>>h1;
            ei.is_hat1=(h1!=0);c.externals.push_back(ei);el--;}
        else if(in&&line.substr(0,5)=="BASE "){std::istringstream ss(line.substr(5));int nb;ss>>c.catalog_id>>c.catalog_type>>c.base_T>>nb;}
        else if(in&&line.substr(0,7)=="ATTACH "){std::istringstream ss(line.substr(7));int tidx,eidx;ss>>tidx>>eidx;
            if(!c.externals.empty()&&c.externals.back().curve_idx<0)c.externals.back().curve_idx=eidx;}
        else if(in&&line.substr(0,8)=="PHYSICS "){std::istringstream ss(line.substr(8));
            ss>>c.T>>c.H_charged>>c.V>>c.H_neutral>>c.det>>c.sig_pos>>c.sig_neg>>c.sig_zero;}
        else if(in&&line.substr(0,3)=="IF "){ifn=std::stoi(line.substr(3));c.IF=Eigen::MatrixXi::Zero(ifn,ifn);ifl=ifn;ifr=0;}
        else if(in&&ifl>0){std::istringstream ss(line);for(int j=0;j<ifn;j++)ss>>c.IF(ifr,j);ifr++;ifl--;}
        else if(in&&line=="END"){es.push_back(std::move(c));in=false;}
    }
    return es;
}

std::vector<Entry> load_cat_dir(const std::string& dir) {
    std::vector<Entry> all;
    for(auto& p:std::filesystem::directory_iterator(dir)){
        if(p.path().extension()!=".cat")continue;
        auto e=load_cat_file(p.path().string());
        for(auto&x:e)all.push_back(std::move(x));
    }
    return all;
}

// ── Quiver renderer (v5/v6 style) ──
// Base quiver: target curves colored, external curve excluded from chain
// External shown as separate annotation: + spec_label
static const char* COLS[]={"red","blue","green!60!black","orange","purple","cyan","magenta","brown!60!red","olive","pink"};

std::string quiver(const Eigen::MatrixXi& IF, const std::vector<Entry::Ext>& exts) {
    int n=IF.rows(); if(n==0)return"\\emptyset";

    // Find external curve indices and their targets
    std::set<int> ext_curves;
    std::map<int,std::vector<int>> target_colors;  // target_curve_idx → list of color indices
    for(int i=0;i<(int)exts.size();i++){
        ext_curves.insert(exts[i].curve_idx);
    }
    for(int i=0;i<(int)exts.size();i++){
        int ei=exts[i].curve_idx;
        for(int j=0;j<n;j++){
            if(j==ei) continue;
            if(IF(ei,j)!=0 && !ext_curves.count(j))
                target_colors[j].push_back(i);
        }
    }

    // Build adjacency excluding external curves
    int n_base=0;
    std::map<int,int> remap;  // original idx → new idx
    std::vector<int> unmap;   // new idx → original idx
    for(int i=0;i<n;i++){
        if(ext_curves.count(i)) continue;
        remap[i]=n_base;
        unmap.push_back(i);
        n_base++;
    }
    if(n_base==0) return"\\emptyset";

    // Build base-only adjacency
    std::vector<std::vector<std::pair<int,int>>>adj(n_base);
    for(int i=0;i<n_base;i++){
        int oi=unmap[i];
        for(int j=i+1;j<n_base;j++){
            int oj=unmap[j];
            if(IF(oi,oj)!=0){adj[i].push_back({j,IF(oi,oj)});adj[j].push_back({i,IF(oj,oi)});}
        }
    }

    // Start from leaf that is NOT a target
    int start=-1;
    for(int i=0;i<n_base;i++)
        if((int)adj[i].size()==1&&!target_colors.count(unmap[i])){start=i;break;}
    if(start<0)for(int i=0;i<n_base;i++)if((int)adj[i].size()==1){start=i;break;}
    if(start<0)start=0;

    std::vector<bool>vis(n_base,false);

    auto lb=[&](int c)->std::string{
        int oc=unmap[c];
        int si=-IF(oc,oc);
        std::string s=std::to_string(si);if(si>=10)s="("+s+")";
        auto it=target_colors.find(oc);
        if(it!=target_colors.end()&&!it->second.empty()){
            // Superscript: colored intersection numbers for each attached external
            std::string sup;
            for(int ei : it->second){
                int ci=ei % 10;
                int ext_curve=exts[ei].curve_idx;
                int int_num=std::abs(IF(oc,ext_curve));
                sup+="\\textcolor{"+std::string(COLS[ci])+"}{"+std::to_string(int_num)+"}";
            }
            return s+"^{"+sup+"}";
        }
        return s;
    };

    std::set<int> target_set;
    for(auto&[k,v]:target_colors) if(remap.count(k)) target_set.insert(remap[k]);

    std::function<std::string(int)>wb;
    std::function<std::string(int,int)>w=[&](int c,int f)->std::string{
        vis[c]=true;std::string l=lb(c);
        std::vector<int>ch;for(auto&[nb,ww]:adj[c])if(!vis[nb])ch.push_back(nb);
        if(ch.empty())return l;if(ch.size()==1)return l+w(ch[0],c);
        int mc=-1;std::vector<int>bc;
        for(int x:ch){if(mc<0&&!target_set.count(x))mc=x;else bc.push_back(x);}
        if(mc<0){mc=ch[0];bc.clear();for(size_t i=1;i<ch.size();i++)bc.push_back(ch[i]);}
        std::vector<std::string>bs;for(int i=(int)bc.size()-1;i>=0;i--)bs.push_back(wb(bc[i]));
        std::string mp=w(mc,c);if(bs.empty())return l+mp;
        if(bs.size()==1)return"\\overset{"+bs[0]+"}{"+l+"}"+mp;
        std::string u;for(size_t i=1;i<bs.size();i++){if(!u.empty())u+=",";u+=bs[i];}
        return"\\overset{"+bs[0]+"}{\\underset{"+u+"}{"+l+"}}"+mp;};
    wb=[&](int c)->std::string{
        vis[c]=true;std::string l=lb(c);
        std::vector<int>ch;for(auto&[nb,ww]:adj[c])if(!vis[nb])ch.push_back(nb);
        if(ch.empty())return l;if(ch.size()==1)return wb(ch[0])+l;
        int mc=-1;std::vector<int>sb;
        for(int x:ch){if(mc<0&&!target_set.count(x))mc=x;else sb.push_back(x);}
        if(mc<0){mc=ch[0];sb.clear();for(size_t i=1;i<ch.size();i++)sb.push_back(ch[i]);}
        std::vector<std::string>ss2;for(int i=(int)sb.size()-1;i>=0;i--)ss2.push_back(wb(sb[i]));
        std::string mp=wb(mc);if(ss2.size()==1)return mp+"\\overset{"+ss2[0]+"}{"+l+"}";
        std::string u;for(size_t i=1;i<ss2.size();i++){if(!u.empty())u+=",";u+=ss2[i];}
        return mp+"\\overset{"+ss2[0]+"}{\\underset{"+u+"}{"+l+"}}";};
    return w(start,-1);
}

// ── NHC cluster: per-curve gauge labels ──
static const std::map<std::string, std::vector<std::string>> nhc_curve_gauges = {
    {"nhc_2_3", {"\\mathfrak{su}_2", "\\mathfrak{g}_2"}},
    {"nhc_2_3_2", {"\\mathfrak{su}_2", "\\mathfrak{so}_7", "\\mathfrak{su}_2"}},
    {"nhc_2_2_3", {"\\mathfrak{g}_2", "\\mathfrak{sp}_1", "(-2)"}},
    {"nhc_2_4", {"\\mathfrak{su}_8", "\\mathfrak{so}_{16}"}},
};

// ── External annotation (shown after quiver) ──
std::string ext_annotation(const std::vector<Entry::Ext>& exts,
                            const std::map<std::string,std::string>& tag2label) {
    std::string s;
    std::set<int> nhc_shown;  // track which NHC cluster externals already displayed

    for(int i=0;i<(int)exts.size();i++){
        auto& e=exts[i];

        // NHC cluster external: show once with per-curve colors
        auto nhc_it = nhc_curve_gauges.find(e.tag);
        if(nhc_it != nhc_curve_gauges.end()){
            if(nhc_shown.count(e.spec_id)) continue;
            nhc_shown.insert(e.spec_id);
            // Collect indices of all curves in this NHC cluster
            std::vector<int> nhc_indices;
            for(int j=i;j<(int)exts.size();j++)
                if(exts[j].tag==e.tag && exts[j].spec_id==e.spec_id)
                    nhc_indices.push_back(j);
            // Build colored cluster label
            s+=" + [";
            auto& gauges = nhc_it->second;
            for(int k=0;k<(int)nhc_indices.size()&&k<(int)gauges.size();k++){
                if(k) s+="\\cdot ";
                int ci=nhc_indices[k] % 10;
                s+="\\textcolor{"+std::string(COLS[ci])+"}{"+gauges[k]+"}";
            }
            s+="]";
            continue;
        }

        std::string lbl;
        if(e.is_hat1){
            if(e.target_si==-1) lbl="\\hat{1}{\\to}1";
            else if(e.target_si==-2) lbl="\\hat{1}{\\to}2";
            else lbl="\\hat{1}";
        } else {
            auto it=tag2label.find(e.tag);
            if(it!=tag2label.end()){
                lbl=it->second;
                if(lbl.size()>=2 && lbl.front()=='$' && lbl.back()=='$')
                    lbl=lbl.substr(1,lbl.size()-2);
            }
            else {
                switch(e.ext_si) {
                    case -2: lbl="\\mathfrak{su}_2"; break;
                    case -3: lbl="\\mathfrak{su}_3"; break;
                    case -4: lbl="\\mathfrak{so}_8"; break;
                    case -5: lbl="\\mathfrak{f}_4"; break;
                    case -6: lbl="\\mathfrak{e}_6"; break;
                    case -7: lbl="\\mathfrak{e}_7'"; break;
                    case -8: lbl="\\mathfrak{e}_7"; break;
                    case -12: lbl="\\mathfrak{e}_8"; break;
                    default: lbl="\\mathfrak{"+e.tag+"}"; break;
                }
            }
        }
        // k is shown in quiver superscripts, not in annotation
        int ci=i % 10;
        s+=" + \\textcolor{"+std::string(COLS[ci])+"}{"+lbl+"}";
    }
    return s;
}

// ── Blowdown ──
Eigen::MatrixXi bd_curve(const Eigen::MatrixXi& IF,int idx){
    int n=IF.rows();Eigen::MatrixXi r=Eigen::MatrixXi::Zero(n-1,n-1);
    for(int i=0,ni=0;i<n;i++){if(i==idx)continue;for(int j=0,nj=0;j<n;j++){if(j==idx)continue;r(ni,nj)=IF(i,j)+IF(i,idx)*IF(idx,j);nj++;}ni++;}return r;}

struct BDResult { Eigen::MatrixXi M; std::string surface; };
BDResult blowdown(const Eigen::MatrixXi& IF, const std::set<int>& h1) {
    int n=IF.rows(),T=0;
    {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>es(IF.cast<double>());for(int i=0;i<n;i++)if(es.eigenvalues()(i)<-1e-8)T++;}
    Eigen::MatrixXi ext=Eigen::MatrixXi::Zero(n+1,n+1);ext.block(0,0,n,n)=IF;
    for(int i=0;i<n;i++){int b=h1.count(i)?-1:(IF(i,i)+2);ext(i,n)=b;ext(n,i)=b;}ext(n,n)=9-T;
    Eigen::MatrixXi cur=ext;std::set<int>ch=h1;int bd=0;
    while(true){int nn=cur.rows(),m=-1;for(int i=0;i<nn-1;i++){if(ch.count(i))continue;if(cur(i,i)==-1){m=i;break;}}
        if(m<0)break;cur=bd_curve(cur,m);std::set<int>nh;for(int h:ch)nh.insert(h>m?h-1:h);ch=nh;bd++;}
    BDResult r;r.M=cur;int nc=cur.rows()-1;
    if(nc==0)r.surface="\\text{pt}";else if(nc==1&&cur(0,0)==1)r.surface="\\mathbb{P}^2";
    else if(nc==2&&cur(0,1)==1){int a=cur(0,0),b=cur(1,1);
        if(a==0||b==0){int fn=(a==0)?std::abs(b):std::abs(a);r.surface="\\mathbb{F}_{"+std::to_string(fn)+"}";}
        else r.surface="\\text{bd}_{"+std::to_string(bd)+"}";}
    else r.surface="\\text{bd}_{"+std::to_string(bd)+"}";return r;}

// ── T_crit ──
struct TCrit { double val; int ival; bool exact; Eigen::MatrixXi bd_M; };

// Helper: build blowdown matrix at given T (extends IF with b0² = 9-T, then
// successively blows down any -1 curves not in the hat1 set).
static inline Eigen::MatrixXi blowdown_at_T(const Eigen::MatrixXi& IF,
                                             const std::set<int>& h1, int T) {
    int n = IF.rows();
    int b0sq = 9 - T;
    Eigen::MatrixXi ext = Eigen::MatrixXi::Zero(n+1, n+1);
    ext.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        int b = h1.count(i) ? -1 : (IF(i,i) + 2);
        ext(i, n) = b; ext(n, i) = b;
    }
    ext(n, n) = b0sq;
    Eigen::MatrixXi cur = ext;
    std::set<int> ch = h1;
    while (true) {
        int nn = cur.rows(), m = -1;
        for (int i = 0; i < nn - 1; i++) {
            if (ch.count(i)) continue;
            if (cur(i, i) == -1) { m = i; break; }
        }
        if (m < 0) break;
        cur = bd_curve(cur, m);
        std::set<int> nh;
        for (int h : ch) nh.insert(h > m ? h - 1 : h);
        ch = nh;
    }
    return cur;
}

// compute_tcrit: try the closed-form det=0 crossing first. If the extended-IF
// determinant is T-independent (A≈0), fall back to scanning T values in
// [sig_neg, T_max] for sig_pos=1 (mirrors the reject_nonuni fallback). The
// scan range T_max is anomaly-feasible by construction (max_t_add_from_fmulti
// + T_sugra), so any T found is both signature- and anomaly-valid.
TCrit compute_tcrit(const Eigen::MatrixXi& IF, const std::set<int>& h1,
                    int sig_neg = -1, int T_max_anom = -1) {
    TCrit res; int n = IF.rows();
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n+1, n+1);
    M.block(0, 0, n, n) = IF.cast<double>();
    for (int i = 0; i < n; i++) {
        double b = h1.count(i) ? -1.0 : (IF(i,i) + 2.0);
        M(i, n) = b; M(n, i) = b;
    }
    M(n, n) = 0; double B = M.determinant();
    M(n, n) = 1; double A = M.determinant() - B;

    if (std::abs(A) < 1e-10) {
        // A=0: extended det is T-constant. Fallback: scan T for sig_pos=1.
        if (sig_neg >= 0 && T_max_anom >= 0) {
            for (int T = sig_neg; T <= T_max_anom; T++) {
                if (check_extended_sig(IF, T, h1)) {
                    res.ival = T; res.val = (double)T; res.exact = true;
                    res.bd_M = blowdown_at_T(IF, h1, T);
                    return res;
                }
            }
        }
        res.val = 1e9; res.exact = false; res.ival = -1; return res;
    }

    res.val = 9.0 - (-B / A);
    int tc = (int)std::round(res.val);
    res.exact = (std::abs(res.val - tc) < 1e-6);
    res.ival = res.exact ? tc : (int)std::ceil(res.val - 1e-9);
    res.bd_M = blowdown_at_T(IF, h1, res.ival);
    return res;
}

// ── Main ──
int main(int argc,char*argv[]){
    if(argc<2){std::cerr<<"Usage: "<<argv[0]<<" <cat_dir> [max_entries=500] [T_max=999]\n";return 1;}
    std::string dir=argv[1]; int max_show=(argc>2)?std::atoi(argv[2]):500;
    int T_filter=(argc>3)?std::atoi(argv[3]):999;
    auto labels=spec_labels();
    std::map<std::string,std::string>tag2label;for(auto&l:labels)tag2label[l.tag]=l.label;

    auto entries=load_cat_dir(dir);
    // T filter
    if(T_filter<999){
        std::vector<Entry> filtered;
        for(auto&e:entries) if(e.T<=T_filter) filtered.push_back(std::move(e));
        entries=std::move(filtered);
    }
    std::cout<<"Loaded "<<entries.size()<<" entries from "<<dir<<" (T<="<<T_filter<<")\n";

    // Sort by combo, T, det
    std::sort(entries.begin(),entries.end(),[](const Entry&a,const Entry&b){
        if(a.combo_key!=b.combo_key)return a.combo_key<b.combo_key;
        if(a.T!=b.T)return a.T<b.T;return std::abs(a.det)<std::abs(b.det);});

    int n_ext=entries.empty()?0:(int)entries[0].externals.size();
    int n_uni=0,n_nonuni=0;
    for(auto&e:entries){if(std::abs(e.det)==1)n_uni++;else n_nonuni++;}

    auto wsm=[](std::ofstream&t,const Eigen::MatrixXi&M){int s=M.rows();t<<"\\left(\\begin{smallmatrix}";
        for(int i=0;i<s;i++){if(i)t<<" \\\\ ";for(int j=0;j<s;j++){if(j)t<<" & ";t<<M(i,j);}}
        t<<"\\end{smallmatrix}\\right)";};

    std::string texfile=dir+"/results.tex";
    std::ofstream tex(texfile);
    tex<<"\\documentclass[10pt]{article}\n\\usepackage[a4paper,margin=0.6in]{geometry}\n"
       <<"\\usepackage{longtable,booktabs,amsmath,amssymb,array}\n\\usepackage[dvipsnames]{xcolor}\n"
       <<"\\setlength{\\parindent}{0pt}\n\\setlength{\\parskip}{2pt}\n\\begin{document}\n\n";
    // Section title: detect NHC external directory
    if(dir.find("nhc") != std::string::npos)
        tex<<"\\section*{NHC w/ more than 2 tensors}\n\n";
    else
        tex<<"\\section*{"<<n_ext<<"-External SUGRA Bases}\n\n";
    tex<<"Colored curves in quiver show attachment points. Matching color in external label: ";
    for(int i=0;i<std::min(n_ext,10);i++)tex<<"\\textcolor{"<<COLS[i]<<"}{ext"<<i+1<<"} ";
    tex<<"\n\nTotal: "<<entries.size()<<" bases ("<<n_uni<<" uni, "<<n_nonuni<<" non-uni).\n\n";

    // T distribution
    std::map<int,int>Tc;for(auto&e:entries)Tc[e.T]++;
    tex<<"\\begin{center}\n\\begin{tabular}{c|c}\n\\toprule\n$T$ & Count \\\\\n\\midrule\n";
    for(auto&[T,c]:Tc)tex<<T<<" & "<<c<<" \\\\\n";
    tex<<"\\bottomrule\n\\end{tabular}\n\\end{center}\n\n";

    // Split entries into unimodular and non-unimodular
    std::vector<Entry> uni_entries, nonuni_entries;
    for(auto&e:entries){
        if(std::abs(e.det)==1) uni_entries.push_back(e);
        else nonuni_entries.push_back(e);
    }

    // Lambda to write a block of entries
    auto write_entries=[&](std::vector<Entry>&ents, int max_n, const std::string& omit_label){
        std::string cur_combo;int cur_T=-1;int shown=0;
        for(auto&e:ents){
            if(shown>=max_n)break;
            if(e.combo_key!=cur_combo){cur_combo=e.combo_key;cur_T=-1;
                // Convert combo key to LaTeX: split by + and map each tag
                auto combo_to_latex = [&](const std::string& combo) -> std::string {
                    std::string result;
                    std::istringstream ss(combo);
                    std::string part;
                    while(std::getline(ss, part, '+')) {
                        if(!result.empty()) result += " + ";
                        if(part=="su2") result+="$\\mathfrak{su}_2$";
                        else if(part=="su3") result+="$\\mathfrak{su}_3$";
                        else if(part=="so8") result+="$\\mathfrak{so}_8$";
                        else if(part=="f4") result+="$\\mathfrak{f}_4$";
                        else if(part=="e6") result+="$\\mathfrak{e}_6$";
                        else if(part=="e7") result+="$\\mathfrak{e}_7$";
                        else if(part=="e7p") result+="$\\mathfrak{e}_7'$";
                        else if(part=="e8") result+="$\\mathfrak{e}_8$";
                        else if(part=="su2n3") result+="$\\mathfrak{su}_2{\\to}(-3)$";
                        else if(part=="su8") result+="$\\mathfrak{su}_8{\\to}\\mathfrak{so}_{16}$";
                        else if(part=="su3n2") result+="$(-3){\\to}(-2)$";
                        else if(part=="so16n2") result+="$\\mathfrak{so}_{16}{\\to}\\mathfrak{su}_8$";
                        else if(part=="hat1m1") result+="$\\hat{1}{\\to}(-1)$";
                        else if(part=="hat1m2") result+="$\\hat{1}{\\to}(-2)$";
                        else if(part=="2" || part=="2mix") result+="$(-2)$";  // no gauge: bare (-2)
                        else if(part.find("nhc_")!=std::string::npos) {
                            if(part=="nhc_2_3") result+="NHC$(-2)(-3)$";
                            else if(part=="nhc_2_3_2") result+="NHC$(-2)(-3)(-2)$";
                            else if(part=="nhc_2_2_3") result+="NHC$(-2)(-2)(-3)$";
                            else result+=part;
                        }
                        else result+=part;
                    }
                    return result;
                };
                tex<<"\n\\subsubsection*{"<<combo_to_latex(e.combo_key)<<"}\n\n";}
            if(e.T!=cur_T){cur_T=e.T;tex<<"\\paragraph{$T="<<cur_T<<"$}\n\n";}

            int delta=e.H_charged-e.V+29*e.T-273;
            int Hn=273+e.V-e.H_charged-29*e.T;
            int Tm=std::min(e.T+max_t_add_from_fmulti(Hn),193);

            // Nominal c_ext: simple sum of central_charges using each external's
            // nominal spec gauge (no NHC merging).
            double c_nom=0.0;
            for(auto&ei:e.externals){
                GaugeInfo g;
                if(ei.is_hat1){
                    g=(ei.target_si==-1)?GAUGE_SU8:GAUGE_SU16;
                } else {
                    g=gauge_from_si(ei.ext_si);
                    if(g.dim==0 && ei.ext_si==-2) g=GAUGE_SU2;
                }
                if(g.dim>0) c_nom+=central_charge(g,ei.int_num);
            }

            // Proper c_ext (NHC-aware) and true c_ext (fiber-weighted, kernel).
            // Match the pipeline: check_nhc → enhance externals (gauge boost for
            // hat1 and m2 attaches) → compute_proper/true.
            std::set<int> ext_set;
            for(auto&ei:e.externals) ext_set.insert(ei.curve_idx);
            NHCResult nhc = check_nhc(e.IF, 8.0);
            double c_proper = -1.0, c_true = -1.0;
            if (nhc.passes) {
                // Enhance each external the way the pipeline does it.
                for (auto& ei : e.externals) {
                    if (ei.is_hat1) {
                        enhance_hat1_gauge(nhc, e.IF, ei.curve_idx, ei.target_si, 8.0);
                    } else if (ei.ext_si == -2 && ei.target_si == -1) {
                        enhance_external_m2_gauge(nhc, e.IF, ei.curve_idx, ei.target_si, 8.0);
                    }
                    if (!nhc.passes) break;
                }
                if (nhc.passes) {
                    c_proper = compute_proper_c_ext(e.IF, nhc, ext_set);
                    c_true   = compute_true_c_ext_fiber(e.IF, nhc, ext_set);
                }
            }

            tex<<"\\texttt{["<<e.entry_id<<"]} $"<<quiver(e.IF,e.externals)<<ext_annotation(e.externals,tag2label)<<"$";
            tex<<" \\quad $T="<<e.T<<",\\;T_{\\max}="<<Tm<<",\\;\\Delta="<<delta
               <<",\\;\\det="<<e.det<<",\\;(n_+,n_-,n_0)=("<<e.sig_pos<<","<<e.sig_neg<<","<<e.sig_zero
               <<"),\\;c_\\mathrm{ext}^{\\mathrm{nom}}="<<std::fixed<<std::setprecision(1)<<c_nom;
            if (c_proper >= 0) tex<<",\\;c_\\mathrm{ext}^{\\mathrm{proper}}="<<std::fixed<<std::setprecision(2)<<c_proper;
            else tex<<",\\;c_\\mathrm{ext}^{\\mathrm{proper}}=\\text{NHC fail}";
            if (c_true >= 0) tex<<",\\;c_\\mathrm{ext}^{\\mathrm{true}}="<<std::fixed<<std::setprecision(2)<<c_true;
            else tex<<",\\;c_\\mathrm{ext}^{\\mathrm{true}}=\\text{n/a}";
            tex<<"$\n\n";

            std::set<int>h1;for(auto&ei:e.externals)if(ei.is_hat1)h1.insert(ei.curve_idx);

            // Unified T_min/T_max + blowdown for both uni and non-uni.
            // Pass sig_neg + T_max_anom so compute_tcrit can use the
            // reject_nonuni-style scan fallback when A=0.
            auto tc=compute_tcrit(e.IF,h1,e.sig_neg,Tm);
            if(tc.ival>=0){
                tex<<"\\hspace{1em}$T_{\\min}=";
                if(tc.exact)tex<<tc.ival;else tex<<std::fixed<<std::setprecision(2)<<tc.val<<"\\to "<<tc.ival;
                tex<<",\\;T_{\\max}="<<Tm<<",\\;T\\in["<<tc.ival<<","<<Tm<<"]$\n\n";
                tex<<"\\hspace{2em}$";wsm(tex,tc.bd_M);tex<<"_{\\!b_0}";
                // Identify surface for unimodular
                if(std::abs(e.det)==1){
                    auto bd=blowdown(e.IF,h1);
                    tex<<"\\to "<<bd.surface;
                }
                tex<<"$\n\n";
            }
            shown++;
        }
        if((int)ents.size()>max_n)
            tex<<"\\bigskip\\ldots "<<ents.size()-max_n<<" more "<<omit_label<<" omitted.\n\n";
    };

    // Unimodular section
    if(!uni_entries.empty()){
        tex<<"\\subsection*{Unimodular ($|\\det|=1$) --- "<<uni_entries.size()<<" bases}\n\n";
        write_entries(uni_entries, max_show, "unimodular");
    }

    // Non-unimodular section
    if(!nonuni_entries.empty()){
        tex<<"\\subsection*{Non-unimodular ($|\\det|>1$) --- "<<nonuni_entries.size()<<" bases}\n\n";
        write_entries(nonuni_entries, max_show, "non-unimodular");
    }

    tex<<"\\end{document}\n";tex.close();
    std::cout<<"Written "<<texfile<<" ("<<uni_entries.size()<<"+"<<nonuni_entries.size()<<" entries)\n";
    return 0;
}
