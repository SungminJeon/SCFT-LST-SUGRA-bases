// glued2tex.cpp — Phase 3 glued → LaTeX with cat2tex-style quivers
// Usage: ./glued2tex <glued_cat_dir> [max_entries=500]

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <map>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <functional>

// ── Reuse Entry struct and quiver from cat2tex ──

struct Entry {
    int entry_id = 0;
    std::string combo_key;
    struct Ext { int curve_idx=-1, spec_id=-1; std::string tag; int ext_si=0, target_si=0, int_num=1; bool is_hat1=false; };
    std::vector<Ext> externals;
    int T=0, H_charged=0, V=0, H_neutral=0, det=0, sig_pos=0, sig_neg=0, sig_zero=0;
    std::string glue_line, base_line;
    Eigen::MatrixXi IF;
    Eigen::MatrixXi BD;
    int bd_sp=0, bd_sn=0, bd_sz=0;
};

// ── Quiver renderer (from cat2tex.cpp) ──
static const char* COLS[]={"red","blue","green!60!black","orange","purple","cyan","magenta","brown!60!red","olive","pink"};

std::string quiver(const Eigen::MatrixXi& IF, const std::vector<Entry::Ext>& exts) {
    int n=IF.rows(); if(n==0)return"\\emptyset";
    std::set<int> ext_curves;
    std::map<int,std::vector<int>> target_colors;
    for(int i=0;i<(int)exts.size();i++) ext_curves.insert(exts[i].curve_idx);
    for(int i=0;i<(int)exts.size();i++){
        int ei=exts[i].curve_idx;
        for(int j=0;j<n;j++){
            if(j==ei) continue;
            if(IF(ei,j)!=0 && !ext_curves.count(j))
                target_colors[j].push_back(i);
        }
    }
    int n_base=0;
    std::map<int,int> remap; std::vector<int> unmap;
    for(int i=0;i<n;i++){
        if(ext_curves.count(i)) continue;
        remap[i]=n_base; unmap.push_back(i); n_base++;
    }
    if(n_base==0) return"\\emptyset";
    std::vector<std::vector<std::pair<int,int>>>adj(n_base);
    for(int i=0;i<n_base;i++){int oi=unmap[i];
        for(int j=i+1;j<n_base;j++){int oj=unmap[j];
            if(IF(oi,oj)!=0){adj[i].push_back({j,IF(oi,oj)});adj[j].push_back({i,IF(oj,oi)});}}}
    int start=-1;
    for(int i=0;i<n_base;i++)
        if((int)adj[i].size()==1&&!target_colors.count(unmap[i])){start=i;break;}
    if(start<0)for(int i=0;i<n_base;i++)if((int)adj[i].size()==1){start=i;break;}
    if(start<0)start=0;
    std::vector<bool>vis(n_base,false);
    auto lb=[&](int c)->std::string{
        int oc=unmap[c]; int si=-IF(oc,oc);
        std::string s=std::to_string(si);if(si>=10)s="("+s+")";
        auto it=target_colors.find(oc);
        if(it!=target_colors.end()&&!it->second.empty()){
            std::string sup;
            for(int ei:it->second){int ci=ei%10;int ext_curve=exts[ei].curve_idx;
                int int_num=std::abs(IF(oc,ext_curve));
                sup+="\\textcolor{"+std::string(COLS[ci])+"}{"+std::to_string(int_num)+"}";}
            return s+"^{"+sup+"}";}
        return s;};
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

// ── Load glued .cat files ──
std::vector<Entry> load_glued(const std::string& dir) {
    std::vector<Entry> entries;
    for (auto& p : std::filesystem::directory_iterator(dir)) {
        if (p.path().extension() != ".cat") continue;
        std::ifstream fin(p.path());
        std::string line;
        Entry cur; bool in_entry=false;
        int ifl=0, bdl=0; bool in_if=false, in_bd=false, in_ext=false;
        int ext_count=0, ifn=0, bdn=0;
        std::vector<std::vector<int>> if_rows, bd_rows;

        while (std::getline(fin, line)) {
            if (line.empty() || line[0]=='#') continue;
            if (line.substr(0,6)=="ENTRY ") {
                if (in_entry && cur.IF.rows()>0) entries.push_back(std::move(cur));
                cur = Entry(); if_rows.clear(); bd_rows.clear();
                cur.entry_id = std::stoi(line.substr(6));
                in_entry=true; in_if=in_bd=in_ext=false;
            }
            else if (in_entry && line.substr(0,6)=="COMBO ") cur.combo_key=line.substr(6);
            else if (in_entry && line.substr(0,5)=="GLUE ") cur.glue_line=line;
            else if (in_entry && line.substr(0,5)=="BASE ") cur.base_line=line;
            else if (in_entry && line.substr(0,8)=="PHYSICS ") {
                std::istringstream ss(line.substr(8));
                ss>>cur.T>>cur.H_charged>>cur.V>>cur.H_neutral>>cur.det>>cur.sig_pos>>cur.sig_neg>>cur.sig_zero;
            }
            else if (in_entry && line.substr(0,10)=="EXTERNALS ") {
                ext_count=std::stoi(line.substr(10)); in_ext=true;
            }
            else if (in_entry && in_ext && ext_count>0) {
                std::istringstream ss(line);
                Entry::Ext ex; ss>>ex.curve_idx>>ex.spec_id>>ex.tag>>ex.ext_si>>ex.target_si>>ex.int_num;
                int h; ss>>h; ex.is_hat1=(h!=0);
                cur.externals.push_back(ex);
                ext_count--; if(ext_count==0) in_ext=false;
            }
            else if (in_entry && line.substr(0,3)=="IF ") {
                ifn=std::stoi(line.substr(3)); in_if=true; ifl=ifn;
            }
            else if (in_entry && in_if && ifl>0) {
                std::istringstream ss(line); std::vector<int> row;
                int v; while(ss>>v) row.push_back(v);
                if_rows.push_back(row); ifl--;
                if(ifl==0) {
                    in_if=false;
                    cur.IF=Eigen::MatrixXi::Zero(ifn,ifn);
                    for(int i=0;i<ifn&&i<(int)if_rows.size();i++)
                        for(int j=0;j<ifn&&j<(int)if_rows[i].size();j++)
                            cur.IF(i,j)=if_rows[i][j];
                }
            }
            else if (in_entry && line.substr(0,9)=="BLOWDOWN ") {
                bdn=std::stoi(line.substr(9)); in_bd=true; bdl=bdn;
            }
            else if (in_entry && in_bd && bdl>0) {
                std::istringstream ss(line); std::vector<int> row;
                int v; while(ss>>v) row.push_back(v);
                bd_rows.push_back(row); bdl--;
                if(bdl==0) {
                    in_bd=false;
                    cur.BD=Eigen::MatrixXi::Zero(bdn,bdn);
                    for(int i=0;i<bdn&&i<(int)bd_rows.size();i++)
                        for(int j=0;j<bdn&&j<(int)bd_rows[i].size();j++)
                            cur.BD(i,j)=bd_rows[i][j];
                }
            }
            else if (in_entry && line.substr(0,6)=="BDSIG ") {
                std::istringstream ss(line.substr(6));
                ss>>cur.bd_sp>>cur.bd_sn>>cur.bd_sz;
            }
            else if (in_entry && line=="END") {
                entries.push_back(std::move(cur));
                cur=Entry(); if_rows.clear(); bd_rows.clear();
                in_entry=false;
            }
        }
        if (in_entry && cur.IF.rows()>0) entries.push_back(std::move(cur));
    }
    return entries;
}

int main(int argc, char* argv[]) {
    if (argc<2) { std::cerr<<"Usage: "<<argv[0]<<" <dir> [max=500]\n"; return 1; }
    std::string dir=argv[1];
    int max_show=(argc>2)?std::atoi(argv[2]):500;

    auto entries=load_glued(dir);
    std::cout<<"Loaded "<<entries.size()<<" entries\n";

    std::sort(entries.begin(),entries.end(),[](const Entry&a,const Entry&b){
        if(a.T!=b.T)return a.T<b.T; return std::abs(a.det)<std::abs(b.det);});

    auto wsm=[](std::ofstream&t,const Eigen::MatrixXi&M){
        int s=M.rows();t<<"\\left(\\begin{smallmatrix}";
        for(int i=0;i<s;i++){if(i)t<<" \\\\ ";for(int j=0;j<s;j++){if(j)t<<" & ";t<<M(i,j);}}
        t<<"\\end{smallmatrix}\\right)";};

    std::string texfile=dir+"/results_glued.tex";
    std::ofstream tex(texfile);
    tex<<"\\documentclass[10pt]{article}\n"
       <<"\\usepackage[a4paper,margin=0.5in]{geometry}\n"
       <<"\\usepackage{amsmath,amssymb,array}\n"
       <<"\\usepackage[dvipsnames]{xcolor}\n"
       <<"\\usepackage{hyperref}\n"
       <<"\\hypersetup{colorlinks=true,linkcolor=blue}\n"
       <<"\\setlength{\\parindent}{0pt}\n\\setlength{\\parskip}{2pt}\n"
       <<"\\begin{document}\n\n"
       <<"\\section*{Phase 3 Glued SUGRA Bases}\n\n"
       <<"Total: "<<entries.size()<<" entries (showing "<<std::min((int)entries.size(),max_show)<<").\n\n"
       <<"$\\mathrm{Base\\;A}\\;\\oplus\\;\\mathrm{Base\\;B}$ via shared external.\n\n";

    int cur_T=-1, shown=0;
    for(auto&e:entries){
        if(shown>=max_show)break;
        if(e.IF.rows()==0)continue;
        shown++;

        int n=e.IF.rows();
        int ext_idx=e.externals.empty()?-1:e.externals[0].curve_idx;

        // Build A sub-IF and B sub-IF
        // A: curves 0..ext_idx (inclusive)
        // B: curves ext_idx..n-1 (inclusive)
        int nA=ext_idx+1, nB=n-ext_idx;

        Eigen::MatrixXi IFA=e.IF.block(0,0,nA,nA);
        Eigen::MatrixXi IFB=Eigen::MatrixXi::Zero(nB,nB);
        for(int i=0;i<nB;i++)for(int j=0;j<nB;j++)
            IFB(i,j)=e.IF(ext_idx+i,ext_idx+j);

        // A's external: last curve (ext_idx) in local coords = nA-1
        std::vector<Entry::Ext> extA={{nA-1, 0, e.combo_key, -2, -1, 1, false}};
        // B's external: first curve (0) in local coords
        std::vector<Entry::Ext> extB={{0, 0, e.combo_key, -2, -1, 1, false}};

        std::string qA=quiver(IFA, extA);
        std::string qB=quiver(IFB, extB);

        // Parse GLUE
        std::string typeA="?",tA="?",typeB="?",tB="?";
        if(!e.glue_line.empty()){
            std::istringstream gs(e.glue_line);std::string d;
            gs>>d>>typeA>>tA>>typeB>>tB;}

        if(e.T!=cur_T){cur_T=e.T;tex<<"\\paragraph{$T="<<cur_T<<"$}\n\n";}

        int Tm=std::min((273+e.V-e.H_charged)/29,193);

        int delta=e.H_charged-e.V+29*e.T-273;

        tex<<"\\texttt{["<<e.entry_id<<"]} $"<<qA
           <<" + \\textcolor{red}{\\mathfrak{su}_2} + "<<qB
           <<"$ \\quad $T="<<e.T<<",\\;T_{\\max}="<<Tm
           <<",\\;\\Delta="<<delta<<",\\;\\det="<<e.det
           <<",\\;(n_+,n_-,n_0)=("<<e.sig_pos<<","<<e.sig_neg<<","<<e.sig_zero<<")$\n\n";

        tex<<"\\hspace{1em}{\\small "<<typeA<<"($T$="<<tA<<") $\\oplus$ "
           <<typeB<<"($T$="<<tB<<")}\n\n";

        if(e.BD.rows()>0){
            int bd_det=(int)std::round(e.BD.cast<double>().determinant());
            tex<<"\\hspace{2em}$";wsm(tex,e.BD);
            tex<<"_{\\!b_0}$ \\quad $\\det="<<bd_det
               <<",\\;\\mathrm{sig}=("<<e.bd_sp<<","<<e.bd_sn<<","<<e.bd_sz<<")$\n\n";
        }
    }

    if((int)entries.size()>max_show)
        tex<<"\\bigskip\\ldots "<<entries.size()-max_show<<" more omitted.\n\n";
    tex<<"\\end{document}\n";tex.close();
    std::cout<<"Written "<<texfile<<" ("<<shown<<" entries)\n";
    return 0;
}
