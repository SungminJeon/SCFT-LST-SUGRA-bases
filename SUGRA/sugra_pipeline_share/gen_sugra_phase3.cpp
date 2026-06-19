// gen_sugra_phase3.cpp
// Phase 3: Glue two Phase 1 entries sharing the same external curve
//
// For each external spec:
//   for each pair (A, B) of Phase 1 entries with that spec:
//     merge A.ext_curve_idx with B.ext_curve_idx
//     → combined IF = LST_A — external — LST_B
//     filter: NHC, anomaly, signature
//
// Compile:
//   g++ -std=c++17 -O2 -o gen_sugra_phase3 gen_sugra_phase3.cpp Tensor.C \
//       TopoLineCompact_enhanced.cpp Topology_enhanced.cpp \
//       -I/usr/include/eigen3 -I.
//
// Usage:
//   ./gen_sugra_phase3 <phase1_dir> [T_sugra_max=193]

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <numeric>
#include <map>
#include <algorithm>
#include <functional>
#include <filesystem>

// === ExternalSpec ===
struct ExternalSpec {
    int id,ext_si,target_si,int_num; bool is_hat1; double cc; std::string label,tag;
};
inline double compute_ext_cc(int es,int in,bool h,int ts){
    if(h)return(ts==-1)?7.0:(ts==-2)?15.0:10.0;
    GaugeInfo g=gauge_from_si(es);if(g.dim>0)return central_charge(g,in);
    if(es==-2&&ts==-1)return central_charge(GAUGE_SU2,in);
    if(es==-2&&ts==-3)return 5.0;if(es==-2&&ts==-4)return 10.0;
    if(es==-3&&ts==-2)return 5.0;if(es==-4&&ts==-2)return 10.0;return 5.0;}
inline std::vector<ExternalSpec> build_all_specs(){
    std::vector<ExternalSpec> sp;int id=0;
    auto a=[&](int es,int ts,int in,bool h,const std::string&l,const std::string&t){
        sp.push_back({id++,es,ts,in,h,compute_ext_cc(es,in,h,ts),l,t});};
    for(int k=1;k<=9;k++)a(-2,-1,k,false,"$\\mathfrak{su}_2$"+(k>1?"$(k\\!=\\!"+std::to_string(k)+")$":std::string("")),k==1?"su2":"su2k"+std::to_string(k));
    a(-2,-3,1,false,"$2{\\to}3$","su2n3");a(-2,-4,2,false,"$\\mathfrak{su}_8{\\to}\\mathfrak{so}_{16}$","su8");
    for(int k=1;k<=9;k++)a(-3,-1,k,false,"$\\mathfrak{su}_3$"+(k>1?"$(k\\!=\\!"+std::to_string(k)+")$":std::string("")),"su3k"+std::to_string(k));
    a(-3,-2,1,false,"$3{\\to}2$","su3n2");a(-4,-1,1,false,"$\\mathfrak{so}_8$","so8");a(-4,-2,2,false,"$4{\\to}2$","so16n2");
    a(-5,-1,1,false,"$\\mathfrak{f}_4$","f4");a(-6,-1,1,false,"$\\mathfrak{e}_6$","e6");
    a(-7,-1,1,false,"$\\mathfrak{e}_7'$","e7p");a(-8,-1,1,false,"$\\mathfrak{e}_7$","e7");
    a(-12,-1,1,false,"$\\mathfrak{e}_8$","e8");
    a(-1,-1,1,true,"$\\hat{1}{\\to}1$","hat1m1");a(-1,-2,1,true,"$\\hat{1}{\\to}2$","hat1m2");
    // 2 / 2mix: gauge-less (-2) attached to LST (-2) that is part of nhc_2_3.
    a(-2,-2,1,false,"$(\\!-\\!2){\\to}(\\!-\\!2)$","2");
    a(-2,-2,1,false,"$(\\!-\\!2)_{mix}$","2mix");
    return sp;}

// === Phase 1 loader ===
struct Phase1Entry {
    int entry_id,spec_id;std::string spec_tag;int catalog_id;std::string catalog_type;
    int base_T,n_base,target_idx,ext_curve_idx;std::vector<int>hat1_multi_targets;
    int T,H_charged,V,H_neutral,sig_pos,sig_neg,sig_zero;long long det;Eigen::MatrixXi final_IF;
    struct RemTarget{int curve_idx,self_int;double available_cc;};std::vector<RemTarget>remaining;};
std::vector<Phase1Entry> load_phase1_cat(const std::string&fn){
    std::vector<Phase1Entry>es;std::ifstream f(fn);if(!f)return es;std::string line;Phase1Entry c;bool in=false;
    int ifl=0,ifr=0,ifn=0,rl=0;
    while(std::getline(f,line)){if(line.empty()||line[0]=='#')continue;
        if(line.substr(0,6)=="ENTRY "){c=Phase1Entry();c.entry_id=std::stoi(line.substr(6));in=true;}
        else if(in&&line.substr(0,5)=="SPEC "){std::istringstream s(line.substr(5));std::string t;int a,b,d,h;s>>c.spec_id>>t>>a>>b>>d>>h;c.spec_tag=t;}
        else if(in&&line.substr(0,5)=="BASE "){std::istringstream s(line.substr(5));s>>c.catalog_id>>c.catalog_type>>c.base_T>>c.n_base;}
        else if(in&&line.substr(0,7)=="ATTACH "){std::istringstream s(line.substr(7));s>>c.target_idx>>c.ext_curve_idx;
            std::string r;if(std::getline(s,r)){auto p=r.find("multi:");if(p!=std::string::npos){std::istringstream m(r.substr(p+6));std::string tk;while(std::getline(m,tk,','))if(!tk.empty())c.hat1_multi_targets.push_back(std::stoi(tk));}}}
        else if(in&&line.substr(0,8)=="PHYSICS "){std::istringstream s(line.substr(8));s>>c.T>>c.H_charged>>c.V>>c.H_neutral>>c.det>>c.sig_pos>>c.sig_neg>>c.sig_zero;}
        else if(in&&line.substr(0,3)=="IF "){ifn=std::stoi(line.substr(3));c.final_IF=Eigen::MatrixXi::Zero(ifn,ifn);ifl=ifn;ifr=0;}
        else if(in&&ifl>0){std::istringstream s(line);for(int j=0;j<ifn;j++)s>>c.final_IF(ifr,j);ifr++;ifl--;}
        else if(in&&line.substr(0,10)=="REMAINING "){rl=std::stoi(line.substr(10));}
        else if(in&&rl>0){Phase1Entry::RemTarget rt;std::istringstream s(line);s>>rt.curve_idx>>rt.self_int>>rt.available_cc;c.remaining.push_back(rt);rl--;}
        else if(in&&line=="END"){es.push_back(std::move(c));in=false;}}return es;}

// === Glue two Phase1 IFs through their external curves ===
// A.ext_curve and B.ext_curve are identified (merged into one curve).
// Result size = n_A + n_B - 1.
// The shared external curve keeps its self-intersection.
// Returned: merged IF, index of shared ext in merged IF, range of B's curves.
struct GlueResult {
    Eigen::MatrixXi IF;
    int shared_ext_idx;      // external curve index in merged IF
    int b_start, b_end;      // range of B's non-ext curves in merged IF
};

GlueResult glue_through_external(const Eigen::MatrixXi& IF_A, int ext_A,
                                   const Eigen::MatrixXi& IF_B, int ext_B) {
    int na = IF_A.rows(), nb = IF_B.rows();
    int nt = na + nb - 1;
    Eigen::MatrixXi G = Eigen::MatrixXi::Zero(nt, nt);

    // A occupies indices [0, na) as-is
    G.block(0, 0, na, na) = IF_A;

    // B's ext_B maps to ext_A; others map to [na, na+nb-2)
    auto map_b = [&](int j) -> int {
        if (j == ext_B) return ext_A;
        return na + (j < ext_B ? j : j - 1);
    };

    // Copy B diagonal (skip ext_B)
    for (int j = 0; j < nb; j++) {
        if (j == ext_B) continue;
        G(map_b(j), map_b(j)) = IF_B(j, j);
    }

    // Copy B off-diagonal (add to existing for the merged ext curve)
    for (int j = 0; j < nb; j++)
        for (int k = j + 1; k < nb; k++)
            if (IF_B(j, k) != 0) {
                G(map_b(j), map_b(k)) += IF_B(j, k);
                G(map_b(k), map_b(j)) += IF_B(k, j);
            }

    GlueResult r;
    r.IF = G;
    r.shared_ext_idx = ext_A;
    r.b_start = na;
    r.b_end = nt;
    return r;
}

// === Phase3 result ===
struct Phase3Result {
    int spec_id; std::string spec_tag;
    int catA_id; std::string catA_type;
    int catB_id; std::string catB_type;
    int shared_ext_idx;
    int b_start, b_end;
    bool ext_is_hat1;
    Eigen::MatrixXi final_IF;
    NHCResult nhc; AnomalyResult anomaly; SigInfo sig;
    std::set<int> hat1_indices;
    int T_min_target = 0;   // the matched T_min
    int T_A = 0, T_B = 0;  // T contributions from each side
};

// === Dedup v2 ===
inline std::string dedup_key_v2(const Eigen::MatrixXi&IF,const std::string&extra){
    int n=IF.rows();Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>es(IF.cast<double>());auto ev=es.eigenvalues();
    std::ostringstream ss;ss<<std::fixed<<std::setprecision(8);ss<<"E:";
    for(int i=0;i<ev.size();i++){double v=ev(i);if(std::abs(v)<1e-10)v=0;if(i)ss<<",";ss<<v;}
    struct VD{int si,deg;std::vector<int>ns;bool operator<(const VD&o)const{return std::tie(si,deg,ns)<std::tie(o.si,o.deg,o.ns);}};
    std::vector<VD>vds(n);for(int i=0;i<n;i++){vds[i].si=IF(i,i);vds[i].deg=0;
        for(int j=0;j<n;j++)if(j!=i&&IF(i,j)!=0){vds[i].deg++;vds[i].ns.push_back(IF(j,j));}
        std::sort(vds[i].ns.begin(),vds[i].ns.end());}std::sort(vds.begin(),vds.end());
    ss<<"|V:";for(auto&v:vds){ss<<"("<<v.si<<","<<v.deg<<",{";for(int k=0;k<(int)v.ns.size();k++){if(k)ss<<",";ss<<v.ns[k];}ss<<"})";}
    std::vector<std::tuple<int,int,int>>edges;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(IF(i,j)!=0){
        int a=IF(i,i),b=IF(j,j);if(a>b)std::swap(a,b);edges.push_back({a,b,IF(i,j)});}
    std::sort(edges.begin(),edges.end());ss<<"|E:";for(auto&[a,b,w]:edges)ss<<"("<<a<<","<<b<<","<<w<<")";
    ss<<"|X:"<<extra;return ss.str();}

// === LaTeX: red=external (shared), blue=side B ===
static const char*CE="red",*CB="blue";
std::string if_to_latex_p3(const Eigen::MatrixXi&IF,int ei,bool eh,int bs,int be){
    int n=IF.rows();if(n==0)return"\\emptyset";if(n==1)return std::to_string(-IF(0,0));
    std::set<int>es,bset,sp;if(ei>=0){es.insert(ei);sp.insert(ei);}
    for(int i=bs;i<be;i++){bset.insert(i);sp.insert(i);}
    std::vector<std::vector<std::pair<int,int>>>adj(n);
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(IF(i,j)!=0){adj[i].push_back({j,IF(i,j)});adj[j].push_back({i,IF(j,i)});}
    int start=-1;for(int i=0;i<n;i++)if((int)adj[i].size()==1&&!sp.count(i)){start=i;break;}
    if(start<0)for(int i=0;i<n;i++)if((int)adj[i].size()==1){start=i;break;}if(start<0)start=0;
    std::vector<bool>vis(n,false);
    auto lb=[&](int c)->std::string{if(c==ei&&eh)return"\\textcolor{"+std::string(CE)+"}{\\hat{1}}";
        int si=-IF(c,c);std::string s=std::to_string(si);if(si>=10)s="("+s+")";
        if(es.count(c))return"\\textcolor{"+std::string(CE)+"}{"+s+"}";
        if(bset.count(c))return"\\textcolor{"+std::string(CB)+"}{"+s+"}";return s;};
    std::function<std::string(int)>wb;
    std::function<std::string(int,int)>w=[&](int c,int f)->std::string{vis[c]=true;std::string l=lb(c);
        std::vector<int>ch;for(auto&[nb,ww]:adj[c])if(!vis[nb])ch.push_back(nb);
        if(ch.empty())return l;if(ch.size()==1)return l+w(ch[0],c);
        int mc=-1;std::vector<int>bc;for(int x:ch){if(mc<0&&!sp.count(x))mc=x;else bc.push_back(x);}
        if(mc<0){mc=ch[0];bc.clear();for(size_t i=1;i<ch.size();i++)bc.push_back(ch[i]);}
        std::vector<std::string>bs;for(int i=(int)bc.size()-1;i>=0;i--)bs.push_back(wb(bc[i]));
        std::string mp=w(mc,c);if(bs.empty())return l+mp;if(bs.size()==1)return"\\overset{"+bs[0]+"}{"+l+"}"+mp;
        std::string u;for(size_t i=1;i<bs.size();i++){if(!u.empty())u+=",";u+=bs[i];}
        return"\\overset{"+bs[0]+"}{\\underset{"+u+"}{"+l+"}}"+mp;};
    wb=[&](int c)->std::string{vis[c]=true;std::string l=lb(c);
        std::vector<int>ch;for(auto&[nb,ww]:adj[c])if(!vis[nb])ch.push_back(nb);
        if(ch.empty())return l;if(ch.size()==1)return wb(ch[0])+l;
        int mc=-1;std::vector<int>sb;for(int x:ch){if(mc<0&&!sp.count(x))mc=x;else sb.push_back(x);}
        if(mc<0){mc=ch[0];sb.clear();for(size_t i=1;i<ch.size();i++)sb.push_back(ch[i]);}
        std::vector<std::string>ss2;for(int i=(int)sb.size()-1;i>=0;i--)ss2.push_back(wb(sb[i]));
        std::string mp=wb(mc);if(ss2.size()==1)return mp+"\\overset{"+ss2[0]+"}{"+l+"}";
        std::string u;for(size_t i=1;i<ss2.size();i++){if(!u.empty())u+=",";u+=ss2[i];}
        return mp+"\\overset{"+ss2[0]+"}{\\underset{"+u+"}{"+l+"}}";};
    return w(start,-1);}

// === Blowdown + CritInfo ===
inline Eigen::MatrixXi blowdown_curve(const Eigen::MatrixXi&IF,int idx){int n=IF.rows();Eigen::MatrixXi r=Eigen::MatrixXi::Zero(n-1,n-1);
    for(int i=0,ni=0;i<n;i++){if(i==idx)continue;for(int j=0,nj=0;j<n;j++){if(j==idx)continue;r(ni,nj)=IF(i,j)+IF(i,idx)*IF(idx,j);nj++;}ni++;}return r;}
inline Eigen::MatrixXi build_b0Q_mixed(const Eigen::MatrixXi&IF,const std::set<int>&h1){int n=IF.rows(),T=0;
    {Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>es(IF.cast<double>());for(int i=0;i<n;i++)if(es.eigenvalues()(i)<-1e-8)T++;}
    Eigen::MatrixXi e=Eigen::MatrixXi::Zero(n+1,n+1);e.block(0,0,n,n)=IF;
    for(int i=0;i<n;i++){int b=h1.count(i)?-1:(IF(i,i)+2);e(i,n)=b;e(n,i)=b;}e(n,n)=9-T;return e;}
struct BlowdownResult{Eigen::MatrixXi extended_IF;std::string surface;int F_n=-1;bool is_P2=false;};
BlowdownResult blowdown_mixed(const Eigen::MatrixXi&IF,const std::set<int>&h1){
    Eigen::MatrixXi cur=build_b0Q_mixed(IF,h1);std::set<int>ch=h1;int bd=0;
    while(true){int nn=cur.rows(),m=-1;for(int i=0;i<nn-1;i++){if(ch.count(i))continue;if(cur(i,i)==-1){m=i;break;}}
        if(m<0)break;cur=blowdown_curve(cur,m);std::set<int>nh;for(int h:ch)nh.insert(h>m?h-1:h);ch=nh;bd++;}
    BlowdownResult r;r.extended_IF=cur;int nc=cur.rows()-1;
    if(nc==0)r.surface="\\text{pt}";else if(nc==1&&cur(0,0)==1){r.is_P2=true;r.surface="\\mathbb{P}^2";}
    else if(nc==2&&cur(0,1)==1){int a=cur(0,0),b=cur(1,1);if(a==0||b==0){r.F_n=(a==0)?std::abs(b):std::abs(a);r.surface="\\mathbb{F}_{"+std::to_string(r.F_n)+"}";}else r.surface="\\text{bd}_{"+std::to_string(bd)+"}";}
    else r.surface="\\text{bd}_{"+std::to_string(bd)+"}";return r;}
struct CritInfo{double T_crit;int T_crit_int;bool exact;int crit_npos,crit_nneg,crit_nzero;BlowdownResult bd_crit;};
CritInfo compute_crit(const Eigen::MatrixXi&IF,const std::set<int>&h1){CritInfo res;int n=IF.rows();
    Eigen::MatrixXd M=Eigen::MatrixXd::Zero(n+1,n+1);M.block(0,0,n,n)=IF.cast<double>();
    for(int i=0;i<n;i++){double b=h1.count(i)?-1.0:(IF(i,i)+2.0);M(i,n)=b;M(n,i)=b;}
    M(n,n)=0;double B=M.determinant();M(n,n)=1;double A=M.determinant()-B;
    if(std::abs(A)<1e-10){res.T_crit=1e9;res.exact=false;res.T_crit_int=-1;return res;}
    double b0c=-B/A;res.T_crit=9.0-b0c;int tc=(int)std::round(res.T_crit);res.exact=(std::abs(res.T_crit-tc)<1e-6);
    res.T_crit_int=res.exact?tc:(int)std::ceil(res.T_crit-1e-9);int b0sq=9-res.T_crit_int;
    M(n,n)=(double)b0sq;Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>es(M);res.crit_npos=res.crit_nneg=res.crit_nzero=0;
    for(int i=0;i<n+1;i++){double ev=es.eigenvalues()(i);if(std::abs(ev)<1e-8)res.crit_nzero++;else if(ev>0)res.crit_npos++;else res.crit_nneg++;}
    Eigen::MatrixXi ext=Eigen::MatrixXi::Zero(n+1,n+1);ext.block(0,0,n,n)=IF;
    for(int i=0;i<n;i++){int b=h1.count(i)?-1:(IF(i,i)+2);ext(i,n)=b;ext(n,i)=b;}ext(n,n)=b0sq;
    Eigen::MatrixXi cur=ext;std::set<int>ch=h1;
    while(true){int nn=cur.rows(),m=-1;for(int i=0;i<nn-1;i++){if(ch.count(i))continue;if(cur(i,i)==-1){m=i;break;}}
        if(m<0)break;cur=blowdown_curve(cur,m);std::set<int>nh;for(int h:ch)nh.insert(h>m?h-1:h);ch=nh;}
    res.bd_crit.extended_IF=cur;int nc=cur.rows()-1;
    if(nc==0)res.bd_crit.surface="\\text{pt}";else if(nc==1&&cur(0,0)==1){res.bd_crit.is_P2=true;res.bd_crit.surface="\\mathbb{P}^2";}
    else if(nc==2&&cur(0,1)==1){int a=cur(0,0),b=cur(1,1);if(a==0||b==0){res.bd_crit.F_n=(a==0)?std::abs(b):std::abs(a);res.bd_crit.surface="\\mathbb{F}_{"+std::to_string(res.bd_crit.F_n)+"}";}else res.bd_crit.surface="\\text{bd}";}
    else res.bd_crit.surface="\\text{bd}";return res;}


// === Main ===
int main(int argc,char*argv[]){
    if(argc<2){std::cerr<<"Usage: "<<argv[0]<<" <phase1_dir> [T_sugra_max=193]\n"
                        <<"       "<<argv[0]<<" --cross <dirA> <dirB> [T_sugra_max=193]\n"
                        <<"       "<<argv[0]<<" --pair-tag <tagA> <tagB> <dir> [T_sugra_max=193]\n";return 1;}

    bool cross_mode = false;
    bool pair_tag_mode = false;
    std::string phase1_dir, dirA, dirB;
    std::string pair_tagA, pair_tagB;
    int T_sugra_max = 193;

    if (std::string(argv[1]) == "--cross") {
        cross_mode = true;
        if (argc < 4) { std::cerr << "Need two dirs for --cross\n"; return 1; }
        dirA = argv[2]; dirB = argv[3];
        if (argc > 4) T_sugra_max = std::atoi(argv[4]);
        std::cout << "=== gen_sugra_phase3: CROSS glue ===\n"
                  << "DirA: " << dirA << "\nDirB: " << dirB << "\n\n";
    } else if (std::string(argv[1]) == "--pair-tag") {
        pair_tag_mode = true;
        if (argc < 5) { std::cerr << "Need <tagA> <tagB> <dir> for --pair-tag\n"; return 1; }
        pair_tagA = argv[2]; pair_tagB = argv[3]; phase1_dir = argv[4];
        if (argc > 5) T_sugra_max = std::atoi(argv[5]);
        std::cout << "=== gen_sugra_phase3: PAIR-TAG glue ===\n"
                  << "TagA: " << pair_tagA << "  TagB: " << pair_tagB
                  << "\nDir: " << phase1_dir << "\n\n";
    } else {
        phase1_dir = argv[1];
        T_sugra_max = (argc > 2) ? std::atoi(argv[2]) : 193;
        std::cout << "=== gen_sugra_phase3: glue through shared external ===\n"
                  << "Phase1: " << phase1_dir << "\n\n";
    }

    auto all_specs=build_all_specs();double cc=8.0;

    // Load ALL Phase 1 entries per spec
    struct P1Entry {
        Phase1Entry entry;
        int T_min;
        const ExternalSpec* spec;
        int source;  // 0=default/dirA, 1=dirB (for cross mode)
    };

    auto load_dir = [&](const std::string& dir, int src, const std::vector<ExternalSpec>& specs,
                         std::map<std::string, std::vector<P1Entry>>& p1s) -> int {
        int loaded = 0;
        for (auto& sp : specs) {
            std::string fn = dir + "/" + sp.tag + ".cat";
            if (!std::filesystem::exists(fn)) continue;
            auto es = load_phase1_cat(fn);
            for (auto& e : es) {
                P1Entry pe;
                pe.entry = std::move(e);
                pe.spec = &sp;
                pe.source = src;
                if (std::abs(pe.entry.det) <= 1) {
                    pe.T_min = pe.entry.T;
                } else {
                    std::set<int> h1;
                    if (sp.is_hat1) h1.insert(pe.entry.ext_curve_idx);
                    auto ci = compute_crit(pe.entry.final_IF, h1);
                    if (ci.T_crit_int < 0) continue;
                    pe.T_min = ci.T_crit_int;
                }
                p1s[sp.tag].push_back(std::move(pe));
                loaded++;
            }
        }
        return loaded;
    };

    std::map<std::string, std::vector<P1Entry>> p1s;
    int total_loaded = 0;

    if (cross_mode) {
        total_loaded += load_dir(dirA, 0, all_specs, p1s);
        total_loaded += load_dir(dirB, 1, all_specs, p1s);
    } else if (pair_tag_mode) {
        // Load only the two requested tags from one dir.
        std::vector<ExternalSpec> only;
        for (auto& s : all_specs) if (s.tag == pair_tagA || s.tag == pair_tagB) only.push_back(s);
        if (only.size() < 2) {
            std::cerr << "Could not find both tags '" << pair_tagA << "' and '"
                      << pair_tagB << "' in spec list\n";
            return 1;
        }
        total_loaded += load_dir(phase1_dir, 0, only, p1s);
    } else {
        total_loaded += load_dir(phase1_dir, 0, all_specs, p1s);
    }
    std::cout << "Loaded: " << total_loaded << " entries across "
              << p1s.size() << " specs\n\n";

    // === Scan pairs ===
    std::vector<Phase3Result> results;
    int tried = 0, passed = 0;
    int rej_sig = 0, rej_nhc = 0, rej_enhance = 0, rej_cext = 0, rej_anom = 0, rej_cc = 0, rej_bdsig = 0;

    // Pair processor: tries to glue A,B through their externals. spec controls
    // gauge enhancement target/hat1 handling and tagging of the result.
    // result_tag overrides spec->tag for the stored result (used in pair-tag mode).
    auto process_pair = [&](P1Entry& A, P1Entry& B, const ExternalSpec* spec,
                             const std::string& result_tag) -> bool {
        // LST-T convention for the glued anomaly: sum of each side's LST
        // tensor count (BASE field 3, base_T) plus 1 for the shared external.
        int T_lst = A.entry.base_T + B.entry.base_T + 1;
        int T_est = A.entry.T + B.entry.T - 1;
        int T_need = std::max(A.T_min, B.T_min);
        if (T_est < T_need) return false;
        if (T_est > T_sugra_max) return false;

        tried++;
        auto gr = glue_through_external(A.entry.final_IF, A.entry.ext_curve_idx,
                                         B.entry.final_IF, B.entry.ext_curve_idx);

        auto sig = compute_sig_fast(gr.IF);
        int T_glued = sig.sig_neg;
        if (sig.sig_pos != 1) { rej_sig++; return false; }
        if (T_glued > T_sugra_max) { rej_sig++; return false; }

        NHCResult nhc = check_nhc(gr.IF, cc);
        if (!nhc.passes) { rej_nhc++; return false; }

        {
            bool in_nhc_cluster = false;
            if (gr.shared_ext_idx < (int)nhc.curve_cluster.size() &&
                nhc.curve_cluster[gr.shared_ext_idx] >= 0) in_nhc_cluster = true;
            if (!in_nhc_cluster) {
                if (spec->is_hat1)
                    enhance_hat1_gauge(nhc, gr.IF, gr.shared_ext_idx, spec->target_si, cc);
                else if (spec->ext_si == -2 && spec->target_si == -1)
                    enhance_external_m2_gauge(nhc, gr.IF, gr.shared_ext_idx, spec->target_si, cc);
            }
        }
        if (!nhc.passes) { rej_enhance++; return false; }

        {
            std::set<int> ext_set = {gr.shared_ext_idx};
            if (compute_proper_c_ext(gr.IF, nhc, ext_set) > g_proper_c_ext_max + 1e-9)
                { rej_cext++; return false; }
        }

        // Set LST-T base T for compute_anomaly (uses thread-local g_current_base_T).
        g_current_base_T = T_lst;
        AnomalyResult anom = compute_anomaly(gr.IF, nhc);
        if (anom.H_neutral < 0) { rej_anom++; return false; }

        bool cc_ok = true;
        int nn = gr.IF.rows();
        for (int k = 0; k < nn && cc_ok; k++) {
            if (gr.IF(k, k) != -1) continue;
            double total_cc = 0.0;
            for (int m = 0; m < nn; m++) {
                if (m == k || gr.IF(k, m) == 0) continue;
                GaugeInfo gm = (m < (int)nhc.curve_gauges.size()) ? nhc.curve_gauges[m] : GAUGE_NONE;
                if (gr.IF(m, m) == -1 && gm.dim == 0) continue;
                int intN = std::abs(gr.IF(k, m));
                GaugeInfo g = (m < (int)nhc.curve_gauges.size() && nhc.curve_gauges[m].dim > 0)
                    ? nhc.curve_gauges[m] : gauge_from_si(gr.IF(m, m));
                if (g.dim > 0) total_cc += central_charge(g, intN);
            }
            if (total_cc > cc + 1e-9) cc_ok = false;
        }
        if (!cc_ok) { rej_cc++; return false; }

        std::set<int> h1_set;
        if (spec->is_hat1) h1_set.insert(gr.shared_ext_idx);
        Eigen::MatrixXi ext_IF = Eigen::MatrixXi::Zero(nn+1, nn+1);
        ext_IF.block(0, 0, nn, nn) = gr.IF;
        for (int k = 0; k < nn; k++) {
            int b0c = h1_set.count(k) ? -1 : (gr.IF(k,k) + 2);
            ext_IF(k, nn) = b0c; ext_IF(nn, k) = b0c;
        }
        // Filter b_0^2 uses sig-derived T (= sig.sig_neg of the partial IF) so
        // ext_IF lands in the (1, T, 0) regime — this matches the bd_sig check.
        ext_IF(nn, nn) = 9 - T_glued;

        Eigen::MatrixXi bd = ext_IF;
        while (true) {
            int sz = bd.rows(), m1 = -1;
            for (int k = 0; k < sz; k++)
                if (bd(k,k) == -1 && bd(k, sz-1) == 1) { m1 = k; break; }
            if (m1 < 0) break;
            bd = blowdown_curve(bd, m1);
        }

        auto bd_sig = compute_sig_fast(bd);
        if (bd_sig.sig_pos != 1) { rej_bdsig++; return false; }

        Phase3Result r;
        r.spec_id = spec->id; r.spec_tag = result_tag;
        r.catA_id = A.entry.catalog_id; r.catA_type = A.entry.catalog_type;
        r.catB_id = B.entry.catalog_id; r.catB_type = B.entry.catalog_type;
        r.shared_ext_idx = gr.shared_ext_idx;
        r.b_start = gr.b_start; r.b_end = gr.b_end;
        r.ext_is_hat1 = spec->is_hat1;
        r.final_IF = gr.IF; r.nhc = nhc; r.anomaly = anom; r.sig = sig;
        r.T_min_target = sig.sig_neg;
        r.T_A = A.entry.T; r.T_B = B.entry.T;
        if (spec->is_hat1) r.hat1_indices.insert(gr.shared_ext_idx);
        results.push_back(std::move(r));
        passed++;
        return true;
    };

    if (pair_tag_mode) {
        auto itA = p1s.find(pair_tagA);
        auto itB = p1s.find(pair_tagB);
        if (itA == p1s.end() || itB == p1s.end()) {
            std::cerr << "Missing entries for tag(s): "
                      << (itA == p1s.end() ? pair_tagA : "") << " "
                      << (itB == p1s.end() ? pair_tagB : "") << "\n";
            return 1;
        }
        auto& entriesA = itA->second;
        auto& entriesB = itB->second;
        const ExternalSpec* specA = entriesA[0].spec;  // A drives enhancement
        std::string combo_tag = pair_tagA + "_x_" + pair_tagB;
        int spec_pass = 0;
        for (auto& A : entriesA) {
            for (auto& B : entriesB) {
                if (process_pair(A, B, specA, combo_tag)) spec_pass++;
            }
        }
        std::cout << "[" << combo_tag << "] " << entriesA.size() << " x "
                  << entriesB.size() << " pairs -> " << spec_pass << " pass\n";
    } else {
        for (auto& [tag, entries] : p1s) {
            const ExternalSpec* spec = entries[0].spec;
            int spec_pass = 0;
            int n = (int)entries.size();
            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    auto& A = entries[i];
                    auto& B = entries[j];
                    if (cross_mode && A.source == B.source) continue;
                    if (process_pair(A, B, spec, tag)) spec_pass++;
                }
            }
            std::cout << "[" << tag << "] " << n << " entries, "
                      << n*(n+1)/2 << " pairs -> " << spec_pass << " pass\n";
        }
    }
    std::cout << "\nTried: " << tried << ", passed: " << passed << "\n";
    std::cout << "Rejected: sig=" << rej_sig << " nhc=" << rej_nhc
              << " enhance=" << rej_enhance << " cext=" << rej_cext
              << " anomaly=" << rej_anom << " cc=" << rej_cc
              << " bdsig=" << rej_bdsig << "\n";

    // Dedup
    std::set<std::string>seen;std::vector<Phase3Result>filtered;
    for(auto&r:results){auto key=dedup_key_v2(r.final_IF,r.spec_tag);if(seen.insert(key).second)filtered.push_back(std::move(r));}
    std::cout<<"After dedup: "<<filtered.size()<<"\n";

    std::sort(filtered.begin(),filtered.end(),[](const Phase3Result&a,const Phase3Result&b){
        if(a.spec_tag!=b.spec_tag)return a.spec_tag<b.spec_tag;if(a.anomaly.T!=b.anomaly.T)return a.anomaly.T<b.anomaly.T;return std::abs(a.sig.det)<std::abs(b.sig.det);});

    std::vector<const Phase3Result*> all_pass;
    for (auto& r : filtered) all_pass.push_back(&r);

    std::map<std::string,int> sc;
    for (auto rp : all_pass) sc[rp->spec_tag]++;
    std::cout << "\nTotal valid: " << all_pass.size() << "\nPer-spec:\n";
    for (auto& [k,v] : sc) std::cout << "  " << k << ": " << v << "\n";

    // === .cat output ===
    std::string cat_dir = cross_mode ? "cat_glued_cross" : "cat_glued";
    std::filesystem::remove_all(cat_dir);
    std::filesystem::create_directories(cat_dir);

    std::map<std::string, std::vector<const Phase3Result*>> by_spec;
    for (auto rp : all_pass) by_spec[rp->spec_tag].push_back(rp);

    int total_cat = 0;
    for (auto& [tag, entries] : by_spec) {
        std::string fname = cat_dir + "/" + tag + ".cat";
        std::ofstream out(fname);
        out << "# Phase3 glued catalog: " << tag << "\n";
        out << "# Total: " << entries.size() << "\n#\n";
        int eid = 0;
        for (auto* rp : entries) {
            auto& r = *rp;
            int n = r.final_IF.rows();
            out << "ENTRY " << eid << "\n";
            out << "COMBO " << r.spec_tag << "\n";
            out << "EXTERNALS 1\n";
            out << r.shared_ext_idx << " " << r.spec_id << " " << r.spec_tag
                << " " << all_specs[r.spec_id].ext_si
                << " " << all_specs[r.spec_id].target_si
                << " " << all_specs[r.spec_id].int_num
                << " " << (r.ext_is_hat1 ? 1 : 0) << "\n";
            out << "BASE " << r.catA_id << " " << r.catA_type
                << "+" << r.catB_type << " " << r.T_A << " " << n << "\n";
            out << "GLUE " << r.catA_type << " " << r.T_A
                << " " << r.catB_type << " " << r.T_B << "\n";
            out << "PHYSICS " << r.anomaly.T << " " << r.anomaly.H_charged
                << " " << r.anomaly.V << " " << r.anomaly.H_neutral
                << " " << r.sig.det << " " << r.sig.sig_pos
                << " " << r.sig.sig_neg << " " << r.sig.sig_zero << "\n";
            out << "IF " << n << "\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (j) out << " ";
                    out << r.final_IF(i, j);
                }
                out << "\n";
            }
            // Blowdown matrix
            std::set<int> h1 = r.hat1_indices;
            Eigen::MatrixXi ext_m = Eigen::MatrixXi::Zero(n+1, n+1);
            ext_m.block(0,0,n,n) = r.final_IF;
            for (int i = 0; i < n; i++) {
                int b0c = h1.count(i) ? -1 : (r.final_IF(i,i) + 2);
                ext_m(i, n) = b0c; ext_m(n, i) = b0c;
            }
            ext_m(n, n) = 9 - r.anomaly.T;
            Eigen::MatrixXi bd = ext_m;
            while (true) {
                int sz = bd.rows(), m1 = -1;
                for (int k = 0; k < sz; k++)
                    if (bd(k,k) == -1 && bd(k, sz-1) == 1) { m1 = k; break; }
                if (m1 < 0) break;
                bd = blowdown_curve(bd, m1);
            }
            auto bds = compute_sig_fast(bd);
            int bdn = bd.rows();
            out << "BLOWDOWN " << bdn << "\n";
            for (int i = 0; i < bdn; i++) {
                for (int j = 0; j < bdn; j++) {
                    if (j) out << " ";
                    out << bd(i, j);
                }
                out << "\n";
            }
            out << "BDSIG " << bds.sig_pos << " " << bds.sig_neg
                << " " << bds.sig_zero << " " << bds.det << "\n";
            out << "END\n\n";
            eid++;
        }
        out.close();
        total_cat += entries.size();
    }

    // Index
    {
        std::ofstream idx(cat_dir + "/INDEX.txt");
        idx << "# Phase3 Glued Catalog\n# Total: " << total_cat << "\n#\n";
        for (auto& [tag, entries] : by_spec)
            idx << tag << " | " << entries.size() << " | " << tag << ".cat\n";
    }
    std::cout << "\nCatalog: " << total_cat << " entries → " << cat_dir << "/\n";

    // === LaTeX ===
    auto wsm=[](std::ofstream&t,const Eigen::MatrixXi&M){int s=M.rows();t<<"\\left(\\begin{smallmatrix}";
        for(int i=0;i<s;i++){if(i)t<<" \\\\ ";for(int j=0;j<s;j++){if(j)t<<" & ";t<<M(i,j);}}t<<"\\end{smallmatrix}\\right)";};
    std::string tf="sugra_phase3.tex";std::ofstream tex(tf);
    tex<<"\\documentclass[10pt]{article}\n\\usepackage[a4paper,margin=0.6in]{geometry}\n\\usepackage{longtable,booktabs,amsmath,amssymb,array}\n\\usepackage[dvipsnames]{xcolor}\n\\setlength{\\parindent}{0pt}\n\\setlength{\\parskip}{2pt}\n\\begin{document}\n\n";

    tex<<"\\section*{Phase 3: LST$_A$ --- External --- LST$_B$}\n\n"
       <<"Total valid: "<<all_pass.size()<<".\n\n";

    std::map<int,int> Tc;
    for (auto rp : all_pass) Tc[rp->anomaly.T]++;
    tex<<"\\begin{center}\n\\begin{tabular}{c|c}\n\\toprule\n$T$ & Count \\\\\n\\midrule\n";
    for(auto&[T,c]:Tc) tex<<T<<" & "<<c<<" \\\\\n";
    tex<<"\\bottomrule\n\\end{tabular}\n\\end{center}\n\n";

    auto write_entry = [&](std::ofstream& tex, const Phase3Result& r) {
        int delta = r.anomaly.H_charged - r.anomaly.V + 29*r.anomaly.T - 273;
        tex << "$" << if_to_latex_p3(r.final_IF, r.shared_ext_idx, r.ext_is_hat1, r.b_start, r.b_end) << "$";
        tex << " \\quad $T=" << r.anomaly.T
            << ",\\; T_A=" << r.T_A << ",\\; T_B=" << r.T_B
            << ",\\; \\Delta=" << delta
            << ",\\; H_n=" << r.anomaly.H_neutral
            << ",\\; \\det=" << r.sig.det
            << ",\\; (n_+,n_-,n_0)=(" << r.sig.sig_pos
            << "," << r.sig.sig_neg << "," << r.sig.sig_zero << ")$\n\n";

        std::set<int> h1 = r.hat1_indices;
        int n = r.final_IF.rows();
        Eigen::MatrixXi ext = Eigen::MatrixXi::Zero(n+1, n+1);
        ext.block(0,0,n,n) = r.final_IF;
        for (int i = 0; i < n; i++) {
            int b0C = h1.count(i) ? -1 : (r.final_IF(i,i) + 2);
            ext(i,n) = b0C; ext(n,i) = b0C;
        }
        ext(n,n) = 9 - r.anomaly.T;

        Eigen::MatrixXi cur = ext;
        int bd_count = 0;
        while (true) {
            int nn2 = cur.rows(), m1 = -1;
            for (int i = 0; i < nn2; i++)
                if (cur(i,i) == -1 && cur(i, nn2-1) == 1) { m1 = i; break; }
            if (m1 < 0) break;
            cur = blowdown_curve(cur, m1);
            bd_count++;
        }
        auto sig_bd = compute_sig_fast(cur);

        tex << "\\hspace{1em}$T=" << r.anomaly.T
            << ",\\; b_0^2=" << 9 - r.anomaly.T
            << ",\\; \\text{bd}=" << bd_count << "$\n\n";
        tex << "\\hspace{1em}blowdown: $";
        wsm(tex, cur);
        tex << "_{\\!b_0}$";
        tex << " \\quad sig$=(" << sig_bd.sig_pos << "," << sig_bd.sig_neg
            << "," << sig_bd.sig_zero << ")$, $\\det=" << sig_bd.det << "$\n\n";
        tex << "\\hspace{1em}{\\scriptsize " << r.catA_type
            << "$(T\\!=\\!" << r.T_A << ")$"
            << " $\\bowtie$ " << r.catB_type
            << "$(T\\!=\\!" << r.T_B << ")$"
            << "}\n\n";
    };

    std::string cs; int cT = -1;
    for (auto rp : all_pass) {
        auto& r = *rp;
        if (r.spec_tag != cs) {
            cs = r.spec_tag; cT = -1;
            tex << "\n\\subsection*{Shared external: " << all_specs[r.spec_id].label
                << " (" << sc[r.spec_tag] << " bases)}\n\n";
        }
        if (r.anomaly.T != cT) {
            cT = r.anomaly.T;
            tex << "\\paragraph{$T=" << cT << "$}\n\n";
        }
        write_entry(tex, r);
    }

    tex<<"\\end{document}\n";tex.close();std::cout<<"\nWritten "<<tf<<"\n";return 0;}
