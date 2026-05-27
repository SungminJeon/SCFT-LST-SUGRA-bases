// Theory_enhanced.h
// 일반화된 porting rule 시스템 - 곡선 인덱스 기반 attachment
// ✨ ENHANCED: TheoryGraph 클래스 추가 with AttachmentPoint support
#pragma once
#include "Tensor.h"
#include <vector>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <map>
#include <sstream>
#include <unordered_map>
#include <set>

// ---- 종류 & 스펙 헬퍼 ----
enum class Kind { SideLink, InteriorLink, Node, External };
struct Spec { Kind kind; int param; };
inline Spec s(int p){ return {Kind::SideLink,     p}; }
inline Spec i(int p){ return {Kind::InteriorLink, p}; }
inline Spec n(int p){ return {Kind::Node,        p}; }
inline Spec e(int p){ return {Kind::External,     p}; }

// ---- 0x0 행렬 안전 출력 ----
template <class M>
inline void PrintMatrixSafe(const M& A, std::ostream& os = std::cout){
    if (A.rows()==0 || A.cols()==0) { os << "[empty 0x0 matrix]\n"; return; }
    os << A << '\n';
}

// ---- Attachment Point 정의 ----
// -1 = Left end (첫 번째 곡선)
// -2 = Right end (마지막 곡선)
// 0, 1, 2, ... = 특정 곡선 인덱스 (0-based)
struct AttachmentPoint {
    int index;  // -1 (Left), -2 (Right), or 0,1,2,... (curve index)
    
    AttachmentPoint() : index(-1) {}
    AttachmentPoint(int idx) : index(idx) {}
    
    bool isLeft() const { return index == -1; }
    bool isRight() const { return index == -2; }
    bool isSpecificCurve() const { return index >= 0; }
    
    // 실제 곡선 인덱스로 변환 (총 곡선 개수 필요)
    int toAbsoluteIndex(int totalCurves) const {
        if (index == -1) return 0;  // Left = first curve
        if (index == -2) return totalCurves - 1;  // Right = last curve
        if (index < 0 || index >= totalCurves) return 0;  // Safety
        return index;  // specific index
    }
};

// ---- Backward compatibility를 위한 Port enum ----
enum class Port { Left, Right, Custom };

inline AttachmentPoint portToAttachment(Port p) {
    return (p == Port::Left) ? AttachmentPoint(-1) : AttachmentPoint(-2);
}

// Port 선택 함수 - AttachmentPoint 기반
inline int pickPortIndex(Kind /*k*/, const Tensor& t, const AttachmentPoint& ap){
    const auto IF = t.GetIntersectionForm();
    const int sz = static_cast<int>(IF.rows());
    if (sz <= 0) return -1;
    return ap.toAbsoluteIndex(sz);
}

// 기존 Port enum 기반 호환 함수
inline int pickPortIndex(Kind /*k*/, const Tensor& t, Port which){
    const auto IF = t.GetIntersectionForm();
    const int sz = static_cast<int>(IF.rows());
    if (sz <= 0) return -1;
    switch(which){
        case Port::Left:   return 0;
        case Port::Right:  return sz - 1;
        case Port::Custom: return (sz >= 2 ? 1 : 0);
    }
    return -1;
}

// ---- Spec -> Tensor (Theory_enhanced 기반의 확장된 build_tensor) ----
inline Tensor build_tensor(const Spec& sp){
    Tensor t;

    switch (sp.kind){
        case Kind::SideLink:
            // s(12) 등 파라미터 내용은 현재 무시: 곡선 2개만 만든다.
            // instantons : notation 88(blowdown induced) 
            if (sp.param == 1) {t.AT(-1);}
			else if (sp.param == 882 ) {t.AT(-2); t.AT(-1);}
			else if (sp.param == 883 ) {t.AT(-2); t.AT(-2); t.AT(-1);}	
			else if (sp.param == 884 ) {t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1);}	
			else if (sp.param == 885 ) {t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1);}	
			else if (sp.param == 886 ) {t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1);}	
			else if (sp.param == 887 ) {t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1);}	
			else if (sp.param == 8881 ) {t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 889 ) {t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 8810 ) {t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 8811 ) {t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-1);}
		

			else if (sp.param == 288 ) {t.AT(-1); t.AT(-2);}
			else if (sp.param == 388 ) {t.AT(-1); t.AT(-2); t.AT(-2);}	
			else if (sp.param == 488 ) {t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-2);}	
			else if (sp.param == 588 ) {t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2);}	
			else if (sp.param == 688 ) {t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2);}	
			else if (sp.param == 788 ) {t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2);}	
			else if (sp.param == 1888 ) {t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2);}
			else if (sp.param == 988 ) {t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2);}
			else if (sp.param == 1088 ) {t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2);}
			else if (sp.param == 1188 ) {t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2); t.AT(-2);}
			

			// interiors
			else if (sp.param == 11 ) { t.AL(1,1); }
			else if (sp.param == 22 ) { t.AL(2,2); }
			else if (sp.param == 33 ) { t.AL(3,3); }
			else if (sp.param == 44 ) { t.AL(4,4); }
			else if (sp.param == 55 ) { t.AL(5,5); }
			else if (sp.param == 331 ) { t.AL(3,3,1); }
			else if (sp.param == 32 ) { t.AL(3,2); }
			else if (sp.param == 23 ) { t.AL(2,3); }
			else if (sp.param == 42 ) { t.AL(4,2); }
			else if (sp.param == 24 ) { t.AL(2,4); }
			else if (sp.param == 43 ) { t.AL(4,3); }
			else if (sp.param == 34 ) { t.AL(3,4); }
			else if (sp.param == 53 ) { t.AL(5,3); }
			else if (sp.param == 35 ) { t.AL(3,5); }
			else if (sp.param == 54 ) { t.AL(5,4); }
			else if (sp.param == 45 ) { t.AL(4,5); }

			// alkali 2 links with no -5 

			else if (sp.param == 991 ) { t.AT(-2); t.ATS(-1,-3); t.AT(-1); }
			else if (sp.param == 9920 ) { t.AT(-1); t.AT(-2); t.ATS(-2,-3); t.AT(-1); }
			else if (sp.param == 9902 ) { t.AT(-1); t.ATS(-2,-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 993 ) { t.AT(-2); t.ATS(-1,-3); t.AT(-2); t.AT(-1); }
			
			// alkali 1 links with no -5

			else if (sp.param == 91) {t.AT(-3); t.ATS(-2,-2); t.AT(-1); }
			else if (sp.param == 92) {t.AT(-2); t.ATS(-2,-3); t.AT(-1); }
			else if (sp.param == 93) {t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 94) {t.AT(-2);  t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1);}  
			else if (sp.param == 95) {t.AT(-2); t.AT(-2);  t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 96) {t.AT(-3);  t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 97) {t.AT(-3);  t.AT(-2); t.AT(-1); }
			else if (sp.param == 98) {t.AT(-2);  t.AT(-3); t.AT(-2); t.AT(-1);} 
			else if (sp.param == 99) {t.AT(-2);  t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1);}
			else if (sp.param == 910) {t.AT(-2); t.AT(-2);  t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1);}
			else if (sp.param == 911) {t.AT(-3);  t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 912) {t.AT(-3);  t.AT(-1);}
			else if (sp.param == 913) {t.AT(-2);  t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1);}
			else if (sp.param == 914) {t.AT(-2); t.AT(-2);  t.AT(-3); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 915) {t.AT(-3);  t.AT(-1); t.AT(-3); t.AT(-1);}
			else if (sp.param == 916) {t.AT(-2);  t.AT(-3); t.AT(-1);}
			else if (sp.param == 917) {t.AT(-2);  t.AT(-2); t.AT(-3); t.AT(-1);}

			// alkali 3 links with one -5 curve
			
			else if (sp.param == 99910) {t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-1);}
			else if (sp.param == 99901) {t.AT(-1);  t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1);}
			else if (sp.param == 99920) {t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1);}
			else if (sp.param == 99902) {t.AT(-1);  t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1);}
			else if (sp.param == 99930) {t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 99903) {t.AT(-1);  t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1);}

			// alkali 2 links with one -5 curve
			
			
			else if (sp.param == 994) {t.AT(-3); t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-1);}
			else if (sp.param == 995) {t.AT(-3); t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1);}
			else if (sp.param == 996) {t.AT(-3); t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 997) {t.AT(-2); t.AT(-3); t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1);}
			else if (sp.param == 998) {t.AT(-2); t.AT(-3); t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 999) {t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1);  t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1);}
			else if (sp.param == 9910) {t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-1);}
			else if (sp.param == 9911) {t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.ATS(-1,-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1);}
			else if (sp.param == 9912) {t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 9913) {t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 9914) {t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); }


			// alkali 1 links with one -5 curve
	
			else if (sp.param == 918) {t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 919) {t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 920) {t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 921) {t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 922) {t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 923) {t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 924) {t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 925) {t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 926) {t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 927) {t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 928) {t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 929) {t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 930) {t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 931) {t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 932) {t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 933) {t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 934) {t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 935) {t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 936) {t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 937) {t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 938) {t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 939) {t.AT(-2); t.AT(-3); t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 940) {t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); }
			else if (sp.param == 941) {t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); }
			else if (sp.param == 942) {t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); }
			else if (sp.param == 943) {t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 944) {t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 945) {t.AT(-2); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }


			// alkali 2 links with two -5 curves
			
			else if (sp.param == 9915) {t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 9916) {t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 9917) {t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }


			// alkali 1 links with two -5 curves

			else if (sp.param == 946) {t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 947) {t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 948) {t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 949) {t.AT(-2); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 950) {t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-2); t.AT(-1); }
			else if (sp.param == 951) {t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 952) {t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 953) {t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 954) {t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-2); t.AT(-1); }
			else if (sp.param == 955) {t.AT(-5); t.AT(-1); t.AT(-2); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 956) {t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }
			else if (sp.param == 957) {t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1); }


			// alkali 1 link with one -5 curve(which is omitted in the appendix of atomic classification paper)

	    		else if (sp.param == 958) {t.AT(-1); t.AT(-5); t.AT(-1); t.AT(-3); t.AT(-1);}



			break;

        case Kind::Node:
            // n(1..12): 각각 g, g-L, ..., g-L^11로 최대 12개 곡선
            t.AT(-sp.param);
            break;

        case Kind::InteriorLink:
            // i(11) = AL(1,1), i(22) = AL(2,2), ...
            // 단순화: param을 적절히 처리
            {
                std::string s = std::to_string(sp.param);
                if (s.length() == 2) {
                    int a = s[0] - '0';
                    int b = s[1] - '0';
                    t.AL(a, b);
                } else if (s.length() == 3) {
                    int a = s[0] - '0';
                    int b = s[1] - '0';
                    int c = s[2] - '0';
                    t.AL(a, b, c);
                }
            }
            break;

        case Kind::External:
            // External: -1 곡선 1개
            t.AT(-sp.param);
            break;
    }
    return t;
}

// ---- getCurveCount helper function ----
inline int getCurveCount(int param, Kind kind) {
    Spec sp;
    sp.kind = kind;
    sp.param = param;
    Tensor t = build_tensor(sp);
    return t.GetIntersectionForm().rows();
}

// ===================== Enhanced TheoryGraph with AttachmentPoint Support =====================

struct EdgeW {
    int  u, v;                    // 노드 id
    AttachmentPoint pu, pv;       // ✨ ENHANCED: AttachmentPoint instead of Port
    int  w;                       // 글루잉 강도(대칭 오프대각에 +w)
};

struct NodeRef { int id; };

class TheoryGraph {
public:
    NodeRef add(Spec sp){
        int id = (int)nodes_.size();
        nodes_.push_back(build_tensor(sp));
        kinds_.push_back(sp.kind);
        params_.push_back(sp.param);
        return NodeRef{id};
    }

    // ✨ ENHANCED: 기본 connect - 자동으로 Right->Left attachment
    void connect(NodeRef a, NodeRef b){
        if (forbidden_(kinds_[a.id], kinds_[b.id]))
            throw std::invalid_argument("Forbidden adjacency s-i");

        // 기본: Right->Left attachment
        AttachmentPoint pa(-2);  // Right
        AttachmentPoint pb(-1);  // Left

        // Porting rule 체크는 기존과 동일하게 작동
        // (AttachmentPoint가 Left/Right인 경우 Port로 변환 가능)
        
        edgesW_.push_back(EdgeW{a.id, b.id, pa, pb, 1});
    }

    // ✨ ENHANCED: AttachmentPoint 기반 연결
    void connect(NodeRef a, const AttachmentPoint& pa, NodeRef b, const AttachmentPoint& pb, int weight=1){
        if (forbidden_(kinds_[a.id], kinds_[b.id]))
            throw std::invalid_argument("Forbidden adjacency s-i");
        if (weight <= 0) throw std::invalid_argument("weight must be positive");

        edgesW_.push_back(EdgeW{a.id, b.id, pa, pb, weight});
    }

    // ✨ BACKWARD COMPATIBILITY: Port 기반 연결도 지원
    void connect(NodeRef a, Port pa, NodeRef b, Port pb, int weight=1){
        connect(a, portToAttachment(pa), b, portToAttachment(pb), weight);
    }

    void print() const {
        std::cout << "TheoryGraph:\n";
        for (auto& e : edgesW_){
            std::cout << "  " << e.u << "(" << e.pu.index << ") --(" << e.w
                      << ")-- " << e.v << "(" << e.pv.index << ")\n";
        }
    }

    auto IF(int node) const { return nodes_.at(node).GetIntersectionForm(); }
    void PrintIF(int node, std::ostream& os = std::cout) const {
        os << "IF[node " << node << "]:\n";
        PrintMatrixSafe(IF(node), os);
    }

    // ✨ ENHANCED: AttachmentPoint를 고려한 ComposeIF_Gluing
    Eigen::MatrixXi ComposeIF_Gluing() const {
        const int N = (int)nodes_.size();
        if (N==0) return Eigen::MatrixXi();

        // 1) 블록 대각합
        std::vector<Eigen::MatrixXi> blocks; blocks.reserve(N);
        std::vector<int> sz; sz.reserve(N);
        for (auto& t : nodes_) { 
            auto M=t.GetIntersectionForm(); 
            blocks.push_back(M); 
            sz.push_back(M.rows()); 
        }
        Eigen::MatrixXi G = BlockDiag_(blocks);

        // 2) prefix offsets
        std::vector<int> off(N+1,0);
        for (int i=0;i<N;++i) off[i+1]=off[i]+sz[i];

        // 3) ✨ ENHANCED: AttachmentPoint를 사용한 간선 처리
        for (const auto& e : edgesW_){
            int iu = pickPortIndex(kinds_[e.u], nodes_[e.u], e.pu);
            int iv = pickPortIndex(kinds_[e.v], nodes_[e.v], e.pv);
            if (iu<0 || iv<0 || iu>=sz[e.u] || iv>=sz[e.v]) continue; // 방어
            int I = off[e.u] + iu;
            int J = off[e.v] + iv;
            G(I,J) += e.w;
            G(J,I) += e.w;
        }
        return G;
    }

    // 호환: 예전 이름 유지
    Eigen::MatrixXi ComposeIF_UnitGluing() const {
        return ComposeIF_Gluing();
    }

    int nodeCount() const { return (int)nodes_.size(); }

private:
    std::vector<Tensor> nodes_;
    std::vector<Kind> kinds_;
    std::vector<int> params_;
    std::vector<EdgeW> edgesW_;

    static bool forbidden_(Kind a, Kind b){
        // s-i can't be glued together
        return ( (a==Kind::SideLink && b==Kind::InteriorLink) ||
                 (a==Kind::InteriorLink && b==Kind::SideLink) );
    }

    static Eigen::MatrixXi BlockDiag_(const std::vector<Eigen::MatrixXi>& blocks){
        if (blocks.empty()) return Eigen::MatrixXi(); // 0x0
        Eigen::Index total = 0;
        for (auto& M : blocks) total += M.rows();
        if (total==0) return Eigen::MatrixXi();
        Eigen::MatrixXi out = Eigen::MatrixXi::Zero(total, total);
        Eigen::Index off = 0;
        for (auto& M : blocks){
            if (M.size()==0) continue;
            out.block(off, off, M.rows(), M.cols()) = M;
            off += M.rows();
        }
        return out.block(0,0,off,off);
    }
};

// ===================== Attachment Rules (from Theory_enhanced.h) =====================

class AttachmentRules {
private:
    // ---- 기존 porting 규칙들을 AttachmentPoint 시스템에 맞게 확장 ----
    
    // SideLink-Node attachment 규칙 (AttachmentPoint 기반)
    static bool isAllowedPortingAttachmentSN(
        int sideParam,
        int nodeParam,
        const AttachmentPoint& sideAttach,
        const AttachmentPoint& nodeAttach
    ) {
        int sideCurveCount = getCurveCount(sideParam, Kind::SideLink);
        int nodeCurveCount = getCurveCount(nodeParam, Kind::Node);
        
        int sideIdx = sideAttach.toAbsoluteIndex(sideCurveCount);
        int nodeIdx = nodeAttach.toAbsoluteIndex(nodeCurveCount);
        
        bool isSideLeft = (sideIdx == 0);
        bool isSideRight = (sideIdx == sideCurveCount - 1);
        bool isNodeLeft = (nodeIdx == 0);
        bool isNodeRight = (nodeIdx == nodeCurveCount - 1);
        
        // 288-1188: Right to Left만 허용
        if ((sideParam == 288 || sideParam == 388 || sideParam == 488 || 
             sideParam == 588 || sideParam == 688 || sideParam == 788 ||
             sideParam == 1888 || sideParam == 988 || sideParam == 1088 || 
             sideParam == 1188) && nodeParam > 0) {
            if (isSideRight && isNodeLeft) return true;
        }
        
        // 882-8811: Left to Right만 허용
        if ((sideParam == 882 || sideParam == 883 || sideParam == 884 || 
             sideParam == 885 || sideParam == 886 || sideParam == 887 ||
             sideParam == 8881 || sideParam == 889 || sideParam == 8810 || 
             sideParam == 8811) && nodeParam > 0) {
            if (isSideLeft && isNodeRight) return true;
        }
        
        // Interior link as side 규칙들
        if (sideParam == 11) return true;
        if (sideParam == 22 || sideParam == 33 || sideParam == 44 || 
            sideParam == 55 || sideParam == 331) {
            return true;
        }
        
        // 비대칭 interior links
        if (sideParam == 32) {
            if (isSideRight && isNodeLeft && nodeParam > 4) return false;
            if (isSideLeft && isNodeRight && nodeParam > 8) return false;
            return true;
        }
        
        if (sideParam == 23) {
            if (isSideLeft && isNodeRight && nodeParam > 4) return false;
            if (isSideRight && isNodeLeft && nodeParam > 8) return false;
            return true;
        }
        
        return true;
    }
    
    // InteriorLink-Node attachment 규칙 (AttachmentPoint 기반)
    static bool checkSpecificAttachmentRulesIN(
        int interiorParam,
        int nodeParam,
        const AttachmentPoint& interiorAttach,
        const AttachmentPoint& nodeAttach
    ) {
        int interiorCurveCount = getCurveCount(interiorParam, Kind::InteriorLink);
        int nodeCurveCount = getCurveCount(nodeParam, Kind::Node);
        
        int interiorIdx = interiorAttach.toAbsoluteIndex(interiorCurveCount);
        int nodeIdx = nodeAttach.toAbsoluteIndex(nodeCurveCount);
        
        bool isInteriorLeft = (interiorIdx == 0);
        bool isInteriorRight = (interiorIdx == interiorCurveCount - 1);
        bool isNodeLeft = (nodeIdx == 0);
        bool isNodeRight = (nodeIdx == nodeCurveCount - 1);
        
        // 기존 규칙들
        if (interiorParam == 32 && nodeParam > 4) {
            if (isInteriorRight && isNodeLeft) return false;
        }
        if (interiorParam == 32 && nodeParam > 8) {
            if (isInteriorLeft && isNodeRight) return false;
        }
        
        if (interiorParam == 23 && nodeParam > 4) {
            if (isInteriorLeft && isNodeRight) return false;
        }
        if (interiorParam == 23 && nodeParam > 8) {
            if (isInteriorRight && isNodeLeft) return false;
        }
        
        // ... (더 많은 규칙들 - Theory_enhanced.h에서 가져옴)
        
        return true;
    }
    
    // 전체 금지 규칙 (포트 무관)
    static bool isBannedPairIN(int interiorParam, int nodeParam) {
        if (interiorParam == 11 && nodeParam > 4) return true;
        if (interiorParam == 22 && nodeParam > 6) return true;
        if (interiorParam == 33 && nodeParam > 8) return true;
        if (interiorParam == 331 && nodeParam > 6) return true;
        if (interiorParam == 44 && (nodeParam > 8 || nodeParam < 6)) return true;
        if (interiorParam == 55 && nodeParam < 6) return true;
        return false;
    }
    
public:
    // ✨ PUBLIC: AttachmentPoint 기반 최종 체크 함수
    static bool isValidAttachment(
        int interiorParam,
        int nodeParam,
        const AttachmentPoint& interiorAttach,
        const AttachmentPoint& nodeAttach
    ) {
        if (isBannedPairIN(interiorParam, nodeParam)) return false;
        return checkSpecificAttachmentRulesIN(interiorParam, nodeParam, interiorAttach, nodeAttach);
    }
};

