#pragma once
#include <Eigen/Dense>
#include <vector>

class Tensor {
	private:
    		//string gauge_alg;				// types of gauge algebra 
								// -> data is needed..?
    		Eigen::MatrixXi  intersection_form;          // intersection form
		int T;
		std::vector<int>  b0_comp;

	public:
    		Tensor();                      
    		~Tensor() = default;
		void Initialize();

    		/* -------- queries -------- */
		Eigen::MatrixXi GetIntersectionForm() const;
    		//string getAnomaly()          const { return anomaly; }
   		double GetDeterminant() const;
	   	int GetExactDet() const;	
		Eigen::VectorXd	GetEigenvalues() const;
		Eigen::VectorXd GetEigenvalues2() const;
	   	int IsUnimodular() const;
		Eigen::VectorXi GetSignature() const;
		int GetT() const;
		int TimeDirection() const;
		int NullDirection() const;
		int SpaceDirection() const; 
		bool IsSUGRA() const;	

    		/* -------- modifiers -------*/
    		void AddTensorMultiplet(int charge);
		void AddT(int charge);		// anomaly += charge
    		//void AddTensorMultiplet(int charge, int anomalyPart); // anomaly += anomalyPart
    		void intersect(int n, int m, int k=1);
		void not_intersect(int n, int m);
		void DeleteTensorMultiplet();
		bool Blowdown5(int n);
		bool Blowdown6(int n);
		void CompleteBlowdown();
		void LSTBlowdown(int ext);
		void FBlowdown();
		void ForcedBlowdown();
		void SetElement(int n, int m, int k);
		void Setb0Q();
		std::vector<int> Getb0Q();
		Eigen::MatrixXi GetIFb0Q();
		// Methods for adding link, node, side link, extra tensor and minimal lst
		void AL(int n, int m, bool b=0);
		void AT(int n);
		void ATS(int n, int m);
		void ATS2(int n, int m, int l);
		void ATS3(int n, int m, int l, int k);	
		void ATE(int n, int m, int l=1);
		void ALSTE(int m, int l=1);
		void SetIF(Eigen::MatrixXi M);
			

	friend std::ostream& operator<<(std::ostream& os, const Tensor& th);
};

