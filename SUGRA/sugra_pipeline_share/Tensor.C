#include "Tensor.h"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>


std::ostream& operator<<(std::ostream& os, const Tensor& th)
{
    /* 원하는 출력 형식을 자유롭게 작성 */
    return os << th.GetIntersectionForm();
}
Tensor::Tensor() {

	Initialize();
}



void Tensor::Initialize() {
	intersection_form.resize(0,0);
	T = 0;
	b0_comp.clear();
}


Eigen::MatrixXi Tensor::GetIntersectionForm() const {

	return intersection_form;
}

double Tensor::GetDeterminant() const {

	return intersection_form.cast<double>().determinant();
}
int Tensor::GetT() const 
{
	return T;
}

Eigen::VectorXd Tensor::GetEigenvalues2() const {

	Eigen::MatrixXd Ad = intersection_form.cast<double>();
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ad);
	Eigen::VectorXd vec = es.eigenvalues();	
	Eigen::MatrixXi A   = intersection_form;	
	int nullity = 0;
	

	// DETERMINE THE NULLITY OF THE INTERSECTION FORM (i.e. number of null directions)
	for (int i = 0; i < T; i++) 
	{

		for (int j = i; j < T; j++ )
		{ 
			if ( j == i && A(j,i) !=0 ) {break;}
			else if ( A(j,i) == 0 ) 
			{ continue; }
			else if ( A(j,i) > 0 || A(j,i) < 0 )
			{ A.row(j).swap(A.row(i)); break;}

		
		}
		
		std::cout << A << std::endl;
		//subtract all the other leading non-zeros
		for (int k = i+1; k < T; k++ )
		{
			if ( A(k,i) > 0 || A(k,i) < 0 )
			{
				A.row(k) = A.row(k)*A(i,i)-A(k,i)*A.row(i);
				std::cout << A << std::endl;

			}
		}
	}

	std::cout << A<< std::endl;
	

	for (int i = 0; i < T; i++)
	{
		if ( A(i,i) == 0 )
		{
			nullity++;
		}
	}
//	std::cout << A << std::endl;
//	std::cout << nullity << std::endl;

	// at the moment, we get the number of null directions.
	// now, what we wanna do is just sorting the eigenvalues in increasing order 
	
	for (int i = 0; i < T; i++)
	{	
		for (int j = i+1; j < T; j++)
		{
			if ( abs(vec(j)) < abs((vec(i))) ) 
			{	
				std::swap(vec(i),vec(j));
				
			}
		}	
	}
	
	// replace the small errors by zeros
	
	for (int k = 0; k < nullity; k++)
	{
		vec(k) = 0;
	}

	return vec;
}
Eigen::VectorXd Tensor::GetEigenvalues() const {
    // ❶ 안전 체크
    if(intersection_form.rows() != intersection_form.cols()) 
        throw std::runtime_error("intersection_form is not square.");

    const int n = intersection_form.rows();

    // ❷ 수치 고윳값
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(intersection_form.cast<double>());
    if(es.info() != Eigen::Success) 
        throw std::runtime_error("Eigen decomposition failed.");

    Eigen::VectorXd vec = es.eigenvalues();

    // ❸ 정수 rank -> nullity
    Eigen::FullPivLU<Eigen::MatrixXd> lu(intersection_form.cast<double>());
    int nullity = n - lu.rank();

    // ❹ 정렬 & 0 덮어쓰기
    std::sort(vec.data(), vec.data()+n, 
              [](double a, double b){ return std::abs(a) < std::abs(b); });
    vec.head(nullity).setZero();

    return vec;
}
int Tensor::IsUnimodular() const
{	
	double det = 1;
	Eigen::VectorXd V = this->GetEigenvalues();

	for ( int i = 0; i < T; i++)
	{	
		if ( V(i) > 0 || V(i) < 0 )
		{
			det *= V(i);
		}
	}

          

	return std::llround(det); 
}
int Tensor::GetExactDet() const
{		
	
	return std::llround(this->GetDeterminant());


}

Eigen::VectorXi Tensor::GetSignature() const
{
	Eigen::VectorXd v = this->GetEigenvalues();
	for (int i = 0; i < T; i++)
	{	
		if ( v(i) > 0 ) { v(i) = 1; }
		if ( v(i) < 0 ) { v(i) = -1;}

	}

	for (int j = 0; j < T; j++)
	{
	
		if ( v(j) > 0 || v(j) < 0 )
		{
			for (int k = j+1; k < T; k++)
			{
				if (v(k) < v(j))
				{
					std::swap(v(k), v(j));
				}
			}
		}
	}

	return v.cast<int>();
}
			

void Tensor::AddTensorMultiplet(int charge)
{	
	T++;
	intersection_form.conservativeResize(T,T);
	
	for (int i = 0; i < T; i++)
	{
		intersection_form(T-1,i) = 0;
		intersection_form(i,T-1) = 0;
	}
	intersection_form(T-1,T-1) = charge;

}
void Tensor::AddT(int charge)
{	
	T++;
	intersection_form.conservativeResize(T,T);
	
	for (int i = 0; i < T; i++)
	{
		intersection_form(T-1,i) = 0;
		intersection_form(i,T-1) = 0;
	}
	intersection_form(T-1,T-1) = charge;

}


void Tensor::intersect(int n, int m, int k)
{
	intersection_form(n-1, m-1) = k;
	intersection_form(m-1, n-1) = k;

}
void Tensor::not_intersect(int n, int m)
{
	intersection_form(n-1,m-1) = 0;
	intersection_form(m-1,n-1) = 0;
}
void Tensor::DeleteTensorMultiplet()
{
	T--;
	intersection_form.conservativeResize(T,T);
}
bool Tensor::IsSUGRA() const
{
	bool c = 0;
	int n = std::llround(std::abs(this->IsUnimodular()));
	int sqrtn = std::llround(std::sqrt((long double)n));
	int timedir = 0;

	bool b = (sqrtn*sqrtn == n && n > 0);

	for (int i = 0; i < (this->GetSignature()).size();i++)
	{
		if ( (this->GetSignature())[i] == 1)
		{
			timedir ++;
		}
	}

	if ( timedir == 1)
	{
		c = 1;
	}
	else if (timedir != 1)
	{
		c = 0;
	}


	return b&&c;
}

void Tensor::CompleteBlowdown()
{
	for (int i = 0; i < T; i++)
	{
		if (intersection_form(i,i) == -1)
		{
			this->Blowdown5(i+1);
			i = -1;
			std::cout << intersection_form << std::endl;
			std::cout << this->GetSignature() << std::endl;
			if (T < 4)
			{
				break;
			}
		}
	}
}

void Tensor::LSTBlowdown(int ext)
{



	bool k;
	if ( ext == 0 )
	{
		//in this case, we are blowing down LST base only//

		this->ForcedBlowdown();
	}

	else if ( ext == 1)
	{
		//in this case, there exists one external curve.. 

		for (int i = 2; i <= T; i++)
		{
			if (intersection_form(0,i-1) == 0 && intersection_form(i-1,i-1) == -1)
			{
				k = this -> Blowdown5(i);
				i = 0;
			}
		}
	}
	else if ( ext == 2)
	{
		for (int i = 2; i <= T-1; i++)
		{
			if (intersection_form(0,i-1) == 0 && intersection_form(i-1,i-1) == -1)
			{
				k = this->Blowdown5(i);
				i = 1;
			}
		}
	}
	else if ( ext == 3)
	{
		for (int i = 2; i <=T-3; i++)
		{
			if (intersection_form(0,i-1) == 0 && intersection_form(i-1, T-1) == 0 && intersection_form(i-1,T-2) == 0 && intersection_form(i-1,i-1) == -1)
			{
				k = this->Blowdown5(i);
				i = 1;
			}
		}
	}
}
void Tensor::FBlowdown()
{

	bool c;
	for( int i = 0; i < T; i++)
	{
		if(intersection_form(i,i) == -1 )

		{
			c = this->Blowdown5(i+1);

			if (c == 1)
			{
				i = -1;
				//std::cout << " Blowdown " << std::endl;
				//std::cout << this->GetIFb0Q() << std::endl;
				//std::cout << this->GetSignature() << std::endl;
			}
			else if (c == 0)
			{
				//std::cout << "NOT SUGRA ANYMORE" << std::endl;

				break;
			}	
		}
		
	}
}
void Tensor::ForcedBlowdown()
{

	bool c;
	for( int i = 0; i < T; i++)
	{
		if(intersection_form(i,i) == -1 && b0_comp[i]==1)

		{
			c = this->Blowdown5(i+1);

			if (c == 1)
			{
				i = -1;
				//std::cout << " Blowdown " << std::endl;
				//std::cout << intersection_form << std::endl;
				//std::cout << this->GetSignature() << std::endl;
				if ( T == 1)
				{
					break;
				}
			}
			else if (c == 0)
			{
				//std::cout << "NOT SUGRA ANYMORE" << std::endl;

			}	
		}
		
	}
}


int Tensor::TimeDirection() const
{
	int timedir = 0;

	for (int i = 0; i < (this -> GetSignature()).size(); i++)
	{

		if ((this -> GetSignature())[i] == 1)
		{
			timedir++;
		}
	}

	return timedir;
}

int Tensor::NullDirection() const
{
	int nulldir = 0;

	for (int i = 0; i <(this->GetSignature()).size(); i++)
	{
		if ((this->GetSignature())[i] == 0)
		{
			nulldir++;
		}
	}
	return nulldir;
}
int Tensor::SpaceDirection() const
{
	int spacedir = 0;

	for (int i = 0; i <(this->GetSignature()).size(); i++)
	{
		if ((this->GetSignature())[i] == -1)
		{
			spacedir++;
		}
	}
	return spacedir;
}
void Tensor::SetElement(int n, int m, int k)
{
	intersection_form(n,m) = k;
}

void Tensor::Setb0Q()
{

	b0_comp.clear();

	for (int i = 0; i < T; i++)
	{
		b0_comp.push_back(intersection_form(i,i)+2);   // so, one must use this method to initial bases. DO NOT USE THIS METHOD in procedure of blowdown
	}

	int t = this->SpaceDirection();
	b0_comp.push_back(9-t);
	
}

std::vector<int> Tensor::Getb0Q()
{
	return b0_comp;
}

Eigen::MatrixXi Tensor::GetIFb0Q()
{
	Eigen::MatrixXi m = intersection_form;
	m.conservativeResize(T+1,T+1);

	for(int i =0; i<T; i++)
	{
		m(i,T) = b0_comp[i];
		m(T,i) = b0_comp[i];
	}

	m(T,T) = b0_comp[T];

	return m;
}



bool Tensor::Blowdown5(int n) 			//THIS METHOD IS FOR BLOWING DOWN b0Q COMPONENT 
{
	

	std::vector<int> v;
	std::vector<int> vcont;
	bool b = 1;

	if(intersection_form(n-1,n-1) != -1)
	{
		std::cout << "THIS CURVE CANNOT BE BLOWN DOWN\n" <<std::endl;

	}
	else
	{


		std::vector<int> vec;
		for ( int a = 0; a < T ;a++ )
		{
			if ( intersection_form(a,n-1) > 0 )
			{
				vec.push_back(a);
			}
		}	
		if ( vec.size() >= 2 )
		{
			/*	
			for (int a =0; a < vec.size(); a++)
			{
				if ( intersection_form(vec[a],vec[a]) >= 0)
				{
					b = 0;
					break;
				}
			}
			*/

			if ( b == 1 )
			{
				T--;
				Eigen::MatrixXi B = Eigen::MatrixXi::Zero(T,T);
				

				for (int i=0; i<T; i++)
				{
					for (int j=0; j<T; j++)
					{
						if ( i < n-1 && j < n-1 )
						{
							B(i,j) = intersection_form(i,j);
						}
						if ( i >= n-1 && j < n-1 )
						{
							B(i,j) = intersection_form(i+1,j);
						}
						if ( i < n-1 && j >= n-1)
						{	
							B(i,j) = intersection_form(i,j+1);
						}
						else if ( i >= n-1 && j >= n-1)
						{
							B(i,j) = intersection_form(i+1,j+1);
						}
					}
				}


				b0_comp[T+1]++;
				b0_comp.erase(b0_comp.begin() + n-1);
				
				for (int k=0; k < T+1; k++)
				{

					if( intersection_form(n-1,k) > 0 )
					{	
						if ( k < n-1 )
						{
							v.push_back(k);
							vcont.push_back(intersection_form(n-1,k));
							b0_comp[k]+=intersection_form(n-1,k);
							B(k,k) += intersection_form(n-1,k) * intersection_form(n-1,k);					

						}
						if ( k > n-1 )
						{
							v.push_back(k-1);
							vcont.push_back(intersection_form(n-1,k));
							b0_comp[k-1]+=intersection_form(n-1,k);
							B(k-1,k-1) += intersection_form(n-1,k) * intersection_form(n-1,k);
						}
					}
					
					
				}


				for (int a2 = 0; a2 < v.size() ;a2++)
				{
					for (int a3 = 0; a3 < a2; a3++)
					{
						B(v[a2],v[a3]) += vcont[a2]*vcont[a3];
						B(v[a3],v[a2]) += vcont[a2]*vcont[a3];
					}
				}	

				intersection_form = B;

			}
		}
		else if ( vec.size() == 1)
		{
			if ( !(intersection_form(vec[0],vec[0]) >= 0) )
			{
				T--;
				Eigen::MatrixXi B = Eigen::MatrixXi::Zero(T,T);
				for (int i=0; i<T; i++)
				{
					for (int j=0; j<T; j++)
					{
						if ( i < n-1 && j < n-1 )
						{
							B(i,j) = intersection_form(i,j);
						}
						if ( i >= n-1 && j < n-1 )
						{
							B(i,j) = intersection_form(i+1,j);
						}
						if ( i < n-1 && j >= n-1)
						{	
							B(i,j) = intersection_form(i,j+1);
						}
						else if ( i >= n-1 && j >= n-1)
						{
							B(i,j) = intersection_form(i+1,j+1);
						}
					}
				}

				
				b0_comp[T+1]++;

				b0_comp.erase(b0_comp.begin() + n-1);
				

				for (int k=0; k < T+1; k++)
				{


					if( intersection_form(n-1,k) > 0)
					{	
						if ( k < n-1 )
						{
							v.push_back(k);
							vcont.push_back(intersection_form(n-1,k));
							b0_comp[k]+=intersection_form(n-1,k);
							B(k,k) += intersection_form(n-1,k) * intersection_form(n-1,k);					

						}
						if ( k > n-1 )
						{
							v.push_back(k-1);
							vcont.push_back(intersection_form(n-1,k));
							b0_comp[k-1]+=intersection_form(n-1,k);
							B(k-1,k-1) +=intersection_form(n-1,k) * intersection_form(n-1,k);
						}
					}
				}


				if (v.size() > 1)
				{
					B(v[0],v[1])++;
					B(v[1],v[0])++;
				}	

				intersection_form = B;

			}
			else 
			{
				b = 0;
			}
		}
		else if (vec.size() == 0)
		{
			b = 0;
			std::cout << "No intersecting curves" << std::endl;
		}
		else if (b == 0)
		{
			std::cout << "Inconsistent base" << std::endl;
		}
				

		return b;

	}
}
bool Tensor::Blowdown6(int n) 			//THIS METHOD IS FOR BLOWING DOWN b0Q COMPONENT 
{
	

	std::vector<int> v;
	std::vector<int> vcont;
	bool b = 1;

	if(intersection_form(n-1,n-1) != -1)
	{
		std::cout << "THIS CURVE CANNOT BE BLOWN DOWN\n" <<std::endl;

	}
	else
	{


		std::vector<int> vec;
		for ( int a = 0; a < T ;a++ )
		{
			if ( intersection_form(a,n-1) > 0 )
			{
				vec.push_back(a);
			}
		}	
		if ( vec.size() >= 2 )
		{
			for (int a =0; a < vec.size(); a++)
			{
				if ( intersection_form(vec[a],vec[a]) >= 0)
				{
					b = 0;
					break;
				}
			}


			if ( b == 1 )
			{
				T--;
				Eigen::MatrixXi B = Eigen::MatrixXi::Zero(T,T);
				

				for (int i=0; i<T; i++)
				{
					for (int j=0; j<T; j++)
					{
						if ( i < n-1 && j < n-1 )
						{
							B(i,j) = intersection_form(i,j);
						}
						if ( i >= n-1 && j < n-1 )
						{
							B(i,j) = intersection_form(i+1,j);
						}
						if ( i < n-1 && j >= n-1)
						{	
							B(i,j) = intersection_form(i,j+1);
						}
						else if ( i >= n-1 && j >= n-1)
						{
							B(i,j) = intersection_form(i+1,j+1);
						}
					}
				}


				b0_comp[T+1]++;
				b0_comp.erase(b0_comp.begin() + n-1);
				
				for (int k=0; k < T+1; k++)
				{

					if( intersection_form(n-1,k) > 0 )
					{	
						if ( k < n-1 )
						{
							v.push_back(k);
							vcont.push_back(intersection_form(n-1,k));
							b0_comp[k]+=intersection_form(n-1,k);
							B(k,k) += intersection_form(n-1,k);					

						}
						if ( k > n-1 )
						{
							v.push_back(k-1);
							vcont.push_back(intersection_form(n-1,k));
							b0_comp[k-1]+=intersection_form(n-1,k);
							B(k-1,k-1) += intersection_form(n-1,k);
						}
					}
					
					
				}


				for (int a2 = 0; a2 < v.size() ;a2++)
				{
					for (int a3 = 0; a3 < a2; a3++)
					{
						B(v[a2],v[a3]) = vcont[a2]*vcont[a3];
						B(v[a3],v[a2]) = vcont[a2]*vcont[a3];
					}
				}	

				intersection_form = B;

			}
		}
		else if ( vec.size() == 1)
		{
			if ( !(intersection_form(vec[0],vec[0]) >= 0) )
			{
				T--;
				Eigen::MatrixXi B = Eigen::MatrixXi::Zero(T,T);
				for (int i=0; i<T; i++)
				{
					for (int j=0; j<T; j++)
					{
						if ( i < n-1 && j < n-1 )
						{
							B(i,j) = intersection_form(i,j);
						}
						if ( i >= n-1 && j < n-1 )
						{
							B(i,j) = intersection_form(i+1,j);
						}
						if ( i < n-1 && j >= n-1)
						{	
							B(i,j) = intersection_form(i,j+1);
						}
						else if ( i >= n-1 && j >= n-1)
						{
							B(i,j) = intersection_form(i+1,j+1);
						}
					}
				}

				
				b0_comp[T+1]++;

				b0_comp.erase(b0_comp.begin() + n-1);
				

				for (int k=0; k < T+1; k++)
				{


					if( intersection_form(n-1,k) > 0)
					{	
						if ( k < n-1 )
						{
							v.push_back(k);
							vcont.push_back(intersection_form(n-1,k));
							b0_comp[k]+=intersection_form(n-1,k);
							B(k,k) += intersection_form(n-1,k);					

						}
						if ( k > n-1 )
						{
							v.push_back(k-1);
							vcont.push_back(intersection_form(n-1,k));
							b0_comp[k-1]+=intersection_form(n-1,k);
							B(k-1,k-1) +=intersection_form(n-1,k);
						}
					}
				}


				if (v.size() > 1)
				{
					B(v[0],v[1])++;
					B(v[1],v[0])++;
				}	

				intersection_form = B;

			}
			else 
			{
				b = 0;
			}
		}
		else if (vec.size() == 0)
		{
			b = 0;
			std::cout << "No intersecting curves" << std::endl;
		}
		else if (b == 0)
		{
			std::cout << "Inconsistent base" << std::endl;
		}
				

		return b;

	}
}


void Tensor::AL(int n, int m, bool b)
{

	///////////////////////////////////////////////////////////////////////////////////
	// Table for nodes and links 
	// Nodes (curve, gauge algebra)
	// 	 (-12, e_8)
	// 	 (-11, e_8')
	// 	 (-10, e_8'')
	// 	 ( -9, e_8''')
	// 	 ( -8, e_7)
	// 	 ( -7, e_7')
	// 	 ( -6, e_6)
	// 	 ( -5, f_4)
	// 	 ( -4, so_8)
	// 	 ( -3, su_2)	
	//
	// Links (n,m) {tensormultiplets}
	// 	 (1,1) {-1}
	// 	 (2,2) {-1,-3,-1}
	// 	 (3,3) {-1,-2,-3,-2,-1}
	// 	 (3,2) {-1,-2,-3,-1}
	// 	 (4,2) {-1,-2,-2,-3,-1}
	// 	 (3,3) {-1,-3,-1,-5,-1,-3,-1} void
	// 	 (3,4) {-1,-3,-1,-5,-1,-3,-2,-1} 
	// 	 (3,5) {-1,-3,-1,-5,-1,-3,-2,-2,-1}
	// 	 (4,4) {-1,-2,-3,-1,-5,-1,-3,-2,-1}
	// 	 (4,5) {-1,-2,-3,-1,-5,-1,-3,-2,-2,-1}
	// 	 (5,5) {-1,-2,-2,-3,-1,-5,-1,-3,-2,-2,-1}
	//
	//
	////////////////////////////////////////////////////////////////////////////////////

	if( T > 0)
	{
		if(b==0)
		{

			if (n==1 && m==1)
			{
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==2 && m==2)
			{
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==3 && m==3)
			{
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==4 && m==4)
			{
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==5 && m==5)
			{
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);


			}
			if (n==3 && m==2)
			{
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==2 && m==3)
			{
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==4 && m==2)
			{
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==2 && m==4)
			{

				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==3 && m==4)
			{

				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==4 && m==3)
			{

				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==3 && m==5)
			{

				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==5 && m==3)
			{

				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}

			if (n==4 && m==5)
			{

				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==5 && m==4)
			{

				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}

		}
		else if ( b==1)
		{
			if( n==3 && m==3)
			{

				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
		}
		else
		{
			std::cout << "\n///////////////    NO SUCH LINK EXISTS  ////////////////\n" <<std::endl;
		}
	}
	else if(T==0)
	{

		if(b==0)
		{

			if (n==1 && m==1)
			{
				AddTensorMultiplet(-1);
			}
			if (n==2 && m==2)
			{
				AddTensorMultiplet(-1);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==3 && m==3)
			{
				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==4 && m==4)
			{
				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==5 && m==5)
			{
				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);


			}
			if (n==3 && m==2)
			{
				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==2 && m==3)
			{
				AddTensorMultiplet(-1);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);

			}
			if (n==4 && m==2)
			{
				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==2 && m==4)
			{

				AddTensorMultiplet(-1);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==3 && m==4)
			{

				AddTensorMultiplet(-1);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==4 && m==3)
			{

				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==3 && m==5)
			{

				AddTensorMultiplet(-1);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==5 && m==3)
			{

				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}

			if (n==4 && m==5)
			{

				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
			if (n==5 && m==4)
			{

				AddTensorMultiplet(-1);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-2);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}

		}
		else if ( b==1)
		{
			if( n==3 && m==3)
			{

				AddTensorMultiplet(-1);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-5);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
				AddTensorMultiplet(-3);
				intersect(T-1,T);
				AddTensorMultiplet(-1);
				intersect(T-1,T);
			}
		}
		else
		{
			std::cout << "\n///////////////    NO SUCH LINK EXISTS  ////////////////\n" <<std::endl;
		}

	}

}

void Tensor::AT(int n)
{
	if (T==0)
	{
		AddT(n);
	}
	else 
	{
		AddT(n);
		intersect(T-1,T);
	}
}

void Tensor::ATS(int n, int m)
{
	if (T==0)
	{
		AddT(n);
		AddT(m);
		intersect(T-1,T);
	}
	else
	{
		AddT(n);
		AddT(m);
		intersect(T-1,T);
		intersect(T-2,T);
	}
}
void Tensor::ATS2(int n, int m, int l)
{
	if (T==0)
	{
		AddT(n);
		AddT(m);
		AddT(l);
		intersect(T-1,T);
		intersect(T-2,T-1);
	}
	else
	{
		AddT(n);
		AddT(m);
		AddT(l);
		intersect(T-1,T);
		intersect(T-2,T-1);
		intersect(T-3,T);
	}
}
void Tensor::ATS3(int n, int m, int l, int k)
{
	if (T==0)
	{
		AddT(n);
		AddT(m);
		AddT(l);
		AddT(k);
		intersect(T-1,T);
		intersect(T-2,T-1);
		intersect(T-3,T-2);
	}
	else
	{
		AddT(n);
		AddT(m);
		AddT(l);
		AddT(k);
		intersect(T-1,T);
		intersect(T-2,T-1);
		intersect(T-3,T-2);
		intersect(T-4,T);
	}
}


void Tensor::ATE(int n, int m, int l)
{
	AddT(n);
	intersect(m,T,l);
}

void Tensor::ALSTE(int m, int l)
{
	AddT(-1);
	AddT(-1);
	intersect(T-1,T);
	intersect(m,T,l);
}

void Tensor::SetIF(Eigen::MatrixXi M)
{
	this -> Initialize();
	int size = M.rows();

	intersection_form = Eigen::MatrixXi::Zero(size,size);
	
	intersection_form = M;

	T = size;
}






		


			

	


		
