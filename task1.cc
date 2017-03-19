#include <iostream>
#include <iomanip>
#include <vector>
#include <valarray>
 #include <sstream>
 #include <cmath>
 #include <fstream>
 #include <string>

using namespace std;
using dvec=vector<double>;
using matrix = vector<dvec>;
using counter = unsigned int;

ostream& operator<< (ostream &os, dvec& x)
{
	stringstream s;
	// s<<'[';
	for(auto i : x) 
	{
		s<<i<<"  ";
	}
	string str = s.str();
	str[str.size()-1] ='\n';
	os<<str;
	return os;
}
ostream& operator<< (ostream &os, matrix& m)
{
	
	for (auto i : m)
	{
		os<<i<<endl;
	}
	
	return os;
}
matrix makeLinearSystem(dvec& u, double t, double h)
{
	dvec a,b,c,d;
	counter s = u.size();

	a.push_back(0);
	b.push_back(1);
	c.push_back(0);
	d.push_back(0);

	for( counter i = 1; i < s-1; i++)
	{
		a.push_back(1/(6*t) - u[i]/(4*h));
		b.push_back(2/(3*t));
		c.push_back(1/(6*t) + u[i]/(4*h));
		d.push_back((u[i+1] + 4*u[i] + u[i-1])/(6*t) + u[i]*(u[i-1] - u[i+1])/(4*h));
	}
	a.push_back(1/(2*t) - u[s-1]/(2*h));
	b.push_back(1/(2*t) + u[s-1]/(2*h));
	c.push_back(0);
	d.push_back((u[s-2] + u[s-1])/(2*t) + u[s-1]*(u[s-2] - u[s-1])/(2*h));

	matrix m;
	m.push_back(a);
	m.push_back(b);
	m.push_back(c);
	m.push_back(d);

	return m;
}

dvec shufle(matrix& m)
{
	dvec a = m[0];
	dvec b = m[1];
	dvec c = m[2];
	dvec d = m[3];
	counter s = b.size();

	dvec et, ks;
	et.push_back(0);
	ks.push_back(0);

	for( counter i = 0; i < s+1; i++)
	{
		ks.push_back(-c[i]/(a[i]*ks[i] + b[i]));
		et.push_back((d[i] - a[i]*et[i])/(a[i]*ks[i] + b[i])); 
	}	
	dvec x(s+1);
	x[s] = 0;
	for ( int i = s-1; i > -1; i--)
	{
		x[i] = ks[i+1]*x[i+1] + et[i+1];
	}
	x.resize(s);
	return x;
}

double makeApp(counter T, counter X)
{
	double t = double(1)/T;
	double h = double(1)/X;
	dvec uu;
	for (counter i = 0; i <= X; i++)
	{
		uu.push_back(i*h);
	}
	dvec yy(uu);
	matrix u,y;
	u.push_back(uu);
	y.push_back(yy);

	for (counter i = 0; i < T; i++)
	{
		matrix&& m = makeLinearSystem(u[i], t, h);
		dvec&& x = shufle(m);
		u.push_back(x);
		dvec yy;
		for (counter j = 0; j <= X; j++)
		{
			yy.push_back(j*h/(1 + (i+1)*t));
		}
		y.push_back(yy);
	}
	matrix e;
	for (counter i = 0; i < u.size(); i++)
	{
		dvec ee;
		for (counter j = 0; j< u[0].size(); j++)
		{
			ee.push_back( abs(y[i][j] - u[i][j]) );
		}
		e.push_back(ee);
	}
	valarray<double> max(e.size());
	for (counter i = 0; i < e.size(); i++)
	{
		valarray<double> val(e[i].data(), e[i].size());
		max[i] = val.max();
	}
	return max.max();
}

int main(int argc, char const *argv[])
{
	ofstream os("out.txt");
	counter N = 5;
	matrix tab;
	const counter A = 5;
	const counter B = 10;
	counter T = pow(2,N);
    dvec TT;
    dvec XX;
	for( counter i = 0; i < A; i++)
	{
		T = 2*T;
		TT.push_back(T);
		dvec e;

		counter X = N;

		for( counter j = 0; j < B; j++)
		{
			X = 2*X;
			if ( i == 0) XX.push_back(X);
			double ee = makeApp(T, X);
			e.push_back(ee);
		}
		tab.push_back(e);
	}
	os<<tab<<endl<<endl;
	os<<"T="<<TT<<endl;
	os<<"X="<<XX<<endl;
	// for (counter i = 0; i < e.size()-1; i++)
	// {
	// 	cout<<"e1="<< e[i]<<endl;
	// 	cout<<"e2="<<e[i+1]<<endl;
	// 	cout<<"p="<<(log2(e[i]/e[i+1]))<<endl;
	// }
	os.close();
	return 0;
}