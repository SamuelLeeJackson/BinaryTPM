#include <iostream>			
#include <iomanip>																									
#include <fstream>																												
#include <sstream>	
using namespace std;

int main(int argc, char *argv[]) {
    string target = argv[1];
	string filename=target+".obj";
	ifstream in(filename.c_str(), ios::in);

	ofstream out(target+".txt");

	string line;
	float x, y, z;
	int a, b, c;
    while (getline(in, line))
    {
        if (line.substr(0,2) == "v ")
        {
            istringstream s(line.substr(2));
            s >> x; s >> y; s >> z;
			out << setprecision(17) << x << "\t" << y << "\t" << z << endl;
        }
		if (line.substr(0,2) == "f ")
        {
            istringstream s(line.substr(2));
            s >> a; s >> b; s >> c;
            a=a-1602;
            b=b-1602;
            c=c-1602;
            a--; b--; c--;
           out << a << "\t" << b << "\t" << c << endl;
		}
	}

};
