#include <fstream>
#include <iostream>
#include <string>
float** read_file(std::string filename,int rows,int cols)
{
	std::fstream file;
	file.open(filename.c_str(), std::ios::in);
	if(!file.is_open()){return 0;}

	float** floats = new float*[cols+1];
	for(int i = 0; i <cols;++i){ floats[i] = new float[rows+1]; }

	//read each row
	for(int i =0;i<rows;++i)
	{
		for(int j =0;j<cols;++j)//push into the col
		{ file >>floats[j][i]; }
	}
	file.close();

	return floats;
}

//int main()
//{
//	int rows = 54;
//	int cols = 2;
//	float** floats;
//	if( !(floats = read_file("energybin.out",rows,cols) ) ){return 0;}

	//write out the data
//	for(int i =0;i<rows;++i)
//	{
//		for(int j =0;j<cols;++j){ std::cout <<floats[j][i] <<"\t"; }
//		std::cout <<"\n";
//	}

//	return 0;
//}
