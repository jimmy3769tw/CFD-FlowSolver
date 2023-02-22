#include <fstream>
#include <sys/types.h>
#include <dirent.h>


bool init_information(){

    if(opendir("Information") == NULL)
    { if (system("mkdir Information") != 0){ return 1; } }

    std::ofstream file;

    std::string fName = "Information/info";
    
    fName += ".csv";
    
    file.open (fName,std::ios::out); // | ios::app

    file 
    
    << "\nnx, ny, nz = " 
    << gA.nx << "," << gA.ny << "," << gA.nz 

    << "\nlx, ly, lz = " 
    << gA.lx << "," << gA.ly << "," << gA.lz

    << "\nRe = "
    << simu.Re 

    << "\nconvectionScheme = "
    << simu.convectionScheme

    << "\ndifussionScheme = "
    << simu.difussionScheme

    << std::endl;
}


bool finish_information(){




}