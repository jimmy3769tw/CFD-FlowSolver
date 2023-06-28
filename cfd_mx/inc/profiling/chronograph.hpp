#pragma once

#include <string>
#include <sstream>

auto get_Chronograph(
    Simulation& simu,
    ClockStruct& chron
){

    std::string A;
    std::string t = ", ";
    std::string n = "\n";

    chron.beginNew.stop();
    chron.beginMainloop.stop();


    static auto s = [&](auto a){ return std::to_string(a);};
    static auto ws = [&](auto a){ std::ostringstream s; s << a; return s.str();};

    
    if (simu.loop == 1)
    {
        A += "MainCount, ";
        A += "ALL, , ";
        A += "SolvePossionEquation, , ";
        A += "#iters, error, ";
        A += "Con & Dif, ";
        A += "T1 to T3, ";
        A += "check SS (L2), ";
        A += "T3 to T0, ";
        A += "BC, ";
        A += "Cd & Cl(IO), ";
        A += "previousRecodingTimer, ";
        A += n;
    }

        A +=  ( s(simu.loop)+t);                                 //MainCount
        A +=  ( s(chron.beginNew.elapsedTime())+t+t);            //ALL
        A +=  ( s(chron.pressure.elapsedTime())+t+t);             //SolvePossionEquation
        A +=  ( s(simu.iters)+t);                                //#iters 
        A +=  ( ws(simu.error)+t);                               //error
        A +=  ( s(chron.convectionDifussion.elapsedTime())+t);   //convection_and_difussion
        A +=  ( s(chron.updateT1toT3.elapsedTime())+t);          //update_UandF_seq
        A +=  ( s(chron.checkL2norm.elapsedTime())+t) ;          //checkL2norm
        A +=  ( s(chron.updateT3toT0.elapsedTime())+t);          //T3toT0
        A +=  ( s(chron.BC.elapsedTime())+t);                    //BC
        if (simu.Dfib_boolT) A +=  ( s(chron.getCdCl.elapsedTime())+t);                    
        A +=  ( s(chron.preIOclock.elapsedTime())+t);                    
        A +=  n;
                                                                // Time_file << chron.beginMainloop.elapsedTime()<< ", ";    //Mainloop
    chron.beginNew.stop();
    chron.beginMainloop.stop();

    return A;
}



auto get_Chronograph_datString(
    Simulation& simu,
    ClockStruct& chron
){

    std::string A;
    std::string t = ", ";
    std::string tab = " ";
    std::string n = "\n";
    // * ---------------------------------------init 
    if (simu.loop == 1)
    {
        std::vector<std::string> variables;
        variables.reserve(3);
        variables.push_back("simulation time");
        variables.push_back("Calculation time (pressure)");
        variables.push_back("#iters");

                    A +=  "TITLE     = \"\"\n";
                    A += "VARIABLES = \"";
                    A += variables.at(0);
                    A += "\",\"";
                    A += variables.at(1);
                    A += "\",\"";
                    A += variables.at(2);
                    A += "\"\n";
                    A += "ZONE T=\"";
                    A += simu.ZONE();
                    A += "\"\n";
    }
    // * ---------------------------------------init 


    static auto s = [&](auto a){ return std::to_string(a);};
    static auto ws = [&](auto a){ std::ostringstream s; s << a; return s.str();};

    A +=  ( ws(simu.getSimuTime())+tab);                        //MainCount
    A +=  ( s(chron.pressure.elapsedTime())+tab);             //SolvePossionEquation
    A +=  ( s(simu.iters)+tab+tab);                           //SolvePossionEquation
    A +=  n;

    return A;
}


auto immediatlyRerecorderTime_csv(
    std::string N,
    Simulation& simu,
    ClockStruct& chron
){
    std::ofstream Time_file;
    std::string fileName = "Information/Chronograph/" + N;
    fileName += std::to_string(simu.TID);
    fileName += ".csv";
    Time_file.open (fileName, std::ios::out|ios::app);
    Time_file << get_Chronograph(simu, chron);
    Time_file.close();
}

auto immediatlyRerecorderTime_dat(
    std::string N,
    Simulation& simu,
    ClockStruct& chron
){

    std::ofstream Time_file;

    std::string fileName = "Information/Chronograph/" + N;

    fileName += std::to_string(simu.TID);

    fileName += ".dat";

    Time_file.open (fileName, std::ios::out|ios::app);

    Time_file << get_Chronograph_datString(simu, chron);

    Time_file.close();

    return true;
}




void recordingTime(
    ClockStruct& sec_chronograph,
    Simulation& simu
)
{

    #if defined (IMMEDIATELY_RECORDING)

    immediatlyRerecorderTime_csv("Time_in_eachLoop",simu, sec_chronograph);
    immediatlyRerecorderTime_dat("Time_in_eachLoop",simu, sec_chronograph);

    if ((simu.loop)%100 == 1){
        immediatlyRerecorderTime_csv("Time_in_100Loop",simu, sec_chronograph);
        immediatlyRerecorderTime_dat("Time_in_100Loop",simu, sec_chronograph);
    }

    if ((simu.loop)%1000 == 1){
        immediatlyRerecorderTime_csv("Time_in_1000Loop",simu, sec_chronograph);
        immediatlyRerecorderTime_dat("Time_in_1000Loop",simu, sec_chronograph);
    }

    if ((simu.loop)%10000 == 1){
        immediatlyRerecorderTime_csv("Time_in_10000Loop",simu, sec_chronograph);
        immediatlyRerecorderTime_dat("Time_in_10000Loop",simu, sec_chronograph);
    }

    #else


    #endif

}