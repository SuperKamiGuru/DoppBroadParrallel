// Created by: Wesley Ford
//This program either takes a single cross-section data file, a directory full of cross-section data files
// or a directory contianing the directories Elastic, Inelastic, Fission and Capture which contian cross-section data files,
// it sorts out which files are valid for doppler broadening, it doppler broadens the data based off the given temperature difference
// and then it outputs the doppler broadened data into the given output filename or a given output directory with the same
// internal structure as the input directory

// many of the ideas used in this file and many of the classes used by this file are copied from the GEANT4 source code.
// This was primarily done to make sure the outputed data is consitent with what the user would get when using the on the fly doppler broadening algorythm.
// A some modifications were made to the copied code however to improve the efficiency and portablitiy of this program.

#define Timer 2

using namespace std;

#include "include/Random/RandomEngine.h"
#include "include/Random/RanecuEngine.h"
//#include <CLHEP/Random/Randomize.h>
#include "include/Random/RandGauss.h"
#include "include/Random/RandGaussQ.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include "zlib.h"
#include <ctime>
#include "include/G4NucleiPropertiesTableAME03.hh"
#include "include/G4NucleiPropertiesTheoreticalTable.hh"
#include <dirent.h>
#include "include/ElementNames.hh"
#include "topc.h"
#include "include/TaskInput.hh"
#include "include/TaskOutput.hh"
#include "include/MarshaledTaskInput.h"
#include "include/MarshaledTaskOutput.h"
#include "include/MarshaledObj.h"
#include <new>

// Error is due to passing of string through buffer, char array is currently not unmarshelling properly, use MyString file in the marshelgen examples
// Error is also due to the accessing of the logfile by the slave, check log creation

// things to do

// X automatically create the temperature directories if they do not exist ie: 1000k
// X recursively search through the directories contianed with in the input file untill the CS data files are found
// X input an overwrite int which will determine whether the program will skip the CS file if the outputfile already exists, skip if the file exists and is complete (check# of data points), or overite all
// X take in a text file of the isotopes to be broadened and what temperatures they will be broadened to
// X take in int to select whether the data will be doppler broadened from the closest existing temperature CS file (See StorkElementData) or from prevTemp
// X take sudo or root password so that it can create files in root protected areas
// X sort the list of input CS files from the macrofile by temperature and isotope so that isotopes with highier temperatures will be broaden from the output of those with lower temp
// X get program to use natural elements when isotopes are unavailable
// X change precision of CS to be within 0.1%
// allow the program to run in parrallel
// get rid of unnecessary data being passed through the task input (the data that the slaves already contain and the data for PrintProgress) and use global variables instead
// create a small code that will filter out points in flat regions and use it before and after DopplerBroad in or to improve the speed of Doppler Broad and G4Stork
// create small code to convert from endf NIST and MCNP to G4NDL
// allow DopplerBroadMacroCreater to take in G4STORK interpolation files so that multiple CS files are created at intervals to cover the temperature interpolation range for all the isotopes in the material
// comment out code
// create manual

// default mass units are m*c^2=MeV, the same as energy
// default momentum units are p*c=MeV, same as energy
// the velocity is kept in terms of v=v/c on purpose for the below calc
// see G4NucleiPropertiesTableAME03 for finding isotope nucleus mass

//use GetDataStream instead
//make the outputfile binary and compressed for speed

//Note the GetDataStream data function used in this program is just a slight modification of
// G4NeutronHPManager::GetDataStream and the doppler broadening algorythm is based off the
// G4NeutronHP[Process]Data::GetCrossSection

const double k_Boltzmann=8.617343*pow(10,-11);
const double nMassC2=939.56563;
const double electron_mass_c2 = 0.510998910;
const double amu_c2 = 931.494028;
const double eV=1.e-6;
const double keV=1.e-3;

void GetClosestTempFileList(string inFileName, string outFileName, stringstream& ss, std::vector<string> &inFileList, std::vector<string> &inFileListNoDep, std::vector<string> &outFileList, std::vector<double> &prevTempList, std::vector<double> &newTempList);
bool GetAllFiles(string inFileName, std::vector<string> &inFileList);
bool FindTemp(string name, double &tempMatch);
bool FindTemp(string name, double tempMatch, int &count);
bool FindProcess(string fileName, int &process);
bool CompleteFile(string fileName);
void SortList(std::vector<string> &inFileList, std::vector<string> &inFileListNoDep, std::vector<string> &outFileList, std::vector<double> &newTempList, std::vector<double> &prevTempList);
void SwapListElem(std::vector<string> &inFileList, std::vector<string> &outFileList, std::vector<double> &newTempList, std::vector<double> &prevTempList, int i, int j);
bool CompareIsotopeNum(string name1, string name2, string comparison);
bool DirectoryExists( const char* pzPath );
bool FindDir(string &inFileName, string &outFileName, double prevTemp, double newTemp);
bool ExtractDouble(string name, double &num);
void GetClosestTempDir(string &inFileName, double &prevTemp, double newTemp);

#if Timer>=1
    void GetFileSize2List(std::vector<string> &fileList, std::vector<double> &fileSize2List, int &totalNumFiles, double &totalFileSize2);
    double GetFileSize2(string fileName);
    void PrintProgress(int index, string fileName, std::vector<double> &fileSize2List, double &totalFileSize2, int &totalNumFiles, double duration, double &sumFileSize2, double &sumDuration, bool success);
    void GetFileSize2List(string inFileName, std::vector<double> &fileSize2List, int &totalNumFiles, double &totalFileSize2);
    void GetDirectoryFileSize2(string inDirName, std::vector<double> &fileSize2List, int &totalNumFiles, double &totalFileSize2);
    void ConvertDirect(string inDirName, string outDirName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile,
    int totalNumFiles, int &index, double &totalFileSize2, double &sumFileSize2, double &sumDuration, std::vector<double> &fileSize2List, bool overWrite);
    bool isApplicable(string fileName);
    void ExtractZA(string fileName, int &Z, int &A);
#else
    void ConvertDirect(string inDirName, string outDirName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile, bool overWrite);
#endif

bool ConvertFile(string inFileName, string outFileName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile, bool overWrite);
bool ConvertFile(string inDirName, string outDirName, string fileName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile, bool overWrite);
void ExtractZA(string fileName, int &Z, int &A, bool log, std::ofstream* logFile);
double GetNuclearMass(double A, double Z, bool log, std::ofstream* logFile);
double NuclearMass(double A, double Z, bool log, std::ofstream* logFile);
double  AtomicMass(double A, double Z);
double  BindingEnergy(double A, double Z);
int DoppBroad(string inFileName, string outFileName, double prevTemp, double newTemp, double isoMassC2, bool ascii, bool log, std::ofstream* logFile);
double findCS(double nKEnerTrans, const double* prevEnVec, const double* prevCSVec, int vecSize);
void GetDataStream( string, std::stringstream&, bool, std::ofstream*);
void SetDataStream( string, std::stringstream&, bool ascii, bool, std::ofstream*);

std::vector<double> theFileSize2List;
std::ofstream *logFileData=NULL;
double theTotalFileSize2=0, theSumFileSize2=0, theSumDuration=0;
int fileIndex=0, fileCount=0, theTotalNumFiles=0;

bool master;
TaskOutput *taskOut;
MarshaledTaskOutput *taskOutMarsh;
MarshaledTaskInput *taskInMarsh;
TaskInput *taskIn;

TOPC_BUF DoBroadening(void *input)
{
    #if Timer>=1
    std::clock_t start;
    start = std::clock();
    #endif

    bool success=0;
    double duration=0;
    string inDirName, outDirName, fileName;

    MarshaledTaskInput inputTaskMarsh(input);
    TaskInput *inputTask = inputTaskMarsh.unmarshal();

    inDirName = (inputTask->inDirName).GetString();
    outDirName = (inputTask->outDirName).GetString();
    fileName = (inputTask->fileName).GetString();

    success=ConvertFile(inDirName, outDirName, fileName, inputTask->prevTemp, inputTask->newTemp, inputTask->ascii, inputTask->log, logFileData, inputTask->overWrite);

    #if Timer>=1
    duration = double(std::clock()-start)/CLOCKS_PER_SEC;
    #endif

    if(taskOut)
        delete taskOut;
    #if Timer>=1
    taskOut = new TaskOutput(inDirName.c_str(), outDirName.c_str(), fileName.c_str(), inputTask->prevTemp, inputTask->newTemp,
                            inputTask->ascii, inputTask->log, inputTask->overWrite, inputTask->index, duration, success);
    #else
    taskOut = new TaskOutput(inDirName.c_str(), outDirName.c_str(), fileName.c_str(), inputTask->prevTemp, inputTask->newTemp,
                            inputTask->ascii, inputTask->log, inputTask->overWrite, duration, success);
    #endif

    if(taskOutMarsh)
        delete taskOutMarsh;
    taskOutMarsh = new MarshaledTaskOutput(taskOut);

    return TOPC_MSG( taskOutMarsh->getBuffer(), taskOutMarsh->getBufferSize());
}

TOPC_ACTION CheckProgress(void *input, void *output)
{
    #if Timer>=1

    MarshaledTaskOutput outputTaskMarsh(output);
    TaskOutput* outputTask = outputTaskMarsh.unmarshal();
    string fileName = (outputTask->fileName).GetString();

    PrintProgress(outputTask->index, fileName, theFileSize2List, theTotalFileSize2, theTotalNumFiles, outputTask->duration, theSumFileSize2, theSumDuration, outputTask->success);
    #endif
    return NO_ACTION;
}

int main(int argc, char **argv)
{
    TOPC_OPT_trace = 0;
	TOPC_OPT_slave_timeout=604800;
    TOPC_init(&argc, &argv);
    master = bool(TOPC_is_master());
    int processID = TOPC_rank();

    CLHEP::HepRandom *CLHEPRand = new CLHEP::HepRandom();

    CLHEP::RanecuEngine *theEngine = new CLHEP::RanecuEngine(2013092304);

    CLHEPRand->setTheEngine(theEngine);

    CLHEPRand->showEngineStatus();
    CLHEPRand->saveEngineStatus();

    ElementNames elementNames;
    elementNames.SetElementNames();

    //double duration;
    int result=0;
    bool macro=false;
    bool closestTemp =false, overWrite =true;
    string closestTempName = "false";
    string macroFileName;
    string inFileName, outFileName;
    string inSubDirName, outSubDirName;
    string outputType="false", createLogFile="false", overWriteFile="true";
    string logFileName;
    double newTemp;
    double prevTemp;
    bool ascii=true, log=false;
    stringstream ss;

    if(argc==9)
    {
        ss << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5] << ' ' << argv[6] << ' ' << argv[7] << ' ' << argv[8];
        ss >> inFileName >> outFileName >> prevTemp >> newTemp >> outputType >> createLogFile >> closestTempName >> overWriteFile;
    }
    else if(argc==8)
    {
        ss << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5] << ' ' << argv[6] << ' ' << argv[7];
        ss >> inFileName >> outFileName >> prevTemp >> newTemp >> outputType >> createLogFile >> closestTempName;
    }
    else if(argc==7)
    {
        ss << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5] << ' ' << argv[6];
        ss >> inFileName >> outFileName >> prevTemp >> newTemp >> outputType >> createLogFile;
    }
    else if(argc==6)
    {
        ss << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4] << ' ' << argv[5];
        ss >> inFileName >> outFileName >> prevTemp >> newTemp >> outputType;
    }
    else if(argc==5)
    {
        ss << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];
        ss >> inFileName >> outFileName >> prevTemp >> newTemp;
    }
    else if(argc==2)
    {
        ss << argv[1];
        ss >> macroFileName;
        if((macroFileName=="help")||(macroFileName=="Help")||(macroFileName=="h")||(macroFileName=="H"))
        {
            cout << "\n### The Potential User Inputs For DopplerBroadData Are:\n" <<
            "inFileName: (file name/directory name) The name of the cross-section data file to be broadened or\n" <<
            " the name of the directory containing the cross-section files to be broadened\n" <<
            "outFileName: (file name/directory name) The name of the cross-section data file for the broadened data to be placed or\n" <<
            " the name of the directory for the broadened cross-section files to be placed\n" <<
            "prevTemp: (number) the temperature of the cross-section data contianed by inFileName, Note: if the\n" <<
            " closestTemp flag is set this input will be ignored\n" <<
            "newTemp: (number) the temperature that the output files will be broadened to, Note if the macroFile\n" <<
            " option is used then this input will be ignored\n" <<
            "outputType: (Default=ascii/compress/...) this determines the format of the output files, if ascii is select then normall\n" <<
            " text files are created if compress is selected then .z compressed files are created\n" <<
            "createLogFile: (Default=false/true) this determines whether a log file is created to contian any errors or signeificant events\n" <<
            "overWriteFile: (Default=true/false) this determines whether existing Doppler broadened files will be regenerated if they have the same name\n" <<
            "See manual for more detials about the inputs and what they do ###" << endl;
            delete CLHEPRand;
            delete theEngine;
            elementNames.ClearStore();
            TOPC_raw_end_master_slave();
            TOPC_finalize();
            return 0;
        }

        else
        {
            ifstream macroFile (macroFileName.c_str(), std::ios::in | std::ios::ate);

            if(macroFile.good())
            {
                int file_size = macroFile.tellg();
                macroFile.seekg( 0 , std::ios::beg );
                char* filedata = new char[ file_size ];
                while ( macroFile )
                {
                    macroFile.read( filedata , file_size );
                }
                macroFile.close();
                string *data = new string ( filedata , file_size );
                delete [] filedata;
                if (data != NULL)
                {
                    ss.str(*data);
                    if(data->back()!='\n')
                        ss << "\n";
                    ss.seekg( 0 , std::ios::beg );
                }
                delete data;

                int numArg;
                ss >> numArg;

                if(numArg==8)
                {
                    ss >> inFileName >> outFileName >> closestTempName >> prevTemp >> outputType >> createLogFile >> overWriteFile;
                }
                else if(numArg==7)
                {
                    ss >> inFileName >> outFileName >> closestTempName >> prevTemp >> outputType >> createLogFile;
                }
                else if(numArg==6)
                {
                    ss >> inFileName >> outFileName >> closestTempName >> prevTemp >> outputType;
                }
                else if(numArg==5)
                {
                    ss >> inFileName >> outFileName >> closestTempName >> prevTemp;
                }

                macro = true;
            }
            else
            {
                cout << "\n### Error Invalid Input, type DoppBroad help to see valid inputs ###" << endl;
            }

        }
    }
    else
    {
        cout << "\n### Error Invalid Number Of Inputs\n" << "The Potential User Inputs For DopplerBroadData Are:\n" <<
            "inFileName: (file name/directory name) The name of the cross-section data file to be broadened or\n" <<
            " the name of the directory containing the cross-section files to be broadened\n" <<
            "outFileName: (file name/directory name) The name of the cross-section data file for the broadened data to be placed or\n" <<
            " the name of the directory for the broadened cross-section files to be placed\n" <<
            "prevTemp: (number) the temperature of the cross-section data contianed by inFileName, Note: if the\n" <<
            " closestTemp flag is set this input will be ignored\n" <<
            "newTemp: (number) the temperature that the output files will be broadened to, Note if the macroFile\n" <<
            " option is used then this input will be ignored\n" <<
            "outputType: (Default=ascii/compress/...) this determines the format of the output files, if ascii is select then normall\n" <<
            " text files are created if compress is selected then .z compressed files are created\n" <<
            "createLogFile: (Default=false/true) this determines whether a log file is created to contian any errors or signeificant events\n" <<
            "overWriteFile: (Default=true/false) this determines whether existing Doppler broadened files will be regenerated if they have the same name\n" <<
            "See manual for more detials about the inputs and what they do ###" << endl;
            delete CLHEPRand;
            delete theEngine;
            elementNames.ClearStore();
            TOPC_raw_end_master_slave();
            TOPC_finalize();
            return 1;
    }

    if(outputType == "compressed"||outputType == "compress"||outputType == "Compressed"||outputType == "Compress"||outputType == "Zipped"||outputType == "Zip"||outputType == "zipped"||outputType == "zip" )
        ascii=false;

    if(createLogFile == "true"||createLogFile == "True"||createLogFile == "1"||createLogFile == "On"||createLogFile == "on")
        log=true;

    if(closestTempName == "true"||closestTempName == "True"||closestTempName == "1"||closestTempName == "On"||closestTempName == "on")
        closestTemp=true;

    if(overWriteFile == "false"||overWriteFile == "False"||overWriteFile == "0"||overWriteFile == "Off"||overWriteFile == "off")
        overWrite=false;


    if((!closestTemp) && (!macro) && (prevTemp>=newTemp))
    {
        cout << "### Invalid Temperature Difference For Doppler Broadening ###\n" << endl;
        delete CLHEPRand;
        delete theEngine;
        elementNames.ClearStore();
        TOPC_raw_end_master_slave();
        TOPC_finalize();
        return 1;
    }

    if(log)
    {
        if (outFileName.back()=='/')
        {
            logFileName=outFileName+"LogDopplerBroadData";
            if(!(DirectoryExists(outFileName.c_str())))
            {
                system( ("mkdir -p -m=666 "+outFileName).c_str());
                if(DirectoryExists(outFileName.c_str()))
                {
                    cout << "### Created Output Directory " << outFileName << " ###" << endl;

                }
                else
                {
                    delete CLHEPRand;
                    delete theEngine;
                    elementNames.ClearStore();
                    TOPC_raw_end_master_slave();
                    TOPC_finalize();
                    return 1;
                }
            }
        }
        else
        {
            size_t pos = outFileName.find_last_of('/');
            if(pos == std::string::npos)
            {
                logFileName="./LogDopplerBroadData";
            }
            else
            {
                logFileName=outFileName.substr(0,pos+1)+"LogDopplerBroadData";
            }

        }

        if(!master)
        {
            stringstream numConv;
            numConv << processID;
            logFileName+="Slave"+numConv.str();
            numConv.clear();
            numConv.str("");
        }
        else
        {
            logFileName+="Master";
        }

        logFileName+=".txt";

        logFileData = new std::ofstream( logFileName.c_str() , std::ios::out | std::ios::trunc );
        if(!logFileData->good())
        {
            cout << "### Error: could not open log file ###" << endl;
            if(logFileData)
                delete logFileData;
            delete CLHEPRand;
            delete theEngine;
            elementNames.ClearStore();
            TOPC_raw_end_master_slave();
            TOPC_finalize();
            return 1;
        }

        if(master)
        {
           logFileData->fill('#');
            (*logFileData) << std::setw(84) << std::left << ""<< endl;
            (*logFileData) << "The Input File Name is " << inFileName << endl;
            (*logFileData) << "The Output File Name is " << outFileName << endl;
            if(macro)
                (*logFileData) << "The Macro File " << macroFileName << " is Being Used" << endl;
            else
                (*logFileData) << "The Temperature Difference is " << (newTemp-prevTemp) << endl;
            (*logFileData) << "The Output Format is " << ((ascii)? "ASCII":"Compressed") << endl;
        }
    }

    TOPC_raw_begin_master_slave( DoBroadening, CheckProgress, NULL );

    #if Timer>=1
        theFileSize2List.reserve(1400);
    #endif
// creates doppler broadened files using the settings in the macrofile
    if(macro)
    {
        if(closestTemp)
        {
            std::vector<double> prevTempList, newTempList;
            std::vector<string> inFileList, inFileListNoDep, outFileList;
            string inDirName, outDirName;
            int temp=1;

            GetClosestTempFileList(inFileName, outFileName, ss, inFileList, inFileListNoDep, outFileList, prevTempList, newTempList);

            cout << "\n" << inFileList.size() << " Files to be Broadened" << endl;

            cout << "\n" << endl;

            cout << "\n" << outFileList.size() << " Files that will be Created are:" << endl;

            cout << "\n" << endl;

            #if Timer>=1
                int index=0;
                GetFileSize2List(inFileListNoDep, theFileSize2List, theTotalNumFiles, theTotalFileSize2);
            #endif

            for(int i=0; i<int(inFileList.size()); i++)
            {

                inDirName = (inFileList[i]).substr(0,(inFileList[i]).find_last_of('/')+1);
                inFileList[i] = (inFileList[i]).substr((inFileList[i]).find_last_of('/')+1, string::npos);

                outDirName = (outFileList[i]).substr(0,(outFileList[i]).find_last_of('/')+1);

                if(!(DirectoryExists(outDirName.c_str())))
                {
                    temp = system( ("mkdir -p -m=666 "+outDirName).c_str());
                    if(DirectoryExists(outDirName.c_str()))
                    {
                        if(log)
                            (*logFileData) << "### Created Output Directory " << outDirName << " ###" << endl;
                        temp=1;
                    }
                }
                if(temp)
                {
                    if(taskIn)
                        delete taskIn;
                    #if Timer>=1
                        taskIn = new TaskInput(inDirName.c_str(), outDirName.c_str(), (inFileList[i]).c_str(), prevTempList[i], newTempList[i], ascii, log, overWrite, fileIndex);
                    #else
                        taskIn = new TaskInput(inDirName, outDirName, inFileList[i], prevTempList[i], newTempList[i], ascii, log, overWrite);
                    #endif

                    if(taskInMarsh)
                        delete taskInMarsh;

                    taskInMarsh = new MarshaledTaskInput(taskIn);
                    TOPC_raw_submit_task_input( TOPC_MSG(taskInMarsh->getBuffer(), taskInMarsh->getBufferSize()) );
                }
                fileIndex++;
            }

            //TOPC_raw_end_master_slave();
        }
        else
        {
            #if Timer>=1
                int index=0;
                GetFileSize2List(inFileName, theFileSize2List, theTotalNumFiles, theTotalFileSize2);
            #endif

            FindDir(inFileName, outFileName, prevTemp, newTemp);

            #if Timer>=1
                ConvertDirect(inFileName, outFileName, prevTemp, newTemp, ascii, log, logFileData, theTotalNumFiles, index, theTotalFileSize2, theSumFileSize2, theSumDuration, theFileSize2List, overWrite);
            #else
                ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, overWrite);
            #endif
        }
    }

    else if (inFileName.back()=='/')
    {
        //if closestTemp is set search in the given input file directory for a temperature directory closest to newTemp
        if(closestTemp)
        {
            GetClosestTempDir(inFileName, prevTemp, newTemp);
        }
        // else if the inFileName doesn't already contain the temperature directory search in the given input file directory for a temperature directory with temperature prevTemp
        else
        {
             FindDir(inFileName, outFileName, prevTemp, newTemp);
        }

        #if Timer>=1
            int index=0;
            GetFileSize2List(inFileName, theFileSize2List, theTotalNumFiles, theTotalFileSize2);
        #endif

        DIR *dir;
        struct dirent *ent;
        if ((dir = opendir (inFileName.c_str())) != NULL)
        {
          /* print all the files and directories within directory */
          while ((ent = readdir (dir)) != NULL)
          {
            if(string(ent->d_name)=="Elastic")
            {
                inSubDirName = inFileName + "Elastic/CrossSection/";
                outSubDirName = outFileName + "Elastic/CrossSection/";

                #if Timer>=1
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, theTotalNumFiles, index, theTotalFileSize2, theSumFileSize2, theSumDuration, theFileSize2List, overWrite);
                #else
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, overWrite);
                #endif

            }
            else if(string(ent->d_name)=="Inelastic")
            {
                inSubDirName = inFileName + "Inelastic/CrossSection/";
                outSubDirName = outFileName + "Inelastic/CrossSection/";

                #if Timer>=1
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, theTotalNumFiles, index, theTotalFileSize2, theSumFileSize2, theSumDuration, theFileSize2List, overWrite);
                #else
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, overWrite);
                #endif
            }
            else if(string(ent->d_name)=="Fission")
            {
                inSubDirName = inFileName + "Fission/CrossSection/";
                outSubDirName = outFileName + "Fission/CrossSection/";

                #if Timer>=1
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, theTotalNumFiles, index, theTotalFileSize2, theSumFileSize2, theSumDuration, theFileSize2List, overWrite);
                #else
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, overWrite);
                #endif
            }
            else if(string(ent->d_name)=="Capture")
            {
                inSubDirName = inFileName + "Capture/CrossSection/";
                outSubDirName = outFileName + "Capture/CrossSection/";

                #if Timer>=1
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, theTotalNumFiles, index, theTotalFileSize2, theSumFileSize2, theSumDuration, theFileSize2List, overWrite);
                #else
                    ConvertDirect(inSubDirName, outSubDirName, prevTemp, newTemp, ascii, log, logFileData, overWrite);
                #endif
            }
            else
            {
                if(taskIn)
                    delete taskIn;
                #if Timer>=1
                    taskIn = new TaskInput(inFileName.c_str(), outFileName.c_str(), ent->d_name, prevTemp, newTemp, ascii, log, overWrite, fileIndex);
                #else
                    taskIn = new TaskInput(inDirName, outDirName, ent->d_name, prevTemp, newTemp, ascii, log, overWrite);
                #endif

                if(taskInMarsh)
                    delete taskInMarsh;

                taskInMarsh = new MarshaledTaskInput(taskIn);
                TOPC_raw_submit_task_input( TOPC_MSG(taskInMarsh->getBuffer(), taskInMarsh->getBufferSize()) );
                fileIndex++;
            }
          }
          closedir(dir);
        }
        else
        {
            cout << "Error: Could not open the given directory" << endl;
            TOPC_raw_end_master_slave();
            if (logFileData)
                delete logFileData;
            if(taskOut)
                delete taskOut;
            if(taskOutMarsh)
                delete taskOutMarsh;
            if(taskIn)
                delete taskIn;
            if(taskInMarsh)
                delete taskInMarsh;
            delete CLHEPRand;
            delete theEngine;
            elementNames.ClearStore();
            TOPC_finalize();
            return 1;
        }

    }
    else
    {
        ConvertFile(inFileName, outFileName, prevTemp, newTemp, ascii, log, logFileData, overWrite);
    }

    ss.str("");
    ss.clear();

    TOPC_raw_end_master_slave();

    elementNames.ClearStore();

    if (logFileData)
        delete logFileData;
    if(taskOut)
        delete taskOut;
    if(taskOutMarsh)
        delete taskOutMarsh;
    if(taskIn)
        delete taskIn;
    if(taskInMarsh)
        delete taskInMarsh;
    delete CLHEPRand;
    delete theEngine;

    TOPC_finalize();

    return result;

}

void GetClosestTempFileList(string inFileName, string outFileName, stringstream& ss, std::vector<string> &inFileList, std::vector<string> &inFileListNoDep, std::vector<string> &outFileList, std::vector<double> &prevTempList, std::vector<double> &newTempList)
{
    int best[4], index, numIso, count;
    string isoName, name, outName;
    double temp, tempMatch;
    double tempBest[4]={0.,0.,0.,0.};
    std::vector<string> allFiles, matchList;
    bool check[4]={false, false, false, false};
    stringstream numConv;

    ss >> numIso;
    GetAllFiles(inFileName, allFiles);

    for(int i=0; i<numIso; i++)
    {
        ss >> isoName >> temp;
        for(int j=0; j<4; j++)
        {
            check[j]=false;
            tempBest[j]=0.;
        }
        for(int j=0; j<int(allFiles.size()); j++)
        {
            name = (allFiles[j]).substr((allFiles[j]).find_last_of('/')+1, std::string::npos);
            name = name.substr(0, name.find_last_of('.'));

            if(CompareIsotopeNum(name, isoName, "=="))
            {
                if(FindTemp(allFiles[j], tempMatch))
                {
                    if(FindProcess(allFiles[j], index))
                    {
                        if(((temp-tempBest[index])>=(temp-tempMatch))&&(tempMatch<temp))
                        {
                            if(CompleteFile(allFiles[j]))
                            {
                                tempBest[index] = tempMatch;
                                best[index] = j;
                                check[index] = true;
                            }
                        }
                    }
                }
            }
        }
        numConv << temp;
        for(int k=0; k<4; k++)
        {
            if(check[k])
            {
                prevTempList.push_back(tempBest[k]);
                newTempList.push_back(temp);
                inFileList.push_back(allFiles[best[k]]);
                // create out file directory from temperature that the file will be doppler broadened too
                FindTemp(allFiles[best[k]], tempBest[k], count);
                outName = outFileName+numConv.str()+'k'+'/'+(allFiles[best[k]]).substr(count+1, std::string::npos);
                outFileList.push_back(outName);
            }
            else if(k==1||k==0||k==3)
            {
                cout << "Error: the input CS file for isotope " << isoName << " could not be found for process " << index << "\n" << endl;
            }
        }
        numConv.clear();
        numConv.str("");
    }

    SortList(inFileList, inFileListNoDep, outFileList, newTempList, prevTempList);
}

bool GetAllFiles(string inFileName, std::vector<string> &inFileList)
{
    DIR *dir;
    struct dirent *ent;
    string name;

    if ((dir = opendir (inFileName.c_str())) != NULL)
    {
      while ((ent = readdir (dir)) != NULL)
      {
        if((string(ent->d_name)!="..")&&(string(ent->d_name)!="."))
        {
            if (GetAllFiles(inFileName+ent->d_name+"/", inFileList))
            {

            }
            else
            {
                inFileList.push_back(inFileName+ent->d_name);
            }
        }
      }
    }
    else
    {
        return false;
    }

    return true;
}

bool FindTemp(string name, double &tempMatch)
{
    bool foundTemp=false;
    stringstream numConv;
    int pos1=name.length(), pos2;
    string partName;

    while((pos1!=0)&&(!foundTemp))
    {
        pos2=name.find_last_of('/', pos1);
        pos1=name.find_last_of('/', pos2-1);
        partName=name.substr(pos1+1, pos2-pos1-1);

        if(((partName[0]>='0')&&(partName[0]<='9'))||(partName[0]=='.'))
        {
            foundTemp=true;
            numConv << partName[0];
            for(int i=1; i<int(partName.length()); i++)
            {
                if(((partName[i]>='0')&&(partName[i]<='9'))||(partName[i]=='.'))
                {
                    numConv << partName[i];
                }
            }
        }

    }

    if(foundTemp)
    {
        numConv >> tempMatch;
    }

    numConv.str("");

    return foundTemp;
}

bool FindTemp(string name, double tempMatch, int &count)
{
    bool foundTemp=false;
    stringstream numConv;
    int pos1=name.length(), pos2;
    string partName;
    count=0;
    double match;

    while((pos1!=0)&&(!foundTemp))
    {
        pos2=name.find_last_of('/', pos1);
        pos1=name.find_last_of('/', pos2-1);
        partName=name.substr(pos1+1, pos2-pos1-1);

        if(((partName[0]>='0')&&(partName[0]<='9'))||(partName[0]=='.'))
        {
            foundTemp=true;
            numConv << partName[0];
            count=pos2;
            for(int i=1; i<int(partName.length()); i++)
            {
                if(((partName[i]>='0')&&(partName[i]<='9'))||(partName[i]=='.'))
                {
                    numConv << partName[i];
                }
            }
        }

    }

    if(foundTemp)
    {
        numConv >> match;
        foundTemp = match==tempMatch;
    }

    numConv.str("");

    return foundTemp;
}

bool FindProcess(string fileName, int &process)
{
    bool found=false;
    int pos1=fileName.length(), pos2;
    string partName;

    while((pos1!=0)&&(!found))
    {
        pos2=fileName.find_last_of('/', pos1);
        pos1=fileName.find_last_of('/', pos2-1);
        partName=fileName.substr(pos1+1, pos2-pos1-1);

        if(partName=="capture"||partName=="Capture")
        {
            process=0;
            found=true;
        }
        else if(partName=="elastic"||partName=="Elastic")
        {
            process=1;
            found=true;
        }
        else if(partName=="fission"||partName=="Fission")
        {
            process=2;
            found=true;
        }
        else if(partName=="inelastic"||partName=="Inelastic"||partName=="inElastic"||partName=="InElastic")
        {
            process=3;
            found=true;
        }

    }

    return found;
}

bool CompleteFile(string fileName)
{
    stringstream stream;
    GetDataStream( fileName, stream, false, NULL);
    int count=0, numPoints;
    double dummy;

    if(stream.str()!="")
    {
        //skips teo dummy variables and gets the number of CS points
        stream >> numPoints;
        stream >> numPoints;
        stream >> numPoints;

        while(stream)
        {
            stream >> dummy >> dummy;
            count++;
        }

        count--;

        if(count==numPoints)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
        return false;
}

void SortList(std::vector<string> &inFileList, std::vector<string> &inFileListNoDep, std::vector<string> &outFileList, std::vector<double> &newTempList, std::vector<double> &prevTempList)
{
    int process1, process2;

    // sort the lists by temperature, isotope and process
    for(int i=0; i<int(inFileList.size()); i++)
    {
        for(int j=i+1; j<int(inFileList.size()); j++)
        {
            if(newTempList[j]<newTempList[i])
            {
                SwapListElem(inFileList, outFileList, newTempList, prevTempList, i, j);
            }
            else if(newTempList[j]==newTempList[i])
            {
                if(CompareIsotopeNum(inFileList[j], inFileList[i], "<"))
                {
                    SwapListElem(inFileList, outFileList, newTempList, prevTempList, i, j);
                }
                else if(CompareIsotopeNum(inFileList[j], inFileList[i], "=="))
                {
                    FindProcess(inFileList[j], process1);
                    FindProcess(inFileList[i], process2);

                    if(process1<process2)
                    {
                        SwapListElem(inFileList, outFileList, newTempList, prevTempList, i, j);
                    }
                    if(process1==process2)
                    {
                        inFileList.erase(inFileList.begin()+j);
                        outFileList.erase(outFileList.begin()+j);
                        newTempList.erase(newTempList.begin()+j);
                        prevTempList.erase(prevTempList.begin()+j);
                    }
                }
            }
        }
    }

    inFileListNoDep=inFileList;

    for(int i=1; i<int(inFileList.size()); i++)
    {
        for(int j=0; j<i; j++)
        {
            if(newTempList[j]>prevTempList[i])
            {
                if(CompareIsotopeNum(inFileList[j], inFileList[i], "=="))
                {
                    FindProcess(inFileList[j], process1);
                    FindProcess(inFileList[i], process2);

                    if(process1==process2)
                    {
                        prevTempList[i]=newTempList[j];
                        inFileList[i]=outFileList[j];
                    }
                }
            }
        }
    }
}

void SwapListElem(std::vector<string> &inFileList, std::vector<string> &outFileList, std::vector<double> &newTempList, std::vector<double> &prevTempList, int i, int j)
{
    double newTemp, prevTemp;
    string inName, outName;

    newTemp=newTempList[i];
    prevTemp=prevTempList[i];
    inName=inFileList[i];
    outName=outFileList[i];

    newTempList[i]=newTempList[j];
    prevTempList[i]=prevTempList[j];
    inFileList[i]=inFileList[j];
    outFileList[i]=outFileList[j];

    newTempList[j]=newTemp;
    prevTempList[j]=prevTemp;
    inFileList[j]=inName;
    outFileList[j]=outName;
}

bool CompareIsotopeNum(string name1, string name2, string comparison)
{
    stringstream numConv;
    int Z1=0, A1=0, Z2=0, A2=0, num1, num2;

    name1 = (name1).substr(name1.find_last_of('/')+1, std::string::npos);
    name1 = name1.substr(0, name1.find_last_of('.'));
    ExtractZA(name1, Z1, A1);

    name2 = (name2).substr(name2.find_last_of('/')+1, std::string::npos);
    name2 = name2.substr(0, name2.find_last_of('.'));
    ExtractZA(name2, Z2, A2);

    if((A1==0)||(A2==0))
    {
        if(comparison=="==")
        {
            if(Z1==Z2)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else if(comparison==">")
        {
            if(Z1>Z2)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else if(comparison=="<")
        {
            if(Z1<Z2)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            cout << "\nError Comparison symbol is not recognized in CompareIsotopeNum\n" << endl;
            return false;
        }
    }

    num1 = Z1*1000+A1;
    num2 = Z2*1000+A2;

    if(comparison=="==")
    {
        if(num1==num2)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else if(comparison==">")
    {
        if(num1>num2)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else if(comparison=="<")
    {
        if(num1<num2)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        cout << "\nError Comparison symbol is not recognized in CompareIsotopeNum\n" << endl;
        return false;
    }
}

bool FindDir(string &inFileName, string &outFileName, double prevTemp, double newTemp)
{
    DIR *dir;
    struct dirent *ent;
    bool check=false;
    stringstream numConv;
    double temp;
    const string originalIn=inFileName;
    string outTemp, fileName;

    numConv << newTemp;
    numConv >> outTemp;

    if ((dir = opendir (inFileName.c_str())) != NULL)
    {
      while ((ent = readdir (dir)) != NULL)
      {
        fileName=string(ent->d_name);
        if((fileName!="..")&&(fileName!="."))
        {
            if ((ExtractDouble( fileName, temp))&&(temp==prevTemp)&&((fileName.back()=='k')||(fileName.back()=='K')))
            {
                inFileName=inFileName+ent->d_name+'/';
                // else if the inFileName doesn't already contain the temperature directory search in the given input file directory for a temperature directory with temperature prevTemp
                double tempEx=-100;
                string partName = outFileName.substr(outFileName.find_last_of('/',outFileName.length()-2), string::npos);
                if(!(ExtractDouble(partName, tempEx)&&(tempEx==newTemp)))
                {
                     outFileName=outFileName+outTemp+'k'+'/';
                }
                check=true;
                break;
            }
            else
            {
                inFileName=inFileName+string(ent->d_name)+'/';
                check=FindDir(inFileName, outFileName, prevTemp, newTemp);
                if(check)
                {
                    break;
                }
                else
                {
                    inFileName=originalIn;
                }
            }
        }
      }
    }
    else
    {
        check= false;
    }

    return check;
}

bool ExtractDouble(string name, double &num)
{
    bool foundDouble = false;
    stringstream numConv;
    for(int i=0; i<int(name.size()); i++)
    {
        if(((name[i]>='0')&&(name[i]<='9'))||(name[i]=='.'))
        {
            foundDouble=true;
            numConv << (name)[i];
        }
    }
    if(foundDouble)
        numConv >> num;
    return foundDouble;
}

void GetClosestTempDir(string &inFileName, double &prevTemp, double newTemp)
{
    string closestName, temp;
    double closest=0, test=0;
    stringstream numConv;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (inFileName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {
        if(((ent->d_name)[0]>='0')&&((ent->d_name)[0]<='9'))
        {
            numConv << ent->d_name;
            numConv >> test;

            if(((newTemp-test)<=(newTemp-closest))&&(test<=newTemp))
            {
                closest=test;
                closestName = inFileName+ent->d_name+'/';
            }

        }
        else
        {
            temp=inFileName+ent->d_name+'/';
            GetClosestTempDir(temp, test, newTemp);
            if(((newTemp-test)<=(newTemp-closest))&&(test<=newTemp))
            {
                closest=test;
                closestName = temp;
            }
        }
      }
      closedir(dir);
    }

    prevTemp=closest;
    inFileName=closestName;
}

#if Timer>=1

void GetFileSize2List(std::vector<string> &fileList, std::vector<double> &fileSize2List, int &totalNumFiles, double &totalFileSize2)
{
    for(int i=0; i<int(fileList.size()); i++)
    {
        fileSize2List.push_back(GetFileSize2(fileList[i]));
        totalFileSize2+=fileSize2List[totalNumFiles];
        totalNumFiles++;
    }

}

double GetFileSize2(string fileName)
{
    std::ifstream in( fileName.c_str() , std::ios::binary | std::ios::ate );
    double file_size2 =0;
    if(in.good())
    {
        file_size2 = (in.tellg());
        file_size2 = file_size2*file_size2;
        in.seekg( 0 , std::ios::beg );
        in.close();
    }
    else
    {
        cout << "### Error: failed to open file " << fileName << " to determine the filesize ###" << endl;
    }

    return file_size2;

}

void PrintProgress(int index, string fileName, std::vector<double> &fileSize2List, double &totalFileSize2, int &totalNumFiles, double duration, double &sumFileSize2, double &sumDuration, bool success)
{
    if(success)
    {
        sumFileSize2+=fileSize2List[index];
        sumDuration+=duration;
        double timeLeft = sumDuration/sumFileSize2*(totalFileSize2-sumFileSize2);
        cout << "Files Converted (" << fileCount+1 << "/" << totalNumFiles << "), Duration Since Start " << sumDuration << "s, Time Remaining " << timeLeft << "s\n"
            << "Total Progress [";
        double sizeRatio = sumFileSize2/totalFileSize2;
        for (int i=0; i<int(60*sizeRatio); i++)
        {
            cout << "#";
        }
        for (int i=0; i<int(60-60*sizeRatio); i++)
        {
            cout << "_";
        }
        cout << "]\n" << endl;

        cout.fill('-');
        cout << std::setw(84) << std::right << "\n\n";

        if(int(sizeRatio)==1)
        {
            cout << "The Total Time Taken Was " << sumDuration << "s, " << (totalNumFiles-index) << " Files Were Not Converted" << endl;
        }
        fileCount++;
    }

    else if(isApplicable(fileName))
    {
        sumDuration+=duration;
        totalFileSize2-=fileSize2List[index];
        fileSize2List.erase(fileSize2List.begin()+index);
        totalNumFiles--;
        fileIndex--;

        if(int(sumFileSize2/totalFileSize2)==1)
        {
            cout << "The Total Time Taken Was " << sumDuration << "s, " << (totalNumFiles-index) << " Files Were Not Converted\n" << endl;
        }
    }

}

void GetFileSize2List(string inFileName, std::vector<double> &fileSize2List, int &totalNumFiles, double &totalFileSize2)
{
    string filename, inSubDirName;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (inFileName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {
        if(string(ent->d_name)=="Elastic")
        {
            inSubDirName = inFileName + "Elastic/CrossSection/";
            GetDirectoryFileSize2(inSubDirName, fileSize2List, totalNumFiles, totalFileSize2);
        }
        else if(string(ent->d_name)=="Inelastic")
        {
            inSubDirName = inFileName + "Inelastic/CrossSection/";
            GetDirectoryFileSize2(inSubDirName, fileSize2List, totalNumFiles, totalFileSize2);
        }
        else if(string(ent->d_name)=="Fission")
        {
            inSubDirName = inFileName + "Fission/CrossSection/";
            GetDirectoryFileSize2(inSubDirName, fileSize2List, totalNumFiles, totalFileSize2);
        }
        else if(string(ent->d_name)=="Capture")
        {
            inSubDirName = inFileName + "Capture/CrossSection/";
            GetDirectoryFileSize2(inSubDirName, fileSize2List, totalNumFiles, totalFileSize2);
        }
        else
        {
            if(isApplicable(ent->d_name))
            {
                filename = inFileName + ent->d_name;
                fileSize2List.push_back(GetFileSize2(filename));
                totalFileSize2+=fileSize2List[totalNumFiles];
                totalNumFiles++;
            }
        }
      }
      closedir(dir);
    }
    else
    {
        cout << "Error: Could not open the given directory" << endl;
    }
}

void GetDirectoryFileSize2(string inDirName, std::vector<double> &fileSize2List, int &totalNumFiles, double &totalFileSize2)
{
    DIR *dir;
    struct dirent *ent;
    string fileName;

    if ((dir = opendir (inDirName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {

        fileName = ent->d_name;
        if(isApplicable(fileName))
        {
            fileName = inDirName + fileName;
            fileSize2List.push_back(GetFileSize2(fileName));
            totalFileSize2+=fileSize2List[totalNumFiles];
            totalNumFiles++;
        }

      }
      closedir(dir);
    }
    else
    {
        cout << "### Error: Could not open directory" << inDirName << " ###" << endl;
    }
}

void ConvertDirect(string inDirName, string outDirName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile,
    int totalNumFiles, int &index, double &totalFileSize2, double &sumFileSize2, double &sumDuration, std::vector<double> &fileSize2List, bool overWrite)
{
    DIR *dir;
    struct dirent *ent;
    string fileName;
    int temp=1;

    if ((dir = opendir (inDirName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {
            if(!(DirectoryExists(outDirName.c_str())))
            {
                temp = system( ("mkdir -p -m=666 "+outDirName).c_str());
                if(DirectoryExists(outDirName.c_str()))
                {
                    if(log)
                        (*logFile) << "### Created Output Directory " << outDirName << " ###" << endl;
                    temp=1;
                }
            }


            if(temp)
            {
                fileName = ent->d_name;
                if(taskIn)
                    delete taskIn;

                taskIn = new TaskInput(inDirName.c_str(), outDirName.c_str(), fileName.c_str(), prevTemp, newTemp, ascii, log, overWrite, fileIndex);

                if(taskInMarsh)
                    delete taskInMarsh;

                taskInMarsh = new MarshaledTaskInput(taskIn);
                TOPC_raw_submit_task_input( TOPC_MSG(taskInMarsh->getBuffer(), taskInMarsh->getBufferSize()) );
                fileIndex++;
            }
            else
            {
                cout << "### Error: Could not create missing directory" << outDirName << " ###" << endl;
                return;
            }

      }
      closedir(dir);
    }
    else
    {
        cout << "### Error: Could not open directory" << inDirName << " ###" << endl;
    }
}

bool isApplicable(string fileName)
{
    int Z=-1, A=-1;
    ExtractZA(fileName, Z, A);
    if((Z==-1)||(A==-1))
        return false;
    else
        return true;
}

void ExtractZA(string fileName, int &Z, int &A)
{
    if(fileName=="")
    {
        Z=A=-1;
        return;
    }
    std::size_t startPos=0;
    stringstream ss;
    ElementNames* elementNames;
    while(startPos!=fileName.length() && (fileName[startPos]<'0' || fileName[startPos]>'9'))
        startPos++;

    if(startPos==fileName.length())
    {
        Z=A=-1;
    }
    else
    {
    ////
        std::size_t found1 = fileName.find_first_of('_', startPos);
        if (found1==std::string::npos)
        {
            Z=A=-1;
        }
        else
        {
            std::size_t found2 = fileName.find_first_of('_', found1+1);
            if (found2==std::string::npos)
            {
                Z=A=-1;
            }
            else
            {

                ss.str(fileName.substr(startPos, found1));
                ss >> Z;
                ss.str("");
                ss.clear();
                if(((found2-found1-1) > 2) && (fileName[found2-2] == 'm'))
                    ss.str(fileName.substr(found1+1, found2-found1-3));
                else
                    ss.str(fileName.substr(found1+1, found2-found1-1));
                ss >> A;
                ss.str("");
                ss.clear();
                ss.str(fileName.substr(found2+1));
                if (!(elementNames->CheckName(ss.str(), Z)))
                {
                    Z=A=-1;
                }
                ss.str("");
                ss.clear();
            }

        }

    }
}

#else

void ConvertDirect(string inDirName, string outDirName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile, overWrite)
{
    DIR *dir;
    struct dirent *ent;
    string fileName;
    int temp=1;
    bool sucess;

    if ((dir = opendir (inDirName.c_str())) != NULL)
    {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL)
      {
            if(!(DirectoryExists(outDirName.c_str())))
                temp = system( ("mkdir -p "+outDirName).c_str());

            if(temp)
            {
                fileName = ent->d_name;
                if(taskIn)
                    delete taskIn;

                taskIn = new TaskInput(inDirName, outDirName, fileName, prevTemp, newTemp, ascii, log, overWrite);

                if(taskInMarsh)
                    delete taskInMarsh;
                taskInMarsh = new MarshaledTaskInput(taskIn);
                TOPC_raw_submit_task_input( TOPC_MSG(taskInMarsh->getBuffer(), taskInMarsh->getBufferSize()) );
            }
            else
            {
                cout << "Could not create missing directory" << outDirName << endl;
            }
            temp=1;

      }
      closedir(dir);
    }
    else
    {
        cout << "Could not open directory" << inDirName << endl;
    }
}

#endif

bool ConvertFile(string inFileName, string outFileName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile, bool overWrite)
{

    size_t pos=std::string::npos;
    string fileName, partFileName;
    bool success=true, notOverWriten=false;
    stringstream stream;

    if(log)
    {
        for(int i=0; i<4; i++)
        {
            pos = inFileName.find_last_of('/', pos);
            if(pos == std::string::npos)
            {
                pos=0;
                break;
            }
            else
                pos--;
        }
        if(pos!=0)
            pos+=2;

        string test = inFileName.substr(inFileName.find_last_of('/', 0), std::string::npos);

        partFileName = inFileName.substr(pos, std::string::npos);

        if((test=="..")||(test=="."))
        {
            return false;
        }

        logFile->fill('-');
        (*logFile) << '\n' << std::setw(84) << std::left << "Start Broadening "+partFileName << '\n' << endl;
    }

    if(CompleteFile(inFileName))
    {
        int Z=-1, A=-1, result=0;
        double isoMassC2;

        pos = inFileName.find_last_of('/');
        if(pos == std::string::npos)
            pos=0;
        else
            pos++;

        fileName = inFileName.substr(pos, std::string::npos);
        ExtractZA(fileName, Z, A, log, logFile);
        isoMassC2 = GetNuclearMass(A,Z, log, logFile);
        if(isoMassC2)
        {
            if(overWrite)
                result = DoppBroad(inFileName, outFileName, prevTemp, newTemp, isoMassC2, ascii, log, logFile);
            else
            {
                if(CompleteFile(outFileName))
                {
                    notOverWriten=true;
                }
                else
                {
                    result = DoppBroad(inFileName, outFileName, prevTemp, newTemp, isoMassC2, ascii, log, logFile);
                }
            }
        }
        else
        {
            if(log)
                (*logFile) << "### Invalid Mass Extract From File " << inFileName << " ###" << endl;
            result = 1;
        }
        if(notOverWriten)
        {
            success=false;
            if(log)
                (*logFile) << "### " << inFileName << " Already Exists And Is Complete And Will Not Be Over Written ###" << endl;
        }
        else if(result)
        {
            success=false;
            if(log)
                (*logFile) << "### Error: Broadening File " << inFileName << " ###" << endl;
        }
    }
    else
    {
        if(log)
        {
            (*logFile) << "Input File " << inFileName << " Is Not Complete, Skipping Task" << endl;
        }
    }


    if(log)
    {
        logFile->fill('-');
        (*logFile) << std::setw(84) << std::left << "End Broadening "+partFileName << endl;
    }

    stream.clear();
    stream.str("");

    return success;
}

bool ConvertFile(string inDirName, string outDirName, string fileName, double prevTemp, double newTemp, bool ascii, bool log, std::ofstream* logFile, bool overWrite)
{
    string partFileName;
    bool success=true, notOverWriten=false;
    stringstream stream;

    if(log)
    {
        size_t pos=std::string::npos;

        for(int i=0; i<4; i++)
        {
            pos = inDirName.find_last_of('/', pos);
            if(pos == std::string::npos)
            {
                pos=0;
                break;
            }
            else
                pos--;
        }
        if(pos!=0)
            pos+=2;

        partFileName = inDirName.substr(pos, std::string::npos)+fileName;

        if((fileName=="..")||(fileName=="."))
        {
            return false;
        }

        logFile->fill('-');
        (*logFile) << '\n' << std::setw(84) << std::left << "Start Broadening "+partFileName << '\n' << endl;
    }

    if(CompleteFile(inDirName+fileName))
    {
        int Z=-1, A=-1, result=0;
        double isoMassC2;
        ExtractZA(fileName, Z, A, log, logFile);
        isoMassC2 = GetNuclearMass(A,Z, log, logFile);
        inDirName += fileName;
        outDirName += fileName;
        if(isoMassC2)
        {
            if(overWrite)
                result = DoppBroad(inDirName, outDirName, prevTemp, newTemp, isoMassC2, ascii, log, logFile);
            else
            {
                if(CompleteFile(outDirName))
                {
                    notOverWriten=true;
                }
                else
                {
                    result = DoppBroad(inDirName, outDirName, prevTemp, newTemp, isoMassC2, ascii, log, logFile);
                }
            }
        }
        else
        {
            if(log)
                (*logFile) << "### Invalid Mass Extract From File " << partFileName << " ###" << endl;
            result = 1;
        }

        if(notOverWriten)
        {
            success=false;
            if(log)
                (*logFile) << "### " << partFileName << " Already Exists And Is Complete And Will Not Be Over Written ###" << endl;
        }
        else if(result)
        {
            success=false;
            if(log)
                (*logFile) << "### Error: Broadening File " << partFileName << " ###" << endl;
        }
    }
    else
    {
        if(log)
        {
            (*logFile) << "Input File " << inDirName+fileName << " Is Not Complete, Skipping Task" << endl;
        }
    }

    if(log)
    {
        logFile->fill('-');
        (*logFile) << std::setw(84) << std::left << "End Broadening "+partFileName << endl;
    }

    return success;
}

bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath);

    if (pDir != NULL)
    {
        bExists = true;
        closedir (pDir);
    }

    return bExists;
}

void ExtractZA(string fileName, int &Z, int &A, bool log, std::ofstream* logFile)
{
        std::size_t startPos=0;
        stringstream ss;
        ElementNames* elementNames;
        while(startPos!=fileName.length() && (fileName[startPos]<'0' || fileName[startPos]>'9'))
            startPos++;

        if(startPos==fileName.length())
        {
            if(log)
                (*logFile) << "### File Name Does Not Contian a Z or an A Value " << fileName << " is Invalid for Broadening ###" << endl;
            Z=A=-1;
        }
        else
        {
        ////
            std::size_t found1 = fileName.find_first_of('_', startPos);
            if (found1==std::string::npos)
            {
                if(log)
                    (*logFile) << "### File Name Does Not Contian a '_', two are needed, one to seperate the Z and A value, \
                and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << endl;
                Z=A=-1;
            }
            else
            {
                std::size_t found2 = fileName.find_first_of('_', found1+1);
                if (found2==std::string::npos)
                {
                    if(log)
                        (*logFile) << "### File Name Does Not Contian a second '_', two are needed, one to seperate the Z and A value, \
                    and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << endl;
                    Z=A=-1;
                }
                else
                {

                    ss.str(fileName.substr(startPos, found1));
                    ss >> Z;
                    ss.str("");
                    ss.clear();
                    if(((found2-found1-1) > 2) && (fileName[found2-2] == 'm'))
                        ss.str(fileName.substr(found1+1, found2-found1-3));
                    else
                        ss.str(fileName.substr(found1+1, found2-found1-1));
                    ss >> A;
                    ss.str("");
                    ss.clear();
                    ss.str(fileName.substr(found2+1));
                    if (!(elementNames->CheckName(ss.str(), Z)))
                    {
                        if(log)
                            (*logFile) << "### " << fileName << " does not include the correct element name at the end ###" << endl;
                        Z=A=-1;
                    }
                    ss.str("");
                    ss.clear();
                }

            }

        }
}

double GetNuclearMass(double A, double Z, bool log, std::ofstream* logFile)
{
    double mass;
    if (Z == -1 || A==-1)
        return 0.;

    if (G4NucleiPropertiesTableAME03::IsInTable(Z,A))
    {
      // AME 03 table
      mass = G4NucleiPropertiesTableAME03::GetNuclearMass(Z,A);
    }
    else if (G4NucleiPropertiesTheoreticalTable::IsInTable(Z,A))
    {
      // Theoretical table
      mass = G4NucleiPropertiesTheoreticalTable::GetNuclearMass(Z,A);
    }
    else
    {
      mass = NuclearMass(double(A),double(Z), log, logFile);
    }
    return mass;
}

double NuclearMass(double A, double Z, bool log, std::ofstream* logFile)
{
  if (A < 1 || Z < 0 || Z > A) {

      if(log)
        (*logFile) << "NuclearMass: Wrong values for A = "
	     << A << " and Z = " << Z << endl;
    return 0.0;
  }

  double mass = AtomicMass(A,Z);
  // atomic mass is converted to nuclear mass according formula in  AME03
  mass -= Z*electron_mass_c2;
  mass += ( 14.4381*std::pow ( Z , 2.39 ) + 1.55468*1e-6*std::pow ( Z , 5.35 ) )*eV;

  return mass;
}

double  AtomicMass(double A, double Z)
{
  const double hydrogen_mass_excess = G4NucleiPropertiesTableAME03::GetMassExcess(1,1);
  const double neutron_mass_excess =  G4NucleiPropertiesTableAME03::GetMassExcess(0,1);

  double mass =
      (A-Z)*neutron_mass_excess + Z*hydrogen_mass_excess - BindingEnergy(A,Z) + A*amu_c2;

  return mass;
}

double  BindingEnergy(double A, double Z)
{
  //
  // Weitzsaecker's Mass formula
  //
  int Npairing = int(A-Z)%2;                  // pairing
  int Zpairing = int(Z)%2;
  double binding =
      - 15.67*A                           // nuclear volume
      + 17.23*std::pow(A,2./3.)                // surface energy
      + 93.15*((A/2.-Z)*(A/2.-Z))/A       // asymmetry
      + 0.6984523*Z*Z*std::pow(A,-1./3.);      // coulomb
  if( Npairing == Zpairing ) binding += (Npairing+Zpairing-1) * 12.0 / std::sqrt(A);  // pairing

  return -binding;
}

int DoppBroad(string inFileName, string outFileName, double prevTemp, double newTemp, double isoMassC2, bool ascii, bool log, std::ofstream* logFile)
{

    stringstream ss;
    #if Timer>=2
        std::clock_t loopStart;
        loopStart = std::clock();
        double loopdur;
    #endif
    #if Timer>=3
        std::clock_t start;
        double duration;
    #endif

    GetDataStream(inFileName, ss, log, logFile);

	// Find the temperature difference between the temperature the cross
	// section was evaluated at versus the temperature of the material
	double tempDiff = newTemp - prevTemp;

	if (tempDiff<=0)
	{
        ss.str("");
        if(log)
            (*logFile) << "the new cross-section temperature is the same as the previous cross-section temperature";
        return 1;
    }


    int tableSize=0;
    ss >> tableSize;

    // Physics Vector
    int vType=0;
    ss >>  vType;

    // contents
    int siz=0;
    ss >> siz;
    if (ss.fail())
    { ss.str(""); return 1; }
    if (siz<=0)
    {
        if(log)
            (*logFile) << " Invalid vector size: " << siz << endl;
        ss.str("");
        return 1;
    }

    double tempEn, tempCS;
    double* prevEnVec = new double[siz];
    double* prevCSVec = new double[siz];

    for(int i = 0; i < siz ; i++)
    {
        tempEn = 0.;
        tempCS = 0.;
        ss >> tempEn >> tempCS;

        if (ss.fail())
        { ss.str(""); return 1; }

        prevEnVec[i] = tempEn;
        prevCSVec[i] = tempCS;
    }

     ss.str("");
     ss.clear();

     if(ascii)
     {
        ss.fill(' ');
        ss << std::setw(14) << std::right << tableSize << "\n" << std::setw(14) << vType << "\n" << std::setw(14) << siz << "\n";
        ss.precision(6);
        ss.setf(std::ios::scientific);
     }

     else
        ss << tableSize << "\n" << vType << "\n" << siz << "\n";


    // Declarations for averaging loop
    int initNumIL = int(std::max(100.,tempDiff/6.)), numIL/*, numMaxOL=3000*initNumIL*/;
    double nKEnerTrans;
    double aXsection;             // random xsec value
    double result;            // final averaged cross section
    double buffer;            // temporary value for previous average
    int counter;              // counts the number averaged xsec sets
    double nPx, nPy=0, nPz=0, isoPx, isoPy, isoPz, isoP2, isoEn /*isoGamma*/, a, x, y, z, p2, nEn, nVMag, sigma=std::sqrt(k_Boltzmann*tempDiff*isoMassC2);
    double isoMassC2Sq = isoMassC2*isoMassC2, nMassC2Sq = nMassC2*nMassC2, invIsoMassC2 = 1/isoMassC2, invNMassC2 = 1/nMassC2, invNVMag=0;
    #if Timer>=3
        double avgFindTime=0, /*prevAvgFindTime=0,*/ avgMathTime=0/*, prevAvgMathTime=0*/;
        double avgNumLoops=0/*, tempMathTime, tempFindTime*/;
    #endif

     #if Timer>=2
            (*logFile) << "Progress of " << inFileName << ":\n [";
            (*logFile).flush();
            int progressCount=1;
    #endif

    for(int i = 0; i < siz ; i++)
    {
        //prevAvgMathTime=avgMathTime;
        //prevAvgFindTime=avgFindTime;

        #if Timer>=3
            start = std::clock();
        #endif
        // Find the Doppler broadened cross section
        result = 0;
        buffer = 0;
        counter = 0;
        numIL = initNumIL;

        nEn = (prevEnVec[i]/1000000)+nMassC2;
        nPx = sqrt(pow(nEn,2)-nMassC2Sq);
        nVMag = nPx*invNMassC2;
        invNVMag = 1/nVMag;

        while(counter == 0 || (/*counter < numMaxOL && */std::abs(buffer-result/std::max(1,counter)) > 0.001*buffer))
        {
            if(counter)
                buffer = result/counter;

            while(counter < numIL)
            {
                isoPx = (CLHEP::RandGaussQ::shoot())*sigma;
                isoPy = (CLHEP::RandGaussQ::shoot())*sigma;
                isoPz = (CLHEP::RandGaussQ::shoot())*sigma;
                isoP2 = isoPx*isoPx+isoPy*isoPy+isoPz*isoPz;

                isoEn = std::sqrt( isoMassC2Sq+isoP2 );

                // this is only needed for the velocity correction when the isoKE is significant compared to the isoMass wich it is not unless temp >10^6K
                //isoGamma = isoEn/(isoMassC2);

            // Find the neutron in the rest frame of the nucleus
                a = ( (nPx*isoPx)/(isoEn+isoMassC2) - nEn )*invIsoMassC2;
                x = nPx+a*isoPx;
                y = nPy+a*isoPy;
                z = nPz+a*isoPz;
                p2 = x*x+y*y+z*z;
                nKEnerTrans = std::sqrt(nMassC2Sq+p2)-nMassC2;

                // Find the cross-section of the element at the energy of the
                // neutron IN THE REST FRAME OF THE NUCLEUS
                #if Timer>=3
                    avgMathTime+=double(std::clock() - start);
                    start = std::clock();
                #endif

                aXsection = findCS(nKEnerTrans, prevEnVec, prevCSVec, siz);

                #if Timer>=3
                    avgFindTime+=double(std::clock() - start);
                    start = std::clock();
                #endif

                // Velocity correction
                aXsection *= sqrt(pow((nVMag-isoPx*invIsoMassC2),2)+pow((isoPy*invIsoMassC2),2)+pow((isoPz*invIsoMassC2),2))*invNVMag;

                // Add cross section to results
                result += aXsection;

                #if Timer>=3
                    avgMathTime+=double(std::clock() - start);
                    start = std::clock();
                #endif

                counter++;
            }

            numIL +=numIL;

        }

        //tempFindTime = double((avgFindTime-prevAvgFindTime)/counter)/CLOCKS_PER_SEC;
        //tempMathTime = double((avgMathTime-prevAvgMathTime)/counter)/CLOCKS_PER_SEC;
        #if Timer>=3
            avgNumLoops += counter;
        #endif

        if (ascii)
        {
            if(!((i)%3) && i!=0 )
                ss << "\n";
            ss << std::setw(14) << prevEnVec[i] << std::setw(14) << (result/counter);
        }
        else
            ss << prevEnVec[i] << ' ' << (result/counter) << ' ';

        #if Timer>=2
            if(siz==1)
            {
                (*logFile).fill('#');
                (*logFile) << std::setw(60) << std::right << "";
                (*logFile).flush();
            }
            else if(progressCount==int(60*i/(siz-1)))
            {
                for(int l=0; l<(int(60*i/(siz-1))-int(60*(abs(i-1))/(siz-1))); l++)
                {
                    (*logFile) << "#";
                }
                (*logFile).flush();
                progressCount++;
            }
            if(i==siz-1)
            {
                loopdur = double(std::clock()-loopStart)/CLOCKS_PER_SEC;
                (*logFile) << "] Time Taken By File " << loopdur << "s\n" << endl;
            }


        #endif
    }

    #if Timer>=3
        avgFindTime /= avgNumLoops;
        avgMathTime /= avgNumLoops;
        avgNumLoops /= siz;

        double findperc, mathperc;
        findperc = avgFindTime*100/(avgFindTime+avgMathTime);
        mathperc = avgMathTime*100/(avgFindTime+avgMathTime);
    #endif

    ss << '\n';

    SetDataStream( outFileName , ss, ascii, log, logFile );

    #if Timer>=3
        if(log)
        {
            duration = (double(std::clock() - timebegin))/CLOCKS_PER_SEC;
            (*logFile) << "\n \n \n ### " << outFileName << " took " << duration <<"s to broaden and it had " << siz << " data points ###" << endl;
            duration = (double(avgMathTime))/CLOCKS_PER_SEC;
            (*logFile) << "### the average time it took to do the math per loop was " << duration << "s, this took " << mathperc << "% of the time ###" << endl;
            duration = (double(avgFindTime))/CLOCKS_PER_SEC;
            (*logFile) << "### the average time it took to find the CS per loop was " << duration << "s, this took " << findperc << "% of the time ###" << endl;
            (*logFile) << "### the average number of loops per data point was " << avgNumLoops << "### \n \n \n" << endl;
        }
    #endif
    ss.str("");
    delete [] prevCSVec;
    delete [] prevEnVec;
    return 0 ;
}

double findCS(double nKEnerTrans, const double* prevEnVec, const double* prevCSVec, int vecSize)
{
    int index=0;
    bool Back, Front;
    int step=int(vecSize/10);

    nKEnerTrans=nKEnerTrans*1000000.;

    if (nKEnerTrans<=prevEnVec[0])
    {
        return  prevCSVec[0];
    }
    else if (nKEnerTrans>=prevEnVec[vecSize-1])
    {
        return  prevCSVec[vecSize-1];
    }

    else
    {
        index = int((vecSize-1)*(nKEnerTrans-prevEnVec[0])
                    /(prevEnVec[vecSize-1]-prevEnVec[0]));

        while(step>1)
        {
            if(index<0)
                index=0;

            if((index+step)>int(vecSize-1))
                index=vecSize-1-step;

            (prevEnVec[index]<=nKEnerTrans)? Back=true : Back=false;
            (prevEnVec[index+step]>nKEnerTrans)? Front=true : Front=false;

            while (!Back)
            {
                index -= step;
                if(index<0)
                    index=0;
                (prevEnVec[index]<=nKEnerTrans)? Back=true : Back=false;
            }

            while (!Front)
            {
                index += step;
                if((index+step)>int(vecSize-1))
                    index=vecSize-1-step;
                (prevEnVec[index+step]>nKEnerTrans)? Front=true : Front=false;
            }

            if(prevEnVec[index]==nKEnerTrans)
                break;

            index = int((step)*(nKEnerTrans-prevEnVec[index])
                    /(prevEnVec[index+step]-prevEnVec[index])+index);

            step/=10;
        }

       while(prevEnVec[index]<nKEnerTrans)
       {
            index++;
       }

       while(prevEnVec[index]>nKEnerTrans)
       {
            index--;
       }

       if (prevEnVec[index]!=nKEnerTrans)
       {
            return ((prevCSVec[index+1]-prevCSVec[index])*(nKEnerTrans-prevEnVec[index])/(prevEnVec[index+1]-prevEnVec[index])+prevCSVec[index]);
       }

       return  prevCSVec[index];

    }

}

void GetDataStream( string filename , std::stringstream& ss, bool log, std::ofstream* logFile)
{
   string* data=NULL;
   std::ifstream* in=NULL;
   //string compfilename(filename);

   if(filename.substr((filename.length()-2),2)==".z")
   {
        in = new std::ifstream ( filename.c_str() , std::ios::binary | std::ios::ate );
   }

   if ( in!=NULL && in->good() )
   {
// Use the compressed file
      uLongf file_size = (uLongf)(in->tellg());
      in->seekg( 0 , std::ios::beg );
      Bytef* compdata = new Bytef[ file_size ];

      while ( *in )
      {
         in->read( (char*)compdata , file_size );
      }

      uLongf complen = (uLongf) ( file_size*4 );
      Bytef* uncompdata = new Bytef[complen];

      while ( Z_OK != uncompress ( uncompdata , &complen , compdata , file_size ) )
      {
         delete[] uncompdata;
         complen *= 2;
         uncompdata = new Bytef[complen];
      }
      delete [] compdata;
      //                                 Now "complen" has uncomplessed size
      data = new string ( (char*)uncompdata , (long)complen );
      delete [] uncompdata;
   }
   else {
// Use regular text file
      std::ifstream thefData( filename.c_str() , std::ios::in | std::ios::ate );
      if ( thefData.good() )
      {
         int file_size = thefData.tellg();
         thefData.seekg( 0 , std::ios::beg );
         char* filedata;
         if(file_size>0)
            filedata = new char[ file_size ];
         else
         {
            if(in!=NULL)
            {
                in->close();
                delete in;
            }
            return;
         }
         while ( thefData )
         {
            thefData.read( filedata , file_size );
         }
         thefData.close();
         data = new string ( filedata , file_size );
         delete [] filedata;
      }
      else
      {
// found no data file
//                 set error bit to the stream
         ss.setstate( std::ios::badbit );
         if(log)
            (*logFile) << endl << "### failed to open ascii file " << filename << " ###" << endl;
      }
   }
   if (data != NULL)
   {
        ss.str(*data);
        if(data->back()!='\n')
            ss << "\n";
        ss.seekg( 0 , std::ios::beg );
    }

   if(in!=NULL)
   {
        in->close();
        delete in;
   }

   delete data;
}


void SetDataStream( string filename , std::stringstream& ss, bool ascii, bool log, std::ofstream* logFile )
{
    //bool cond=true;
   if (!ascii)
   {
        string compfilename(filename);

        if(compfilename.back()!='z')
            compfilename += ".z";

       std::ofstream* out = new std::ofstream ( compfilename.c_str() , std::ios::binary | std::ios::trunc);
       if ( ss.good() )
       {
       //
    // Create the compressed file
          ss.seekg( 0 , std::ios::end );
          uLongf file_size = (uLongf)(ss.tellg());
          ss.seekg( 0 , std::ios::beg );
          Bytef* uncompdata = new Bytef[ file_size ];

          while ( ss ) {
              ss.read( (char*)uncompdata , file_size );
          }

          uLongf complen = compressBound(file_size);

          Bytef* compdata = new Bytef[complen];

          if ( Z_OK == compress ( compdata , &complen , uncompdata , file_size ) )
          {
            out->write((char*)compdata, (long)complen);
            if (out->fail())
            {
                if(log)
                    (*logFile) << endl << "writing the compressed data to the output file " << compfilename << " failed" << endl
                    << " may not have permission to delete an older version of the file" << endl;
            }
          }
          else
          {
            if(log)
                (*logFile) << endl << "compressing the data failed" << endl;
          }

          delete [] uncompdata;
          delete [] compdata;
       }
       else
       {
            if(log)
                (*logFile) << endl << "### failed to write to binary file ###" << endl;
       }

       out->close(); delete out;
   }
   else
   {
// Use regular text file
    string compfilename(filename);

    if(compfilename.substr((compfilename.length()-2),2)==".z")
    {
        compfilename.pop_back();
        compfilename.pop_back();
    }

      std::ofstream out( compfilename.c_str() , std::ios::out | std::ios::trunc );
      if ( ss.good() )
      {
         ss.seekg( 0 , std::ios::end );
         int file_size = ss.tellg();
         ss.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( ss ) {
            ss.read( filedata , file_size );
            if(!file_size)
            {
                if(log)
                    (*logFile) << "\n #### Error the size of the stringstream is invalid ###" << endl;
                break;
            }
         }
         out.write(filedata, file_size);
         if (out.fail())
        {
            if(log)
                (*logFile) << endl << "writing the ascii data to the output file " << compfilename << " failed" << endl
                 << " may not have permission to delete an older version of the file" << endl;
        }
         out.close();
         delete [] filedata;
      }
      else
      {
// found no data file
//                 set error bit to the stream
         ss.setstate( std::ios::badbit );

         if(log)
            (*logFile) << endl << "### failed to write to ascii file " << compfilename << " ###" << endl;
      }
   }
   ss.str("");
}
