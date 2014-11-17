#ifndef TASKINPUT_HH
#define TASKINPUT_HH

//MSH_BEGIN
class TaskInput
{
    public:
        TaskInput()
        {
            inDirName="";
            outDirName="";
            fileName="";
            prevTemp=0.;
            newTemp=0.;
            totalFileSize2=0.;
            sumFileSize2=0.;
            sumDuration=0.;
            ascii=0;
            log=0;
            overWrite=0;
            index=0;
            totalNumFiles=0;
        }

        TaskInput(string in1, string in2, string in3, double in4, double in5, double in6, double in7, double in8,
                    bool in9, bool in10, bool in11, int in12, int in13)
        {
            inDirName=in1;
            outDirName=in2;
            fileName=in3;
            prevTemp=in4;
            newTemp=in5;
            totalFileSize2=in6;
            sumFileSize2=in7;
            sumDuration=in8;
            ascii=in9;
            log=in10;
            overWrite=in11;
            index=in12;
            totalNumFiles=in13;
        }

        TaskInput(string in1, string in2, string in3, double in4, double in5, bool in6, bool in7, bool in8)
        : totalFileSize2(0), sumFileSize2(0), sumDuration(0), index(0), totalNumFiles(0)
        {
            inDirName=in1;
            outDirName=in2;
            fileName=in3;
            prevTemp=in4;
            newTemp=in5;
            ascii=in6;
            log=in7;
            overWrite=in8;
        }

        TaskInput& operator=(const TaskInput& other)
        {
            inDirName=other.inDirName;
            outDirName=other.outDirName;
            fileName=other.fileName;
            prevTemp=other.prevTemp;
            newTemp=other.newTemp;
            totalFileSize2=other.totalFileSize2;
            sumFileSize2=other.sumFileSize2;
            sumDuration=other.sumDuration;
            ascii=other.ascii;
            log=other.log;
            overWrite=other.overWrite;
            index=other.index;
            totalNumFiles=other.totalNumFiles;

            return *this;
        }

        string inDirName; //MSH: primitive
        string outDirName; //MSH: primitive
        string fileName; //MSH: primitive
        double prevTemp; //MSH: primitive
        double newTemp; //MSH: primitive
        double totalFileSize2; //MSH: primitive
        double sumFileSize2; //MSH: primitive
        double sumDuration; //MSH: primitive
        bool ascii; //MSH: primitive
        bool log; //MSH: primitive
        bool overWrite; //MSH: primitive
        int index; //MSH: primitive
        int totalNumFiles; //MSH: primitive
};
//MSH_END

#endif // TASK_INPUT_HH
