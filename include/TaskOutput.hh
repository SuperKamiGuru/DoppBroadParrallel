#ifndef TASKOUTPUT_HH
#define TASKOUTPUT_HH

#include "MyString.h"
#include <string.h>

//MSH_include_begin
#include "MarshaledMyString.h"
//MSH_include_end

//MSH_BEGIN
class TaskOutput
{
    public:
        TaskOutput()
        {
            duration=0.;
            success=0;
            inDirSize=0;
            outDirSize=0;
            fileSize=0;
        }
        virtual ~TaskOutput()
        {

        }
        TaskOutput(const char *in1, const char *in2, const char *in3, double in4, double in5, double in6, double in7, double in8,
                    bool in9, bool in10, bool in11, int in12, int in13, double in14, bool in15)
        {
            inDirName = *(new MyString(in1));
            outDirName = *(new MyString(in2));
            fileName = *(new MyString(in3));

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
            duration=in14;
            success=in15;
        }

        TaskOutput(const char *in1, const char *in2, const char *in3, double in4, double in5, bool in6, bool in7, bool in8, double in9, bool in10)
        : totalFileSize2(0), sumFileSize2(0), sumDuration(0), index(0), totalNumFiles(0)
        {
            inDirName = *(new MyString(in1));
            outDirName = *(new MyString(in2));
            fileName = *(new MyString(in3));

            prevTemp=in4;
            newTemp=in5;
            ascii=in6;
            log=in7;
            overWrite=in8;
            duration=in9;
            success=in10;
        }

        TaskOutput& operator=(const TaskOutput& other)
        {
            inDirName = *(new MyString(other.inDirName));
            outDirName = *(new MyString(other.outDirName));
            fileName = *(new MyString(other.fileName));

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
            duration=other.duration;
            success=other.success;
            return *this;
        }

        MyString inDirName; /* MSH: predefined */
        MyString outDirName; /* MSH: predefined */
        MyString fileName; /* MSH: predefined */

        double prevTemp; //MSH: primitive
        double newTemp; //MSH: primitive
        double totalFileSize2; //MSH: primitive
        double sumFileSize2; //MSH: primitive
        double sumDuration; //MSH: primitive
        bool ascii; //MSH: primitive
        bool log; //MSH: primitive
        bool overWrite; //MSH: primitive
        int inDirSize; //MSH: primitive
        int outDirSize; //MSH: primitive
        int fileSize; //MSH: primitive
        int index; //MSH: primitive
        int totalNumFiles; //MSH: primitive
        double duration;//MSH: primitive
        bool success;//MSH: primitive
};
//MSH_END
#endif // TASK_OUTPUT_HH

